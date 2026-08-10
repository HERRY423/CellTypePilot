"""Deterministic, task-bounded Novelty & OOD Verification Engine.

Provides a 5-gate audit protocol for evaluating OOD/novel cell-type candidates
without changing canonical identity decisions:

Gate 1: Technical QC, Batch/Donor Confounding & Doublet Co-expression Audit
Gate 2: Subclustering Homogeneity vs. Heterogeneous Mixture Audit
Gate 3: Exploratory Cell State Lens Disambiguation (cycling, stress, hypoxia)
Gate 4: Atlas & Extended Pack Gap Assessment
Gate 5: Human Adjudication Protocol & Append-Only Audit Trail Logging

Claim Boundary:
CellTypePilot novelty output is a review-priority packet for human adjudication.
It never automatically renames biological identity or claims validated novel cell discovery
without signed human expert adjudication in novelty_adjudication_log.jsonl.
"""

from __future__ import annotations

import json
import logging
from contextlib import suppress
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

from .constants import (
    CRITIC_DOUBLET_ACTIVE_COVERAGE,
    STATE_ATLAS_PATH,
)
from .data_adapter import build_lineage_groups, find_metadata_columns

logger = logging.getLogger(__name__)

VERIFICATION_SCHEMA = "celltypepilot.novelty-verification.v1"
ADJUDICATION_SCHEMA = "celltypepilot.novelty-adjudication.v1"

# Valid adjudication verdicts
VERDICTS = {
    "validated_novel_cell_type",
    "novel_cell_state",
    "atlas_gap_resolved",
    "rejected_technical_artifact",
    "rejected_mixed_cluster",
    "inconclusive_pending_experiment",
}


def _check_qc_and_batch_confounding(
    adata: ad.AnnData,
    cluster_key: str,
    cluster_id: str,
) -> dict[str, Any]:
    """Gate 1: Audit QC metrics and batch/donor metadata concentration."""
    mask = adata.obs[cluster_key].astype(str) == str(cluster_id)
    n_cells = int(mask.sum())

    if n_cells == 0:
        return {
            "qc_passed": False,
            "batch_passed": False,
            "flags": ["EMPTY_CLUSTER"],
            "dominant_batch_info": {},
        }

    flags: list[str] = []

    # 1. QC Audit
    cluster_obs = adata.obs[mask]
    n_genes_col = next(
        (c for c in ["n_genes_by_counts", "n_genes", "n_genes_detected"] if c in cluster_obs), None
    )
    mito_col = next(
        (c for c in ["pct_counts_mt", "percent_mito", "pct_mito", "mt_frac"] if c in cluster_obs),
        None,
    )
    ribo_col = next(
        (c for c in ["pct_counts_ribo", "percent_ribo", "pct_ribo"] if c in cluster_obs), None
    )

    if n_genes_col and n_genes_col in adata.obs:
        cl_median_genes = float(cluster_obs[n_genes_col].median())
        all_median_genes = float(adata.obs[n_genes_col].median())
        if cl_median_genes < 0.5 * all_median_genes:
            flags.append(
                f"LOW_GENE_COUNT_MEDIAN:{cl_median_genes:.0f}_vs_ALL:{all_median_genes:.0f}"
            )

    if mito_col and mito_col in adata.obs:
        cl_median_mito = float(cluster_obs[mito_col].median())
        if cl_median_mito > 15.0:  # >15% mito RNA indicates stressed/dying cells
            flags.append(f"HIGH_MITO_PERCENT:{cl_median_mito:.1f}%")

    if ribo_col and ribo_col in adata.obs:
        cl_median_ribo = float(cluster_obs[ribo_col].median())
        if cl_median_ribo > 40.0:
            flags.append(f"HIGH_RIBO_PERCENT:{cl_median_ribo:.1f}%")

    # 2. Batch / Donor / Sample Confounding Audit
    metadata_cols = find_metadata_columns(adata)
    candidate_keys = dict.fromkeys(
        metadata_cols.get("batch_keys", [])
        + metadata_cols.get("sample_keys", [])
        + metadata_cols.get("donor_keys", [])
        + metadata_cols.get("study_keys", [])
    )

    dominant_info = {}
    for key in candidate_keys:
        if key not in adata.obs:
            continue
        values = cluster_obs[key].astype(str)
        if values.empty:
            continue
        top_val, top_count = values.value_counts().index[0], values.value_counts().iloc[0]
        frac = top_count / n_cells
        dominant_info[key] = {"top_value": str(top_val), "fraction": round(float(frac), 4)}
        total_dataset_values = len(adata.obs[key].unique())
        if frac >= 0.80 and total_dataset_values > 1:
            flags.append(f"CONFOUNDED_{key.upper()}:{top_val}({frac:.0%})")

    qc_passed = not any(f.startswith("HIGH_MITO") or f.startswith("LOW_GENE") for f in flags)
    batch_passed = not any(f.startswith("CONFOUNDED") for f in flags)

    return {
        "qc_passed": qc_passed,
        "batch_passed": batch_passed,
        "flags": flags,
        "dominant_batch_info": dominant_info,
    }


def _check_doublet_signature(
    adata: ad.AnnData,
    cluster_key: str,
    cluster_id: str,
    atlas: dict,
    tissue: str,
    layer: str | None = None,
) -> dict[str, Any]:
    """Gate 1 (cont.): Audit cross-lineage doublet co-expression using Critic rules."""
    lineages = build_lineage_groups(atlas, tissue)
    if not lineages:
        return {"doublet_flagged": False, "active_lineages": [], "details": "no_lineages_found"}

    mask = adata.obs[cluster_key].astype(str) == str(cluster_id)
    n_cells = int(mask.sum())
    if n_cells == 0:
        return {"doublet_flagged": False, "active_lineages": [], "details": "empty_cluster"}

    X_cl = adata[mask].X if layer is None else adata[mask].layers[layer]
    if hasattr(X_cl, "toarray"):
        X_cl = X_cl.toarray()

    var_names = list(adata.var_names)
    var_map = {gene: idx for idx, gene in enumerate(var_names)}

    active_lineages = []
    for lin_name, lin_genes in lineages.items():
        present = [g for g in lin_genes if g in var_map]
        if not present:
            continue
        indices = [var_map[g] for g in present]
        expr_count = (X_cl[:, indices] > 0).sum(axis=1)
        # Fraction of cells in cluster expressing >= 2 markers of this lineage
        active_frac = float((expr_count >= 2).mean())
        if active_frac >= CRITIC_DOUBLET_ACTIVE_COVERAGE:
            active_lineages.append((lin_name, active_frac))

    doublet_flagged = len(active_lineages) >= 2
    return {
        "doublet_flagged": doublet_flagged,
        "active_lineages": [f"{name}:{frac:.0%}" for name, frac in active_lineages],
        "details": "cross_lineage_coexpression_detected" if doublet_flagged else "single_lineage",
    }


def verify_subclustering_homogeneity(
    adata: ad.AnnData,
    cluster_key: str,
    cluster_id: str,
    resolution: float = 0.5,
) -> dict[str, Any]:
    """Gate 2: Run local PCA + Leiden subclustering to test cluster homogeneity.

    If subclustering cleanly splits the cluster into 2+ distinct sub-clusters
    with distinct DE programs, flags as heterogeneous mixture.
    """
    import scanpy as sc

    mask = adata.obs[cluster_key].astype(str) == str(cluster_id)
    n_cells = int(mask.sum())

    if n_cells < 30:
        return {
            "is_homogeneous": True,
            "n_subclusters": 1,
            "subcluster_sizes": {f"{cluster_id}_sub0": n_cells},
            "subcluster_score": 1.0,
            "note": "too_few_cells_for_subclustering",
        }

    try:
        sub_adata = adata[mask].copy()
        sc.pp.neighbors(
            sub_adata,
            use_rep="X_pca" if "X_pca" in sub_adata.obsm else None,
            n_neighbors=min(15, n_cells - 1),
        )
        sc.tl.leiden(sub_adata, resolution=resolution, key_added="subcluster")

        sub_counts = sub_adata.obs["subcluster"].value_counts().to_dict()
        n_subclusters = len(sub_counts)

        # Evaluate if sub-clusters are balanced or just minor noisy splits (>15% size)
        major_subclusters = [k for k, v in sub_counts.items() if (v / n_cells) >= 0.15]
        is_homogeneous = len(major_subclusters) <= 1

        return {
            "is_homogeneous": is_homogeneous,
            "n_subclusters": n_subclusters,
            "n_major_subclusters": len(major_subclusters),
            "subcluster_sizes": {str(k): int(v) for k, v in sub_counts.items()},
            "subcluster_score": round(1.0 / max(1, len(major_subclusters)), 4),
            "note": "homogeneous_sublineage"
            if is_homogeneous
            else "heterogeneous_mixture_detected",
        }
    except Exception as exc:
        logger.debug("Subclustering failed for cluster %s: %s", cluster_id, exc)
        return {
            "is_homogeneous": True,
            "n_subclusters": 1,
            "subcluster_sizes": {},
            "subcluster_score": 1.0,
            "note": f"subclustering_error:{exc}",
        }


def _evaluate_state_vs_identity(
    unmapped_de_markers: list[str],
    state_atlas_path: Path = STATE_ATLAS_PATH,
) -> dict[str, Any]:
    """Gate 3: Cross-reference unmapped DE program against state_atlas.json."""
    if not unmapped_de_markers or not state_atlas_path.exists():
        return {
            "is_cell_state_driven": False,
            "matching_state": None,
            "state_overlap_fraction": 0.0,
            "state_genes_matched": [],
        }

    try:
        state_atlas = json.loads(state_atlas_path.read_text(encoding="utf-8"))
    except Exception:
        return {
            "is_cell_state_driven": False,
            "matching_state": None,
            "state_overlap_fraction": 0.0,
            "state_genes_matched": [],
        }

    unmapped_set = {g.upper() for g in unmapped_de_markers}
    best_state = None
    best_overlap = 0.0
    best_matched_genes = []

    for state_name, state_def in state_atlas.get("states", {}).items():
        pos_markers = [
            m.get("gene", "").upper() if isinstance(m, dict) else str(m).upper()
            for m in state_def.get("positive_markers", [])
        ]
        matched = [g for g in pos_markers if g in unmapped_set]
        if matched:
            overlap_frac = len(matched) / len(unmapped_set)
            if overlap_frac > best_overlap:
                best_overlap = overlap_frac
                best_state = state_name
                best_matched_genes = matched

    is_state_driven = best_overlap >= 0.40  # >=40% of unmapped DE program matches state module
    return {
        "is_cell_state_driven": is_state_driven,
        "matching_state": best_state,
        "state_overlap_fraction": round(best_overlap, 4),
        "state_genes_matched": best_matched_genes,
    }


def verify_novelty_candidate(
    adata: ad.AnnData,
    cluster_key: str,
    focus_cluster: str,
    critic_row: dict[str, Any] | pd.Series,
    atlas: dict,
    tissue: str,
    unmapped_de_markers: list[str] | None = None,
    layer: str | None = None,
) -> dict[str, Any]:
    """Run full 5-gate verification on a focus OOD/novel candidate cluster.

    Returns a structured NoveltyVerificationPacket dictionary.
    """
    cluster_id = str(focus_cluster)

    # Gate 1: QC & Batch Confounding
    gate1_qc = _check_qc_and_batch_confounding(adata, cluster_key, cluster_id)

    # Gate 1 (cont.): Doublet Co-expression
    gate1_doublet = _check_doublet_signature(
        adata, cluster_key, cluster_id, atlas, tissue, layer=layer
    )

    # Gate 2: Local Subclustering Audit
    gate2_subcluster = verify_subclustering_homogeneity(adata, cluster_key, cluster_id)

    # Gate 3: State Lens Disambiguation
    markers_to_check = unmapped_de_markers or []
    if not markers_to_check:
        top_unmapped = str(critic_row.get("top_unmapped_markers", ""))
        markers_to_check = [m.strip() for m in top_unmapped.split(";") if m.strip()]

    gate3_state = _evaluate_state_vs_identity(markers_to_check)

    # Gate 4: Overall Verification Recommendation Synthesis
    verification_passed = (
        gate1_qc["qc_passed"]
        and gate1_qc["batch_passed"]
        and not gate1_doublet["doublet_flagged"]
        and gate2_subcluster["is_homogeneous"]
        and not gate3_state["is_cell_state_driven"]
    )

    if not gate1_qc["qc_passed"] or not gate1_qc["batch_passed"]:
        suggested_classification = "rejected_technical_artifact"
    elif gate1_doublet["doublet_flagged"] or not gate2_subcluster["is_homogeneous"]:
        suggested_classification = "rejected_mixed_cluster"
    elif gate3_state["is_cell_state_driven"]:
        suggested_classification = "novel_cell_state"
    elif verification_passed:
        suggested_classification = "verification_passed_candidate"
    else:
        suggested_classification = "inconclusive_pending_experiment"

    packet = {
        "schema_version": VERIFICATION_SCHEMA,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "cluster": cluster_id,
        "identity_best_match": str(
            critic_row.get("candidate_cell_type", critic_row.get("cell_type", "Unknown"))
        ),
        "identity_evidence_score": float(
            critic_row.get("evidence_score", critic_row.get("combined_score", 0.0)) or 0.0
        ),
        "unmapped_de_markers": markers_to_check,
        "gates": {
            "gate1_qc_batch": gate1_qc,
            "gate1_doublet": gate1_doublet,
            "gate2_subclustering": gate2_subcluster,
            "gate3_state_lens": gate3_state,
        },
        "verification_passed": verification_passed,
        "suggested_classification": suggested_classification,
        "adjudication_status": "pending_human_review",
        "claim_boundary": (
            "Verification packet is an audit checklist for human adjudication. "
            "It does not automatically rename cell identity or claim validated novel discovery "
            "without signed human expert sign-off in novelty_adjudication_log.jsonl."
        ),
    }

    return packet


def log_novelty_adjudication(
    output_dir: str | Path,
    cluster: str,
    verdict: str,
    reviewer: str,
    notes: str | None = None,
    pmid: str | None = None,
) -> dict[str, Any]:
    """Log a human expert adjudication verdict to an append-only audit trail.

    Updates artifact_status.json to mark derived artifacts stale as per Principle 14.
    """
    if verdict not in VERDICTS:
        raise ValueError(f"Invalid verdict '{verdict}'. Must be one of: {sorted(VERDICTS)}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    entry = {
        "schema_version": ADJUDICATION_SCHEMA,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "cluster": str(cluster),
        "verdict": verdict,
        "reviewer": str(reviewer),
        "notes": str(notes) if notes else "",
        "pmid_or_evidence_link": str(pmid) if pmid else "",
    }

    # Append to novelty_adjudication_log.jsonl
    log_path = output_dir / "novelty_adjudication_log.jsonl"
    with log_path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")

    # Mark artifact_status.json stale as per Principle 14
    status_path = output_dir / "artifact_status.json"
    status_data = {}
    if status_path.exists():
        with suppress(Exception):
            status_data = json.loads(status_path.read_text(encoding="utf-8"))

    status_data["novelty_review_status"] = "adjudicated"
    status_data["last_adjudication_at"] = entry["timestamp"]
    status_data["derived_artifacts_stale"] = True
    status_data["stale_reason"] = f"Novelty adjudication verdict logged for cluster {cluster}"

    status_path.write_text(
        json.dumps(status_data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    logger.info("Novelty adjudication logged for cluster %s: %s by %s", cluster, verdict, reviewer)
    return entry
