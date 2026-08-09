"""Fail-closed OOD / novel-cell-type candidate detection.

The detector does not name new biology. It adds a third, independent review
axis that asks whether a cluster looks insufficiently supported by the current
atlas/reference while still carrying a distinctive DE marker program worth
human review.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from .constants import MARKER_FC_THRESHOLD, MARKER_FDR_THRESHOLD
from .data_adapter import find_metadata_columns

NOVELTY_SCHEMA = "celltypepilot.novelty-ood.v1"
NOVELTY_SCORE_SEMANTICS = "review_priority_score_not_discovery_probability"
NOVELTY_CLAIM_BOUNDARY = (
    "Novelty/OOD output is an automated review-priority signal. It is not a validated new "
    "cell-type discovery, not a cell ontology term assignment, and not a batch-robustness claim."
)

KNOWN_SUPPORTED = "known_supported"
ATLAS_GAP_CANDIDATE = "atlas_gap_candidate"
OOD_NOVEL_CANDIDATE = "ood_novel_candidate"
REVIEW_ARTIFACT_OR_MIXED = "review_artifact_or_mixed"
INSUFFICIENT_SIGNAL = "insufficient_signal"


def _atlas_marker_genes(marker_definitions: Mapping[str, Mapping[str, Any]]) -> set[str]:
    genes: set[str] = set()
    for definition in marker_definitions.values():
        for key in ("positive_markers", "negative_markers"):
            genes.update(str(gene) for gene in definition.get(key, []) if gene)
    return genes


def _get_de_results(
    adata: ad.AnnData,
    cluster_key: str,
    de_results: dict[str, pd.DataFrame] | None = None,
) -> dict[str, pd.DataFrame]:
    if de_results is not None:
        return {str(key): value.copy() for key, value in de_results.items()}
    try:
        from .marker_scorer import _extract_de_results

        return _extract_de_results(adata, cluster_key)
    except Exception:
        return {}


def _top_unmapped_program(
    cluster_de: pd.DataFrame,
    atlas_genes: set[str],
    max_markers: int = 12,
) -> dict[str, Any]:
    if cluster_de.empty:
        return {
            "n_unmapped_de_markers": 0,
            "top_unmapped_markers": "",
            "marker_program_strength": 0.0,
        }
    frame = cluster_de.copy()
    frame["gene"] = frame["gene"].astype(str)
    frame["logfoldchange"] = pd.to_numeric(frame["logfoldchange"], errors="coerce").fillna(0.0)
    frame["pval_adj"] = pd.to_numeric(frame["pval_adj"], errors="coerce").fillna(1.0)
    supported = frame[
        (frame["logfoldchange"] >= MARKER_FC_THRESHOLD)
        & (frame["pval_adj"] <= MARKER_FDR_THRESHOLD)
        & (~frame["gene"].isin(atlas_genes))
    ].sort_values(["logfoldchange", "pval_adj"], ascending=[False, True])
    if supported.empty:
        return {
            "n_unmapped_de_markers": 0,
            "top_unmapped_markers": "",
            "marker_program_strength": 0.0,
        }
    top = supported.head(max_markers)
    strength = float(np.clip(top["logfoldchange"].mean() / 2.0, 0.0, 1.0))
    return {
        "n_unmapped_de_markers": int(len(supported)),
        "top_unmapped_markers": ";".join(top["gene"].tolist()),
        "marker_program_strength": round(strength, 4),
    }


def _reference_summary(ref_scores: pd.DataFrame | None, cluster: str) -> dict[str, Any]:
    if ref_scores is None or ref_scores.empty:
        return {
            "reference_top_match": "",
            "reference_top_score": np.nan,
            "reference_margin": np.nan,
            "reference_entropy": np.nan,
        }
    rows = ref_scores[ref_scores["cluster"].astype(str) == str(cluster)].sort_values("ref_rank")
    if rows.empty:
        return {
            "reference_top_match": "",
            "reference_top_score": np.nan,
            "reference_margin": np.nan,
            "reference_entropy": np.nan,
        }
    scores = pd.to_numeric(rows["ref_score"], errors="coerce").fillna(0.0).to_numpy(dtype=float)
    top_score = float(scores[0]) if len(scores) else np.nan
    second = float(scores[1]) if len(scores) > 1 else 0.0
    top5 = scores[: min(5, len(scores))]
    norm = top5 / (top5.sum() + 1e-12)
    entropy = float(-np.sum(norm * np.log2(norm + 1e-12))) if len(norm) else np.nan
    return {
        "reference_top_match": str(rows.iloc[0]["cell_type"]),
        "reference_top_score": round(top_score, 4),
        "reference_margin": round(top_score - second, 4),
        "reference_entropy": round(entropy, 4),
    }


def _cluster_metadata_warnings(adata: ad.AnnData, cluster_key: str, cluster: str) -> list[str]:
    metadata = find_metadata_columns(adata)
    candidate_keys = (
        metadata.get("batch_keys", [])
        + metadata.get("sample_keys", [])
        + metadata.get("donor_keys", [])
        + metadata.get("study_keys", [])
    )
    cluster_mask = adata.obs[cluster_key].astype(str) == str(cluster)
    warnings: list[str] = []
    for key in dict.fromkeys(candidate_keys):
        if key not in adata.obs:
            continue
        values = adata.obs.loc[cluster_mask, key].astype(str)
        if values.empty:
            continue
        dominant_fraction = float(values.value_counts(normalize=True).iloc[0])
        if dominant_fraction >= 0.8:
            warnings.append(f"{key}_enriched:{dominant_fraction:.0%}")
    if not warnings:
        warnings.append("batch_sample_enrichment_not_flagged")
    return warnings


def score_novelty_candidates(
    adata: ad.AnnData,
    cluster_key: str,
    critic_results: pd.DataFrame,
    marker_definitions: Mapping[str, Mapping[str, Any]],
    ref_scores: pd.DataFrame | None = None,
    de_results: dict[str, pd.DataFrame] | None = None,
    min_unmapped_markers: int = 3,
) -> pd.DataFrame:
    """Score clusters for OOD/novel review candidacy without changing identity calls."""
    atlas_genes = _atlas_marker_genes(marker_definitions)
    de_by_cluster = _get_de_results(adata, cluster_key, de_results)
    rows: list[dict[str, Any]] = []

    for _, row in critic_results.iterrows():
        cluster = str(row.get("cluster", ""))
        flags = str(row.get("critic_flags", ""))
        decision = str(row.get("decision", "accepted"))
        best_score = float(row.get("evidence_score", row.get("combined_score", 0.0)) or 0.0)
        best_overlap = float(row.get("pct_overlap", 0.0) or 0.0)
        neg_conflict = float(row.get("neg_conflict", 0.0) or 0.0)
        cluster_de = de_by_cluster.get(cluster, pd.DataFrame())
        program = _top_unmapped_program(cluster_de, atlas_genes)
        ref = _reference_summary(ref_scores, cluster)
        ref_entropy = float(ref["reference_entropy"]) if pd.notna(ref["reference_entropy"]) else 0.0
        ref_margin = float(ref["reference_margin"]) if pd.notna(ref["reference_margin"]) else 0.0

        low_identity_support = (
            decision == "abstain"
            or best_score < 0.45
            or best_overlap < 0.30
            or any(token in flags for token in ("LOW_EVIDENCE", "LOW_DE_SUPPORT", "PARTIAL"))
        )
        disagreement = any(token in flags for token in ("ENSEMBLE_DISAGREEMENT", "WEAK_REFERENCE"))
        diffuse_reference = bool(ref_entropy >= 1.5 or (ref["reference_top_match"] and ref_margin < 0.15))
        artifact_or_mixed = any(
            token in flags
            for token in ("POSSIBLE_DOUBLET", "NEG_MARKER_CONFLICT", "ONTOLOGY_MISMATCH")
        ) or neg_conflict >= 0.2
        has_distinctive_unmapped_program = (
            program["n_unmapped_de_markers"] >= min_unmapped_markers
            and program["marker_program_strength"] >= 0.25
        )

        identity_gap = max(0.0, 1.0 - best_score)
        support_gap = max(0.0, 1.0 - best_overlap)
        novelty_score = float(
            np.clip(
                0.35 * identity_gap
                + 0.25 * support_gap
                + 0.25 * program["marker_program_strength"]
                + 0.10 * float(disagreement or diffuse_reference)
                + 0.05 * float(decision == "abstain"),
                0.0,
                1.0,
            )
        )

        alternatives = []
        if artifact_or_mixed:
            alternatives.append("doublet_mixed_or_marker_conflict_possible")
        if disagreement:
            alternatives.append("marker_reference_disagreement")
        if diffuse_reference:
            alternatives.append("diffuse_or_low_margin_reference_match")
        alternatives.extend(_cluster_metadata_warnings(adata, cluster_key, cluster))

        if artifact_or_mixed and low_identity_support:
            novelty_decision = REVIEW_ARTIFACT_OR_MIXED
            candidate_label = "Mixed/artifact review candidate"
        elif low_identity_support and has_distinctive_unmapped_program:
            novelty_decision = OOD_NOVEL_CANDIDATE
            candidate_label = "Novel/OOD candidate"
        elif has_distinctive_unmapped_program and (disagreement or diffuse_reference):
            novelty_decision = ATLAS_GAP_CANDIDATE
            candidate_label = "Atlas-gap candidate"
        elif not has_distinctive_unmapped_program and low_identity_support:
            novelty_decision = INSUFFICIENT_SIGNAL
            candidate_label = "Insufficient signal"
        else:
            novelty_decision = KNOWN_SUPPORTED
            candidate_label = "Known-supported or low-priority"

        rows.append(
            {
                "cluster": cluster,
                "novelty_decision": novelty_decision,
                "novelty_candidate_label": candidate_label,
                "novelty_score": round(novelty_score, 4),
                "novelty_score_semantics": NOVELTY_SCORE_SEMANTICS,
                "atlas_best_match": row.get("candidate_cell_type", row.get("cell_type", "")),
                "atlas_best_score": round(best_score, 4),
                "atlas_best_overlap": round(best_overlap, 4),
                "identity_decision": decision,
                "critic_flags": flags,
                **program,
                **ref,
                "alternative_explanations": ";".join(dict.fromkeys(alternatives)),
                "recommended_next_actions": _recommended_actions(novelty_decision),
                "claim_boundary": NOVELTY_CLAIM_BOUNDARY,
                "novelty_schema": NOVELTY_SCHEMA,
            }
        )

    return pd.DataFrame(rows)


def _recommended_actions(decision: str) -> str:
    if decision == OOD_NOVEL_CANDIDATE:
        return (
            "review_top_unmapped_markers; check_batch_sample_distribution; compare_external_atlas; "
            "subcluster; literature_or_experimental_validation_before_naming"
        )
    if decision == ATLAS_GAP_CANDIDATE:
        return (
            "review_marker_reference_disagreement; evaluate_as_atlas_curation_gap; "
            "add_reviewed_context_or_extension_pack_only_after evidence review"
        )
    if decision == REVIEW_ARTIFACT_OR_MIXED:
        return "check_doublets_qc_batch; subcluster; do_not_name_as_novel_until_artifact_excluded"
    if decision == INSUFFICIENT_SIGNAL:
        return "increase_resolution_or_cells; inspect_qc; keep_identity_unknown"
    return "routine_review"


def attach_novelty_results(
    critic_results: pd.DataFrame,
    novelty_results: pd.DataFrame,
) -> pd.DataFrame:
    """Attach novelty review columns without modifying canonical identity columns."""
    if novelty_results.empty:
        return critic_results.copy()
    columns = [
        "cluster",
        "novelty_decision",
        "novelty_candidate_label",
        "novelty_score",
        "novelty_score_semantics",
        "top_unmapped_markers",
        "n_unmapped_de_markers",
        "marker_program_strength",
        "alternative_explanations",
        "recommended_next_actions",
        "claim_boundary",
        "novelty_schema",
    ]
    output = critic_results.merge(
        novelty_results[[column for column in columns if column in novelty_results]],
        on="cluster",
        how="left",
        validate="one_to_one",
    )
    output["ood_novelty_signal"] = output["novelty_decision"].fillna("not_assessed")
    return output


def build_novelty_manifest(novelty_results: pd.DataFrame) -> dict[str, Any]:
    counts = (
        novelty_results["novelty_decision"].value_counts().to_dict()
        if not novelty_results.empty and "novelty_decision" in novelty_results
        else {}
    )
    return {
        "schema_version": NOVELTY_SCHEMA,
        "enabled": True,
        "score_semantics": NOVELTY_SCORE_SEMANTICS,
        "decision_counts": counts,
        "claim_boundary": NOVELTY_CLAIM_BOUNDARY,
        "identity_invariant": True,
        "required_validation_before_naming": [
            "artifact/doublet/QC review",
            "batch/sample distribution review",
            "external atlas or held-out reference comparison when available",
            "human expert review of top unmapped markers",
            "literature or experimental validation before declaring a new cell type",
        ],
    }
