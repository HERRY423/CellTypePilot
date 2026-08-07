"""Annotation Critic — independent review of cell type annotations.

The Critic is CellTypePilot's trust layer. It doesn't just report confidence scores;
it performs a skeptical, evidence-based review of each annotation.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .constants import (
    CONFIDENCE_HIGH,
    CONFIDENCE_LOW,
    CONFIDENCE_MEDIUM,
    CONFIDENCE_REVIEW,
    CRITIC_DOUBLET_COEXPR_THRESHOLD,
    CRITIC_LOW_COVERAGE_THRESHOLD,
    CRITIC_NEG_MARKER_PCT_THRESHOLD,
    ENSEMBLE_AGREEMENT_THRESHOLD,
    MARKER_PCT_THRESHOLD,
)
from .data_adapter import get_all_markers_for_tissue


def run_critic(
    adata: ad.AnnData,
    cluster_key: str,
    annotations: pd.DataFrame,
    atlas: dict,
    tissue: str,
    ensemble_info: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Run the full critic pipeline on annotations.

    For each cluster annotation, performs:
    1. Evidence sufficiency check
    2. Negative marker conflict check
    3. Doublet/mixed signal heuristic
    4. Ontology consistency check
    5. Ensemble agreement check (if ensemble_info provided)
    6. Structured confidence recalibration

    Args:
        adata: AnnData object
        cluster_key: Column in obs with cluster labels
        annotations: Per-cluster annotations (from marker_scorer summary)
        atlas: Marker atlas dict
        tissue: Tissue context
        ensemble_info: Optional ensemble scores DataFrame from ensemble_scorer

    Returns the annotations DataFrame with added critic columns:
        critic_flags, critic_evidence, critic_confidence, critic_notes
    """
    markers = get_all_markers_for_tissue(atlas, tissue)
    results = annotations.copy()

    # Build ensemble lookup if provided
    ensemble_lookup = {}
    if ensemble_info is not None and not ensemble_info.empty:
        for _, row in ensemble_info.iterrows():
            cl = str(row.get("cluster", ""))
            if cl not in ensemble_lookup:
                ensemble_lookup[cl] = {
                    "marker_score": row.get("marker_score", 0),
                    "ref_score": row.get("ref_score", 0),
                    "agreement": row.get("agreement", True),
                    "source": row.get("source", ""),
                }

    critic_flags_list = []
    critic_evidence_list = []
    critic_confidence_list = []
    critic_notes_list = []

    for _, row in results.iterrows():
        cluster = row["cluster"]
        cell_type = row["cell_type"]
        flags = []
        evidence_parts = []
        notes = []

        ct_info = markers.get(cell_type, {})
        pos_markers = ct_info.get("positive_markers", [])
        neg_markers = ct_info.get("negative_markers", [])

        # 1. Evidence sufficiency check
        sufficiency_result = _check_evidence_sufficiency(adata, cluster, cluster_key, pos_markers)
        if sufficiency_result["flag"]:
            flags.append(sufficiency_result["flag"])
        evidence_parts.append(sufficiency_result["evidence"])
        if sufficiency_result["note"]:
            notes.append(sufficiency_result["note"])

        # 2. Negative marker conflict check
        neg_result = _check_negative_markers(adata, cluster, cluster_key, neg_markers)
        if neg_result["flag"]:
            flags.append(neg_result["flag"])
        evidence_parts.append(neg_result["evidence"])
        if neg_result["note"]:
            notes.append(neg_result["note"])

        # 3. Doublet / mixed signal heuristic
        doublet_result = _check_doublet_signal(adata, cluster, cluster_key, markers)
        if doublet_result["flag"]:
            flags.append(doublet_result["flag"])
        evidence_parts.append(doublet_result["evidence"])
        if doublet_result["note"]:
            notes.append(doublet_result["note"])

        # 4. Ontology consistency check
        onto_result = _check_ontology_consistency(cell_type, row.get("cl_id", ""))
        if onto_result["flag"]:
            flags.append(onto_result["flag"])
        if onto_result["note"]:
            notes.append(onto_result["note"])

        # 5. Ensemble agreement check
        ens_result = _check_ensemble_agreement(cluster, cell_type, ensemble_lookup)
        if ens_result["flag"]:
            flags.append(ens_result["flag"])
        evidence_parts.append(ens_result["evidence"])
        if ens_result["note"]:
            notes.append(ens_result["note"])

        # 6. Recalibrate confidence
        new_confidence = _recalibrate_confidence(
            original_confidence=row.get("confidence", CONFIDENCE_LOW),
            flags=flags,
            sufficiency=sufficiency_result,
            neg_conflict=neg_result,
        )

        critic_flags_list.append("; ".join(flags) if flags else "PASS")
        critic_evidence_list.append(" | ".join(evidence_parts))
        critic_confidence_list.append(new_confidence)
        critic_notes_list.append("; ".join(notes) if notes else "")

    results["critic_flags"] = critic_flags_list
    results["critic_evidence"] = critic_evidence_list
    results["critic_confidence"] = critic_confidence_list
    results["critic_notes"] = critic_notes_list

    return results


def _check_evidence_sufficiency(
    adata: ad.AnnData,
    cluster: str,
    cluster_key: str,
    pos_markers: list[str],
) -> dict:
    """Check if positive markers provide sufficient evidence for the annotation."""
    if not pos_markers:
        return {
            "flag": "NO_MARKERS",
            "evidence": "No positive markers defined for this cell type",
            "note": "Cannot validate without marker definitions",
        }

    detected = [g for g in pos_markers if g in adata.var_names]
    mask = adata.obs[cluster_key] == cluster
    subset = adata[mask]

    expressed_count = 0
    total_expr_details = []

    for gene in detected:
        idx = list(adata.var_names).index(gene)
        expr = subset.X[:, idx]
        expr = expr.toarray().flatten() if sparse.issparse(expr) else np.asarray(expr).flatten()

        pct = np.mean(expr > 0)
        mean_expr = np.mean(expr)
        if pct >= MARKER_PCT_THRESHOLD:
            expressed_count += 1
        total_expr_details.append(f"{gene}: {pct:.0%} cells, mean={mean_expr:.2f}")

    coverage = expressed_count / max(len(detected), 1)
    evidence = (
        f"Coverage: {expressed_count}/{len(detected)} ({coverage:.0%}) positive markers expressed"
    )

    if coverage < CRITIC_LOW_COVERAGE_THRESHOLD:
        return {
            "flag": "LOW_EVIDENCE",
            "evidence": evidence,
            "note": f"Only {coverage:.0%} of expected markers detected — annotation may be unreliable",
        }
    elif coverage < 0.5:
        return {
            "flag": "PARTIAL_EVIDENCE",
            "evidence": evidence,
            "note": f"Moderate marker coverage ({coverage:.0%}) — consider manual review",
        }
    else:
        return {"flag": "", "evidence": evidence, "note": ""}


def _check_negative_markers(
    adata: ad.AnnData,
    cluster: str,
    cluster_key: str,
    neg_markers: list[str],
) -> dict:
    """Check for negative marker conflicts — markers that should NOT be expressed."""
    if not neg_markers:
        return {"flag": "", "evidence": "No negative markers defined", "note": ""}

    detected = [g for g in neg_markers if g in adata.var_names]
    mask = adata.obs[cluster_key] == cluster
    subset = adata[mask]

    conflicts = []
    for gene in detected:
        idx = list(adata.var_names).index(gene)
        expr = subset.X[:, idx]
        expr = expr.toarray().flatten() if sparse.issparse(expr) else np.asarray(expr).flatten()

        pct = np.mean(expr > 0)
        if pct > CRITIC_NEG_MARKER_PCT_THRESHOLD:
            conflicts.append(f"{gene} ({pct:.0%})")

    if conflicts:
        return {
            "flag": "NEG_MARKER_CONFLICT",
            "evidence": f"Negative markers expressed: {', '.join(conflicts)}",
            "note": f"{len(conflicts)} negative marker(s) unexpectedly expressed — possible misannotation or doublet",
        }
    else:
        return {
            "flag": "",
            "evidence": f"All {len(detected)} negative markers appropriately absent",
            "note": "",
        }


def _check_doublet_signal(
    adata: ad.AnnData,
    cluster: str,
    cluster_key: str,
    all_markers: dict[str, dict],
) -> dict:
    """Check for doublet/mixed signal — two mutually exclusive lineage markers co-expressed."""
    mask = adata.obs[cluster_key] == cluster
    subset = adata[mask]

    # Find which cell types have strong marker expression in this cluster
    active_types = []
    for ct_name, ct_info in all_markers.items():
        pos = ct_info.get("positive_markers", [])
        detected = [g for g in pos if g in adata.var_names]
        if not detected:
            continue

        n_expressed = 0
        for gene in detected:
            idx = list(adata.var_names).index(gene)
            expr = subset.X[:, idx]
            expr = expr.toarray().flatten() if sparse.issparse(expr) else np.asarray(expr).flatten()
            pct = np.mean(expr > 0)
            if pct >= MARKER_PCT_THRESHOLD:
                n_expressed += 1

        coverage = n_expressed / max(len(detected), 1)
        if coverage >= 0.4:
            active_types.append((ct_name, coverage))

    active_types.sort(key=lambda x: -x[1])

    # Check for mutually exclusive lineages co-expressed
    if len(active_types) >= 2:
        top1, top1_cov = active_types[0]
        top2, top2_cov = active_types[1]

        # Check if these are from different lineages (not subtypes)
        if (
            top1_cov >= CRITIC_DOUBLET_COEXPR_THRESHOLD
            and top2_cov >= CRITIC_DOUBLET_COEXPR_THRESHOLD
        ):
            return {
                "flag": "POSSIBLE_DOUBLET",
                "evidence": f"Co-expression of {top1} ({top1_cov:.0%}) and {top2} ({top2_cov:.0%}) markers",
                "note": "Two distinct lineage signatures co-expressed — possible doublet or transitional state. Consider sub-clustering.",
            }

    return {"flag": "", "evidence": "No doublet signal detected", "note": ""}


def _check_ontology_consistency(cell_type: str, cl_id: str) -> dict:
    """Basic ontology consistency check."""
    if not cl_id:
        return {
            "flag": "NO_CL_ID",
            "note": "No Cell Ontology ID assigned — cannot verify term validity",
        }

    # Check CL ID format
    if not cl_id.startswith("CL:"):
        return {
            "flag": "INVALID_CL_FORMAT",
            "note": f"CL ID '{cl_id}' does not match expected format (CL:XXXXXXX)",
        }

    return {"flag": "", "note": ""}


def _recalibrate_confidence(
    original_confidence: str,
    flags: list[str],
    sufficiency: dict,
    neg_conflict: dict,
) -> str:
    """Recalibrate confidence based on critic findings."""
    confidence_order = {
        CONFIDENCE_HIGH: 4,
        CONFIDENCE_MEDIUM: 3,
        CONFIDENCE_LOW: 2,
        CONFIDENCE_REVIEW: 1,
    }
    reverse_order = {v: k for k, v in confidence_order.items()}

    level = confidence_order.get(original_confidence, 2)

    # Downgrade based on flags
    for flag in flags:
        if flag == "NO_MARKERS" or flag == "LOW_EVIDENCE":
            level = min(level, 1)
        elif flag == "PARTIAL_EVIDENCE":
            level = min(level, 2)
        elif flag == "NEG_MARKER_CONFLICT" or flag == "POSSIBLE_DOUBLET":
            level = min(level, 1)

    return reverse_order.get(level, CONFIDENCE_REVIEW)


def _check_ensemble_agreement(
    cluster: str,
    cell_type: str,
    ensemble_lookup: dict,
) -> dict:
    """Check if ensemble scoring agrees with the marker-based annotation."""
    if not ensemble_lookup:
        return {"flag": "", "evidence": "Ensemble: not available", "note": ""}

    ens = ensemble_lookup.get(str(cluster))
    if not ens:
        return {"flag": "", "evidence": "Ensemble: no data for this cluster", "note": ""}

    m_score = ens.get("marker_score", 0)
    r_score = ens.get("ref_score", 0)
    agreement = ens.get("agreement", True)
    source = ens.get("source", "")

    evidence = f"Ensemble: marker={m_score:.2f}, ref={r_score:.2f}, source={source}"

    # Flag if methods strongly disagree
    if not agreement:
        gap = abs(m_score - r_score)
        if gap > ENSEMBLE_AGREEMENT_THRESHOLD * 2:
            return {
                "flag": "ENSEMBLE_DISAGREEMENT",
                "evidence": evidence,
                "note": (
                    f"Marker ({m_score:.2f}) and reference ({r_score:.2f}) "
                    f"strongly disagree — possible transitional state"
                ),
            }
        else:
            return {
                "flag": "ENSEMBLE_MILD_DISAGREEMENT",
                "evidence": evidence,
                "note": f"Mild disagreement between methods (gap={gap:.2f})",
            }

    # Flag if single-source only (lower confidence)
    if source == "reference" and r_score < 0.4:
        return {
            "flag": "WEAK_REFERENCE_ONLY",
            "evidence": evidence,
            "note": "Only reference embedding supports this annotation with low confidence",
        }

    return {"flag": "", "evidence": evidence, "note": ""}


def generate_critic_summary(critic_results: pd.DataFrame) -> dict:
    """Generate summary statistics from critic results."""
    summary = {
        "total_clusters": len(critic_results),
        "pass": 0,
        "flagged": 0,
        "flag_types": {},
        "confidence_distribution": {},
        "clusters_needing_review": [],
    }

    for _, row in critic_results.iterrows():
        flags = row.get("critic_flags", "PASS")
        conf = row.get("critic_confidence", CONFIDENCE_REVIEW)

        # Confidence distribution
        summary["confidence_distribution"][conf] = (
            summary["confidence_distribution"].get(conf, 0) + 1
        )

        if flags == "PASS":
            summary["pass"] += 1
        else:
            summary["flagged"] += 1
            summary["clusters_needing_review"].append(row["cluster"])
            for flag in flags.split("; "):
                summary["flag_types"][flag] = summary["flag_types"].get(flag, 0) + 1

    return summary
