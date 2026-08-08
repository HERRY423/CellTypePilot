"""Annotation Critic — rules-based review of cell type annotations.

The Critic is CellTypePilot's trust layer. It doesn't just report confidence scores;
it performs a skeptical, evidence-based review of each annotation. It is not an
independent biological replicate or an independently trained annotator.
"""

from __future__ import annotations

import re
from contextlib import suppress

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .constants import (
    CONFIDENCE_HIGH,
    CONFIDENCE_LOW,
    CONFIDENCE_MEDIUM,
    CONFIDENCE_REVIEW,
    CRITIC_DOUBLET_ACTIVE_COVERAGE,
    CRITIC_DOUBLET_COEXPR_THRESHOLD,
    CRITIC_DOUBLET_OVERLAP_JACCARD,
    CRITIC_LOW_COVERAGE_THRESHOLD,
    CRITIC_NEG_MARKER_PCT_THRESHOLD,
    ENSEMBLE_AGREEMENT_THRESHOLD,
    MARKER_PCT_THRESHOLD,
)
from .data_adapter import build_lineage_groups, get_all_markers_for_tissue


def run_critic(
    adata: ad.AnnData,
    cluster_key: str,
    annotations: pd.DataFrame,
    atlas: dict,
    tissue: str,
    ensemble_info: pd.DataFrame | None = None,
    layer: str | None = None,
    evidence_policy: str = "database",
    marker_definitions: dict[str, dict] | None = None,
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
    markers = (
        marker_definitions
        if marker_definitions is not None
        else get_all_markers_for_tissue(atlas, tissue, evidence_policy=evidence_policy)
    )
    lineage_groups = build_lineage_groups(atlas, tissue)
    gene_idx = _gene_index_map(adata)
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
    sufficiency_details = []

    for _, row in results.iterrows():
        cluster = row["cluster"]
        cell_type = row["cell_type"]
        flags = []
        evidence_parts = []
        notes = []

        ct_info = markers.get(cell_type, {})
        pos_markers = ct_info.get("positive_markers", [])
        neg_markers = ct_info.get("negative_markers", [])

        status_counts = ct_info.get("evidence_status_counts", {})
        edge_verified = sum(
            count
            for status, count in status_counts.items()
            if status in {"database_record_verified", "primary_source_verified"}
        )
        if status_counts and edge_verified == 0:
            flags.append("AGGREGATE_PROVENANCE_ONLY")
            notes.append(
                "Marker relations cite database-level sources but have not been verified "
                "against a stable database record or primary-source evidence locator"
            )

        # 1. Evidence sufficiency check
        sufficiency_result = _check_evidence_sufficiency(
            adata, cluster, cluster_key, pos_markers, gene_idx, layer=layer
        )
        if sufficiency_result["flag"]:
            flags.append(sufficiency_result["flag"])
        evidence_parts.append(sufficiency_result["evidence"])
        if sufficiency_result["note"]:
            notes.append(sufficiency_result["note"])

        # Expression alone is not directional evidence. Require the candidate's
        # expected panel to contain markers that also pass the DE gates.
        de_support = float(row.get("pct_overlap", 0.0) or 0.0)
        if de_support < CRITIC_LOW_COVERAGE_THRESHOLD:
            flags.append("LOW_DE_SUPPORT")
            evidence_parts.append(f"DE support: {de_support:.0%} of expected markers")
            notes.append(
                "Too few expected markers pass direction, logFC, FDR, and expression gates"
            )
        elif de_support < 0.5:
            flags.append("PARTIAL_DE_SUPPORT")
            evidence_parts.append(f"DE support: {de_support:.0%} of expected markers")
            notes.append("Partial directional DE support; candidate requires review")

        # 2. Negative marker conflict check
        neg_result = _check_negative_markers(
            adata, cluster, cluster_key, neg_markers, gene_idx, layer=layer
        )
        if neg_result["flag"]:
            flags.append(neg_result["flag"])
        evidence_parts.append(neg_result["evidence"])
        if neg_result["note"]:
            notes.append(neg_result["note"])

        # 3. Doublet / mixed signal heuristic
        doublet_result = _check_doublet_signal(
            adata, cluster, cluster_key, markers, lineage_groups, gene_idx, layer=layer
        )
        if doublet_result["flag"]:
            flags.append(doublet_result["flag"])
        evidence_parts.append(doublet_result["evidence"])
        if doublet_result["note"]:
            notes.append(doublet_result["note"])

        # 4. Ontology consistency check
        onto_result = _check_ontology_consistency(
            cell_type,
            row.get("cl_id", ""),
            expected_cl_id=ct_info.get("cl_id", ""),
        )
        if onto_result["flag"]:
            flags.append(onto_result["flag"])
        if onto_result["note"]:
            notes.append(onto_result["note"])

        # User context may expand candidates, but unreviewed context-only support
        # is never allowed to upgrade a draft hypothesis into an accepted identity.
        if bool(row.get("context_only_support", False)):
            if row.get("context_review_status") == "reviewed":
                flags.append("REVIEWED_CONTEXT_SUPPORT")
                notes.append("Identity support depends only on an explicitly reviewed context pack")
            else:
                flags.append("UNREVIEWED_CONTEXT_ONLY")
                notes.append("Identity support depends only on unreviewed user context")

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

        if new_confidence == CONFIDENCE_REVIEW and not flags:
            flags.append("LOW_MODEL_CONFIDENCE")

        critic_flags_list.append("; ".join(flags) if flags else "PASS")
        critic_evidence_list.append(" | ".join(evidence_parts))
        critic_confidence_list.append(new_confidence)
        critic_notes_list.append("; ".join(notes) if notes else "")
        sufficiency_details.append(sufficiency_result)

    results["critic_flags"] = critic_flags_list
    results["critic_evidence"] = critic_evidence_list
    results["critic_confidence"] = critic_confidence_list
    results["critic_notes"] = critic_notes_list
    results["critic_method"] = "rules_based_same_run_review"
    results["critic_independence"] = "not_independent"

    for column in (
        "n_expected_markers",
        "n_present_markers",
        "n_expressed_markers",
        "n_missing_markers",
        "n_silent_markers",
        "expected_marker_coverage",
        "present_markers",
        "expressed_markers",
        "missing_markers",
        "silent_markers",
    ):
        results[column] = [
            detail.get(column, 0 if column.startswith("n_") else "")
            for detail in sufficiency_details
        ]

    # Fail closed: preserve the best candidate, but publish Unknown for any
    # annotation that lacks adequate, non-conflicting evidence.
    abstain_flags = {
        "NO_MARKERS",
        "LOW_EVIDENCE",
        "PARTIAL_EVIDENCE",
        "LOW_DE_SUPPORT",
        "PARTIAL_DE_SUPPORT",
        "NEG_MARKER_CONFLICT",
        "POSSIBLE_DOUBLET",
        "ENSEMBLE_DISAGREEMENT",
        "WEAK_REFERENCE_ONLY",
        "LOW_MODEL_CONFIDENCE",
        "ONTOLOGY_MISMATCH",
        "UNKNOWN_ATLAS_LABEL",
        "INVALID_CL_FORMAT",
        "NO_CL_ID",
        "UNREVIEWED_CONTEXT_ONLY",
    }
    results["candidate_cell_type"] = results["cell_type"]
    results["candidate_cl_id"] = results.get("cl_id", "")
    decisions = []
    reasons = []
    for _, row in results.iterrows():
        row_flags = set(str(row["critic_flags"]).split("; "))
        active = sorted(row_flags & abstain_flags)
        abstain = bool(active) or row["critic_confidence"] == CONFIDENCE_REVIEW
        decisions.append("abstain" if abstain else "accepted")
        reasons.append("; ".join(active) if abstain else "")
    results["decision"] = decisions
    results["abstain_reason"] = reasons
    abstain_mask = results["decision"] == "abstain"
    results.loc[abstain_mask, "cell_type"] = "Unknown"
    if "cl_id" in results.columns:
        results.loc[abstain_mask, "cl_id"] = ""

    # Layer-1 evidence summary: one glanceable verdict line per cluster
    results["evidence_summary"] = [format_evidence_summary(row) for _, row in results.iterrows()]

    return results


def _gene_index_map(adata: ad.AnnData) -> dict[str, int]:
    """Precompute gene name → column index once per run."""
    return {g: i for i, g in enumerate(adata.var_names)}


def _expression_pct(
    adata: ad.AnnData,
    mask,
    gene: str,
    gene_idx: dict[str, int],
    layer: str | None = None,
) -> float:
    """Fraction of masked cells with detectable expression of a gene."""
    idx = gene_idx[gene]
    matrix = adata.layers[layer] if layer is not None else adata.X
    expr = matrix[mask][:, idx]
    expr = expr.toarray().flatten() if sparse.issparse(expr) else np.asarray(expr).flatten()
    return float(np.mean(expr > 0))


def _check_evidence_sufficiency(
    adata: ad.AnnData,
    cluster: str,
    cluster_key: str,
    pos_markers: list[str],
    gene_idx: dict[str, int] | None = None,
    layer: str | None = None,
) -> dict:
    """Check if positive markers provide sufficient evidence for the annotation."""
    if not pos_markers:
        return {
            "flag": "NO_MARKERS",
            "evidence": "No positive markers defined for this cell type",
            "note": "Cannot validate without marker definitions",
            "n_expected_markers": 0,
            "n_present_markers": 0,
            "n_expressed_markers": 0,
            "n_missing_markers": 0,
            "n_silent_markers": 0,
            "expected_marker_coverage": 0.0,
            "present_markers": "",
            "expressed_markers": "",
            "missing_markers": "",
            "silent_markers": "",
        }

    if gene_idx is None:
        gene_idx = _gene_index_map(adata)

    present = [g for g in pos_markers if g in gene_idx]
    missing = [g for g in pos_markers if g not in gene_idx]
    mask = (adata.obs[cluster_key].astype(str) == str(cluster)).values

    expressed = []
    for gene in present:
        pct = _expression_pct(adata, mask, gene, gene_idx, layer=layer)
        if pct >= MARKER_PCT_THRESHOLD:
            expressed.append(gene)
    silent = [g for g in present if g not in expressed]

    coverage = len(expressed) / len(pos_markers)
    evidence = (
        f"Coverage: {len(expressed)}/{len(pos_markers)} ({coverage:.0%}) expected markers expressed; "
        f"{len(missing)} missing from matrix; {len(silent)} present but silent"
    )
    details = {
        "n_expected_markers": len(pos_markers),
        "n_present_markers": len(present),
        "n_expressed_markers": len(expressed),
        "n_missing_markers": len(missing),
        "n_silent_markers": len(silent),
        "expected_marker_coverage": round(coverage, 4),
        "present_markers": ";".join(present),
        "expressed_markers": ";".join(expressed),
        "missing_markers": ";".join(missing),
        "silent_markers": ";".join(silent),
    }

    if coverage < CRITIC_LOW_COVERAGE_THRESHOLD:
        return {
            "flag": "LOW_EVIDENCE",
            "evidence": evidence,
            **details,
            "note": f"Only {coverage:.0%} of expected markers detected — annotation may be unreliable",
        }
    elif coverage < 0.5:
        return {
            "flag": "PARTIAL_EVIDENCE",
            "evidence": evidence,
            **details,
            "note": f"Moderate marker coverage ({coverage:.0%}) — consider manual review",
        }
    else:
        return {"flag": "", "evidence": evidence, "note": "", **details}


def _check_negative_markers(
    adata: ad.AnnData,
    cluster: str,
    cluster_key: str,
    neg_markers: list[str],
    gene_idx: dict[str, int] | None = None,
    layer: str | None = None,
) -> dict:
    """Check for negative marker conflicts — markers that should NOT be expressed."""
    if not neg_markers:
        return {"flag": "", "evidence": "No negative markers defined", "note": ""}

    if gene_idx is None:
        gene_idx = _gene_index_map(adata)

    detected = [g for g in neg_markers if g in gene_idx]
    mask = (adata.obs[cluster_key].astype(str) == str(cluster)).values

    conflicts = []
    for gene in detected:
        pct = _expression_pct(adata, mask, gene, gene_idx, layer=layer)
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
    lineage_groups: dict[str, str] | None = None,
    gene_idx: dict[str, int] | None = None,
    layer: str | None = None,
) -> dict:
    """Check for doublet/mixed signal — two mutually exclusive lineages co-expressed.

    Calibrated to avoid the three main false-positive modes:
    1. Same-lineage co-expression (e.g. "T cell" + "CD4+ T cell") is subtype
       refinement, not a doublet — only pairs from different root lineages count.
    2. Heavily overlapping marker sets are redundant signatures of one biology,
       not independent lineages (Jaccard guard).
    3. Genes shared across panels (e.g. the cytotoxic program GNLY/PRF1 shared
       by NK and CD8+ T cells) are weak lineage evidence on their own, so they
       are excluded from panel coverage (fractional weighting as fallback).
    """
    lineage_groups = lineage_groups or {}
    if gene_idx is None:
        gene_idx = _gene_index_map(adata)

    mask = (adata.obs[cluster_key].astype(str) == str(cluster)).values

    # Gene specificity: count how many panels list each gene. Genes shared
    # across panels are weak lineage evidence and are excluded from coverage;
    # a purely shared panel falls back to fractional weighting (1/n_panels).
    panel_count: dict[str, int] = {}
    for ct_info in all_markers.values():
        for g in set(ct_info.get("positive_markers", [])):
            panel_count[g] = panel_count.get(g, 0) + 1

    # Find which cell types have strong marker expression in this cluster
    active_types: list[tuple[str, float]] = []
    for ct_name, ct_info in all_markers.items():
        pos = ct_info.get("positive_markers", [])
        detected = [g for g in pos if g in gene_idx]
        if not detected:
            continue

        specific = [g for g in pos if panel_count.get(g, 1) == 1]
        if specific:
            n_expressed = sum(
                1
                for gene in specific
                if gene in gene_idx
                and _expression_pct(adata, mask, gene, gene_idx, layer=layer)
                >= MARKER_PCT_THRESHOLD
            )
            coverage = n_expressed / len(specific)
        else:
            weighted = sum(
                1.0 / panel_count[gene]
                for gene in detected
                if _expression_pct(adata, mask, gene, gene_idx, layer=layer) >= MARKER_PCT_THRESHOLD
            )
            coverage = weighted / len(detected)
        if coverage >= CRITIC_DOUBLET_ACTIVE_COVERAGE:
            active_types.append((ct_name, coverage))

    active_types.sort(key=lambda x: -x[1])

    if len(active_types) < 2:
        return {"flag": "", "evidence": "No doublet signal detected", "note": ""}

    # Look for a strong co-active signature from a DIFFERENT root lineage.
    # Same-lineage pairs (parent + subtype) are expected biology, not doublets.
    top1, top1_cov = active_types[0]
    top1_markers = set(all_markers[top1].get("positive_markers", []))

    for top2, top2_cov in active_types[1:]:
        if lineage_groups.get(top1, top1) == lineage_groups.get(top2, top2):
            continue

        overlap = top1_markers & set(all_markers[top2].get("positive_markers", []))
        union = top1_markers | set(all_markers[top2].get("positive_markers", []))
        jaccard = len(overlap) / len(union) if union else 0.0
        if jaccard >= CRITIC_DOUBLET_OVERLAP_JACCARD:
            continue  # redundant signatures, not independent lineages

        if (
            top1_cov >= CRITIC_DOUBLET_COEXPR_THRESHOLD
            and top2_cov >= CRITIC_DOUBLET_COEXPR_THRESHOLD
        ):
            return {
                "flag": "POSSIBLE_DOUBLET",
                "evidence": (
                    f"Cross-lineage co-expression of {top1} ({top1_cov:.0%}) "
                    f"and {top2} ({top2_cov:.0%}) markers"
                ),
                "note": "Two distinct lineage signatures co-expressed — possible doublet or transitional state. Consider sub-clustering.",
            }

    return {"flag": "", "evidence": "No doublet signal detected", "note": ""}


def _check_ontology_consistency(
    cell_type: str,
    cl_id: str,
    expected_cl_id: str | None = None,
) -> dict:
    """Validate identifier syntax and the atlas-declared label-to-CL mapping.

    This deliberately does not claim live Cell Ontology validation. The bundled
    atlas mapping is checked exactly; live term existence/version validation is a
    separate provenance task.
    """
    if expected_cl_id is None:
        expected_cl_id = cl_id
    elif not expected_cl_id:
        return {
            "flag": "UNKNOWN_ATLAS_LABEL",
            "note": f"'{cell_type}' has no declared Cell Ontology mapping in the active atlas",
        }
    if not cl_id:
        return {
            "flag": "NO_CL_ID",
            "note": "No Cell Ontology ID assigned — cannot verify term validity",
        }

    if re.fullmatch(r"CL:\d{7}", str(cl_id)) is None:
        return {
            "flag": "INVALID_CL_FORMAT",
            "note": f"CL ID '{cl_id}' does not match expected format (CL:XXXXXXX)",
        }

    if str(cl_id) != str(expected_cl_id):
        return {
            "flag": "ONTOLOGY_MISMATCH",
            "note": (
                f"Label '{cell_type}' is mapped to {expected_cl_id} in the active atlas, "
                f"but the annotation supplied {cl_id}"
            ),
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
        if flag in {"NO_MARKERS", "LOW_EVIDENCE", "LOW_DE_SUPPORT"}:
            level = min(level, 1)
        elif flag in {"PARTIAL_EVIDENCE", "PARTIAL_DE_SUPPORT"}:
            level = min(level, 2)
        elif flag in {"AGGREGATE_PROVENANCE_ONLY", "REVIEWED_CONTEXT_SUPPORT"}:
            level = min(level, 3)
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

    # Narrative one-liner for reports / CLI / agent presentation
    summary["narrative"] = format_run_narrative(summary)

    return summary


# ──────────────────────────────────────────────
# Evidence summary layer
# ──────────────────────────────────────────────

# Fixed action guidance per flag — turns a flag into a next step.
FLAG_ACTIONS = {
    "PASS": "Accept annotation.",
    "LOW_EVIDENCE": "Manual review needed; consider sub-clustering.",
    "PARTIAL_EVIDENCE": "Review; may be correct for rare/transitional states.",
    "LOW_DE_SUPPORT": "Insufficient directional DE support; keep as Unknown.",
    "PARTIAL_DE_SUPPORT": "Partial directional DE support; review before labeling.",
    "NEG_MARKER_CONFLICT": "Likely misannotation or doublet; verify markers.",
    "POSSIBLE_DOUBLET": "Sub-cluster or mark as doublet.",
    "NO_MARKERS": "Add marker definitions for this cell type.",
    "NO_CL_ID": "Assign a Cell Ontology ID manually if needed.",
    "INVALID_CL_FORMAT": "Correct the Cell Ontology ID format.",
    "ONTOLOGY_MISMATCH": "Correct the label-to-CL mapping before accepting the annotation.",
    "UNKNOWN_ATLAS_LABEL": "Add a versioned label-to-CL mapping before acceptance.",
    "AGGREGATE_PROVENANCE_ONLY": (
        "Treat as a draft; verify marker edges against stable records or primary sources."
    ),
    "REVIEWED_CONTEXT_SUPPORT": (
        "Retain the medium-confidence cap and independently review the custom marker panel."
    ),
    "UNREVIEWED_CONTEXT_ONLY": "Keep as Unknown until the custom marker panel is reviewed.",
    "CALIBRATED_LOW_CONFIDENCE": "Keep as Unknown under the locked calibration policy.",
    "ENSEMBLE_DISAGREEMENT": "Review; scoring methods strongly disagree.",
    "ENSEMBLE_MILD_DISAGREEMENT": "Minor method disagreement; usually acceptable.",
    "WEAK_REFERENCE_ONLY": "Weak reference support; seek marker evidence.",
}


def format_evidence_summary(row) -> str:
    """Render one glanceable evidence-summary line for a cluster.

    Bridges the gap between the conclusion layer (confidence/flag) and the
    full evidence layer (evidence_table.csv): verdict + key evidence + action.
    """
    cluster = row.get("cluster", "?")
    cell_type = row.get("cell_type", "?")
    candidate = row.get("candidate_cell_type", cell_type)
    conf = row.get("critic_confidence", "unknown")
    flags = row.get("critic_flags", "PASS")
    overlap = row.get("pct_overlap", None)
    score = row.get("combined_score", None)

    verdict = "PASS" if flags == "PASS" else f"FLAGGED ({flags})"

    label = cell_type if cell_type == candidate else f"{cell_type} (candidate: {candidate})"
    parts = [f"Cluster {cluster} → {label} [{conf.upper()}] {verdict}"]

    evidence_bits = []
    if overlap is not None:
        with suppress(TypeError, ValueError):
            evidence_bits.append(f"marker overlap {float(overlap):.0%}")
    if score is not None:
        with suppress(TypeError, ValueError):
            evidence_bits.append(f"score {float(score):.2f}")
    if evidence_bits:
        parts.append("Evidence: " + ", ".join(evidence_bits))

    if flags == "PASS":
        parts.append("Action: " + FLAG_ACTIONS["PASS"])
    else:
        actions = []
        for flag in flags.split("; "):
            action = FLAG_ACTIONS.get(flag)
            if action and action not in actions:
                actions.append(action)
        if actions:
            parts.append("Action: " + " ".join(actions))

    return " | ".join(parts)


def format_run_narrative(summary: dict) -> str:
    """Render a one-paragraph narrative of a full critic run."""
    total = summary.get("total_clusters", 0)
    passed = summary.get("pass", 0)
    flagged = summary.get("flagged", 0)
    conf_dist = summary.get("confidence_distribution", {})
    flag_types = summary.get("flag_types", {})

    if total == 0:
        return "No clusters were reviewed."

    conf_parts = [
        f"{count} {level}" for level, count in sorted(conf_dist.items(), key=lambda kv: -kv[1])
    ]
    narrative = (
        f"{total} cluster(s) reviewed: {passed} passed critic checks, "
        f"{flagged} flagged for review. Confidence distribution: " + ", ".join(conf_parts) + "."
    )
    if flag_types:
        flag_desc = ", ".join(f"{flag} ×{count}" for flag, count in flag_types.items())
        narrative += f" Flag breakdown: {flag_desc}."
        narrative += (
            " Flagged clusters need human adjudication before publication; "
            "unflagged clusters carry full marker evidence in evidence_table.csv."
        )
    else:
        narrative += " All annotations carry full marker evidence in evidence_table.csv."
    return narrative
