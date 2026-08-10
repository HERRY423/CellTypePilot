"""Actionable Unknown diagnoses and contrastive candidate evidence.

This module explains existing deterministic results. It never rescales,
reranks, accepts, or rescues an identity decision.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

CONTRAST_SCHEMA = "celltypepilot.contrastive-evidence.v1"
GAP_SCHEMA = "celltypepilot.actionable-evidence-gaps.v1"
CONTRAST_FILE = "contrastive_evidence.csv"
GAP_FILE = "evidence_gaps.json"

FORBIDDEN_ACTIONS = [
    "do_not_convert_ranking_margin_to_probability",
    "do_not_accept_a_candidate_only_because_it_is_ranked_first",
    "do_not_use_free_text_or_literature_search_as_scoring_evidence",
    "do_not_use_state_or_novelty_to_rescue_identity",
    "do_not_override_unknown_without_an_explicit_human_decision",
]

GAP_RULES: dict[str, dict[str, Any]] = {
    "atlas_marker_definition_missing": {
        "why": "The candidate has no governed positive-marker definition in scope.",
        "actions": ["review_compatible_data_only_pack", "keep_unknown"],
    },
    "marker_not_measured": {
        "why": "Expected markers are absent from the active expression feature space.",
        "actions": [
            "inspect_gene_identity_and_assay_coverage",
            "request_a_compatible_expression_matrix",
            "keep_unknown",
        ],
    },
    "marker_present_but_silent": {
        "why": "Expected markers are measured but do not meet the expression-fraction gate.",
        "actions": [
            "review_expression_layer_and_cluster_qc",
            "inspect_marker_expression",
            "keep_unknown",
        ],
    },
    "directional_support_gap": {
        "why": "Too few expected markers pass direction, logFC, FDR, and expression gates.",
        "actions": ["inspect_cluster_level_de", "review_clustering_resolution", "keep_unknown"],
    },
    "negative_marker_conflict": {
        "why": "Markers expected to contradict the candidate are expressed.",
        "actions": [
            "inspect_negative_marker_expression",
            "assess_ambient_rna_or_mixed_identity",
            "keep_unknown",
        ],
    },
    "mixed_lineage_or_doublet": {
        "why": "More than one lineage program is active in the cluster.",
        "actions": ["run_or_import_doublet_diagnostics", "review_subclustering", "keep_unknown"],
    },
    "reference_disagreement": {
        "why": "Marker and reference evidence do not support the same identity.",
        "actions": [
            "compare_marker_and_reference_evidence",
            "review_reference_scope_and_provenance",
            "keep_unknown",
        ],
    },
    "candidate_separation_unresolved": {
        "why": "Top candidates lack candidate-specific directional marker support.",
        "actions": [
            "review_contrastive_evidence",
            "inspect_discriminating_markers",
            "keep_unknown",
        ],
    },
    "ontology_or_atlas_scope_gap": {
        "why": "The candidate identity or ontology mapping is not valid in the active Atlas scope.",
        "actions": [
            "review_atlas_identity_mapping",
            "review_compatible_data_only_pack",
            "keep_unknown",
        ],
    },
    "context_review_gap": {
        "why": "The candidate depends on unreviewed user-supplied context.",
        "actions": ["request_human_review_of_context_pack", "keep_unknown"],
    },
    "calibration_downgrade": {
        "why": "A locked calibration policy downgraded this identity.",
        "actions": ["inspect_calibration_artifact", "keep_unknown"],
    },
    "aggregate_provenance_gap": {
        "why": "Marker relationships have aggregate provenance but no locked edge-level record.",
        "actions": ["queue_marker_edges_for_human_curation", "keep_unknown"],
    },
    "unresolved_uncertainty": {
        "why": "The current deterministic evidence does not authorize an identity.",
        "actions": ["review_full_evidence_table", "keep_unknown"],
    },
}


def _split_markers(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    if isinstance(value, (list, tuple, set)):
        return sorted({str(item).strip() for item in value if str(item).strip()})
    return sorted({item.strip() for item in str(value).split(";") if item.strip()})


def _as_float(value: Any, default: float = 0.0) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return default if pd.isna(parsed) else parsed


def _as_int(value: Any) -> int:
    return int(_as_float(value))


def _as_text(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return str(value)


def _candidate_details(marker_scores: pd.DataFrame, cluster: str, cell_type: str) -> dict:
    if marker_scores.empty:
        return {}
    rows = marker_scores[
        (marker_scores["cluster"].astype(str) == str(cluster))
        & (marker_scores["cell_type"].astype(str) == str(cell_type))
    ]
    return rows.iloc[0].to_dict() if not rows.empty else {}


def _join(items: set[str] | list[str]) -> str:
    return ";".join(sorted(set(items)))


def build_contrastive_evidence(
    marker_scores: pd.DataFrame,
    ensemble_scores: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Describe top-1 versus top-2 evidence without changing their ranking."""
    use_ensemble = ensemble_scores is not None and not ensemble_scores.empty
    ranking = ensemble_scores.copy() if use_ensemble else marker_scores.copy()
    if ranking.empty:
        return pd.DataFrame()
    score_column = "ensemble_score" if use_ensemble else "combined_score"
    rank_column = "rank"
    records: list[dict[str, Any]] = []
    for cluster, group in ranking.groupby("cluster", sort=True):
        ordered = group.sort_values(rank_column, kind="stable").head(2)
        if ordered.empty:
            continue
        selected = ordered.iloc[0]
        alternative = ordered.iloc[1] if len(ordered) > 1 else None
        selected_type = str(selected["cell_type"])
        alternative_type = str(alternative["cell_type"]) if alternative is not None else ""
        selected_marker = _candidate_details(marker_scores, str(cluster), selected_type)
        alternative_marker = _candidate_details(marker_scores, str(cluster), alternative_type)
        selected_support = set(_split_markers(selected_marker.get("pos_supporting_markers")))
        alternative_support = set(_split_markers(alternative_marker.get("pos_supporting_markers")))
        selected_only = selected_support - alternative_support
        alternative_only = alternative_support - selected_support
        shared = selected_support & alternative_support
        selected_score = _as_float(selected.get(score_column, 0.0))
        alternative_score = (
            _as_float(alternative.get(score_column, 0.0)) if alternative is not None else 0.0
        )
        if alternative is None:
            interpretation = "No second candidate was available in the active scoring scope."
        elif selected_only:
            interpretation = (
                "Selected candidate has directional supporting markers not shared by the "
                "alternative; review negative conflicts and missing markers before acceptance."
            )
        elif alternative_only:
            interpretation = (
                "Alternative candidate has candidate-specific directional support while the "
                "selected candidate does not; ranking requires human review."
            )
        else:
            interpretation = (
                "No candidate-specific directional supporting markers separate the top two "
                "candidates; the ranking signal alone is insufficient evidence."
            )
        records.append(
            {
                "schema_version": CONTRAST_SCHEMA,
                "cluster": str(cluster),
                "ranking_source": "ensemble" if use_ensemble else "marker",
                "ranking_semantics": "relative_evidence_signal_not_probability",
                "selected_candidate": selected_type,
                "selected_cl_id": _as_text(selected_marker.get("cl_id", "")),
                "selected_score": round(selected_score, 4),
                "alternative_candidate": alternative_type,
                "alternative_cl_id": _as_text(alternative_marker.get("cl_id", "")),
                "alternative_score": round(alternative_score, 4),
                "score_margin": round(selected_score - alternative_score, 4),
                "shared_supporting_markers": _join(shared),
                "selected_only_supporting_markers": _join(selected_only),
                "alternative_only_supporting_markers": _join(alternative_only),
                "selected_missing_markers": _as_text(
                    selected_marker.get("pos_missing_markers", "")
                ),
                "alternative_missing_markers": _as_text(
                    alternative_marker.get("pos_missing_markers", "")
                ),
                "selected_silent_markers": _as_text(selected_marker.get("pos_silent_markers", "")),
                "alternative_silent_markers": _as_text(
                    alternative_marker.get("pos_silent_markers", "")
                ),
                "selected_negative_conflicts": _as_text(
                    selected_marker.get("neg_expressed_markers", "")
                ),
                "alternative_negative_conflicts": _as_text(
                    alternative_marker.get("neg_expressed_markers", "")
                ),
                "selected_marker_provenance_status": _as_text(
                    selected_marker.get("marker_provenance_status", "")
                ),
                "alternative_marker_provenance_status": _as_text(
                    alternative_marker.get("marker_provenance_status", "")
                ),
                "interpretation": interpretation,
                "claim_boundary": (
                    "Contrastive evidence explains an existing ranking; it does not validate "
                    "either identity or authorize an override."
                ),
            }
        )
    return pd.DataFrame(records)


def attach_contrastive_evidence(
    critic_results: pd.DataFrame,
    contrastive: pd.DataFrame,
) -> pd.DataFrame:
    """Attach compact contrast fields to the evidence table."""
    if critic_results.empty or contrastive.empty:
        return critic_results
    columns = [
        "cluster",
        "alternative_candidate",
        "alternative_cl_id",
        "alternative_score",
        "score_margin",
        "shared_supporting_markers",
        "selected_only_supporting_markers",
        "alternative_only_supporting_markers",
        "interpretation",
    ]
    compact = contrastive[columns].rename(
        columns={
            "interpretation": "contrastive_evidence_summary",
            "score_margin": "candidate_score_margin",
        }
    )
    result = critic_results.copy()
    result["cluster"] = result["cluster"].astype(str)
    compact = compact.copy()
    compact["cluster"] = compact["cluster"].astype(str)
    return result.merge(compact, on="cluster", how="left")


def _gap(gap_type: str, observed: str) -> dict[str, Any]:
    rule = GAP_RULES[gap_type]
    return {
        "gap_type": gap_type,
        "observed": observed,
        "why_it_matters": rule["why"],
        "allowed_next_actions": list(rule["actions"]),
    }


def build_actionable_evidence_gaps(critic_results: pd.DataFrame) -> dict[str, Any]:
    """Convert each abstention into bounded, non-label-producing next actions."""
    clusters: list[dict[str, Any]] = []
    counts: Counter[str] = Counter()
    for row in critic_results.to_dict(orient="records"):
        decision = str(row.get("decision", "")).strip().casefold()
        if decision != "abstain" and str(row.get("cell_type", "")) != "Unknown":
            continue
        flags = set(_split_markers(str(row.get("critic_flags", "")).replace("; ", ";")))
        gaps: list[dict[str, Any]] = []
        expected = _as_int(row.get("n_expected_markers", 0))
        missing = _as_int(row.get("n_missing_markers", 0))
        silent = _as_int(row.get("n_silent_markers", 0))
        if expected == 0 or "NO_MARKERS" in flags:
            gaps.append(_gap("atlas_marker_definition_missing", "no governed markers in scope"))
        if missing:
            gaps.append(
                _gap("marker_not_measured", f"{missing}/{expected} expected markers absent")
            )
        if silent:
            gaps.append(
                _gap("marker_present_but_silent", f"{silent}/{expected} expected markers silent")
            )
        if flags & {"LOW_DE_SUPPORT", "PARTIAL_DE_SUPPORT", "LOW_EVIDENCE", "PARTIAL_EVIDENCE"}:
            gaps.append(
                _gap(
                    "directional_support_gap",
                    f"pct_overlap={row.get('pct_overlap', '')}; flags={';'.join(sorted(flags))}",
                )
            )
        if "NEG_MARKER_CONFLICT" in flags:
            gaps.append(
                _gap(
                    "negative_marker_conflict",
                    str(row.get("neg_expressed_markers", "") or "negative markers expressed"),
                )
            )
        if "POSSIBLE_DOUBLET" in flags:
            gaps.append(_gap("mixed_lineage_or_doublet", "POSSIBLE_DOUBLET critic flag"))
        if flags & {"ENSEMBLE_DISAGREEMENT", "WEAK_REFERENCE_ONLY"}:
            gaps.append(_gap("reference_disagreement", ";".join(sorted(flags))))
        if flags & {"ONTOLOGY_MISMATCH", "UNKNOWN_ATLAS_LABEL", "INVALID_CL_FORMAT", "NO_CL_ID"}:
            gaps.append(_gap("ontology_or_atlas_scope_gap", ";".join(sorted(flags))))
        if "UNREVIEWED_CONTEXT_ONLY" in flags:
            gaps.append(_gap("context_review_gap", "UNREVIEWED_CONTEXT_ONLY"))
        if "CALIBRATED_LOW_CONFIDENCE" in flags:
            gaps.append(_gap("calibration_downgrade", "CALIBRATED_LOW_CONFIDENCE"))
        if "AGGREGATE_PROVENANCE_ONLY" in flags:
            gaps.append(_gap("aggregate_provenance_gap", "AGGREGATE_PROVENANCE_ONLY"))
        selected_only = _split_markers(row.get("selected_only_supporting_markers"))
        if row.get("alternative_candidate") and not selected_only:
            gaps.append(
                _gap(
                    "candidate_separation_unresolved",
                    f"candidate={row.get('candidate_cell_type', '')}; "
                    f"alternative={row.get('alternative_candidate', '')}; "
                    "no selected-only directional markers",
                )
            )
        if not gaps:
            gaps.append(_gap("unresolved_uncertainty", str(row.get("abstain_reason", ""))))
        allowed: list[str] = []
        for item in gaps:
            counts[item["gap_type"]] += 1
            for action in item["allowed_next_actions"]:
                if action not in allowed:
                    allowed.append(action)
        clusters.append(
            {
                "cluster": str(row.get("cluster", "")),
                "published_identity": str(row.get("cell_type", "Unknown")),
                "candidate_cell_type": str(row.get("candidate_cell_type", "")),
                "candidate_cl_id": str(row.get("candidate_cl_id", "")),
                "alternative_candidate": str(row.get("alternative_candidate", "")),
                "decision": str(row.get("decision", "")),
                "critic_confidence": str(row.get("critic_confidence", "")),
                "critic_flags": str(row.get("critic_flags", "")),
                "gaps": gaps,
                "allowed_next_actions": allowed,
                "forbidden_actions": list(FORBIDDEN_ACTIONS),
                "human_action_required": True,
            }
        )
    return {
        "schema_version": GAP_SCHEMA,
        "status": "action_required" if clusters else "no_unknown_clusters",
        "n_unknown_clusters": len(clusters),
        "gap_type_counts": dict(sorted(counts.items())),
        "clusters": clusters,
        "claim_boundary": (
            "Evidence gaps route investigation and review. They do not select a replacement "
            "identity, weaken abstention, or authorize autonomous evidence promotion."
        ),
    }


def write_agent_evidence_artifacts(
    output_dir: str | Path,
    contrastive: pd.DataFrame,
    gaps: dict[str, Any],
) -> dict[str, Path]:
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    contrast_path = root / CONTRAST_FILE
    gap_path = root / GAP_FILE
    contrastive.to_csv(contrast_path, index=False)
    gap_path.write_text(json.dumps(gaps, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return {"contrastive_evidence": contrast_path, "evidence_gaps": gap_path}


def load_agent_evidence_indexes(output_dir: str | Path) -> tuple[dict[str, dict], dict[str, dict]]:
    """Load artifacts for the review queue; tolerate legacy runs without them."""
    root = Path(output_dir)
    contrast_index: dict[str, dict] = {}
    gap_index: dict[str, dict] = {}
    contrast_path = root / CONTRAST_FILE
    if contrast_path.is_file():
        frame = pd.read_csv(contrast_path, dtype=str).fillna("")
        contrast_index = {str(row["cluster"]): row for row in frame.to_dict(orient="records")}
    gap_path = root / GAP_FILE
    if gap_path.is_file():
        payload = json.loads(gap_path.read_text(encoding="utf-8"))
        gap_index = {str(row.get("cluster", "")): row for row in payload.get("clusters", [])}
    return contrast_index, gap_index
