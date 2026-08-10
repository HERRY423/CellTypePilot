"""Human review panel payload: Identity × State × Novelty side-by-side.

Builds a single JSON document for the Web Inspector / Agent hosts with:
supporting/opposing markers, neighbor candidates, donor/batch strata when
present, and literature hooks. Never mutates annotations.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

REVIEW_PANEL_SCHEMA = "celltypepilot.review-panel.v1"


def _split_markers(value: Any) -> list[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return []
    if "," in text:
        return [part.strip() for part in text.split(",") if part.strip()]
    if ";" in text:
        return [part.strip() for part in text.split(";") if part.strip()]
    return [text]


def _row_for_cluster(evidence: pd.DataFrame, cluster: str) -> dict[str, Any]:
    if evidence.empty or "cluster" not in evidence.columns:
        return {}
    rows = evidence[evidence["cluster"].astype(str) == str(cluster)]
    if rows.empty:
        return {}
    return rows.iloc[0].to_dict()


def _neighbor_candidates(evidence: pd.DataFrame, cluster: str, limit: int = 5) -> list[dict]:
    """Pull runner-up candidates from evidence columns when present."""
    row = _row_for_cluster(evidence, cluster)
    if not row:
        return []
    candidates: list[dict] = []
    # Common column patterns from scorer / ensemble.
    for rank in range(1, limit + 1):
        for prefix in ("candidate", "alt", "rank"):
            name_key = f"{prefix}_{rank}_type"
            score_key = f"{prefix}_{rank}_score"
            if name_key in row and str(row.get(name_key, "")).strip():
                candidates.append(
                    {
                        "rank": rank,
                        "cell_type": str(row[name_key]),
                        "score": float(row.get(score_key) or 0.0) if score_key in row else None,
                    }
                )
    # Fallback: second-best fields
    if not candidates:
        for key, score_key in (
            ("second_cell_type", "second_score"),
            ("runner_up", "runner_up_score"),
            ("candidate_cell_type", "combined_score"),
        ):
            if key in row and str(row.get(key, "")).strip():
                candidates.append(
                    {
                        "rank": len(candidates) + 1,
                        "cell_type": str(row[key]),
                        "score": float(row.get(score_key) or 0.0) if score_key in row else None,
                    }
                )
    # De-dupe against primary identity
    primary = str(row.get("cell_type", ""))
    out = []
    seen = {primary}
    for item in candidates:
        if item["cell_type"] in seen:
            continue
        seen.add(item["cell_type"])
        out.append(item)
    return out[:limit]


def _donor_batch_strata(adata, cluster: str, cluster_key: str) -> dict[str, Any]:
    if adata is None or cluster_key not in adata.obs.columns:
        return {
            "status": "not_assessed_missing_cluster_key",
            "donors": [],
            "batches": [],
            "samples": [],
        }
    mask = adata.obs[cluster_key].astype(str) == str(cluster)
    sub = adata.obs.loc[mask]
    if sub.empty:
        return {"status": "not_assessed_empty_cluster", "donors": [], "batches": [], "samples": []}

    def _counts(aliases: tuple[str, ...]) -> dict[str, Any]:
        for col in aliases:
            if col in sub.columns:
                vc = sub[col].astype(str).value_counts()
                return {
                    "status": "assessed",
                    "column": col,
                    "levels": [
                        {"label": str(idx), "n_cells": int(count)}
                        for idx, count in vc.head(12).items()
                    ],
                    "n_levels": int(vc.shape[0]),
                }
        return {
            "status": "not_assessed_missing_metadata",
            "column": None,
            "levels": [],
            "n_levels": 0,
        }

    return {
        "status": "assessed",
        "n_cells": int(len(sub)),
        "donors": _counts(("donor", "donor_id", "Donor", "subject_id", "patient")),
        "batches": _counts(
            ("batch", "batch_id", "assay", "platform", "channel", "ProcessingMethod")
        ),
        "samples": _counts(("sample", "sample_id", "biosample_id", "library_id", "Source")),
    }


def _literature_hooks(row: dict[str, Any]) -> dict[str, Any]:
    """Surface literature-related fields if present; never invent PMIDs."""
    pmids = _split_markers(row.get("literature_pmids") or row.get("pmids"))
    sources = _split_markers(row.get("marker_provenance_sources") or row.get("sources"))
    status = row.get("marker_provenance_status") or row.get("literature_status") or "not_assessed"
    if not pmids and not sources:
        return {
            "status": "not_assessed",
            "pmids": [],
            "sources": [],
            "note": "No literature fields on this evidence row; run celltypepilot literature for ad-hoc checks.",
        }
    return {
        "status": str(status),
        "pmids": pmids,
        "sources": sources,
        "note": "Provenance/literature is display-only; human must verify primary sources.",
    }


def build_cluster_review_panel(
    *,
    cluster: str,
    evidence: pd.DataFrame,
    state_results: pd.DataFrame | None = None,
    novelty_results: pd.DataFrame | None = None,
    adata=None,
    cluster_key: str = "ctp_cl_id",
    pending_override: dict | None = None,
    audit_events: list[dict] | None = None,
) -> dict[str, Any]:
    """Assemble Identity × State × Novelty panel for one cluster."""
    row = _row_for_cluster(evidence, cluster)
    state_row = _row_for_cluster(state_results, cluster) if state_results is not None else {}
    novelty_row = _row_for_cluster(novelty_results, cluster) if novelty_results is not None else {}

    supporting = _split_markers(
        row.get("pos_supporting_markers")
        or row.get("supporting_markers")
        or row.get("pos_present_markers")
    )
    opposing = _split_markers(
        row.get("neg_expressed_markers")
        or row.get("neg_conflict_markers")
        or row.get("conflicting_markers")
    )
    silent = _split_markers(row.get("pos_silent_markers"))
    missing = _split_markers(row.get("pos_missing_markers"))

    identity = {
        "axis": "identity",
        "cell_type": row.get("cell_type", "Unknown"),
        "candidate_cell_type": row.get("candidate_cell_type", row.get("cell_type")),
        "decision": row.get("decision", "unknown"),
        "abstain_reason": row.get("abstain_reason", ""),
        "cl_id": row.get("cl_id", ""),
        "evidence_score": float(row.get("combined_score") or row.get("evidence_score") or 0.0),
        "critic_confidence": row.get("critic_confidence", "unknown"),
        "critic_flags": row.get("critic_flags", "PASS"),
        "critic_notes": row.get("critic_notes", ""),
        "supporting_markers": supporting,
        "opposing_markers": opposing,
        "silent_markers": silent,
        "missing_markers": missing,
        "pct_overlap": float(row.get("pct_overlap") or 0.0),
        "neighbor_candidates": _neighbor_candidates(evidence, cluster),
    }

    state = {
        "axis": "state",
        "state_candidate": state_row.get("state_candidate")
        or row.get("cell_state_candidate")
        or "Unknown",
        "state_decision": state_row.get("state_decision")
        or row.get("state_decision")
        or "not_assessed",
        "state_score": float(state_row.get("state_score") or row.get("state_score") or 0.0),
        "state_confidence": state_row.get("state_confidence")
        or row.get("state_confidence")
        or "needs_review",
        "state_evidence": state_row.get("state_evidence") or row.get("state_evidence") or "",
        "supporting_markers": _split_markers(
            state_row.get("pos_supporting_markers") or state_row.get("supporting_markers")
        ),
        "missing_markers": _split_markers(state_row.get("pos_missing_markers")),
        "note": "State is independent; cannot overwrite or rescue identity.",
    }

    novelty = {
        "axis": "novelty",
        "novelty_decision": novelty_row.get("novelty_decision")
        or row.get("novelty_decision")
        or "not_assessed",
        "novelty_score": float(novelty_row.get("novelty_score") or row.get("novelty_score") or 0.0),
        "top_unmapped_markers": _split_markers(
            novelty_row.get("top_unmapped_markers") or row.get("top_unmapped_markers")
        ),
        "alternative_explanations": novelty_row.get("alternative_explanations")
        or row.get("alternative_explanations")
        or "",
        "recommended_next_actions": novelty_row.get("recommended_next_actions")
        or row.get("recommended_next_actions")
        or "",
        "note": "Novelty/OOD is a review priority axis, not a validated discovery.",
    }

    return {
        "schema_version": REVIEW_PANEL_SCHEMA,
        "cluster": str(cluster),
        "axes": {
            "identity": identity,
            "state": state,
            "novelty": novelty,
        },
        "donor_batch_strata": _donor_batch_strata(adata, cluster, cluster_key),
        "literature": _literature_hooks(row),
        "pending_override": pending_override,
        "audit_events": audit_events or [],
        "edit_policy": {
            "append_only_audit": True,
            "derived_artifacts_stale_after_apply": True,
            "resign_after_regenerate": True,
        },
    }


def load_optional_csv(path: Path) -> pd.DataFrame:
    if path.is_file():
        return pd.read_csv(path)
    return pd.DataFrame()
