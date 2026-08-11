"""Backend-neutral candidate contract for identity generation.

Prediction systems generate candidates; they do not publish CellTypePilot's
final identity.  This module normalizes heterogeneous backend output into one
auditable table.  Marker scoring is represented in the same table, but its
role is permanently ``evidence_only``.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd

from .identity_contract import resolve_identity_label

CANDIDATE_SCHEMA = "celltypepilot.backend-candidates.v1"


class CandidateContractError(ValueError):
    """Raised when a candidate artifact cannot be normalized safely."""


BACKEND_ALIASES = {
    "celltypist": "celltypist",
    "popv": "popv",
    "singler": "singler",
    "single_r": "singler",
    "scanvi": "scanvi",
    "scanvi/scvi": "scanvi",
    "custom_reference": "custom_reference",
    "reference": "custom_reference",
    "knn": "knn",
    "correlation": "correlation",
    "llm": "llm",
    "marker": "marker_evidence",
    "marker_scorer": "marker_evidence",
    "marker_evidence": "marker_evidence",
}

BACKEND_POLICY = {
    "celltypist": ("reference_linear", "decision_candidate"),
    "popv": ("reference_ensemble", "decision_candidate"),
    "singler": ("reference_correlation", "decision_candidate"),
    "scanvi": ("reference_latent", "decision_candidate"),
    "custom_reference": ("custom_reference", "decision_candidate"),
    "knn": ("custom_reference", "decision_candidate"),
    "correlation": ("custom_reference", "decision_candidate"),
    "llm": ("llm_hypothesis", "hypothesis_only"),
    "marker_evidence": ("marker_evidence", "evidence_only"),
}

OUTPUT_COLUMNS = [
    "schema_version",
    "cluster",
    "backend",
    "backend_family",
    "independence_group",
    "decision_role",
    "raw_cell_type",
    "canonical_cell_type",
    "raw_cl_id",
    "canonical_cl_id",
    "identity_resolution",
    "identity_resolved",
    "raw_score",
    "score_semantics",
    "rank",
    "source_artifact",
    "source_version",
    "claim_boundary",
]


def _backend_name(value: Any) -> str:
    key = str(value or "").strip().casefold().replace("-", "_").replace(" ", "_")
    backend = BACKEND_ALIASES.get(key)
    if backend is None:
        allowed = ", ".join(sorted(BACKEND_POLICY))
        raise CandidateContractError(f"Unsupported candidate backend {value!r}; allowed: {allowed}")
    return backend


def _read_candidate_artifact(path: str | Path) -> pd.DataFrame:
    candidate_path = Path(path)
    if not candidate_path.exists():
        raise CandidateContractError(f"Candidate artifact not found: {candidate_path}")
    if candidate_path.suffix.casefold() == ".csv":
        return pd.read_csv(candidate_path)
    if candidate_path.suffix.casefold() == ".json":
        try:
            payload = json.loads(candidate_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise CandidateContractError(f"Invalid candidate JSON: {exc}") from exc
        if isinstance(payload, dict):
            payload = payload.get("candidates", payload.get("predictions"))
        if not isinstance(payload, list):
            raise CandidateContractError(
                "Candidate JSON must be a list or contain a candidates/predictions list"
            )
        return pd.DataFrame(payload)
    raise CandidateContractError("Candidate artifacts must be CSV or JSON")


def aggregate_cell_candidates(frame: pd.DataFrame, cluster_assignments: pd.Series) -> pd.DataFrame:
    """Aggregate one top-1 cell prediction per backend to cluster vote fractions."""
    cell_col = _pick_column(frame, ("cell_id", "barcode", "obs_name"), required=True)
    label_col = _pick_column(
        frame,
        ("cell_type", "candidate_cell_type", "predicted_label", "prediction", "label"),
        required=True,
    )
    backend_col = _pick_column(frame, ("backend", "method", "model_backend"), required=True)
    rank_col = _pick_column(frame, ("rank", "ref_rank", "candidate_rank"), required=False)
    cl_col = _pick_column(frame, ("cl_id", "candidate_cl_id", "ontology_id"), required=False)
    work = frame.copy()
    if rank_col:
        ranks = pd.to_numeric(work[rank_col], errors="coerce")
        work = work[ranks == 1].copy()
    work[cell_col] = work[cell_col].astype(str)
    work[backend_col] = work[backend_col].astype(str)
    if work.duplicated([cell_col, backend_col]).any():
        raise CandidateContractError(
            "Cell-level candidate artifacts require exactly one top-1 row per cell/backend"
        )
    assignments = pd.Series(
        cluster_assignments.astype(str).to_numpy(),
        index=cluster_assignments.index.astype(str),
    )
    unknown = sorted(set(work[cell_col]) - set(assignments.index))
    if unknown:
        raise CandidateContractError(
            f"Candidate artifact contains {len(unknown)} cell IDs absent from query obs; "
            f"examples: {unknown[:3]}"
        )
    work["cluster"] = work[cell_col].map(assignments)
    group_cols = ["cluster", backend_col, label_col]
    if cl_col:
        group_cols.append(cl_col)
    counts = work.groupby(group_cols, dropna=False).size().rename("n_cells").reset_index()
    totals = counts.groupby(["cluster", backend_col])["n_cells"].transform("sum")
    counts["score"] = counts["n_cells"] / totals
    counts["score_semantics"] = "within_backend_cluster_top1_vote_fraction_not_probability"
    counts["rank"] = (
        counts.groupby(["cluster", backend_col])["score"]
        .rank(ascending=False, method="first")
        .astype(int)
    )
    return counts


def _pick_column(frame: pd.DataFrame, names: tuple[str, ...], *, required: bool) -> str | None:
    match = next((name for name in names if name in frame.columns), None)
    if required and match is None:
        raise CandidateContractError(f"Candidate table requires one of columns: {', '.join(names)}")
    return match


def normalize_candidate_table(
    frame: pd.DataFrame,
    resolver: dict,
    *,
    source_artifact: str = "in_memory",
    default_backend: str | None = None,
    source_version: str = "",
) -> pd.DataFrame:
    """Normalize one candidate table without interpreting backend scores.

    Scores are retained with their declared semantics.  Ranking and agreement,
    rather than cross-backend score arithmetic, drive the downstream selector.
    """
    if frame.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    cluster_col = _pick_column(frame, ("cluster", "cluster_id"), required=True)
    label_col = _pick_column(
        frame,
        ("cell_type", "candidate_cell_type", "predicted_label", "prediction", "label"),
        required=True,
    )
    backend_col = _pick_column(frame, ("backend", "method", "model_backend"), required=False)
    if backend_col is None and not default_backend:
        raise CandidateContractError("Candidate table requires backend/method or a default backend")
    score_col = _pick_column(
        frame,
        ("raw_score", "score", "probability", "confidence", "ref_score", "combined_score"),
        required=False,
    )
    rank_col = _pick_column(frame, ("rank", "ref_rank", "candidate_rank"), required=False)
    cl_col = _pick_column(frame, ("cl_id", "candidate_cl_id", "ontology_id"), required=False)
    semantics_col = _pick_column(frame, ("score_semantics",), required=False)
    version_col = _pick_column(
        frame, ("source_version", "backend_version", "model_version"), required=False
    )

    work = pd.DataFrame(
        {
            "cluster": frame[cluster_col].astype(str),
            "raw_cell_type": frame[label_col].fillna("").astype(str),
            "raw_cl_id": frame[cl_col].fillna("").astype(str) if cl_col else "",
            "raw_score": pd.to_numeric(frame[score_col], errors="coerce")
            if score_col
            else float("nan"),
        }
    )
    raw_backends = (
        frame[backend_col] if backend_col else pd.Series(default_backend, index=frame.index)
    )
    work["backend"] = [_backend_name(value) for value in raw_backends]
    work["backend_family"] = [BACKEND_POLICY[value][0] for value in work["backend"]]
    work["independence_group"] = work["backend_family"]
    work["decision_role"] = [BACKEND_POLICY[value][1] for value in work["backend"]]
    work["score_semantics"] = (
        frame[semantics_col].fillna("backend_declared_uninterpreted").astype(str).to_numpy()
        if semantics_col
        else "backend_score_not_cross_backend_comparable"
    )
    work["source_version"] = (
        frame[version_col].fillna("").astype(str).to_numpy() if version_col else source_version
    )

    resolutions = [
        resolve_identity_label(label, resolver, cl_id)
        for label, cl_id in zip(work["raw_cell_type"], work["raw_cl_id"], strict=True)
    ]
    work["canonical_cell_type"] = [item["canonical_label"] for item in resolutions]
    work["identity_resolution"] = [item["resolution"] for item in resolutions]
    work["identity_resolved"] = [
        item["resolution"] not in {"unresolved", "ambiguous_alias_unresolved"}
        and item["canonical_label"] != "Unknown"
        for item in resolutions
    ]
    cl_by_name = resolver.get("cl_by_name", {})
    work["canonical_cl_id"] = [
        cl_by_name.get(item["canonical_label"], item.get("cl_id", "")) if resolved else ""
        for item, resolved in zip(resolutions, work["identity_resolved"], strict=True)
    ]

    if rank_col:
        ranks = pd.to_numeric(frame[rank_col], errors="coerce")
        if ranks.isna().any() or (ranks < 1).any():
            raise CandidateContractError("Candidate ranks must be positive integers")
        work["rank"] = ranks.astype(int).to_numpy()
    else:
        order = work.assign(_row_order=range(len(work)))
        order = order.sort_values(
            ["cluster", "backend", "raw_score", "_row_order"],
            ascending=[True, True, False, True],
            na_position="last",
        )
        order["rank"] = order.groupby(["cluster", "backend"]).cumcount() + 1
        work["rank"] = order.sort_index()["rank"].astype(int)

    duplicate_top = work[(work["rank"] == 1)].duplicated(["cluster", "backend"], keep=False)
    if duplicate_top.any():
        raise CandidateContractError(
            "Each backend may provide only one rank-1 candidate per cluster"
        )

    work["schema_version"] = CANDIDATE_SCHEMA
    work["source_artifact"] = source_artifact
    work["claim_boundary"] = (
        "backend_candidate_not_final_identity; scores_not_cross_backend_probabilities"
    )
    return work[OUTPUT_COLUMNS].sort_values(["cluster", "backend", "rank"]).reset_index(drop=True)


def load_candidate_artifacts(
    paths: Iterable[str | Path],
    resolver: dict,
    cluster_assignments: pd.Series | None = None,
) -> pd.DataFrame:
    """Load and concatenate governed candidate artifacts."""
    tables = []
    for path in paths:
        frame = _read_candidate_artifact(path)
        if "cluster" not in frame and "cluster_id" not in frame:
            if cluster_assignments is None:
                raise CandidateContractError(
                    "Cell-level candidate artifact requires query cluster assignments"
                )
            frame = aggregate_cell_candidates(frame, cluster_assignments)
        tables.append(
            normalize_candidate_table(frame, resolver, source_artifact=str(Path(path).resolve()))
        )
    return concatenate_candidates(tables)


def marker_scores_as_evidence(scores: pd.DataFrame, resolver: dict) -> pd.DataFrame:
    """Expose marker rankings in the contract without granting a voting role."""
    if scores.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    frame = scores.rename(columns={"combined_score": "raw_score"}).copy()
    frame["backend"] = "marker_evidence"
    frame["score_semantics"] = "deterministic_marker_evidence_ranking_not_probability"
    return normalize_candidate_table(frame, resolver, source_artifact="marker_scores.csv")


def reference_scores_as_candidates(scores: pd.DataFrame, resolver: dict) -> pd.DataFrame:
    """Normalize an in-process CellTypist/scANVI/custom-reference result."""
    if scores.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    backend = str(scores.attrs.get("reference_backend", "custom_reference"))
    frame = scores.rename(columns={"ref_score": "raw_score", "ref_rank": "rank"}).copy()
    frame["backend"] = backend
    frame["score_semantics"] = "backend_cluster_support_not_cross_backend_probability"
    return normalize_candidate_table(frame, resolver, source_artifact="reference_scores.csv")


def concatenate_candidates(tables: Iterable[pd.DataFrame]) -> pd.DataFrame:
    nonempty = [table for table in tables if table is not None and not table.empty]
    if not nonempty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    combined = pd.concat(nonempty, ignore_index=True)
    return (
        combined[OUTPUT_COLUMNS].sort_values(["cluster", "backend", "rank"]).reset_index(drop=True)
    )
