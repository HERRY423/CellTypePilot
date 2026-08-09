"""Composable QC diagnostic axes that cannot rescue identity labels.

Diagnostic axes (low RNA, high mito, doublet, ambient-RNA, sample enrichment,
batch sensitivity) are independent review signals.  Missing metadata or tool
output yields ``not_assessed_*`` statuses — never ``clean`` / ``pass`` as a claim
that an artifact is absent.

External doublet / ambient tools may contribute scores via artifact contracts.
Their presence can escalate review priority; they never overwrite, upgrade, or
rescue identity decisions produced by the annotation critic.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

QC_SCHEMA = "celltypepilot.qc-diagnostics.v1"
IDENTITY_INVARIANT = (
    "QC diagnostics are a review axis only. They must not change identity "
    "labels, decisions, CL IDs, or abstain reasons. They cannot rescue Unknown."
)

# Canonical diagnostic names for composable contracts.
QC_AXES = (
    "low_rna",
    "high_mito",
    "doublet",
    "ambient_rna",
    "sample_enrichment",
    "batch_sensitivity",
)

# Forbidden outcome tokens: absence of assessment must never be sold as clean.
FORBIDDEN_CLEAN_TOKENS = frozenset(
    {
        "clean",
        "qc_clean",
        "no_artifact",
        "artifact_free",
        "doublet_free",
        "ambient_free",
        "batch_free",
    }
)

# Common obs column aliases (detection only; missing → not_assessed).
LOW_RNA_ALIASES = (
    "n_genes_by_counts",
    "n_genes",
    "nGene",
    "n_genes_detected",
    "total_features_by_counts",
)
HIGH_MITO_ALIASES = (
    "pct_counts_mt",
    "percent_mito",
    "pct_mito",
    "mt_frac",
    "percent.mt",
    "mito_percent",
    "qc_mitoFraction",
)
DOUBLET_ALIASES = (
    "doublet_score",
    "doublet_scores",
    "scrublet_score",
    "predicted_doublet",
    "is_doublet",
    "scDblFinder_score",
    "scdblfinder_class",
)
AMBIENT_ALIASES = (
    "ambient_score",
    "soup_fraction",
    "soup_contamination",
    "ambient_rna_fraction",
    "decontX_contamination",
    "contamination",
)
SAMPLE_ALIASES = (
    "sample",
    "sample_id",
    "Sample",
    "biosample_id",
    "library_id",
    "channel",
    "Source",
)
BATCH_ALIASES = (
    "batch",
    "batch_id",
    "Batch",
    "assay",
    "platform",
    "chemistry",
    "ProcessingMethod",
)


class QCDiagnosticError(ValueError):
    """Raised for invalid QC contracts or identity-invariant violations."""


@dataclass(frozen=True)
class DiagnosticContract:
    """One composable diagnostic axis evaluation."""

    axis: str
    status: str
    flag: str
    n_cells_assessed: int
    n_cells_flagged: int
    n_cells_missing: int
    source: str
    details: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return {
            "axis": self.axis,
            "status": self.status,
            "flag": self.flag,
            "n_cells_assessed": self.n_cells_assessed,
            "n_cells_flagged": self.n_cells_flagged,
            "n_cells_missing": self.n_cells_missing,
            "source": self.source,
            "details": self.details,
            "identity_effect": "none",
            "can_rescue_identity": False,
        }


def _first_present(columns: Any, candidates: tuple[str, ...]) -> str | None:
    for name in candidates:
        if name in columns:
            return name
    return None


def _not_assessed(
    axis: str,
    *,
    reason: str,
    source: str = "none",
    details: dict[str, Any] | None = None,
) -> DiagnosticContract:
    status = reason if reason.startswith("not_assessed") else f"not_assessed_{reason}"
    return DiagnosticContract(
        axis=axis,
        status=status,
        flag="NOT_ASSESSED",
        n_cells_assessed=0,
        n_cells_flagged=0,
        n_cells_missing=0,
        source=source,
        details=details or {"reason": reason},
    )


def _validate_no_clean_claim(status: str, flag: str) -> None:
    tokens = {status.lower(), flag.lower()}
    if tokens & FORBIDDEN_CLEAN_TOKENS:
        raise QCDiagnosticError(
            f"Forbidden clean-claim token in QC status/flag: {status!r}/{flag!r}. "
            "Missing evidence must use not_assessed_*."
        )
    if "clean" in status.lower() or "clean" in flag.lower():
        raise QCDiagnosticError(
            f"QC diagnostics forbid 'clean' language ({status!r}/{flag!r}); use not_assessed_*."
        )


def load_external_tool_table(
    path: str | Path,
    *,
    axis: str,
    cell_id_column: str = "cell_id",
    score_column: str | None = None,
    flag_column: str | None = None,
    threshold: float | None = None,
    higher_is_worse: bool = True,
) -> pd.DataFrame:
    """Load an external doublet/ambient tool CSV into a cell-level diagnostic table.

    Returns columns: cell_id, score (optional), flagged (bool|NA), axis, source.
    Missing score/flag cells remain NA so strata become missing_diagnostic / not_assessed.
    """
    if axis not in {"doublet", "ambient_rna"}:
        raise QCDiagnosticError(
            f"External tool tables are supported for doublet/ambient_rna only; got {axis!r}"
        )
    table_path = Path(path)
    if not table_path.is_file():
        raise QCDiagnosticError(f"External {axis} table not found: {table_path}")
    frame = pd.read_csv(table_path, dtype={cell_id_column: str})
    if cell_id_column not in frame.columns:
        raise QCDiagnosticError(
            f"External {axis} table missing cell id column {cell_id_column!r}"
        )
    out = pd.DataFrame({"cell_id": frame[cell_id_column].astype(str)})
    out["axis"] = axis
    out["source"] = f"external_file:{table_path.name}"

    score = None
    if score_column and score_column in frame.columns:
        score = pd.to_numeric(frame[score_column], errors="coerce")
        out["score"] = score
    elif score_column:
        raise QCDiagnosticError(f"score_column {score_column!r} not in external table")

    flagged = pd.Series(pd.NA, index=out.index, dtype="boolean")
    if flag_column and flag_column in frame.columns:
        raw = frame[flag_column]
        positives = {"true", "1", "yes", "doublet", "doublet_found", "high", "flagged"}
        flagged = raw.astype(str).str.strip().str.lower().isin(positives)
        flagged = flagged.mask(raw.isna())
    elif score is not None and threshold is not None:
        if higher_is_worse:
            flagged = score.ge(threshold)
        else:
            flagged = score.le(threshold)
        flagged = flagged.mask(score.isna())
    elif score is not None:
        # Score present but no threshold: assessed_descriptive only at cell level.
        out["flagged"] = pd.Series(pd.NA, index=out.index, dtype="boolean")
        return out
    else:
        raise QCDiagnosticError(
            f"External {axis} table needs flag_column or score_column(+threshold)"
        )

    out["flagged"] = flagged.astype("boolean")
    return out


def _threshold_flag(
    values: pd.Series,
    *,
    threshold: float,
    direction: str,
) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    if direction == "le":
        flag = numeric.le(threshold)
    else:
        flag = numeric.ge(threshold)
    return flag.mask(numeric.isna())


def _count_true(flag: pd.Series) -> int:
    """Count True without treating NA as False via deprecated fillna downcast."""
    return int(np.sum(np.asarray(flag) == True))  # noqa: E712


def assess_low_rna(
    obs: pd.DataFrame,
    *,
    column: str | None = None,
    threshold: float = 200.0,
) -> DiagnosticContract:
    key = column or _first_present(obs.columns, LOW_RNA_ALIASES)
    if not key:
        return _not_assessed("low_rna", reason="missing_metadata", source="obs")
    flag = _threshold_flag(obs[key], threshold=threshold, direction="le")
    assessed = int(flag.notna().sum())
    flagged = _count_true(flag)
    missing = int(flag.isna().sum())
    if assessed == 0:
        return _not_assessed(
            "low_rna",
            reason="missing_cell_values",
            source=f"obs:{key}",
            details={"column": key, "threshold": threshold},
        )
    status = "assessed_with_flags" if flagged else "assessed_no_cells_flagged"
    flag_label = "LOW_RNA_PRESENT" if flagged else "NO_CELLS_FLAGGED"
    _validate_no_clean_claim(status, flag_label)
    return DiagnosticContract(
        axis="low_rna",
        status=status,
        flag=flag_label,
        n_cells_assessed=assessed,
        n_cells_flagged=flagged,
        n_cells_missing=missing,
        source=f"obs:{key}",
        details={
            "column": key,
            "threshold": threshold,
            "direction": "le",
            "note": "NO_CELLS_FLAGGED is not a claim that low-RNA cells are absent from biology.",
        },
    )


def assess_high_mito(
    obs: pd.DataFrame,
    *,
    column: str | None = None,
    threshold: float = 0.2,
) -> DiagnosticContract:
    key = column or _first_present(obs.columns, HIGH_MITO_ALIASES)
    if not key:
        return _not_assessed("high_mito", reason="missing_metadata", source="obs")
    # Heuristic: values > 1 are treated as percent (0-100).
    values = pd.to_numeric(obs[key], errors="coerce")
    thr = threshold * 100.0 if float(values.dropna().median() or 0) > 1.5 else threshold
    flag = _threshold_flag(obs[key], threshold=thr, direction="ge")
    assessed = int(flag.notna().sum())
    flagged = _count_true(flag)
    missing = int(flag.isna().sum())
    if assessed == 0:
        return _not_assessed(
            "high_mito",
            reason="missing_cell_values",
            source=f"obs:{key}",
            details={"column": key, "threshold": thr},
        )
    status = "assessed_with_flags" if flagged else "assessed_no_cells_flagged"
    flag_label = "HIGH_MITO_PRESENT" if flagged else "NO_CELLS_FLAGGED"
    _validate_no_clean_claim(status, flag_label)
    return DiagnosticContract(
        axis="high_mito",
        status=status,
        flag=flag_label,
        n_cells_assessed=assessed,
        n_cells_flagged=flagged,
        n_cells_missing=missing,
        source=f"obs:{key}",
        details={
            "column": key,
            "threshold": thr,
            "direction": "ge",
            "note": "NO_CELLS_FLAGGED is not a claim that high-mito cells are absent.",
        },
    )


def assess_tool_axis(
    obs: pd.DataFrame,
    *,
    axis: str,
    external: pd.DataFrame | None = None,
    obs_aliases: tuple[str, ...] = (),
    threshold: float | None = None,
) -> DiagnosticContract:
    """Assess doublet or ambient axis from external tool table and/or obs columns."""
    if axis not in {"doublet", "ambient_rna"}:
        raise QCDiagnosticError(f"assess_tool_axis expects doublet/ambient_rna, got {axis}")

    n = len(obs)
    index = obs.index.astype(str)
    flagged = pd.Series(pd.NA, index=index, dtype="boolean")
    source_parts: list[str] = []

    # Prefer external tool output when provided.
    if external is not None and not external.empty:
        if "cell_id" not in external.columns or "flagged" not in external.columns:
            return _not_assessed(
                axis,
                reason="invalid_external_contract",
                source="external",
                details={"required_columns": ["cell_id", "flagged"]},
            )
        ext = external.copy()
        ext["cell_id"] = ext["cell_id"].astype(str)
        mapped = ext.drop_duplicates("cell_id").set_index("cell_id")["flagged"]
        aligned = mapped.reindex(index)
        flagged = aligned.astype("boolean")
        source_parts.append(str(ext["source"].iloc[0]) if "source" in ext else "external")

    key = _first_present(obs.columns, obs_aliases)
    if key and flagged.isna().all():
        values = obs[key]
        if values.dtype == bool or set(pd.Series(values).dropna().unique()) <= {0, 1, True, False}:
            f = values.astype("boolean")
        else:
            if threshold is None:
                return DiagnosticContract(
                    axis=axis,
                    status="assessed_descriptive_scores_only",
                    flag="SCORES_PRESENT_NO_THRESHOLD",
                    n_cells_assessed=int(pd.to_numeric(values, errors="coerce").notna().sum()),
                    n_cells_flagged=0,
                    n_cells_missing=int(pd.to_numeric(values, errors="coerce").isna().sum()),
                    source=f"obs:{key}",
                    details={
                        "column": key,
                        "note": "Scores recorded without a predeclared threshold; not a cleanliness claim.",
                    },
                )
            f = _threshold_flag(values, threshold=threshold, direction="ge").astype("boolean")
        flagged = f
        source_parts.append(f"obs:{key}")

    if flagged.isna().all() and not source_parts:
        return _not_assessed(axis, reason="missing_tool_output_and_metadata", source="none")

    assessed = int(flagged.notna().sum())
    n_flagged = _count_true(flagged)
    missing = int(flagged.isna().sum())
    if assessed == 0:
        return _not_assessed(
            axis,
            reason="missing_cell_values",
            source="+".join(source_parts) or "none",
        )
    status = "assessed_with_flags" if n_flagged else "assessed_no_cells_flagged"
    flag_label = f"{axis.upper()}_PRESENT" if n_flagged else "NO_CELLS_FLAGGED"
    # Normalize ambient_rna flag token
    if axis == "ambient_rna" and n_flagged:
        flag_label = "AMBIENT_RNA_PRESENT"
    if axis == "doublet" and n_flagged:
        flag_label = "DOUBLET_PRESENT"
    _validate_no_clean_claim(status, flag_label)
    return DiagnosticContract(
        axis=axis,
        status=status,
        flag=flag_label,
        n_cells_assessed=assessed,
        n_cells_flagged=n_flagged,
        n_cells_missing=missing,
        source="+".join(source_parts),
        details={
            "n_cells": n,
            "note": (
                "NO_CELLS_FLAGGED means no positive flags among assessed cells; "
                "it is not a claim that doublets/ambient RNA are biologically absent."
            ),
        },
    )


def assess_sample_enrichment(
    obs: pd.DataFrame,
    *,
    cluster_key: str | None,
    sample_key: str | None = None,
    dominant_fraction_threshold: float = 0.8,
) -> DiagnosticContract:
    from .robustness import sample_enrichment_diagnostics

    sample_col = sample_key or _first_present(obs.columns, SAMPLE_ALIASES)
    if not cluster_key or cluster_key not in obs.columns:
        return _not_assessed(
            "sample_enrichment",
            reason="missing_cluster_metadata",
            source="obs",
        )
    if not sample_col:
        return _not_assessed(
            "sample_enrichment",
            reason="missing_sample_metadata",
            source="obs",
        )
    table = sample_enrichment_diagnostics(
        obs,
        cluster_key=cluster_key,
        sample_key=sample_col,
        dominant_fraction_threshold=dominant_fraction_threshold,
    )
    # Normalize robustness module's missing row to not_assessed_*
    if len(table) == 1 and str(table.iloc[0].get("status", "")).startswith("missing"):
        return _not_assessed(
            "sample_enrichment",
            reason="missing_cluster_or_sample_metadata",
            source="obs",
        )
    n_flagged = int((table["flag"].astype(str) == "SAMPLE_ENRICHED").sum()) if not table.empty else 0
    status = "assessed_with_flags" if n_flagged else "assessed_descriptive"
    flag_label = "SAMPLE_ENRICHMENT_PRESENT" if n_flagged else "NO_CLUSTER_SAMPLE_ENRICHED"
    _validate_no_clean_claim(status, flag_label)
    return DiagnosticContract(
        axis="sample_enrichment",
        status=status,
        flag=flag_label,
        n_cells_assessed=int(table["n_cells"].sum()) if "n_cells" in table else 0,
        n_cells_flagged=n_flagged,
        n_cells_missing=0,
        source=f"obs:{sample_col}",
        details={
            "cluster_key": cluster_key,
            "sample_key": sample_col,
            "dominant_fraction_threshold": dominant_fraction_threshold,
            "n_clusters": int(len(table)),
            "n_clusters_flagged": n_flagged,
            "note": "Descriptive cluster–sample dominance; not batch correction.",
        },
    )


def assess_batch_sensitivity_axis(
    obs: pd.DataFrame,
    *,
    batch_key: str | None = None,
) -> DiagnosticContract:
    """Presence of batch metadata is required for sensitivity audits.

    Full donor-level performance ranges require locked benchmark predictions and
    are computed in the release path. Annotation-time this axis only reports
    whether batch/platform metadata exists to enable such audits.
    """
    key = batch_key or _first_present(obs.columns, BATCH_ALIASES)
    if not key:
        return _not_assessed(
            "batch_sensitivity",
            reason="missing_batch_metadata",
            source="obs",
            details={
                "note": (
                    "No batch/platform column; batch sensitivity is not_assessed. "
                    "This is not evidence of batch-free data."
                )
            },
        )
    n_levels = int(obs[key].astype(str).nunique(dropna=True))
    return DiagnosticContract(
        axis="batch_sensitivity",
        status="metadata_present_performance_not_estimated",
        flag="BATCH_METADATA_PRESENT",
        n_cells_assessed=int(obs[key].notna().sum()),
        n_cells_flagged=0,
        n_cells_missing=int(obs[key].isna().sum()),
        source=f"obs:{key}",
        details={
            "column": key,
            "n_levels": n_levels,
            "note": (
                "Metadata enables later batch-sensitivity estimation on locked "
                "benchmark predictions; annotation-time does not claim robustness."
            ),
        },
    )


def assert_identity_invariant(
    identity_before: Mapping[str, Any],
    identity_after: Mapping[str, Any],
) -> None:
    """Fail closed if QC assembly mutated identity fields."""
    keys = (
        "cell_type",
        "decision",
        "cl_id",
        "abstain_reason",
        "ctp_cell_type",
        "ctp_decision",
        "ctp_cl_id",
    )
    for key in keys:
        if key in identity_before or key in identity_after:
            if identity_before.get(key) != identity_after.get(key):
                raise QCDiagnosticError(
                    f"Identity invariant violated for {key!r}: QC diagnostics "
                    "must not change identity labels or decisions."
                )


def assemble_qc_diagnostics(
    adata,
    *,
    cluster_key: str | None = None,
    doublet_table: pd.DataFrame | None = None,
    ambient_table: pd.DataFrame | None = None,
    low_rna_threshold: float = 200.0,
    high_mito_threshold: float = 0.2,
    doublet_threshold: float | None = 0.25,
    ambient_threshold: float | None = 0.2,
    sample_key: str | None = None,
    batch_key: str | None = None,
    identity_snapshot: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a composable, fail-closed QC diagnostic artifact for one AnnData.

    Parameters
    ----------
    identity_snapshot:
        Optional mapping of identity fields taken *before* QC assembly. When
        provided, they are compared to the same fields after assembly (must match).
    """
    obs = adata.obs
    contracts = [
        assess_low_rna(obs, threshold=low_rna_threshold),
        assess_high_mito(obs, threshold=high_mito_threshold),
        assess_tool_axis(
            obs,
            axis="doublet",
            external=doublet_table,
            obs_aliases=DOUBLET_ALIASES,
            threshold=doublet_threshold,
        ),
        assess_tool_axis(
            obs,
            axis="ambient_rna",
            external=ambient_table,
            obs_aliases=AMBIENT_ALIASES,
            threshold=ambient_threshold,
        ),
        assess_sample_enrichment(obs, cluster_key=cluster_key, sample_key=sample_key),
        assess_batch_sensitivity_axis(obs, batch_key=batch_key),
    ]

    for contract in contracts:
        _validate_no_clean_claim(contract.status, contract.flag)

    axes = {c.axis: c.to_dict() for c in contracts}
    # Ensure every canonical axis is present even if assembly list drifts.
    for axis in QC_AXES:
        if axis not in axes:
            axes[axis] = _not_assessed(axis, reason="not_in_assembly").to_dict()

    any_flagged = any(
        c.n_cells_flagged > 0 or c.flag.endswith("_PRESENT") for c in contracts
    )
    any_assessed = any(not c.status.startswith("not_assessed") for c in contracts)

    if not any_assessed:
        rollup = "not_assessed"
        rollup_flag = "NOT_ASSESSED"
    elif any_flagged:
        rollup = "diagnostics_flagged"
        rollup_flag = "QC_FLAGS_PRESENT"
    else:
        # Explicitly avoid "clean"
        rollup = "assessed_no_flags_among_available_axes"
        rollup_flag = "NO_FLAGS_AMONG_ASSESSED_AXES"

    _validate_no_clean_claim(rollup, rollup_flag)

    if identity_snapshot is not None:
        # Re-read identity-like columns from obs if present.
        after = {
            key: (
                list(obs[key].astype(str).head(5))
                if key in obs.columns
                else identity_snapshot.get(key)
            )
            for key in identity_snapshot
        }
        # Only compare non-obs snapshot keys (cluster-level summary strings).
        for key, value in identity_snapshot.items():
            if key not in obs.columns and after.get(key) != value:
                raise QCDiagnosticError(
                    f"Identity invariant violated for snapshot key {key!r}"
                )

    report = {
        "schema_version": QC_SCHEMA,
        "identity_invariant": IDENTITY_INVARIANT,
        "can_rescue_identity": False,
        "identity_effect": "none",
        "forbidden_clean_tokens": sorted(FORBIDDEN_CLEAN_TOKENS),
        "missing_metadata_policy": "not_assessed_never_clean",
        "rollup_status": rollup,
        "rollup_flag": rollup_flag,
        "axes": axes,
        "axis_order": list(QC_AXES),
        "agent_guidance": (
            "Present QC as a diagnostic review axis. Do not upgrade Unknown identity "
            "because QC looks good. Missing axes are not_assessed, not clean."
        ),
    }
    return report


def write_qc_diagnostics(
    report: dict[str, Any],
    output_dir: str | Path,
) -> dict[str, Path]:
    """Write JSON + long CSV for QC diagnostics."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "qc_diagnostics.json"
    csv_path = output / "qc_diagnostics.csv"
    json_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    rows = []
    for axis in report.get("axis_order", QC_AXES):
        payload = report["axes"].get(axis, {})
        rows.append(
            {
                "axis": axis,
                "status": payload.get("status"),
                "flag": payload.get("flag"),
                "n_cells_assessed": payload.get("n_cells_assessed"),
                "n_cells_flagged": payload.get("n_cells_flagged"),
                "n_cells_missing": payload.get("n_cells_missing"),
                "source": payload.get("source"),
                "can_rescue_identity": False,
                "identity_effect": "none",
            }
        )
    pd.DataFrame(rows).to_csv(csv_path, index=False, lineterminator="\n")
    return {"qc_diagnostics_json": json_path, "qc_diagnostics_csv": csv_path}


def build_benchmark_diagnostic_specs(
    diagnostics: Mapping[str, Mapping[str, object]] | None,
) -> dict[str, dict[str, object]]:
    """Normalize registry cohort diagnostics into composable axis specs.

    Missing axes are left empty so ``qc_stratified_performance`` returns
    ``not_assessed_missing_predeclared_column`` rather than assuming clean cells.
    """
    base = {
        "low_quality": {},
        "doublet": {},
        "ambient_rna": {},
    }
    if not diagnostics:
        return base
    out = dict(base)
    for key in base:
        if key in diagnostics and diagnostics[key]:
            out[key] = dict(diagnostics[key])
    # Allow registry "low_quality" mito specs; also accept high_mito alias.
    if "high_mito" in diagnostics and diagnostics["high_mito"] and not out["low_quality"]:
        out["low_quality"] = dict(diagnostics["high_mito"])
    return out
