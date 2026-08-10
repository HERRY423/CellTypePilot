"""Donor-aware robustness diagnostics for public benchmark releases.

The functions in this module operate on locked out-of-fold predictions.  Cells
are retained for classification metrics, but uncertainty and method comparison
use donors as the independent unit.  Diagnostic absence is represented as an
explicit ``not_assessed`` row rather than silently interpreted as a clean result.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd

from .benchmark import BenchmarkValidationError, _classification_metrics
from .bootstrap import grouped_bootstrap_metric_ci

REPORT_METRICS = (
    "accuracy",
    "macro_f1",
    "balanced_accuracy",
    "coverage",
    "selective_accuracy",
)


def merge_prediction_metadata(
    truth: pd.Series,
    metadata: pd.DataFrame,
    predictions: pd.DataFrame,
) -> pd.DataFrame:
    """Join predictions to truth/metadata by cell identifier with fail-closed checks."""
    required = {"cell_id", "method", "predicted_label"}
    missing = required - set(predictions)
    if missing:
        raise BenchmarkValidationError(f"Prediction table missing columns: {sorted(missing)}")
    if predictions.duplicated(["cell_id", "method"]).any():
        raise BenchmarkValidationError("Each method may provide only one prediction per cell")

    obs = metadata.copy()
    obs.index = obs.index.astype(str)
    if obs.index.has_duplicates:
        raise BenchmarkValidationError("Metadata cell identifiers must be unique")
    truth_map = truth.copy()
    truth_map.index = truth_map.index.astype(str)
    obs["__truth__"] = truth_map.reindex(obs.index)
    if obs["__truth__"].isna().any():
        raise BenchmarkValidationError("Ground truth is missing for metadata cells")

    frame = predictions.copy()
    frame["cell_id"] = frame["cell_id"].astype(str)
    unknown = set(frame["cell_id"]) - set(obs.index)
    if unknown:
        raise BenchmarkValidationError(
            f"Predictions contain {len(unknown)} cells absent from benchmark metadata"
        )
    frame = frame.merge(
        obs, left_on="cell_id", right_index=True, how="left", validate="many_to_one"
    )
    frame["predicted_label"] = frame["predicted_label"].fillna("Unknown").astype(str)
    if "confidence" in frame:
        confidence = pd.to_numeric(frame["confidence"], errors="coerce")
        if confidence.isna().any() or ((confidence < 0) | (confidence > 1)).any():
            raise BenchmarkValidationError("confidence must be numeric and within [0, 1]")
        frame["confidence"] = confidence
    return frame


def _metrics_for_frame(frame: pd.DataFrame) -> dict[str, float | int]:
    confidence = frame["confidence"].to_numpy() if "confidence" in frame else None
    return _classification_metrics(
        frame["__truth__"].astype(str).to_numpy(),
        frame["predicted_label"].astype(str).to_numpy(),
        confidence,
    )


def evaluate_by_independent_unit(
    merged: pd.DataFrame,
    *,
    study_key: str,
    donor_key: str,
    cohort_key: str = "__cohort__",
) -> pd.DataFrame:
    """Calculate per-donor metrics; these rows are the inferential observations."""
    missing = [key for key in (study_key, donor_key) if key not in merged]
    if missing:
        raise BenchmarkValidationError(f"Missing independence columns: {missing}")
    if merged[[study_key, donor_key]].isna().any().any():
        raise BenchmarkValidationError("Study/donor identifiers must not be missing")

    frame = merged.copy()
    frame["study_id"] = frame[study_key].astype(str)
    frame["donor_id"] = frame[donor_key].astype(str)
    frame["donor_unit"] = frame["study_id"] + "::" + frame["donor_id"]
    cohort_present = cohort_key in frame
    rows = []
    for (method, donor_unit), group in frame.groupby(["method", "donor_unit"], sort=True):
        studies = group["study_id"].unique()
        if len(studies) != 1:
            raise BenchmarkValidationError(f"Donor unit {donor_unit!r} spans multiple studies")
        row = {
            "method": str(method),
            "donor_unit": str(donor_unit),
            "study_id": str(studies[0]),
            "donor_id": str(group["donor_id"].iloc[0]),
            **_metrics_for_frame(group),
        }
        if cohort_present:
            cohorts = group[cohort_key].astype(str).unique()
            row["cohort_id"] = cohorts[0] if len(cohorts) == 1 else "__mixed__"
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_independent_units(
    donor_metrics: pd.DataFrame,
    *,
    metrics: Iterable[str] = REPORT_METRICS,
    n_boot: int = 2000,
    seed: int = 42,
) -> pd.DataFrame:
    """Report equally weighted donor means and hierarchical bootstrap intervals."""
    rows = []
    for method, group in donor_metrics.groupby("method", sort=True):
        for metric in metrics:
            if metric not in group:
                continue
            finite = group.loc[np.isfinite(pd.to_numeric(group[metric], errors="coerce"))].copy()
            if finite.empty:
                rows.append(
                    {
                        "method": method,
                        "metric": metric,
                        "estimate": np.nan,
                        "ci_lower": np.nan,
                        "ci_upper": np.nan,
                        "se": np.nan,
                        "n_donors": 0,
                        "n_studies": 0,
                        "status": "not_estimable",
                    }
                )
                continue
            result = grouped_bootstrap_metric_ci(
                finite[metric].to_numpy(dtype=float),
                finite["donor_unit"].to_numpy(dtype=str),
                strata=(
                    finite["study_id"].to_numpy(dtype=str)
                    if finite["study_id"].nunique() > 1
                    else None
                ),
                n_boot=n_boot,
                seed=seed,
            )
            rows.append(
                {
                    "method": method,
                    "metric": metric,
                    "estimate": result.point_estimate,
                    "ci_lower": result.ci_lower,
                    "ci_upper": result.ci_upper,
                    "se": result.se,
                    "n_donors": int(finite["donor_unit"].nunique()),
                    "n_studies": int(finite["study_id"].nunique()),
                    "status": "estimated" if len(finite) >= 3 else "descriptive_only_lt3_donors",
                }
            )
    return pd.DataFrame(rows)


def _paired_sign_flip_pvalue(differences: np.ndarray, seed: int, n_resamples: int) -> float:
    differences = np.asarray(differences, dtype=float)
    differences = differences[np.isfinite(differences)]
    if len(differences) == 0:
        return np.nan
    observed = abs(float(differences.mean()))
    if np.allclose(differences, 0):
        return 1.0
    if len(differences) <= 16:
        signs = np.asarray(list(itertools.product((-1.0, 1.0), repeat=len(differences))))
        null = np.abs((signs * differences).mean(axis=1))
        return float(np.mean(null >= observed - 1e-15))
    else:
        rng = np.random.default_rng(seed)
        signs = rng.choice((-1.0, 1.0), size=(n_resamples, len(differences)))
        null = np.abs((signs * differences).mean(axis=1))
    return float((np.sum(null >= observed - 1e-15) + 1) / (len(null) + 1))


def _benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    output = pd.Series(np.nan, index=p_values.index, dtype=float)
    finite = p_values.dropna().astype(float)
    if finite.empty:
        return output
    ordered = finite.sort_values()
    ranks = np.arange(1, len(ordered) + 1)
    adjusted = ordered.to_numpy() * len(ordered) / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output.loc[ordered.index] = np.minimum(adjusted, 1.0)
    return output


def paired_method_comparisons(
    donor_metrics: pd.DataFrame,
    *,
    metrics: Iterable[str] = ("macro_f1", "balanced_accuracy", "coverage"),
    reference_method: str = "celltypepilot",
    n_resamples: int = 20000,
    seed: int = 42,
) -> pd.DataFrame:
    """Compare methods on paired donor units with multiplicity correction."""
    methods = sorted(donor_metrics["method"].astype(str).unique())
    if reference_method in methods:
        pairs = [(reference_method, method) for method in methods if method != reference_method]
    else:
        pairs = list(itertools.combinations(methods, 2))
    rows = []
    for metric in metrics:
        if metric not in donor_metrics:
            continue
        pivot = donor_metrics.pivot(index="donor_unit", columns="method", values=metric)
        for method_a, method_b in pairs:
            if method_a not in pivot or method_b not in pivot:
                continue
            paired = pivot[[method_a, method_b]].dropna()
            differences = (paired[method_a] - paired[method_b]).to_numpy(dtype=float)
            if len(differences) < 3:
                rows.append(
                    {
                        "method_a": method_a,
                        "method_b": method_b,
                        "metric": metric,
                        "n_paired_donors": len(differences),
                        "mean_difference": np.nan if not len(differences) else differences.mean(),
                        "median_difference": np.nan
                        if not len(differences)
                        else np.median(differences),
                        "win_fraction": np.nan
                        if not len(differences)
                        else np.mean(differences > 0),
                        "p_value": np.nan,
                        "status": "underpowered_lt3_paired_donors",
                    }
                )
                continue
            rng = np.random.default_rng(seed)
            boot = np.asarray(
                [
                    rng.choice(differences, len(differences), replace=True).mean()
                    for _ in range(2000)
                ]
            )
            rows.append(
                {
                    "method_a": method_a,
                    "method_b": method_b,
                    "metric": metric,
                    "n_paired_donors": len(differences),
                    "mean_difference": float(differences.mean()),
                    "median_difference": float(np.median(differences)),
                    "difference_ci_lower": float(np.quantile(boot, 0.025)),
                    "difference_ci_upper": float(np.quantile(boot, 0.975)),
                    "win_fraction": float(np.mean(differences > 0)),
                    "p_value": _paired_sign_flip_pvalue(differences, seed, n_resamples),
                    "status": "estimated",
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["p_value_bh"] = _benjamini_hochberg(result["p_value"])
        result["reject_fdr_005"] = result["p_value_bh"].lt(0.05).fillna(False)
    return result


def batch_sensitivity(
    merged: pd.DataFrame,
    *,
    study_key: str,
    donor_key: str,
    axes: Mapping[str, str | None],
    metric: str = "macro_f1",
    min_donors_per_level: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Describe donor-level performance ranges across platform/batch/condition axes."""
    level_rows = []
    summary_rows = []
    for axis_name, column in axes.items():
        if not column or column not in merged:
            summary_rows.append(
                {
                    "axis": axis_name,
                    "status": "not_assessed_missing_metadata",
                    "metric": metric,
                }
            )
            continue
        working = merged.dropna(subset=[column]).copy()
        working["__donor_unit__"] = (
            working[study_key].astype(str) + "::" + working[donor_key].astype(str)
        )
        rows = []
        for (method, level, donor), group in working.groupby(
            ["method", column, "__donor_unit__"], sort=True, observed=True
        ):
            metrics = _metrics_for_frame(group)
            rows.append(
                {
                    "axis": axis_name,
                    "level": str(level),
                    "method": str(method),
                    "donor_unit": str(donor),
                    **metrics,
                }
            )
        donor_level = pd.DataFrame(rows)
        if donor_level.empty:
            summary_rows.append({"axis": axis_name, "status": "not_estimable", "metric": metric})
            continue
        for (method, level), group in donor_level.groupby(
            ["method", "level"], sort=True, observed=True
        ):
            level_rows.append(
                {
                    "axis": axis_name,
                    "level": level,
                    "method": method,
                    "metric": metric,
                    "estimate": float(group[metric].mean()),
                    "sd_across_donors": (
                        float(group[metric].std(ddof=1)) if len(group) > 1 else np.nan
                    ),
                    "n_donors": int(group["donor_unit"].nunique()),
                    "status": (
                        "estimated"
                        if group["donor_unit"].nunique() >= min_donors_per_level
                        else "descriptive_only_low_donor_count"
                    ),
                }
            )
        level_frame = pd.DataFrame(level_rows)
        axis_frame = level_frame[level_frame["axis"] == axis_name]
        for method, group in axis_frame.groupby("method", sort=True):
            valid = group[group["n_donors"] >= min_donors_per_level]
            if len(valid) < 2:
                summary_rows.append(
                    {
                        "axis": axis_name,
                        "method": method,
                        "metric": metric,
                        "status": "not_estimable_lt2_supported_levels",
                        "n_supported_levels": len(valid),
                    }
                )
                continue
            worst = valid.loc[valid["estimate"].idxmin()]
            best = valid.loc[valid["estimate"].idxmax()]
            summary_rows.append(
                {
                    "axis": axis_name,
                    "method": method,
                    "metric": metric,
                    "status": "estimated_descriptive",
                    "n_supported_levels": len(valid),
                    "worst_level": worst["level"],
                    "worst_estimate": worst["estimate"],
                    "best_level": best["level"],
                    "best_estimate": best["estimate"],
                    "max_minus_min": float(best["estimate"] - worst["estimate"]),
                }
            )
    return pd.DataFrame(level_rows), pd.DataFrame(summary_rows)


def sample_enrichment_diagnostics(
    metadata: pd.DataFrame,
    *,
    cluster_key: str | None,
    sample_key: str | None,
    dominant_fraction_threshold: float = 0.8,
) -> pd.DataFrame:
    """Flag clusters dominated by one biological sample; never call this batch correction."""
    columns = [
        "cluster",
        "n_cells",
        "n_samples",
        "dominant_sample",
        "dominant_sample_fraction",
        "normalized_sample_entropy",
        "flag",
        "status",
    ]
    if (
        not cluster_key
        or not sample_key
        or cluster_key not in metadata
        or sample_key not in metadata
    ):
        return pd.DataFrame(
            [
                {
                    "cluster": "__all__",
                    "flag": "NOT_ASSESSED",
                    # Fail closed: missing metadata is not_assessed, never "clean".
                    "status": "not_assessed_missing_cluster_or_sample_metadata",
                }
            ],
            columns=columns,
        )
    rows = []
    frame = metadata.dropna(subset=[cluster_key, sample_key]).copy()
    for cluster, group in frame.groupby(cluster_key, sort=True):
        proportions = group[sample_key].astype(str).value_counts(normalize=True)
        n_samples = len(proportions)
        entropy = -float(np.sum(proportions * np.log(proportions)))
        normalized_entropy = entropy / math.log(n_samples) if n_samples > 1 else 0.0
        dominant_fraction = float(proportions.iloc[0])
        rows.append(
            {
                "cluster": str(cluster),
                "n_cells": len(group),
                "n_samples": n_samples,
                "dominant_sample": str(proportions.index[0]),
                "dominant_sample_fraction": dominant_fraction,
                "normalized_sample_entropy": normalized_entropy,
                "flag": "SAMPLE_ENRICHED"
                if dominant_fraction >= dominant_fraction_threshold
                else "PASS",
                "status": "estimated_descriptive",
            }
        )
    return pd.DataFrame(rows, columns=columns)


def _flag_series(frame: pd.DataFrame, spec: Mapping[str, object]) -> pd.Series | None:
    key = str(spec.get("key", ""))
    if not key or key not in frame:
        return None
    values = frame[key]
    if "threshold" in spec:
        numeric = pd.to_numeric(values, errors="coerce")
        threshold = float(spec["threshold"])
        direction = str(spec.get("direction", "ge"))
        flag = numeric.le(threshold) if direction == "le" else numeric.ge(threshold)
        return flag.mask(numeric.isna())
    positives = {str(value).strip().lower() for value in spec.get("positive_values", [True, 1])}
    flag = values.astype(str).str.strip().str.lower().isin(positives)
    return flag.mask(values.isna())


def qc_stratified_performance(
    merged: pd.DataFrame,
    *,
    study_key: str,
    donor_key: str,
    diagnostics: Mapping[str, Mapping[str, object]],
    metric: str = "macro_f1",
) -> pd.DataFrame:
    """Stratify performance by predeclared low-quality/doublet/ambient flags.

    Missing diagnostic columns or cell values yield ``not_assessed_*`` rows.
    Absence of flags is never rewritten as ``clean``.
    """
    rows = []
    # Keep low_quality as registry alias (mito/low-RNA style thresholds) plus
    # explicit doublet / ambient_rna axes for external tool contracts.
    for diagnostic in ("low_quality", "doublet", "ambient_rna", "low_rna", "high_mito"):
        spec = diagnostics.get(diagnostic, {})
        if not spec:
            rows.append(
                {
                    "diagnostic": diagnostic,
                    "status": "not_assessed_missing_predeclared_column",
                    "metric": metric,
                    "stratum": "not_assessed",
                    "identity_effect": "none",
                }
            )
            continue
        flag = _flag_series(merged, spec)
        if flag is None:
            rows.append(
                {
                    "diagnostic": diagnostic,
                    "status": "not_assessed_missing_predeclared_column",
                    "metric": metric,
                    "stratum": "not_assessed",
                    "identity_effect": "none",
                }
            )
            continue
        working = merged.copy()
        working["__stratum__"] = np.select(
            [flag.eq(True), flag.eq(False)],
            ["flagged", "not_flagged"],
            default="missing_diagnostic",
        )
        working["__donor_unit__"] = (
            working[study_key].astype(str) + "::" + working[donor_key].astype(str)
        )
        for (method, stratum, donor), group in working.groupby(
            ["method", "__stratum__", "__donor_unit__"], sort=True, observed=True
        ):
            missing_values = str(stratum) == "missing_diagnostic"
            rows.append(
                {
                    "diagnostic": diagnostic,
                    "method": str(method),
                    "stratum": str(stratum),
                    "donor_unit": str(donor),
                    "metric": metric,
                    "estimate": np.nan if missing_values else _metrics_for_frame(group)[metric],
                    "n_cells": len(group),
                    "status": (
                        "not_assessed_missing_cell_values"
                        if missing_values
                        else "estimated_donor_stratum"
                    ),
                    "identity_effect": "none",
                }
            )
    return pd.DataFrame(rows)
