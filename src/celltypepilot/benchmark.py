"""Leakage-resistant study/donor holdout benchmark utilities.

The module evaluates predictions; it never silently trains a comparator on the
test study or fabricates predictions for an unavailable external method.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .bootstrap import grouped_bootstrap_metric_ci
from .calibration import CalibrationError, calibration_diagnostics

COMPARATOR_METHODS = ("celltypepilot", "celltypist", "singler", "azimuth", "popv")
ABSTAIN_LABELS = {"unknown", "abstain", "unassigned", "na", "nan", ""}


class BenchmarkValidationError(ValueError):
    """Raised when a benchmark would violate its declared independence design."""


def apply_truth_label_map(truth: pd.Series, label_map: pd.DataFrame | None) -> pd.Series:
    """Apply an exhaustive, predeclared truth map to the evaluation endpoint."""
    if label_map is None:
        return truth.astype(str)
    required = {"method", "raw_label", "canonical_label"}
    missing = required - set(label_map)
    if missing:
        raise BenchmarkValidationError(f"Label map missing columns: {sorted(missing)}")
    truth_rows = label_map[label_map["method"].astype(str) == "__truth__"]
    if truth_rows.empty:
        raise BenchmarkValidationError("Label map must include __truth__ rows")
    if truth_rows["raw_label"].astype(str).duplicated().any():
        raise BenchmarkValidationError("Truth label map contains duplicate raw labels")
    mapping = dict(
        zip(
            truth_rows["raw_label"].astype(str),
            truth_rows["canonical_label"].astype(str),
            strict=True,
        )
    )
    observed = set(truth.astype(str))
    missing_truth = observed - set(mapping)
    if missing_truth:
        raise BenchmarkValidationError(
            f"Truth label map is not exhaustive (examples: {sorted(missing_truth)[:5]})"
        )
    mapped = truth.astype(str).map(mapping)
    mapped.index = truth.index
    return mapped


def build_holdout_assignments(
    metadata: pd.DataFrame,
    study_key: str,
    donor_key: str,
    strategy: str = "study",
) -> pd.DataFrame:
    """Create exhaustive leave-one-study or leave-one-donor test assignments."""
    missing = [key for key in (study_key, donor_key) if key not in metadata.columns]
    if missing:
        raise BenchmarkValidationError(f"Missing metadata columns: {missing}")
    if metadata.index.has_duplicates:
        raise BenchmarkValidationError("Cell identifiers must be unique")
    if metadata[[study_key, donor_key]].isna().any().any():
        raise BenchmarkValidationError("Study and donor metadata must not contain missing values")
    if strategy not in {"study", "donor"}:
        raise BenchmarkValidationError("strategy must be 'study' or 'donor'")

    frame = metadata[[study_key, donor_key]].astype(str).copy()
    donor_study_counts = frame.groupby(donor_key)[study_key].nunique()
    ambiguous_donors = donor_study_counts[donor_study_counts > 1].index.tolist()
    if ambiguous_donors:
        raise BenchmarkValidationError(
            "Donor identifiers occur in multiple studies; supply a globally unique donor key "
            f"before benchmarking (examples: {ambiguous_donors[:3]})"
        )
    frame["cell_id"] = frame.index.astype(str)
    frame["donor_unit"] = frame[study_key] + "::" + frame[donor_key]
    held_out = frame[study_key] if strategy == "study" else frame["donor_unit"]
    frame["fold_id"] = strategy + "=" + held_out
    frame["held_out_study"] = frame[study_key]
    frame["held_out_donor"] = frame["donor_unit"]
    frame["role"] = "test"
    return frame[
        [
            "cell_id",
            "fold_id",
            "role",
            study_key,
            donor_key,
            "held_out_study",
            "held_out_donor",
        ]
    ].reset_index(drop=True)


def validate_out_of_fold_predictions(
    assignments: pd.DataFrame,
    predictions: pd.DataFrame,
) -> None:
    """Require every scored prediction to match its predeclared test fold."""
    required = {"cell_id", "fold_id", "method", "predicted_label"}
    missing = required - set(predictions.columns)
    if missing:
        raise BenchmarkValidationError(f"Prediction table missing columns: {sorted(missing)}")
    if predictions.duplicated(["cell_id", "method"]).any():
        raise BenchmarkValidationError("Each method may provide only one prediction per cell")
    if "confidence" in predictions.columns:
        confidence = pd.to_numeric(predictions["confidence"], errors="coerce")
        if confidence.isna().any() or ((confidence < 0) | (confidence > 1)).any():
            raise BenchmarkValidationError("confidence must be numeric and within [0, 1]")

    expected = assignments.set_index("cell_id")["fold_id"].astype(str)
    unknown_cells = set(predictions["cell_id"].astype(str)) - set(expected.index.astype(str))
    if unknown_cells:
        raise BenchmarkValidationError(
            f"Predictions contain {len(unknown_cells)} cells absent from the holdout plan"
        )
    for row in predictions.itertuples(index=False):
        cell_id = str(row.cell_id)
        if str(row.fold_id) != str(expected.loc[cell_id]):
            raise BenchmarkValidationError(
                f"Prediction for cell {cell_id} is assigned to {row.fold_id}, "
                f"expected test fold {expected.loc[cell_id]}"
            )


def build_cluster_level_track(
    truth: pd.Series,
    assignments: pd.DataFrame,
    predictions: pd.DataFrame,
    cluster_labels: pd.Series,
    *,
    min_truth_purity: float = 0.0,
    min_prediction_support: float = 0.5,
) -> tuple[pd.Series, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Aggregate every method to the same predeclared cluster evaluation unit.

    CellTypePilot is a cluster annotator. Comparing its repeated cluster call
    with cell-level calls from reference tools creates an asymmetric endpoint.
    This function creates a separate cluster track: majority truth and every
    method prediction are aggregated by the same locked rule. Cluster purity
    and prediction support are retained as diagnostics rather than hidden.
    """
    if not 0 <= min_truth_purity <= 1:
        raise BenchmarkValidationError("min_truth_purity must be within [0, 1]")
    if not 0 < min_prediction_support <= 1:
        raise BenchmarkValidationError("min_prediction_support must be within (0, 1]")
    validate_out_of_fold_predictions(assignments, predictions)

    assignment = assignments.copy()
    assignment["cell_id"] = assignment["cell_id"].astype(str)
    assignment = assignment.set_index("cell_id", drop=False)
    clusters = cluster_labels.copy()
    clusters.index = clusters.index.astype(str)
    missing = assignment.index.difference(clusters.index)
    extra = clusters.index.difference(assignment.index)
    if len(missing) or len(extra):
        raise BenchmarkValidationError(
            f"cluster labels must be exhaustive: missing={len(missing)}, extra={len(extra)}"
        )
    assignment["cluster"] = clusters.reindex(assignment.index).astype(str)
    assignment["cluster_unit"] = (
        assignment["fold_id"].astype(str) + "::cluster=" + assignment["cluster"]
    )

    truth_map = truth.copy()
    truth_map.index = truth_map.index.astype(str)
    if truth_map.reindex(assignment.index).isna().any():
        raise BenchmarkValidationError("Ground truth is missing for cluster-track cells")
    assignment["truth"] = truth_map.reindex(assignment.index).astype(str)

    cluster_rows = []
    cluster_truth: dict[str, str] = {}
    kept_units: set[str] = set()
    for unit, group in assignment.groupby("cluster_unit", sort=True):
        counts = group["truth"].value_counts()
        top_count = int(counts.iloc[0])
        tied = sorted(counts[counts == top_count].index.astype(str))
        majority = tied[0] if len(tied) == 1 else "__ambiguous_truth__"
        purity = top_count / len(group)
        status = "eligible"
        if len(tied) > 1:
            status = "ambiguous_truth_tie"
        elif purity < min_truth_purity:
            status = "below_predeclared_truth_purity"
        if status == "eligible":
            kept_units.add(str(unit))
            cluster_truth[str(unit)] = majority
        cluster_rows.append(
            {
                "cluster_unit": str(unit),
                "fold_id": str(group["fold_id"].iloc[0]),
                "cluster": str(group["cluster"].iloc[0]),
                "n_cells": int(len(group)),
                "majority_truth": majority,
                "truth_purity": float(purity),
                "n_truth_labels": int(len(counts)),
                "truth_distribution": ";".join(
                    f"{label}:{count}" for label, count in counts.items()
                ),
                "status": status,
            }
        )

    diagnostics = pd.DataFrame(cluster_rows)
    cluster_assignments = (
        assignment[assignment["cluster_unit"].isin(kept_units)]
        .groupby("cluster_unit", sort=True)
        .first()
        .reset_index()
    )
    cluster_assignments["cell_id"] = cluster_assignments["cluster_unit"].astype(str)
    cluster_assignments = cluster_assignments[
        [
            "cell_id",
            "fold_id",
            "role",
            "held_out_study",
            "held_out_donor",
        ]
    ]
    cluster_truth_series = pd.Series(cluster_truth, name="truth", dtype=str)

    prediction = predictions.copy()
    prediction["cell_id"] = prediction["cell_id"].astype(str)
    prediction = prediction.merge(
        assignment[["cell_id", "cluster_unit"]].reset_index(drop=True),
        on="cell_id",
        how="left",
        validate="many_to_one",
    )
    prediction = prediction[prediction["cluster_unit"].isin(kept_units)]
    prediction_rows = []
    support_rows = []
    for (method, fold_id, unit), group in prediction.groupby(
        ["method", "fold_id", "cluster_unit"], sort=True
    ):
        labels = group["predicted_label"].fillna("Unknown").astype(str)
        covered = ~labels.str.strip().str.casefold().isin(ABSTAIN_LABELS)
        covered_labels = labels[covered]
        support = float(len(covered_labels) / len(group)) if len(group) else 0.0
        decision = "Unknown"
        aggregation_status = "insufficient_prediction_support"
        if support >= min_prediction_support and not covered_labels.empty:
            counts = covered_labels.value_counts()
            top_count = int(counts.iloc[0])
            winners = sorted(counts[counts == top_count].index.astype(str))
            if len(winners) == 1:
                decision = winners[0]
                aggregation_status = "majority"
            else:
                aggregation_status = "prediction_tie_abstain"
        confidence = np.nan
        if "confidence" in group and decision != "Unknown":
            winning = group.loc[labels == decision, "confidence"]
            confidence = float(pd.to_numeric(winning, errors="coerce").mean())
        row = {
            "cell_id": str(unit),
            "fold_id": str(fold_id),
            "method": str(method),
            "predicted_label": decision,
        }
        if "confidence" in group:
            row["confidence"] = 0.0 if np.isnan(confidence) else confidence
        prediction_rows.append(row)
        support_rows.append(
            {
                "cluster_unit": str(unit),
                "method": str(method),
                "prediction_support_fraction": support,
                "aggregation_status": aggregation_status,
            }
        )
    cluster_predictions = pd.DataFrame(prediction_rows)
    if support_rows:
        diagnostics = diagnostics.merge(pd.DataFrame(support_rows), on="cluster_unit", how="left")
    return cluster_truth_series, cluster_assignments, cluster_predictions, diagnostics


def _calc_acc(y_true, y_pred):
    return float(np.mean(y_true == y_pred))


def _calc_cov(y_true, y_pred):
    abstained = np.array([str(value).strip().lower() in ABSTAIN_LABELS for value in y_pred])
    return float(np.mean(~abstained))


def evaluate_holdout_predictions(
    truth: pd.Series,
    assignments: pd.DataFrame,
    predictions: pd.DataFrame,
    expected_methods: tuple[str, ...] = COMPARATOR_METHODS,
    bootstrap_ci: bool = True,
    n_boot: int = 1000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute abstention-aware aggregate and per-fold classification metrics."""
    validate_out_of_fold_predictions(assignments, predictions)
    truth_map = truth.copy()
    truth_map.index = truth_map.index.astype(str)
    fold_rows = []

    for (method, fold_id), group in predictions.groupby(["method", "fold_id"], sort=True):
        y_true = truth_map.reindex(group["cell_id"].astype(str)).astype(str)
        if y_true.isna().any():
            raise BenchmarkValidationError("Ground truth is missing for predicted cells")
        y_pred = group["predicted_label"].fillna("Unknown").astype(str)
        confidence = group["confidence"].to_numpy() if "confidence" in group else None
        metrics = _classification_metrics(y_true.to_numpy(), y_pred.to_numpy(), confidence)
        fold_rows.append({"method": method, "fold_id": fold_id, **metrics})

    per_fold = pd.DataFrame(fold_rows)
    aggregate_rows = []

    # Calculate mean and SE across folds
    fold_means = per_fold.groupby("method").mean(numeric_only=True)
    fold_se = per_fold.groupby("method").sem(numeric_only=True)

    for method in expected_methods:
        method_predictions = predictions[predictions["method"] == method]
        if method_predictions.empty:
            aggregate_rows.append(
                {
                    "method": method,
                    "status": "not_provided",
                    "n_cells": 0,
                    "n_folds": 0,
                }
            )
            continue

        y_true = truth_map.reindex(method_predictions["cell_id"].astype(str)).astype(str).to_numpy()
        y_pred = method_predictions["predicted_label"].fillna("Unknown").astype(str).to_numpy()
        conf = (
            method_predictions["confidence"].to_numpy()
            if "confidence" in method_predictions
            else None
        )

        base_metrics = _classification_metrics(y_true, y_pred, conf)

        agg = {
            "method": method,
            "status": "evaluated",
            "n_folds": int(method_predictions["fold_id"].nunique()),
            **base_metrics,
        }

        if bootstrap_ci:
            assignment_metadata = assignments.set_index(assignments["cell_id"].astype(str))
            method_metadata = assignment_metadata.reindex(method_predictions["cell_id"].astype(str))
            unit_rows = []
            for donor_unit, indices in method_metadata.groupby("held_out_donor").groups.items():
                positions = method_metadata.index.get_indexer(indices)
                unit_truth = y_true[positions]
                unit_prediction = y_pred[positions]
                unit_rows.append(
                    {
                        "donor_unit": str(donor_unit),
                        "study": str(method_metadata.loc[indices[0], "held_out_study"]),
                        "accuracy": _calc_acc(unit_truth, unit_prediction),
                        "coverage": _calc_cov(unit_truth, unit_prediction),
                    }
                )
            unit_frame = pd.DataFrame(unit_rows)
            strata = unit_frame["study"].to_numpy() if unit_frame["study"].nunique() > 1 else None
            acc_res = grouped_bootstrap_metric_ci(
                unit_frame["accuracy"].to_numpy(),
                unit_frame["donor_unit"].to_numpy(),
                strata=strata,
                n_boot=n_boot,
            )
            agg["accuracy_ci_lower"] = acc_res.ci_lower
            agg["accuracy_ci_upper"] = acc_res.ci_upper
            cov_res = grouped_bootstrap_metric_ci(
                unit_frame["coverage"].to_numpy(),
                unit_frame["donor_unit"].to_numpy(),
                strata=strata,
                n_boot=n_boot,
            )
            agg["coverage_ci_lower"] = cov_res.ci_lower
            agg["coverage_ci_upper"] = cov_res.ci_upper
            agg["ci_independent_unit"] = "donor"

        # Add cross-fold stats
        if method in fold_means.index:
            for col in fold_means.columns:
                if col not in ["n_cells"]:
                    agg[f"mean_cv_{col}"] = fold_means.loc[method, col]
                    agg[f"se_cv_{col}"] = fold_se.loc[method, col]

        aggregate_rows.append(agg)

    return pd.DataFrame(aggregate_rows), per_fold


def compare_methods_significance(
    per_fold_df: pd.DataFrame, method_a: str, method_b: str, metric: str
) -> dict:
    """Compare two methods using paired Wilcoxon signed-rank test on fold-level metrics."""
    fold_col = None
    for candidate in ["fold_id", "fold", "test_fold"]:
        if candidate in per_fold_df.columns:
            fold_col = candidate
            break
    if fold_col is None:
        return {"p_value": np.nan, "effect_size": np.nan, "significant_at_005": False}

    a_df = per_fold_df[per_fold_df["method"] == method_a].set_index(fold_col)
    b_df = per_fold_df[per_fold_df["method"] == method_b].set_index(fold_col)

    common_folds = a_df.index.intersection(b_df.index)
    if len(common_folds) < 3:
        return {"p_value": np.nan, "effect_size": np.nan, "significant_at_005": False}

    a_vals = a_df.loc[common_folds, metric].values
    b_vals = b_df.loc[common_folds, metric].values

    diffs = a_vals - b_vals
    if np.all(diffs == 0):
        return {"p_value": 1.0, "effect_size": 0.0, "significant_at_005": False}

    try:
        stat, p_val = stats.wilcoxon(diffs)
    except Exception:
        p_val = np.nan

    mean_diff = np.mean(diffs)
    effect_size = mean_diff / np.std(diffs) if np.std(diffs) > 0 else 0.0

    return {
        "p_value": float(p_val),
        "effect_size": float(effect_size),
        "significant_at_005": bool(p_val < 0.05),
        "method_a": method_a,
        "method_b": method_b,
        "metric": metric,
        "mean_difference": float(mean_diff),
        "n_folds": len(common_folds),
    }


def _classification_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    confidence: np.ndarray | None = None,
) -> dict:
    abstained = np.array([str(value).strip().lower() in ABSTAIN_LABELS for value in y_pred])
    covered = ~abstained
    correct = y_true == y_pred
    labels = sorted(set(y_true))
    recalls = []
    f1_scores = []
    for label in labels:
        true_label = y_true == label
        pred_label = y_pred == label
        tp = int(np.sum(true_label & pred_label))
        fn = int(np.sum(true_label & ~pred_label))
        fp = int(np.sum(~true_label & pred_label))
        recall = tp / (tp + fn) if tp + fn else 0.0
        precision = tp / (tp + fp) if tp + fp else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        recalls.append(recall)
        f1_scores.append(f1)

    metrics = {
        "n_cells": int(len(y_true)),
        "accuracy": float(np.mean(correct)) if len(y_true) else 0.0,
        "macro_f1": float(np.mean(f1_scores)) if f1_scores else 0.0,
        "balanced_accuracy": float(np.mean(recalls)) if recalls else 0.0,
        "coverage": float(np.mean(covered)) if len(y_true) else 0.0,
        "abstain_rate": float(np.mean(abstained)) if len(y_true) else 0.0,
        "selective_accuracy": float(np.mean(correct[covered])) if np.any(covered) else np.nan,
    }
    if confidence is not None:
        try:
            diagnostics, _, _, _ = calibration_diagnostics(y_true, y_pred, confidence)
        except CalibrationError as exc:
            raise BenchmarkValidationError(str(exc)) from exc
        metrics.update(
            {
                "top_label_brier": diagnostics["top_label_brier"],
                "ece": diagnostics["ece"],
                "aurc": diagnostics["aurc"],
            }
        )
    return metrics


def build_calibration_artifacts(
    truth: pd.Series,
    predictions: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build per-method calibration bins and risk-coverage curves."""
    if "confidence" not in predictions.columns:
        raise BenchmarkValidationError("Predictions need a confidence column for calibration")
    truth_map = truth.copy()
    truth_map.index = truth_map.index.astype(str)
    bin_frames = []
    risk_frames = []
    for method, frame in predictions.groupby("method", sort=True):
        y_true = truth_map.reindex(frame["cell_id"].astype(str))
        if y_true.isna().any():
            raise BenchmarkValidationError("Ground truth is missing for predicted cells")
        _, bins, risk, _ = calibration_diagnostics(
            y_true.astype(str).to_numpy(),
            frame["predicted_label"].fillna("Unknown").astype(str).to_numpy(),
            frame["confidence"].astype(float).to_numpy(),
        )
        bins.insert(0, "method", method)
        risk.insert(0, "method", method)
        bin_frames.append(bins)
        risk_frames.append(risk)
    return (
        pd.concat(bin_frames, ignore_index=True) if bin_frames else pd.DataFrame(),
        pd.concat(risk_frames, ignore_index=True) if risk_frames else pd.DataFrame(),
    )


def save_benchmark_plan(
    assignments: pd.DataFrame,
    output_dir: str | Path,
    study_key: str,
    donor_key: str,
    strategy: str,
) -> dict[str, Path]:
    """Persist the locked split plan and its machine-readable design metadata."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    assignments_path = output / "holdout_assignments.csv"
    assignments.to_csv(assignments_path, index=False)
    assignments_sha256 = hashlib.sha256(assignments_path.read_bytes()).hexdigest()
    manifest = {
        "schema_version": "celltypepilot.benchmark.v1",
        "strategy": strategy,
        "study_key": study_key,
        "donor_key": donor_key,
        "n_cells": len(assignments),
        "n_folds": int(assignments["fold_id"].nunique()),
        "assignments_sha256": assignments_sha256,
        "comparators": list(COMPARATOR_METHODS),
        "prediction_policy": "out_of_fold_only",
        "missing_method_policy": "report_not_provided_do_not_impute",
        "validation_scope": {
            "schema_version": "celltypepilot.validation-scope.v1",
            "run_role": "locked_study_donor_holdout_benchmark",
            "batch_robustness_claim": "evaluated_only_after_predictions_are_scored",
            "complex_sample_robustness_claim": "limited_to_declared_benchmark_distribution",
            "statistical_independence_claim": "study_or_donor_fold_isolated",
            "claim_boundary": (
                "Benchmark results support only the predeclared species, tissue, labels, "
                "studies, donors, comparators, and calibration policy recorded in this manifest."
            ),
            "agent_plugin_guidance": (
                "Plugin hosts may use benchmark_results.csv and benchmark_results_by_fold.csv "
                "as robustness evidence only after all requested methods report completed or "
                "not_provided status without fold-validation errors."
            ),
        },
    }
    manifest_path = output / "benchmark_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return {"assignments": assignments_path, "manifest": manifest_path}
