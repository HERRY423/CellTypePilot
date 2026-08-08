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

from .calibration import CalibrationError, calibration_diagnostics

COMPARATOR_METHODS = ("celltypepilot", "celltypist", "singler", "azimuth", "popv")
ABSTAIN_LABELS = {"unknown", "abstain", "unassigned", "na", "nan", ""}


class BenchmarkValidationError(ValueError):
    """Raised when a benchmark would violate its declared independence design."""


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


def evaluate_holdout_predictions(
    truth: pd.Series,
    assignments: pd.DataFrame,
    predictions: pd.DataFrame,
    expected_methods: tuple[str, ...] = COMPARATOR_METHODS,
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
        y_true = truth_map.reindex(method_predictions["cell_id"].astype(str)).astype(str)
        y_pred = method_predictions["predicted_label"].fillna("Unknown").astype(str)
        aggregate_rows.append(
            {
                "method": method,
                "status": "evaluated",
                "n_folds": int(method_predictions["fold_id"].nunique()),
                **_classification_metrics(
                    y_true.to_numpy(),
                    y_pred.to_numpy(),
                    method_predictions["confidence"].to_numpy()
                    if "confidence" in method_predictions
                    else None,
                ),
            }
        )
    return pd.DataFrame(aggregate_rows), per_fold


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
            diagnostics, _, _ = calibration_diagnostics(y_true, y_pred, confidence)
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
        _, bins, risk = calibration_diagnostics(
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
    }
    manifest_path = output / "benchmark_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return {"assignments": assignments_path, "manifest": manifest_path}
