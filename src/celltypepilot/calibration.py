"""Leakage-aware confidence diagnostics and abstention policy utilities."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ABSTAIN_LABELS = {"unknown", "abstain", "unassigned", "na", "nan", ""}


class CalibrationError(ValueError):
    """Raised when calibration inputs are invalid or would mix evaluation roles."""


def _is_abstained(values: np.ndarray) -> np.ndarray:
    return np.array([str(value).strip().lower() in ABSTAIN_LABELS for value in values])


def calibration_diagnostics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    confidence: np.ndarray,
    n_bins: int = 10,
) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """Compute top-label calibration bins and a risk-coverage curve.

    ``confidence`` is interpreted as confidence in the predicted top label, so
    the reported Brier score is binary correctness Brier, not multiclass Brier.
    """
    y_true = np.asarray(y_true).astype(str)
    y_pred = np.asarray(y_pred).astype(str)
    confidence = np.asarray(confidence, dtype=float)
    if not (len(y_true) == len(y_pred) == len(confidence)):
        raise CalibrationError("truth, predictions, and confidence must have equal length")
    if len(confidence) == 0 or np.any(~np.isfinite(confidence)):
        raise CalibrationError("confidence must be non-empty and finite")
    if np.any((confidence < 0) | (confidence > 1)):
        raise CalibrationError("confidence values must be in [0, 1]")

    abstained = _is_abstained(y_pred)
    eligible = ~abstained
    correct = (y_true == y_pred).astype(float)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bin_rows = []
    ece = 0.0
    n_eligible = int(np.sum(eligible))
    for index in range(n_bins):
        lower, upper = edges[index], edges[index + 1]
        in_bin = eligible & (confidence >= lower)
        in_bin &= confidence <= upper if index == n_bins - 1 else confidence < upper
        count = int(np.sum(in_bin))
        if count:
            mean_confidence = float(np.mean(confidence[in_bin]))
            empirical_accuracy = float(np.mean(correct[in_bin]))
            gap = abs(mean_confidence - empirical_accuracy)
            ece += (count / n_eligible) * gap if n_eligible else 0.0
        else:
            mean_confidence = np.nan
            empirical_accuracy = np.nan
            gap = np.nan
        bin_rows.append(
            {
                "bin": index,
                "lower": lower,
                "upper": upper,
                "n": count,
                "mean_confidence": mean_confidence,
                "empirical_accuracy": empirical_accuracy,
                "absolute_gap": gap,
            }
        )

    risk_rows = []
    if n_eligible:
        thresholds = np.unique(confidence[eligible])[::-1]
        for threshold in thresholds:
            retained = eligible & (confidence >= threshold)
            retained_count = int(np.sum(retained))
            coverage = retained_count / len(y_true)
            risk = 1.0 - float(np.mean(correct[retained]))
            risk_rows.append(
                {
                    "n_retained": retained_count,
                    "threshold": float(threshold),
                    "coverage": coverage,
                    "selective_risk": risk,
                    "selective_accuracy": 1.0 - risk,
                }
            )
    risk_curve = pd.DataFrame(risk_rows)
    aurc = (
        float(np.trapezoid(risk_curve["selective_risk"], risk_curve["coverage"]))
        if len(risk_curve) > 1
        else np.nan
    )
    diagnostics = {
        "n_predictions": int(len(y_true)),
        "n_scored": n_eligible,
        "top_label_brier": float(np.mean((confidence[eligible] - correct[eligible]) ** 2))
        if n_eligible
        else np.nan,
        "ece": float(ece) if n_eligible else np.nan,
        "aurc": aurc,
        "confidence_semantics": "predicted_top_label_correctness",
    }
    return diagnostics, pd.DataFrame(bin_rows), risk_curve


def fit_abstention_policy(
    truth: pd.Series,
    predictions: pd.DataFrame,
    method: str,
    max_selective_error: float = 0.1,
    min_coverage: float = 0.2,
    dataset_role: str = "calibration",
) -> dict:
    """Fit a threshold only on an explicitly designated calibration dataset."""
    if dataset_role != "calibration":
        raise CalibrationError("Abstention thresholds may only be fit on role='calibration'")
    if not 0 <= max_selective_error < 1 or not 0 < min_coverage <= 1:
        raise CalibrationError("Invalid max_selective_error or min_coverage")
    required = {"cell_id", "method", "predicted_label", "confidence"}
    missing = required - set(predictions.columns)
    if missing:
        raise CalibrationError(f"Calibration predictions missing columns: {sorted(missing)}")
    frame = predictions[predictions["method"] == method].copy()
    if frame.empty:
        raise CalibrationError(f"No calibration predictions for method {method!r}")
    truth_map = truth.copy()
    truth_map.index = truth_map.index.astype(str)
    y_true = truth_map.reindex(frame["cell_id"].astype(str))
    if y_true.isna().any():
        raise CalibrationError("Calibration truth is missing for predicted cells")
    confidence = frame["confidence"].astype(float).to_numpy()
    y_pred = frame["predicted_label"].fillna("Unknown").astype(str).to_numpy()
    diagnostics, bins, risk_curve = calibration_diagnostics(
        y_true.astype(str).to_numpy(), y_pred, confidence
    )
    candidates = risk_curve[
        (risk_curve["selective_risk"] <= max_selective_error)
        & (risk_curve["coverage"] >= min_coverage)
    ]
    if candidates.empty:
        raise CalibrationError(
            "No confidence threshold satisfies the requested risk and coverage constraints"
        )
    chosen = candidates.sort_values(["coverage", "threshold"], ascending=[False, False]).iloc[0]
    source_bytes = frame.sort_values("cell_id").to_csv(index=False).encode("utf-8")
    return {
        "schema_version": "celltypepilot.abstention-policy.v1",
        "method": method,
        "dataset_role": dataset_role,
        "threshold": float(chosen["threshold"]),
        "max_selective_error": max_selective_error,
        "min_coverage": min_coverage,
        "empirical_selective_error": float(chosen["selective_risk"]),
        "empirical_coverage": float(chosen["coverage"]),
        "calibration_predictions_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "diagnostics": diagnostics,
        "calibration_bins": bins.to_dict(orient="records"),
    }


def save_abstention_policy(policy: dict, path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(_json_safe(policy), indent=2, allow_nan=False),
        encoding="utf-8",
    )
    return output


def _json_safe(value):
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def apply_policy_to_annotations(results: pd.DataFrame, policy: dict) -> pd.DataFrame:
    """Only downgrade plugin annotations; a policy can never upgrade a critic abstention."""
    if policy.get("schema_version") != "celltypepilot.abstention-policy.v1":
        raise CalibrationError("Unsupported abstention policy schema")
    if policy.get("method") != "celltypepilot":
        raise CalibrationError("Annotation pipeline accepts only a celltypepilot policy")
    threshold = float(policy["threshold"])
    output = results.copy()
    if "combined_score" in output:
        scores = pd.to_numeric(output["combined_score"], errors="coerce").fillna(0.0)
    else:
        scores = pd.Series(0.0, index=output.index)
    forced = scores < threshold
    for index in output.index[forced]:
        if output.at[index, "decision"] != "abstain":
            output.at[index, "candidate_cell_type"] = output.at[index, "cell_type"]
        output.at[index, "decision"] = "abstain"
        output.at[index, "cell_type"] = "Unknown"
        if "cl_id" in output:
            output.at[index, "cl_id"] = ""
        prior_reason = str(output.at[index, "abstain_reason"] or "").strip("; ")
        output.at[index, "abstain_reason"] = "; ".join(
            value for value in (prior_reason, "CALIBRATED_LOW_CONFIDENCE") if value
        )
        flags = str(output.at[index, "critic_flags"] or "").strip()
        output.at[index, "critic_flags"] = (
            "CALIBRATED_LOW_CONFIDENCE"
            if flags == "PASS" or not flags
            else flags + "; CALIBRATED_LOW_CONFIDENCE"
        )
        output.at[index, "critic_confidence"] = "needs_review"
    output["calibration_threshold"] = threshold
    output["calibration_policy_applied"] = True
    return output
