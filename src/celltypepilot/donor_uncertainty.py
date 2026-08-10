"""Donor-aware uncertainty assessment.

Provides donor-stratified calibration diagnostics and leave-one-donor-out
stability analysis for cluster annotations. Donor-awareness is critical
because cells from the same donor are not independent observations.
"""

from __future__ import annotations

import logging
from collections.abc import Callable

import numpy as np
import pandas as pd

from .calibration import calibration_diagnostics

logger = logging.getLogger(__name__)


def donor_stratified_ece(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    confidence: np.ndarray,
    donor_ids: np.ndarray,
    n_bins: int = 10,
) -> dict:
    """Compute per-donor ECE and global aggregate statistics.

    Parameters
    ----------
    y_true : true labels
    y_pred : predicted labels
    confidence : prediction confidence scores
    donor_ids : donor identifiers for each cell
    n_bins : number of calibration bins

    Returns
    -------
    dict with: global_ece, per_donor_ece (DataFrame), ece_cv, outlier_donors
    """
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    confidence = np.asarray(confidence)
    donor_ids = np.asarray(donor_ids)

    # Global ECE
    try:
        global_diag, _, _, _ = calibration_diagnostics(y_true, y_pred, confidence, n_bins)
        global_ece = global_diag.get("ece", np.nan)
    except Exception:
        global_ece = np.nan

    # Per-donor ECE
    df = pd.DataFrame(
        {
            "y_true": y_true,
            "y_pred": y_pred,
            "confidence": confidence,
            "donor": donor_ids,
        }
    )

    donor_eces = []
    for d, group in df.groupby("donor"):
        if len(group) < 10:
            # Too few cells for meaningful ECE
            donor_eces.append({"donor": d, "ece": np.nan, "n_cells": len(group)})
            continue
        try:
            diag, _, _, _ = calibration_diagnostics(
                group["y_true"].values,
                group["y_pred"].values,
                group["confidence"].values,
                n_bins,
            )
            ece = diag.get("ece", np.nan)
        except Exception:
            ece = np.nan
        donor_eces.append({"donor": d, "ece": ece, "n_cells": len(group)})

    per_donor_df = pd.DataFrame(donor_eces)

    if per_donor_df.empty or per_donor_df["ece"].isna().all():
        return {
            "global_ece": global_ece,
            "per_donor_ece": per_donor_df,
            "ece_cv": np.nan,
            "outlier_donors": [],
        }

    eces = per_donor_df["ece"].dropna().values
    if len(eces) < 2:
        ece_cv = np.nan
        outliers = []
    else:
        mean_ece = np.mean(eces)
        ece_cv = float(np.std(eces, ddof=1) / mean_ece) if mean_ece > 0 else np.nan
        median_ece = np.median(eces)
        outliers = per_donor_df[per_donor_df["ece"] > 2 * median_ece]["donor"].tolist()

    return {
        "global_ece": float(global_ece),
        "per_donor_ece": per_donor_df,
        "ece_cv": float(ece_cv) if not np.isnan(ece_cv) else np.nan,
        "outlier_donors": outliers,
    }


def donor_stability_assessment(
    results_df: pd.DataFrame,
    cluster_key: str,
    donor_key: str,
    scorer_fn: Callable | None = None,
) -> pd.DataFrame:
    """Leave-one-donor-out stability assessment for cluster annotations.

    For each donor, removes that donor's cells and checks whether the
    top annotation for each cluster changes. This reveals whether any
    single donor is disproportionately driving a cluster's identity.

    Parameters
    ----------
    results_df : DataFrame with at minimum cluster_key, donor_key,
        'cell_type' (annotation), and 'combined_score' columns
    cluster_key : column identifying clusters
    donor_key : column identifying donors
    scorer_fn : optional re-scoring function (not used in simple mode)

    Returns
    -------
    DataFrame with per-cluster: annotation_stable, unstable_donors,
        max_score_delta, n_donors
    """
    if donor_key not in results_df.columns:
        logger.warning("Donor key '%s' not found in results DataFrame.", donor_key)
        return pd.DataFrame()

    if cluster_key not in results_df.columns:
        logger.warning("Cluster key '%s' not found in results DataFrame.", cluster_key)
        return pd.DataFrame()

    donors = results_df[donor_key].unique()
    if len(donors) < 3:
        logger.warning(
            "Fewer than 3 donors (%d) available; "
            "leave-one-donor-out stability requires at least 3.",
            len(donors),
        )
        return pd.DataFrame()

    # Determine annotation column
    annot_col = "cell_type"
    if annot_col not in results_df.columns:
        for candidate in ["decision", "annotation", "predicted_label"]:
            if candidate in results_df.columns:
                annot_col = candidate
                break

    score_col = "combined_score"
    if score_col not in results_df.columns:
        score_col = None

    clusters = results_df[cluster_key].unique()
    results = []

    for cl in clusters:
        cl_data = results_df[results_df[cluster_key] == cl]

        if annot_col not in cl_data.columns:
            results.append(
                {
                    "cluster": cl,
                    "annotation_stable": True,
                    "unstable_donors": [],
                    "max_score_delta": 0.0,
                    "n_donors": len(donors),
                }
            )
            continue

        # Baseline annotation: most frequent cell type across all donors
        baseline_annot = cl_data[annot_col].mode()
        baseline_annot = baseline_annot.iloc[0] if len(baseline_annot) > 0 else "Unknown"

        # Baseline score
        baseline_score = (
            float(cl_data[score_col].mean()) if score_col and score_col in cl_data.columns else 0.0
        )

        unstable_donors = []
        max_delta = 0.0

        for donor in donors:
            # Remove this donor's cells
            loo_data = cl_data[cl_data[donor_key] != donor]

            if len(loo_data) == 0:
                continue

            # Check if annotation changes
            loo_annot = loo_data[annot_col].mode()
            loo_annot = loo_annot.iloc[0] if len(loo_annot) > 0 else "Unknown"

            if loo_annot != baseline_annot:
                unstable_donors.append(donor)

            # Score delta
            if score_col and score_col in loo_data.columns:
                loo_score = float(loo_data[score_col].mean())
                delta = abs(loo_score - baseline_score)
                max_delta = max(max_delta, delta)

        results.append(
            {
                "cluster": cl,
                "annotation_stable": len(unstable_donors) == 0,
                "unstable_donors": unstable_donors,
                "max_score_delta": float(max_delta),
                "n_donors": int(len(donors)),
            }
        )

    return pd.DataFrame(results)


def donor_uncertainty_summary(
    stability_df: pd.DataFrame,
    ece_result: dict,
) -> dict:
    """Aggregate stability and ECE results into a summary.

    Returns
    -------
    dict with: pct_stable_clusters, mean_donor_ece, ece_cv,
        n_outlier_donors, overall_assessment (stable/caution/unstable)
    """
    pct_stable = 1.0
    if not stability_df.empty and "annotation_stable" in stability_df.columns:
        pct_stable = float(stability_df["annotation_stable"].mean())

    per_donor_ece = ece_result.get("per_donor_ece", pd.DataFrame())
    mean_ece = (
        float(per_donor_ece["ece"].mean())
        if not per_donor_ece.empty and "ece" in per_donor_ece.columns
        else np.nan
    )
    ece_cv = ece_result.get("ece_cv", np.nan)
    n_outliers = len(ece_result.get("outlier_donors", []))

    # Assessment logic
    if pct_stable > 0.9 and (np.isnan(ece_cv) or ece_cv < 0.5) and n_outliers == 0:
        assessment = "stable"
    elif pct_stable < 0.7 or n_outliers > 2:
        assessment = "unstable"
    else:
        assessment = "caution"

    return {
        "schema_version": "celltypepilot.donor-uncertainty.v1",
        "pct_stable_clusters": float(pct_stable),
        "mean_donor_ece": float(mean_ece) if not np.isnan(mean_ece) else None,
        "ece_cv": float(ece_cv) if not np.isnan(ece_cv) else None,
        "n_outlier_donors": n_outliers,
        "overall_assessment": assessment,
    }
