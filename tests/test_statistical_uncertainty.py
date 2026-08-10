"""Unit and integration tests for Statistical Uncertainty Assessment (Phase B)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from celltypepilot.benchmark import compare_methods_significance
from celltypepilot.bootstrap import (
    BootstrapResult,
    bootstrap_cluster_stability,
    bootstrap_metric_ci,
)
from celltypepilot.calibration import (
    calibration_diagnostics,
    conformal_risk_bound,
)
from celltypepilot.calibration_transforms import (
    IsotonicCalibrator,
    PlattCalibrator,
    TemperatureCalibrator,
    auto_select_calibrator,
)
from celltypepilot.donor_uncertainty import (
    donor_stability_assessment,
    donor_stratified_ece,
)


def test_bootstrap_metric_ci():
    y_true = np.array(["A", "A", "A", "A", "B", "B", "B", "B"])
    y_pred = np.array(["A", "A", "A", "B", "B", "B", "B", "A"])

    def acc_fn(t, p):
        return float(np.mean(t == p))

    res = bootstrap_metric_ci(y_true, y_pred, acc_fn, n_boot=100, seed=42)
    assert isinstance(res, BootstrapResult)
    assert res.point_estimate == 0.75
    assert 0.0 <= res.ci_lower <= res.ci_upper <= 1.0
    assert res.se >= 0.0


def test_bootstrap_cluster_stability():
    import anndata as ad

    rng = np.random.default_rng(42)
    X = rng.normal(size=(100, 10))
    # Two clear clusters
    X[:50] += 5.0
    obs = pd.DataFrame({"cluster": ["0"] * 50 + ["1"] * 50})
    var = pd.DataFrame(index=[f"g{i}" for i in range(10)])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    df = bootstrap_cluster_stability(adata, "cluster", n_boot=20, subsample_frac=0.8)
    assert isinstance(df, pd.DataFrame)
    assert "cluster" in df.columns
    assert "stability_score" in df.columns


def test_donor_stratified_ece():
    y_true = np.array(["A", "A", "B", "B", "A", "A", "B", "B"])
    y_pred = np.array(["A", "A", "B", "A", "A", "B", "B", "B"])
    conf = np.array([0.9, 0.8, 0.85, 0.9, 0.7, 0.6, 0.95, 0.8])
    donors = np.array(["d1", "d1", "d1", "d1", "d2", "d2", "d2", "d2"])

    res = donor_stratified_ece(y_true, y_pred, conf, donors)
    assert "global_ece" in res
    assert "per_donor_ece" in res
    assert isinstance(res["per_donor_ece"], pd.DataFrame)


def test_donor_stability_assessment():
    df = pd.DataFrame(
        {
            "cluster": ["0", "0", "0", "1", "1", "1"],
            "donor": ["d1", "d2", "d3", "d1", "d2", "d3"],
            "cell_type": ["T cell", "T cell", "T cell", "B cell", "B cell", "B cell"],
            "combined_score": [0.8, 0.85, 0.79, 0.9, 0.92, 0.88],
        }
    )

    stab = donor_stability_assessment(df, "cluster", "donor")
    assert isinstance(stab, pd.DataFrame)
    assert "annotation_stable" in stab.columns
    assert stab["annotation_stable"].all()


def test_calibration_transforms_fit_transform():
    scores = np.linspace(0.1, 0.9, 50)
    # Binary correctness correlates with score
    y_correct = (scores + np.random.uniform(-0.1, 0.1, 50) > 0.5).astype(int)

    iso = IsotonicCalibrator().fit(scores, y_correct)
    calib_iso = iso.transform(scores)
    assert len(calib_iso) == 50
    assert (calib_iso >= 0).all() and (calib_iso <= 1).all()

    platt = PlattCalibrator().fit(scores, y_correct)
    calib_platt = platt.transform(scores)
    assert len(calib_platt) == 50

    temp = TemperatureCalibrator().fit(scores, y_correct)
    calib_temp = temp.transform(scores)
    assert len(calib_temp) == 50


def test_auto_select_calibrator():
    scores = np.linspace(0.1, 0.9, 50)
    y_correct = (scores > 0.5).astype(int)

    best_model, report = auto_select_calibrator(scores, y_correct, n_cv=3)
    assert best_model is not None
    assert "selected_method" in report


def test_conformal_risk_bound():
    y_true = np.array(["A"] * 20 + ["B"] * 20)
    y_pred = np.array(["A"] * 20 + ["B"] * 20)
    conf = np.linspace(0.5, 0.99, 40)

    _, _, risk_curve, _ = calibration_diagnostics(y_true, y_pred, conf)
    bound = conformal_risk_bound(risk_curve, alpha=0.2, delta=0.05)
    assert isinstance(bound, dict)


def test_compare_methods_significance():
    per_fold_df = pd.DataFrame(
        {
            "fold": [0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6],
            "method": ["celltypepilot"] * 7 + ["celltypist"] * 7,
            "macro_f1": [
                0.85,
                0.88,
                0.86,
                0.89,
                0.87,
                0.88,
                0.86,
                0.70,
                0.72,
                0.71,
                0.73,
                0.69,
                0.71,
                0.70,
            ],
        }
    )

    res = compare_methods_significance(per_fold_df, "celltypepilot", "celltypist", "macro_f1")
    assert "p_value" in res
    assert "significant_at_005" in res
    assert res["significant_at_005"] is True
