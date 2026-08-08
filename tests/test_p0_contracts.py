import anndata as ad
import numpy as np
import pandas as pd

from celltypepilot.benchmark import build_holdout_assignments, evaluate_holdout_predictions
from celltypepilot.benchmark_runner import materialize_fold
from celltypepilot.calibration import CalibrationError, fit_abstention_policy
from celltypepilot.data_adapter import (
    get_all_markers_for_tissue,
    load_marker_atlas,
    summarize_atlas_evidence,
    validate_atlas_provenance,
)
from celltypepilot.reference_registry import (
    ReferenceContractError,
    select_registered_celltypist_model,
    validate_reference_adata,
)


def test_holdout_and_calibration_contracts_are_distinct():
    metadata = pd.DataFrame(
        {"study": ["s1", "s1", "s2", "s2"], "donor": ["d1", "d1", "d2", "d2"]},
        index=["c1", "c2", "c3", "c4"],
    )
    assignments = build_holdout_assignments(metadata, "study", "donor", "study")
    predictions = assignments[["cell_id", "fold_id"]].copy()
    predictions["method"] = "celltypepilot"
    predictions["predicted_label"] = ["T", "B", "B", "B"]
    predictions["confidence"] = [0.95, 0.2, 0.9, 0.8]
    truth = pd.Series(["T", "T", "B", "B"], index=metadata.index)
    aggregate, _ = evaluate_holdout_predictions(truth, assignments, predictions)
    assert aggregate.loc[aggregate["method"] == "celltypepilot", "ece"].notna().all()

    policy = fit_abstention_policy(
        truth,
        predictions,
        "celltypepilot",
        max_selective_error=0.0,
        min_coverage=0.5,
        dataset_role="calibration",
    )
    assert policy["threshold"] == 0.8
    try:
        fit_abstention_policy(truth, predictions, "celltypepilot", dataset_role="test")
    except CalibrationError:
        pass
    else:
        raise AssertionError("test predictions must not fit a released threshold")


def test_reference_contract_fails_closed_for_missing_or_wrong_scope():
    reference = ad.AnnData(
        X=np.ones((2, 2)),
        obs=pd.DataFrame({"cell_type": ["T", "B"]}, index=["r1", "r2"]),
        var=pd.DataFrame(index=["A", "B"]),
    )
    try:
        validate_reference_adata(reference, "human", "blood", "cell_type")
    except ReferenceContractError:
        pass
    else:
        raise AssertionError("missing reference contract must fail")

    reference.uns["celltypepilot_reference"] = {
        "species": "human",
        "tissues": ["blood"],
        "source": "test",
        "version": "v1",
        "label_ontology": "CL",
        "training_studies": ["s1"],
    }
    assert validate_reference_adata(reference, "human", "blood", "cell_type")["status"] == (
        "verified"
    )
    try:
        validate_reference_adata(reference, "human", "lung", "cell_type")
    except ReferenceContractError:
        pass
    else:
        raise AssertionError("wrong-tissue reference must fail")


def test_celltypist_registry_does_not_guess_outside_declared_scope():
    assert select_registered_celltypist_model("human", "blood")["model_name"] == (
        "Immune_All_Low.pkl"
    )
    try:
        select_registered_celltypist_model("human", "lung")
    except ReferenceContractError:
        pass
    else:
        raise AssertionError("immune model must not silently become a lung model")


def test_atlas_evidence_levels_do_not_upgrade_aggregate_citations():
    atlas = load_marker_atlas("human")
    assert validate_atlas_provenance(atlas) == []
    summary = summarize_atlas_evidence(atlas)
    assert summary["total_relationships"] == 599
    assert summary["edge_verified_fraction"] == 0.0
    exploratory = get_all_markers_for_tissue(atlas, "blood", "database")
    strict = get_all_markers_for_tissue(atlas, "blood", "edge_verified")
    assert exploratory["T cell"]["positive_markers"]
    assert strict["T cell"]["positive_markers"] == []


def test_benchmark_fold_never_exposes_test_truth(tmp_path):
    obs = pd.DataFrame(
        {
            "cluster": ["0", "0", "1", "1"],
            "truth": ["T", "T", "B", "B"],
            "study": ["s1", "s1", "s2", "s2"],
            "donor": ["d1", "d1", "d2", "d2"],
        },
        index=["c1", "c2", "c3", "c4"],
    )
    adata = ad.AnnData(X=np.ones((4, 2)), obs=obs, var=pd.DataFrame(index=["A", "B"]))
    assignments = build_holdout_assignments(obs, "study", "donor", "study")
    paths = materialize_fold(
        adata,
        assignments,
        "study=s1",
        "truth",
        "cluster",
        tmp_path,
        "human",
        "blood",
    )
    train = ad.read_h5ad(paths["train"])
    test = ad.read_h5ad(paths["test"])
    assert "truth" not in test.obs
    assert "cell_type" not in test.obs
    assert "cell_type" in train.obs
    assert train.uns["celltypepilot_reference"]["held_out_fold"] == "study=s1"
