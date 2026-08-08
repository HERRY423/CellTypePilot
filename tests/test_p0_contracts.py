import json

import anndata as ad
import numpy as np
import pandas as pd

from celltypepilot.benchmark import (
    build_holdout_assignments,
    evaluate_holdout_predictions,
    save_benchmark_plan,
)
from celltypepilot.benchmark_runner import materialize_fold
from celltypepilot.calibration import CalibrationError, fit_abstention_policy
from celltypepilot.data_adapter import (
    get_all_markers_for_tissue,
    load_marker_atlas,
    summarize_atlas_evidence,
    validate_atlas_provenance,
)
from celltypepilot.orchestrator import assess_annotation_validation_scope
from celltypepilot.reference_registry import (
    ReferenceContractError,
    select_registered_celltypist_model,
    validate_reference_adata,
)
from celltypepilot.uncertainty import (
    EVIDENCE_SCORE_SEMANTICS,
    attach_uncertainty_language,
    build_uncertainty_language_manifest,
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


def test_atlas_evidence_levels_reflect_literature_sweep():
    atlas = load_marker_atlas("human")
    assert validate_atlas_provenance(atlas) == []
    summary = summarize_atlas_evidence(atlas)
    assert summary["total_relationships"] == 599
    # After literature sweep: 280 edges at literature_cooccurrence_supported (rank 1)
    # edge_verified_fraction counts rank >= 1 (includes literature + database + primary)
    assert 0.4 < summary["edge_verified_fraction"] < 0.5
    # literature_cooccurrence edges are present
    assert summary["verification_counts"].get("literature_cooccurrence_supported", 0) == 280
    # edge_verified policy (rank >= 2) still excludes literature-only edges
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


def test_benchmark_manifest_is_the_robustness_evidence_entrypoint(tmp_path):
    metadata = pd.DataFrame(
        {"study": ["s1", "s1", "s2", "s2"], "donor": ["d1", "d1", "d2", "d2"]},
        index=["c1", "c2", "c3", "c4"],
    )
    assignments = build_holdout_assignments(metadata, "study", "donor", "study")
    paths = save_benchmark_plan(assignments, tmp_path, "study", "donor", "study")
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    validation_scope = manifest["validation_scope"]
    assert validation_scope["run_role"] == "locked_study_donor_holdout_benchmark"
    assert validation_scope["statistical_independence_claim"] == "study_or_donor_fold_isolated"
    assert "claim_boundary" in validation_scope


def test_annotation_validation_scope_never_claims_batch_robustness():
    adata = ad.AnnData(
        X=np.ones((4, 2)),
        obs=pd.DataFrame(
            {
                "study": ["s1", "s1", "s2", "s2"],
                "donor": ["d1", "d1", "d2", "d2"],
                "batch": ["b1", "b1", "b2", "b2"],
            },
            index=["c1", "c2", "c3", "c4"],
        ),
        var=pd.DataFrame(index=["A", "B"]),
    )
    scope = assess_annotation_validation_scope(adata)
    assert scope["run_role"] == "draft_annotation_for_human_review"
    assert scope["batch_robustness_claim"] == "not_assessed"
    assert scope["complex_sample_robustness_claim"] == "not_assessed"
    assert scope["metadata_candidates"]["study_keys"] == ["study"]
    assert scope["metadata_candidates"]["donor_keys"] == ["donor"]
    assert "benchmark" in scope["required_for_robustness_claims"]


def test_uncertainty_language_never_rebrands_scores_as_probabilities():
    annotations = pd.DataFrame(
        {
            "cluster": ["0", "1"],
            "cell_type": ["T cell", "Unknown"],
            "combined_score": [0.82, 0.21],
            "critic_confidence": ["high", "needs_review"],
            "decision": ["accepted", "abstain"],
        }
    )
    enriched = attach_uncertainty_language(annotations)
    assert enriched["evidence_score"].tolist() == [0.82, 0.21]
    assert set(enriched["evidence_score_semantics"]) == {EVIDENCE_SCORE_SEMANTICS}
    assert set(enriched["critic_confidence_semantics"]) == {
        "rule_based_review_category_not_probability"
    }
    assert enriched["calibrated_probability"].isna().all()
    assert set(enriched["calibrated_probability_semantics"]) == {
        "not_available_from_annotation_run"
    }
    assert set(enriched["ood_novelty_signal"]) == {"not_assessed_in_annotation_run"}
    assert enriched.loc[1, "unknown_label_semantics"] == "safety_abstention_not_biological_class"

    manifest_block = build_uncertainty_language_manifest()
    assert manifest_block["probability_columns"]["calibrated_probability"] == (
        "not_available_from_annotation_run"
    )
    assert "not calibrated probabilities" in manifest_block["product_claim_boundary"]
