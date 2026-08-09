import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from celltypepilot import marker_scorer
from celltypepilot.benchmark import (
    BenchmarkValidationError,
    build_holdout_assignments,
    evaluate_holdout_predictions,
)
from celltypepilot.benchmark_runner import materialize_fold
from celltypepilot.calibration import (
    CalibrationError,
    apply_policy_to_annotations,
    calibration_diagnostics,
    fit_abstention_policy,
)
from celltypepilot.critic import run_critic
from celltypepilot.data_adapter import (
    get_all_markers_for_tissue,
    load_marker_atlas,
    summarize_atlas_evidence,
    validate_atlas_provenance,
)
from celltypepilot.orchestrator import run_annotation_pipeline
from celltypepilot.reference_registry import (
    ReferenceContractError,
    validate_reference_adata,
)


def _small_adata() -> ad.AnnData:
    obs = pd.DataFrame(
        {"cluster": ["0"] * 4 + ["1"] * 4},
        index=[f"cell-{i}" for i in range(8)],
    )
    # A and B are expressed in cluster 0; C is present but silent; D is absent.
    X = np.array(
        [[2.0, 2.0, 0.0]] * 4 + [[0.0, 1.0, 0.0]] * 4,
        dtype=float,
    )
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=["A", "B", "C"]))


def test_de_support_requires_fdr_logfc_expression_and_expected_denominator(monkeypatch):
    adata = _small_adata()
    monkeypatch.setattr(marker_scorer.sc.tl, "rank_genes_groups", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        marker_scorer,
        "_extract_de_results",
        lambda *_: {
            "0": pd.DataFrame(
                {
                    "gene": ["A", "B", "C"],
                    "logfoldchange": [1.0, 1.0, -1.0],
                    "pval": [0.001, 0.2, 0.001],
                    "pval_adj": [0.01, 0.2, 0.01],
                }
            ),
            "1": pd.DataFrame(
                {
                    "gene": ["A"],
                    "logfoldchange": [-1.0],
                    "pval": [0.01],
                    "pval_adj": [0.01],
                }
            ),
        },
    )
    scores = marker_scorer.compute_marker_scores(
        adata,
        "cluster",
        {"candidate": {"positive_markers": ["A", "B", "C", "D"], "negative_markers": []}},
    )
    row = scores[scores["cluster"] == "0"].iloc[0]
    assert row["pos_supporting_markers"] == "A"
    assert row["pct_overlap"] == 0.25
    assert row["n_pos_missing"] == 1
    assert row["pos_missing_markers"] == "D"
    assert row["n_pos_silent"] == 1
    assert row["pos_silent_markers"] == "C"


def test_scanpy_de_extraction_preserves_adjusted_p_values():
    obs = pd.DataFrame(
        {"cluster": ["0"] * 20 + ["1"] * 20},
        index=[f"cell-{i}" for i in range(40)],
    )
    X = np.zeros((40, 2), dtype=float)
    X[:20, 0] = np.log1p(5)
    X[20:, 1] = np.log1p(5)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=["A", "B"]))
    scores = marker_scorer.compute_marker_scores(
        adata,
        "cluster",
        {
            "type_a": {"positive_markers": ["A"], "negative_markers": ["B"]},
            "type_b": {"positive_markers": ["B"], "negative_markers": ["A"]},
        },
    )
    top = scores[scores["rank"] == 1].set_index("cluster")
    assert top.loc["0", "cell_type"] == "type_a"
    assert top.loc["0", "n_pos_de"] == 1
    assert top.loc["1", "cell_type"] == "type_b"
    assert top.loc["1", "n_pos_de"] == 1


def test_partial_evidence_abstains_and_preserves_candidate():
    adata = _small_adata()
    atlas = {
        "tissues": {
            "general": {
                "cell_types": {
                    "candidate": {
                        "cl_id": "CL:0000001",
                        "positive_markers": ["A", "C", "D", "E", "F"],
                        "negative_markers": [],
                    }
                }
            }
        }
    }
    annotations = pd.DataFrame(
        {
            "cluster": ["0"],
            "cell_type": ["candidate"],
            "cl_id": ["CL:0000001"],
            "combined_score": [0.8],
            "confidence": ["high"],
            "pct_overlap": [0.2],
        }
    )
    result = run_critic(adata, "cluster", annotations, atlas, "general")
    row = result.iloc[0]
    assert "PARTIAL_EVIDENCE" in row["critic_flags"]
    assert "PARTIAL_DE_SUPPORT" in row["critic_flags"]
    assert row["decision"] == "abstain"
    assert row["cell_type"] == "Unknown"
    assert row["candidate_cell_type"] == "candidate"
    assert row["n_missing_markers"] == 3
    assert row["n_silent_markers"] == 1


def test_study_holdout_evaluation_is_abstention_aware_and_reports_missing_methods():
    metadata = pd.DataFrame(
        {"study": ["s1", "s1", "s2", "s2"], "donor": ["d1", "d1", "d2", "d2"]},
        index=["c1", "c2", "c3", "c4"],
    )
    assignments = build_holdout_assignments(metadata, "study", "donor", "study")
    predictions = assignments[["cell_id", "fold_id"]].copy()
    predictions["method"] = "celltypepilot"
    predictions["predicted_label"] = ["T", "Unknown", "B", "T"]
    truth = pd.Series(["T", "T", "B", "B"], index=metadata.index)
    aggregate, per_fold = evaluate_holdout_predictions(truth, assignments, predictions)
    ctp = aggregate[aggregate["method"] == "celltypepilot"].iloc[0]
    assert ctp["status"] == "evaluated"
    assert ctp["coverage"] == 0.75
    assert ctp["selective_accuracy"] == 2 / 3
    assert aggregate[aggregate["method"] == "celltypist"].iloc[0]["status"] == "not_provided"
    assert len(per_fold) == 2

    bad = predictions.copy()
    bad.loc[0, "fold_id"] = "study=wrong"
    try:
        evaluate_holdout_predictions(truth, assignments, bad)
    except BenchmarkValidationError:
        pass
    else:
        raise AssertionError("fold mismatch must fail closed")

    ambiguous = metadata.copy()
    ambiguous.loc["c3", "donor"] = "d1"
    try:
        build_holdout_assignments(ambiguous, "study", "donor", "study")
    except BenchmarkValidationError:
        pass
    else:
        raise AssertionError("cross-study donor reuse must fail closed")


def test_every_bundled_marker_relationship_has_structured_provenance():
    assert validate_atlas_provenance(load_marker_atlas("human")) == []
    premium_path = (
        Path(__file__).parents[1]
        / "src/celltypepilot/data/packs/premium/marker_atlas.json"
    )
    premium = json.loads(premium_path.read_text(encoding="utf-8"))
    assert validate_atlas_provenance(premium) == []

    evidence = summarize_atlas_evidence(load_marker_atlas("human"))
    assert evidence["total_relationships"] == 599
    # After literature sweep: 280 edges at literature_cooccurrence_supported
    assert 0.4 < evidence["edge_verified_fraction"] < 0.5
    assert evidence["primary_verified_fraction"] == 0.0


def test_strict_marker_policy_excludes_unverified_edges_without_upgrading_them():
    atlas = load_marker_atlas("human")
    exploratory = get_all_markers_for_tissue(atlas, "blood", evidence_policy="database")
    strict = get_all_markers_for_tissue(atlas, "blood", evidence_policy="edge_verified")
    assert exploratory["T cell"]["positive_markers"]
    assert strict["T cell"]["positive_markers"] == []


def test_reference_contract_is_required_and_scope_checked():
    ref = _small_adata()
    ref.obs["cell_type"] = "T cell"
    try:
        validate_reference_adata(ref, "human", "blood", "cell_type")
    except ReferenceContractError:
        pass
    else:
        raise AssertionError("missing custom reference contract must fail closed")

    ref.uns["celltypepilot_reference"] = {
        "species": "human",
        "tissues": ["blood"],
        "source": "unit-test",
        "version": "v1",
        "label_ontology": "CL",
        "training_studies": ["train-study"],
    }
    contract = validate_reference_adata(ref, "human", "blood", "cell_type")
    assert contract["status"] == "verified"
    try:
        validate_reference_adata(ref, "human", "lung", "cell_type")
    except ReferenceContractError:
        pass
    else:
        raise AssertionError("tissue-mismatched reference must fail closed")


def test_calibration_is_risk_coverage_aware_and_can_only_downgrade():
    truth = pd.Series(["T", "T", "B", "B"], index=["c1", "c2", "c3", "c4"])
    predictions = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "method": ["celltypepilot"] * 4,
            "predicted_label": ["T", "B", "B", "B"],
            "confidence": [0.95, 0.2, 0.9, 0.8],
        }
    )
    diagnostics, bins, risk, _ = calibration_diagnostics(
        truth.to_numpy(),
        predictions["predicted_label"].to_numpy(),
        predictions["confidence"].to_numpy(),
    )
    assert 0 <= diagnostics["ece"] <= 1
    assert not bins.empty
    assert not risk.empty
    policy = fit_abstention_policy(
        truth,
        predictions,
        "celltypepilot",
        max_selective_error=0.0,
        min_coverage=0.5,
        dataset_role="calibration",
    )
    assert policy["threshold"] >= 0.8
    try:
        fit_abstention_policy(truth, predictions, "celltypepilot", dataset_role="test")
    except CalibrationError:
        pass
    else:
        raise AssertionError("test data must never fit an abstention policy")

    annotations = pd.DataFrame(
        {
            "cell_type": ["T", "B"],
            "candidate_cell_type": ["T", "B"],
            "cl_id": ["CL:0000084", "CL:0000236"],
            "combined_score": [0.95, 0.2],
            "decision": ["accepted", "accepted"],
            "abstain_reason": ["", ""],
            "critic_flags": ["PASS", "PASS"],
            "critic_confidence": ["high", "high"],
        }
    )
    calibrated = apply_policy_to_annotations(annotations, policy)
    assert calibrated.loc[0, "cell_type"] == "T"
    assert calibrated.loc[1, "cell_type"] == "Unknown"
    assert calibrated.loc[1, "decision"] == "abstain"


def test_fold_materialization_strips_test_truth_and_declares_training_reference(tmp_path):
    adata = _small_adata()
    adata.obs["truth"] = ["T"] * 4 + ["B"] * 4
    adata.obs["study"] = ["s1"] * 4 + ["s2"] * 4
    adata.obs["donor"] = ["d1"] * 4 + ["d2"] * 4
    assignments = build_holdout_assignments(adata.obs, "study", "donor", "study")
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


def test_reference_ensemble_critic_writeback_report_and_manifest_are_one_pipeline(tmp_path):
    markers = [
        "CD3D",
        "CD3E",
        "CD2",
        "TRAC",
        "CD7",
        "CD19",
        "MS4A1",
        "CD79A",
        "CD79B",
        "PAX5",
        "CD14",
        "NCAM1",
    ]
    genes = markers + [f"G{i}" for i in range(108)]
    query_X = np.zeros((40, 120))
    reference_X = np.zeros((40, 120))
    for matrix in (query_X, reference_X):
        matrix[:20, :5] = 5
        matrix[20:, 5:10] = 5
        matrix[:20, 12:62] = 2
        matrix[20:, 62:112] = 2
    query = ad.AnnData(
        X=np.log1p(query_X),
        obs=pd.DataFrame(
            {"leiden": ["0"] * 20 + ["1"] * 20},
            index=[f"q{i}" for i in range(40)],
        ),
        var=pd.DataFrame(index=genes),
    )
    reference = ad.AnnData(
        X=reference_X,
        obs=pd.DataFrame(
            {"cell_type": ["T cell"] * 20 + ["B cell"] * 20},
            index=[f"r{i}" for i in range(40)],
        ),
        var=pd.DataFrame(index=genes),
    )
    reference.uns["celltypepilot_reference"] = {
        "species": "human",
        "tissues": ["blood"],
        "source": "synthetic-unit-test",
        "version": "v1",
        "label_ontology": "test-labels",
        "training_studies": ["synthetic-reference"],
    }
    query_path = tmp_path / "query.h5ad"
    reference_path = tmp_path / "reference.h5ad"
    query.write(query_path)
    reference.write(reference_path)

    result = run_annotation_pipeline(
        query_path,
        "leiden",
        tmp_path / "out",
        species="human",
        tissue="blood",
        reference_path=reference_path,
        reference_backend="correlation",
        no_figures=True,
    )
    assert set(result["critic_results"]["cell_type"]) == {"T cell", "B cell"}
    assert set(result["critic_results"]["decision"]) == {"accepted"}
    assert "reference_scores" in result["paths"]
    assert "ensemble_scores" in result["paths"]
    assert "data.annotated.h5ad" in result["manifest"]["outputs"]
    assert result["manifest"]["parameters"]["pipeline_stages"] == [
        "context",
        "marker",
        "reference",
        "ensemble",
        "critic",
        "state",
        "novelty_ood",
        "writeback",
        "report",
        "manifest",
    ]
    validation_scope = result["manifest"]["parameters"]["validation_scope"]
    assert validation_scope["run_role"] == "draft_annotation_for_human_review"
    assert validation_scope["batch_robustness_claim"] == "not_assessed"
    assert validation_scope["complex_sample_robustness_claim"] == "not_assessed"
    assert "benchmark" in validation_scope["required_for_robustness_claims"]
    annotated = ad.read_h5ad(result["paths"]["annotated"])
    assert {"ctp_cell_type", "ctp_candidate_cell_type", "ctp_decision"} <= set(
        annotated.obs.columns
    )
