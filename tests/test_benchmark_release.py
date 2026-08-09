from __future__ import annotations

import hashlib
import json

import anndata as ad
import numpy as np
import pandas as pd

from celltypepilot.benchmark import BenchmarkValidationError, build_holdout_assignments
from celltypepilot.benchmark_release import (
    build_public_benchmark_release,
    validate_benchmark_registry,
)
from celltypepilot.benchmark_runner import apply_locked_label_map
from celltypepilot.bootstrap import grouped_bootstrap_metric_ci
from celltypepilot.robustness import (
    evaluate_by_independent_unit,
    merge_prediction_metadata,
    paired_method_comparisons,
    qc_stratified_performance,
    sample_enrichment_diagnostics,
    summarize_independent_units,
)


def test_grouped_bootstrap_uses_biological_units_not_cells():
    result = grouped_bootstrap_metric_ci(
        np.array([1.0, 0.0]),
        np.array(["large_donor", "small_donor"]),
        n_boot=200,
        seed=7,
    )
    assert result.point_estimate == 0.5
    assert result.ci_lower <= 0.5 <= result.ci_upper


def test_donor_summary_is_not_cell_count_weighted():
    truth = pd.Series(["A"] * 101, index=[f"c{i}" for i in range(101)])
    metadata = pd.DataFrame(
        {
            "study": ["s"] * 101,
            "donor": ["large"] * 100 + ["small"],
        },
        index=truth.index,
    )
    predictions = pd.DataFrame(
        {
            "cell_id": truth.index,
            "method": "celltypepilot",
            "predicted_label": ["A"] * 100 + ["B"],
        }
    )
    merged = merge_prediction_metadata(truth, metadata, predictions)
    donor = evaluate_by_independent_unit(merged, study_key="study", donor_key="donor")
    summary = summarize_independent_units(donor, metrics=("accuracy",), n_boot=200)
    assert np.isclose(summary.iloc[0]["estimate"], 0.5)
    assert np.isclose((merged["__truth__"] == merged["predicted_label"]).mean(), 100 / 101)


def test_paired_comparison_requires_common_donors_and_adjusts_p_values():
    frame = pd.DataFrame(
        {
            "method": ["celltypepilot"] * 4 + ["celltypist"] * 4,
            "donor_unit": ["d1", "d2", "d3", "d4"] * 2,
            "study_id": ["s"] * 8,
            "macro_f1": [0.9, 0.8, 0.7, 0.6, 0.6, 0.5, 0.4, 0.3],
            "balanced_accuracy": [0.9, 0.8, 0.7, 0.6, 0.6, 0.5, 0.4, 0.3],
            "coverage": [0.8] * 8,
        }
    )
    result = paired_method_comparisons(frame, n_resamples=1000)
    assert set(result["n_paired_donors"]) == {4}
    assert {"p_value", "p_value_bh", "mean_difference", "difference_ci_lower"} <= set(result)
    # Four same-direction paired differences have two equally extreme sign
    # assignments out of all 2**4 exact permutations.
    nonzero = result[result["metric"] != "coverage"]
    assert set(nonzero["p_value"]) == {0.125}


def test_sample_enrichment_has_explicit_not_assessed_state():
    missing = sample_enrichment_diagnostics(
        pd.DataFrame({"cluster": ["0"]}), cluster_key="cluster", sample_key="sample"
    )
    assert missing.iloc[0]["flag"] == "NOT_ASSESSED"

    frame = pd.DataFrame({"cluster": ["0"] * 9 + ["0"], "sample": ["s1"] * 9 + ["s2"]})
    result = sample_enrichment_diagnostics(
        frame, cluster_key="cluster", sample_key="sample", dominant_fraction_threshold=0.8
    )
    assert result.iloc[0]["flag"] == "SAMPLE_ENRICHED"


def test_comparator_label_map_is_exhaustive_and_idempotent():
    predictions = pd.DataFrame(
        {
            "method": ["m", "m"],
            "predicted_label": ["raw_a", "raw_b"],
        }
    )
    incomplete = pd.DataFrame({"method": ["m"], "raw_label": ["raw_a"], "canonical_label": ["A"]})
    with np.testing.assert_raises(BenchmarkValidationError):
        apply_locked_label_map(predictions, incomplete)

    complete = pd.DataFrame(
        {
            "method": ["m", "m"],
            "raw_label": ["raw_a", "raw_b"],
            "canonical_label": ["A", "B"],
        }
    )
    mapped = apply_locked_label_map(predictions, complete)
    assert mapped["predicted_label"].tolist() == ["A", "B"]
    pd.testing.assert_frame_equal(apply_locked_label_map(mapped, complete), mapped)


def test_comparator_label_map_allows_only_unambiguous_casefold_fallback():
    predictions = pd.DataFrame(
        {"method": ["celltypepilot"], "predicted_label": ["Endothelial cell"]}
    )
    locked = pd.DataFrame(
        {
            "method": ["celltypepilot"],
            "raw_label": ["endothelial cell"],
            "canonical_label": ["endothelial"],
        }
    )
    mapped = apply_locked_label_map(predictions, locked)
    assert mapped["predicted_label"].tolist() == ["endothelial"]
    assert mapped["raw_predicted_label"].tolist() == ["Endothelial cell"]

    ambiguous = pd.concat(
        [
            locked,
            pd.DataFrame(
                {
                    "method": ["celltypepilot"],
                    "raw_label": ["Endothelial Cell"],
                    "canonical_label": ["other"],
                }
            ),
        ],
        ignore_index=True,
    )
    with np.testing.assert_raises(BenchmarkValidationError):
        apply_locked_label_map(predictions, ambiguous)


def test_missing_qc_values_are_not_treated_as_clean_cells():
    merged = pd.DataFrame(
        {
            "method": ["celltypepilot"] * 3,
            "study": ["s"] * 3,
            "donor": ["d"] * 3,
            "__truth__": ["A"] * 3,
            "predicted_label": ["A"] * 3,
            "mito": [0.3, 0.1, np.nan],
        }
    )
    result = qc_stratified_performance(
        merged,
        study_key="study",
        donor_key="donor",
        diagnostics={"low_quality": {"key": "mito", "threshold": 0.2}},
    )
    missing = result[
        (result["diagnostic"] == "low_quality") & (result["stratum"] == "missing_diagnostic")
    ].iloc[0]
    assert missing["status"] == "not_assessed_missing_cell_values"
    assert np.isnan(missing["estimate"])


def _write_cohort(tmp_path, cohort_id: str, start: int) -> dict:
    cells = [f"{cohort_id}_c{i}" for i in range(8)]
    obs = pd.DataFrame(
        {
            "truth": ["A", "A", "B", "B"] * 2,
            "donor": [f"d{start}"] * 4 + [f"d{start + 1}"] * 4,
            "assay": ["p1"] * 4 + ["p2"] * 4,
            "sample": [f"x{start}"] * 4 + [f"x{start + 1}"] * 4,
            "condition": ["case"] * 4 + ["control"] * 4,
            "cluster": ["0", "0", "1", "1"] * 2,
        },
        index=cells,
    )
    data_path = tmp_path / f"{cohort_id}.h5ad"
    ad.AnnData(X=np.ones((8, 2)), obs=obs, var=pd.DataFrame(index=["G1", "G2"])).write_h5ad(
        data_path
    )
    assignments = build_holdout_assignments(obs.assign(study=cohort_id), "study", "donor", "donor")
    assignment_path = tmp_path / f"{cohort_id}.assignments.csv"
    assignments.to_csv(assignment_path, index=False)
    predictions = []
    for method in ("celltypepilot", "celltypist"):
        frame = assignments[["cell_id", "fold_id"]].copy()
        frame["method"] = method
        frame["predicted_label"] = obs.loc[frame["cell_id"], "truth"].to_numpy()
        if method == "celltypist":
            frame.loc[frame.index[-1], "predicted_label"] = "A"
        frame["confidence"] = 0.8
        predictions.append(frame)
    prediction_path = tmp_path / f"{cohort_id}.predictions.csv"
    pd.concat(predictions, ignore_index=True).to_csv(prediction_path, index=False)
    status_path = tmp_path / f"{cohort_id}.status.csv"
    pd.DataFrame(
        [
            {"method": method, "fold_id": fold, "status": "completed"}
            for method in ("celltypepilot", "celltypist")
            for fold in assignments["fold_id"].unique()
        ]
    ).to_csv(status_path, index=False)
    label_map_path = tmp_path / f"{cohort_id}.label_map.csv"
    pd.DataFrame(
        {
            "method": [
                "__truth__",
                "__truth__",
                "celltypepilot",
                "celltypepilot",
                "celltypist",
                "celltypist",
            ],
            "raw_label": ["A", "B", "A", "B", "A", "B"],
            "canonical_label": ["A", "B", "A", "B", "A", "B"],
        }
    ).to_csv(label_map_path, index=False)
    return {
        "cohort_id": cohort_id,
        "title": cohort_id,
        "species": "human",
        "tissue": "general",
        "collection_url": "https://example.org/collection",
        "dataset_url": "https://example.org/data.h5ad",
        "dataset_version_id": f"version-{cohort_id}",
        "citation_doi": "10.0000/example",
        "truth_provenance": "synthetic test truth",
        "constant_study_id": cohort_id,
        "local_path": data_path.name,
        "expected_bytes": data_path.stat().st_size,
        "expected_cells": 8,
        "expected_sha256": hashlib.sha256(data_path.read_bytes()).hexdigest(),
        "assignments_path": assignment_path.name,
        "predictions_path": prediction_path.name,
        "comparator_status_path": status_path.name,
        "label_map_path": label_map_path.name,
        "label_map_sha256": hashlib.sha256(label_map_path.read_bytes()).hexdigest(),
        "metadata": {
            "truth_key": "truth",
            "study_key": None,
            "donor_key": "donor",
            "platform_key": "assay",
            "sample_key": "sample",
            "condition_key": "condition",
            "batch_key": None,
            "cluster_key": "cluster",
        },
        "diagnostics": {},
    }


def test_build_public_release_writes_claim_bounded_artifacts(tmp_path):
    cohorts = [_write_cohort(tmp_path, "c1", 1), _write_cohort(tmp_path, "c2", 3)]
    registry = {
        "schema_version": "celltypepilot.public-benchmark-registry.v1",
        "release_id": "synthetic-v1",
        "required_methods": ["celltypepilot", "celltypist"],
        "analysis_plan": {
            "primary_metric": "macro_f1",
            "independent_unit": "donor",
            "holdout_policy": "leave-one-donor-out",
            "multiplicity_correction": "Benjamini-Hochberg",
            "negative_result_policy": "retain all",
        },
        "cohorts": cohorts,
    }
    validate_benchmark_registry(registry)
    registry_path = tmp_path / "registry.json"
    registry_path.write_text(json.dumps(registry), encoding="utf-8")
    output = tmp_path / "release"
    paths = build_public_benchmark_release(registry_path, output, n_boot=100)
    manifest = json.loads(paths["release_manifest"].read_text(encoding="utf-8"))
    assert manifest["readiness"]["status"] == "claim_ready"
    assert manifest["readiness"]["n_independent_donors"] == 4
    assert manifest["validation_scope"]["cell_weighted_inference_prohibited"] is True
    assert (output / "negative_results.csv").exists()
    assert (output / "batch_sensitivity_summary.csv").exists()
    assert (output / "benchmark_report.md").exists()


def test_incomplete_release_retains_comparator_failure_without_predictions(tmp_path):
    cohorts = [_write_cohort(tmp_path, "c1", 1), _write_cohort(tmp_path, "c2", 3)]
    for cohort in cohorts:
        (tmp_path / cohort["predictions_path"]).unlink()
        status_path = tmp_path / cohort["comparator_status_path"]
        pd.DataFrame(
            [
                {
                    "method": "celltypepilot",
                    "fold_id": "preflight",
                    "status": "failed_runtime_limit",
                    "detail": "exceeded locked limit",
                },
                {
                    "method": "celltypist",
                    "fold_id": "preflight",
                    "status": "dependency_unavailable",
                    "detail": "not installed",
                },
            ]
        ).to_csv(status_path, index=False)
    registry = {
        "schema_version": "celltypepilot.public-benchmark-registry.v1",
        "release_id": "synthetic-incomplete-v1",
        "required_methods": ["celltypepilot", "celltypist"],
        "analysis_plan": {
            "primary_metric": "macro_f1",
            "independent_unit": "donor",
            "holdout_policy": "leave-one-donor-out",
            "multiplicity_correction": "Benjamini-Hochberg",
            "negative_result_policy": "retain all",
        },
        "cohorts": cohorts,
    }
    registry_path = tmp_path / "registry.json"
    registry_path.write_text(json.dumps(registry), encoding="utf-8")
    paths = build_public_benchmark_release(
        registry_path,
        tmp_path / "release",
        allow_incomplete=True,
        n_boot=100,
    )
    manifest = json.loads(paths["release_manifest"].read_text(encoding="utf-8"))
    negative = pd.read_csv(paths["negative_results"])
    assert manifest["readiness"]["status"] == "incomplete_not_claim_ready"
    assert "celltypepilot:failed_runtime_limit" in set(negative["finding"])
    assert "celltypist:dependency_unavailable" in set(negative["finding"])
    with np.testing.assert_raises(BenchmarkValidationError):
        build_public_benchmark_release(
            registry_path,
            tmp_path / "release_strict",
            n_boot=100,
        )
