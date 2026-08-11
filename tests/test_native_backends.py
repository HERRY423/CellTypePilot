"""Contracts for native annotate backends and depth-domain planning."""

from __future__ import annotations

import json
import sys
import types
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from celltypepilot.domain_validation_pipeline import (
    build_domain_validation_plan,
    execute_domain_validation_plan,
)
from celltypepilot.native_backend_config import (
    NativeBackendConfigError,
    hash_native_backend_dependencies,
    load_native_backend_config,
)
from celltypepilot.native_backends import (
    NativeBackendRunError,
    NativeBackendUnavailableError,
    _run_reference_backend,
    run_fold_native_backend,
    run_native_backends,
)


def _query() -> ad.AnnData:
    obs = pd.DataFrame(
        {"leiden": ["0", "0", "1", "1"]},
        index=["cell-1", "cell-2", "cell-3", "cell-4"],
    )
    var = pd.DataFrame(index=[f"G{i}" for i in range(60)])
    return ad.AnnData(np.ones((4, 60)), obs=obs, var=var)


def _reference(path: Path) -> Path:
    obs = pd.DataFrame({"cell_type": ["T cell", "B cell"]}, index=["ref-1", "ref-2"])
    var = pd.DataFrame(index=[f"G{i}" for i in range(60)])
    ad.AnnData(np.ones((2, 60)), obs=obs, var=var).write_h5ad(path)
    return path


def _resolver() -> dict:
    return {
        "by_label": {"t cell": "T cell", "b cell": "B cell"},
        "by_cl": {},
        "cl_by_name": {"T cell": "CL:0000084", "B cell": "CL:0000236"},
        "safe_parent_fallbacks": {},
    }


def _write_config(path: Path, reference: Path, backend: str = "popv") -> Path:
    path.write_text(
        json.dumps(
            {
                "schema_version": "celltypepilot.native-backends.v1",
                "continue_on_failure": True,
                "resume": True,
                "backends": [
                    {
                        "backend": backend,
                        "reference_path": str(reference),
                        "label_key": "cell_type",
                        **({"mode": "retrain"} if backend == "popv" else {}),
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    return path


def test_native_config_resolves_and_hashes_reference(tmp_path):
    reference = _reference(tmp_path / "reference.h5ad")
    config = load_native_backend_config(_write_config(tmp_path / "native.json", reference))

    assert config["backends"][0]["reference_path"] == str(reference.resolve())
    assert hash_native_backend_dependencies(config)[str(reference.resolve())]
    assert len(config["config_sha256"]) == 64


def test_llm_config_requires_explicit_network_authorization(tmp_path):
    path = tmp_path / "native.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": "celltypepilot.native-backends.v1",
                "backends": [{"backend": "llm", "model": "bounded-model"}],
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(NativeBackendConfigError, match="allow_network=true"):
        load_native_backend_config(path)


def test_native_backend_cell_outputs_are_aggregated_and_resumed(tmp_path, monkeypatch):
    reference = _reference(tmp_path / "reference.h5ad")
    config = load_native_backend_config(_write_config(tmp_path / "native.json", reference))
    calls = []

    def fake_popv(query, cluster_key, entry, run_dir):
        calls.append(entry["backend"])
        return (
            pd.DataFrame(
                {
                    "cell_id": query.obs_names.astype(str),
                    "predicted_label": ["T cell", "T cell", "B cell", "B cell"],
                    "confidence": [0.8, 0.9, 0.7, 0.6],
                    "backend": "popv",
                }
            ),
            {"backend_version": "test"},
        )

    monkeypatch.setattr("celltypepilot.native_backends._run_popv", fake_popv)
    first = run_native_backends(
        _query(),
        "leiden",
        pd.DataFrame(),
        _resolver(),
        tmp_path / "output",
        config,
        species="human",
        tissue="general",
        input_sha256="input-hash",
    )
    second = run_native_backends(
        _query(),
        "leiden",
        pd.DataFrame(),
        _resolver(),
        tmp_path / "output",
        config,
        species="human",
        tissue="general",
        input_sha256="input-hash",
    )

    assert calls == ["popv"]
    assert set(first["candidates"]["canonical_cell_type"]) == {"T cell", "B cell"}
    assert set(first["candidates"]["decision_role"]) == {"decision_candidate"}
    assert (tmp_path / "output/native_backends/popv/raw_candidates.csv").is_file()
    assert second["status"].iloc[0]["status"] == "completed_from_checkpoint"


def test_unavailable_backend_is_retained_without_candidate_fallback(tmp_path, monkeypatch):
    reference = _reference(tmp_path / "reference.h5ad")
    config = load_native_backend_config(_write_config(tmp_path / "native.json", reference))

    def unavailable(*args, **kwargs):
        raise NativeBackendUnavailableError("popV unavailable in test")

    monkeypatch.setattr("celltypepilot.native_backends._run_popv", unavailable)
    result = run_native_backends(
        _query(),
        "leiden",
        pd.DataFrame(),
        _resolver(),
        tmp_path / "output",
        config,
        species="human",
        tissue="general",
        input_sha256="input-hash",
    )

    assert result["candidates"].empty
    assert result["status"].iloc[0]["status"] == "unavailable"
    assert result["status"].iloc[0]["error_type"] == "NativeBackendUnavailableError"


def test_llm_runner_is_strictly_hypothesis_only(tmp_path, monkeypatch):
    config_path = tmp_path / "llm.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": "celltypepilot.native-backends.v1",
                "backends": [
                    {
                        "backend": "llm",
                        "provider": "openai",
                        "model": "bounded-test-model",
                        "allow_network": True,
                        "api_key_env": "TEST_CTP_API_KEY",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    class FakeResponses:
        def create(self, **kwargs):
            assert kwargs["text"]["format"]["strict"] is True
            return types.SimpleNamespace(
                id="response-test",
                output_text=json.dumps(
                    {
                        "decisions": [
                            {
                                "cluster": "0",
                                "selected_label": "T cell",
                                "abstain": False,
                                "rationale": "bounded test",
                            },
                            {
                                "cluster": "1",
                                "selected_label": "",
                                "abstain": True,
                                "rationale": "insufficient evidence",
                            },
                        ]
                    }
                ),
                usage=None,
            )

    class FakeOpenAI:
        def __init__(self, api_key):
            assert api_key == "secret"
            self.responses = FakeResponses()

    monkeypatch.setitem(sys.modules, "openai", types.SimpleNamespace(OpenAI=FakeOpenAI))
    monkeypatch.setenv("TEST_CTP_API_KEY", "secret")
    marker_scores = pd.DataFrame(
        {
            "cluster": ["0", "1"],
            "cell_type": ["T cell", "B cell"],
            "combined_score": [0.8, 0.7],
        }
    )
    result = run_native_backends(
        _query(),
        "leiden",
        marker_scores,
        _resolver(),
        tmp_path / "output",
        load_native_backend_config(config_path),
        species="human",
        tissue="general",
        input_sha256="input-hash",
    )

    assert len(result["candidates"]) == 1
    assert result["candidates"].iloc[0]["decision_role"] == "hypothesis_only"
    assert result["candidates"].iloc[0]["independence_group"] == "llm_hypothesis"


def test_fold_scanvi_dispatch_does_not_silently_use_correlation(tmp_path, monkeypatch):
    train = _reference(tmp_path / "train.h5ad")
    query_path = tmp_path / "query.h5ad"
    _query().write_h5ad(query_path)
    observed = {}

    def fake_reference(query, cluster_key, entry, species, tissue):
        observed.update(entry)
        return (
            pd.DataFrame(
                {
                    "cluster": ["0", "1"],
                    "cell_type": ["T cell", "B cell"],
                    "score": [0.8, 0.7],
                    "rank": [1, 1],
                }
            ),
            {"backend_version": "test"},
        )

    monkeypatch.setattr("celltypepilot.native_backends._run_reference_backend", fake_reference)
    predictions, _ = run_fold_native_backend(
        train,
        query_path,
        "leiden",
        "scanvi",
        tmp_path / "fold",
        species="human",
        tissue="general",
    )

    assert observed["backend"] == "scanvi"
    assert "method" not in observed
    assert len(predictions) == 4


def test_scanvi_native_runner_rejects_non_count_expression(tmp_path):
    reference = _reference(tmp_path / "reference.h5ad")
    query = _query()
    query.X = np.full(query.shape, 0.5)

    with pytest.raises(NativeBackendRunError, match="integer raw counts"):
        _run_reference_backend(
            query,
            "leiden",
            {
                "backend": "scanvi",
                "reference_path": str(reference),
                "label_key": "cell_type",
            },
            "human",
            "general",
        )


def test_fold_runner_reuses_product_raw_candidates(tmp_path, monkeypatch):
    train = _reference(tmp_path / "train.h5ad")
    query_path = tmp_path / "query.h5ad"
    _query().write_h5ad(query_path)
    raw_path = tmp_path / "product/popv/raw_candidates.csv"
    raw_path.parent.mkdir(parents=True)
    pd.DataFrame(
        {
            "cell_id": ["cell-1", "cell-2", "cell-3", "cell-4"],
            "predicted_label": ["T cell", "T cell", "B cell", "B cell"],
            "confidence": [0.8, 0.8, 0.7, 0.7],
        }
    ).to_csv(raw_path, index=False)

    def must_not_run(*args, **kwargs):
        raise AssertionError("popV was re-executed instead of reusing the product artifact")

    monkeypatch.setattr("celltypepilot.native_backends._run_popv", must_not_run)
    predictions, provenance = run_fold_native_backend(
        train,
        query_path,
        "leiden",
        "popv",
        tmp_path / "fold",
        species="human",
        tissue="general",
        entry_overrides={"reuse_raw_path": str(raw_path)},
    )

    assert len(predictions) == 4
    assert provenance["reused_product_native_artifact"] == str(raw_path)


def test_custom_reference_native_runner_executes_without_optional_runtime(tmp_path):
    genes = [f"G{i}" for i in range(250)]
    reference_x = np.vstack(
        [
            np.r_[np.full(125, 8.0), np.ones(125)],
            np.r_[np.full(125, 7.0), np.ones(125)],
            np.r_[np.ones(125), np.full(125, 8.0)],
            np.r_[np.ones(125), np.full(125, 7.0)],
        ]
    )
    reference_path = tmp_path / "reference.h5ad"
    ad.AnnData(
        reference_x,
        obs=pd.DataFrame(
            {"cell_type": ["T cell", "T cell", "B cell", "B cell"]},
            index=["r1", "r2", "r3", "r4"],
        ),
        var=pd.DataFrame(index=genes),
    ).write_h5ad(reference_path)
    query = ad.AnnData(
        reference_x.copy(),
        obs=pd.DataFrame(
            {"leiden": ["0", "0", "1", "1"]},
            index=["cell-1", "cell-2", "cell-3", "cell-4"],
        ),
        var=pd.DataFrame(index=genes),
    )
    config_path = tmp_path / "native.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": "celltypepilot.native-backends.v1",
                "backends": [
                    {
                        "backend": "custom_reference",
                        "method": "correlation",
                        "reference_path": str(reference_path),
                        "label_key": "cell_type",
                        "allow_unverified_reference": True,
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    result = run_native_backends(
        query,
        "leiden",
        pd.DataFrame(),
        _resolver(),
        tmp_path / "output",
        load_native_backend_config(config_path),
        species="human",
        tissue="general",
        input_sha256="input-hash",
    )

    assert result["status"].iloc[0]["status"] == "completed"
    assert set(result["candidates"].query("rank == 1")["canonical_cell_type"]) == {
        "T cell",
        "B cell",
    }


def test_domain_plan_exposes_current_evidence_deficits(tmp_path):
    registry = tmp_path / "registry.json"
    registry.write_text(
        json.dumps(
            {
                "schema_version": "celltypepilot.public-benchmark-registry.v1",
                "cohorts": [
                    {
                        "cohort_id": "lung_missing",
                        "title": "Lung cohort",
                        "species": "human",
                        "tissue": "lung",
                        "constant_study_id": "study-1",
                        "local_path": "missing.h5ad",
                        "label_map_path": "missing.csv",
                        "cluster_map_path": "missing-clusters.csv",
                        "metadata": {
                            "truth_key": "cell_type",
                            "study_key": None,
                            "donor_key": "donor",
                            "platform_key": "assay",
                            "cluster_key": "cluster",
                        },
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    result = build_domain_validation_plan(registry, tmp_path / "out")
    lung = result["plan"]["domains"]["lung"]

    assert lung["registered_cohorts"] == 1
    assert lung["ready_cohorts"] == 0
    assert lung["claim_ready"] is False
    assert "INSUFFICIENT_PUBLIC_COHORTS" in lung["blockers"]
    assert "COHORT_ASSETS_OR_LABEL_MAPS_INCOMPLETE" in lung["blockers"]


def test_domain_run_requires_every_method_fold_checkpoint(tmp_path, monkeypatch):
    data_path = tmp_path / "lung.h5ad"
    obs = pd.DataFrame(
        {
            "cell_type": ["A", "A", "B", "B"],
            "donor": ["d1", "d1", "d2", "d2"],
            "assay": ["10x"] * 4,
            "cluster": ["0", "0", "1", "1"],
        },
        index=["c1", "c2", "c3", "c4"],
    )
    ad.AnnData(np.ones((4, 2)), obs=obs, var=pd.DataFrame(index=["G1", "G2"])).write_h5ad(data_path)
    label_map_path = tmp_path / "labels.csv"
    rows = []
    for method in ["__truth__", "celltypepilot", "celltypist", "popv", "singler"]:
        rows.extend(
            {"method": method, "raw_label": label, "canonical_label": label} for label in ("A", "B")
        )
    pd.DataFrame(rows).to_csv(label_map_path, index=False)
    registry = tmp_path / "registry.json"
    registry.write_text(
        json.dumps(
            {
                "schema_version": "celltypepilot.public-benchmark-registry.v1",
                "cohorts": [
                    {
                        "cohort_id": "lung_ready",
                        "title": "Lung cohort",
                        "species": "human",
                        "tissue": "lung",
                        "constant_study_id": "study-1",
                        "donor_namespace": "study-1",
                        "local_path": str(data_path),
                        "label_map_path": str(label_map_path),
                        "metadata": {
                            "truth_key": "cell_type",
                            "study_key": None,
                            "donor_key": "donor",
                            "platform_key": "assay",
                            "cluster_key": "cluster",
                        },
                    }
                ],
            }
        ),
        encoding="utf-8",
    )
    planned = build_domain_validation_plan(registry, tmp_path / "out")
    assert planned["plan"]["cohorts"][0]["execution_status"] == "ready"

    def fake_comparators(adata, assignments, *args, methods, **kwargs):
        folds = list(assignments["fold_id"].drop_duplicates().astype(str))
        status = pd.DataFrame(
            [
                {"method": method, "fold_id": fold, "status": "completed"}
                for method in methods
                for fold in folds
                if not (method == "scanvi" and fold == folds[-1])
            ]
        )
        return pd.DataFrame(), status

    monkeypatch.setattr(
        "celltypepilot.benchmark_runner.run_benchmark_comparators", fake_comparators
    )
    result = execute_domain_validation_plan(planned["plan_path"])

    assert result["status"].iloc[0]["status"] == "incomplete"
    assert result["manifest"]["claim_ready"] is False
