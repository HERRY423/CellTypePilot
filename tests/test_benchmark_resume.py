from __future__ import annotations

import json
import tempfile

import anndata as ad
import numpy as np
import pandas as pd

from celltypepilot import benchmark_runner
from celltypepilot.benchmark import build_holdout_assignments


def _tiny_benchmark():
    obs = pd.DataFrame(
        {
            "truth": ["A", "A", "B", "B"],
            "study": ["s"] * 4,
            "donor": ["d1", "d1", "d2", "d2"],
            "cluster": ["0", "0", "1", "1"],
        },
        index=["c1", "c2", "c3", "c4"],
    )
    dataset = ad.AnnData(
        X=np.ones((4, 2)),
        obs=obs,
        var=pd.DataFrame(index=["G1", "G2"]),
    )
    assignments = build_holdout_assignments(obs, "study", "donor", "donor")
    return dataset, assignments


def test_benchmark_runtime_uses_output_owned_writable_cache(tmp_path):
    configured = benchmark_runner.configure_benchmark_runtime(tmp_path)
    assert tempfile.gettempdir() == configured["temp_dir"]
    assert (tmp_path / "_runtime_cache" / "tmp").is_dir()
    assert (tmp_path / "_runtime_cache" / "numba").is_dir()


def test_completed_fold_checkpoints_resume_without_reexecution(tmp_path, monkeypatch):
    dataset, assignments = _tiny_benchmark()
    calls: list[str] = []

    def fake_run(paths, cluster_key, species, tissue):
        test = ad.read_h5ad(paths["test"])
        calls.append(str(test.obs_names[0]))
        return (
            pd.DataFrame(
                {
                    "cell_id": test.obs_names.astype(str),
                    "predicted_label": "A",
                    "confidence": 0.5,
                }
            ),
            {
                "implementation": "test",
                "version": "1",
                "reference_policy": "fold_train_only",
                "confidence_semantics": "test_score",
            },
        )

    monkeypatch.setattr(benchmark_runner, "run_celltypepilot_fold", fake_run)
    first_predictions, first_status = benchmark_runner.run_benchmark_comparators(
        dataset,
        assignments,
        "truth",
        "cluster",
        tmp_path,
        "human",
        "general",
        methods=("celltypepilot",),
    )
    assert len(calls) == 2
    assert len(first_predictions) == 4
    assert set(first_status["status"]) == {"completed"}
    assert (tmp_path / "out_of_fold_predictions.csv").exists()
    assert (tmp_path / "comparator_status.csv").exists()

    def must_not_run(*args, **kwargs):
        raise AssertionError("completed checkpoint was executed again")

    monkeypatch.setattr(benchmark_runner, "run_celltypepilot_fold", must_not_run)
    resumed_predictions, resumed_status = benchmark_runner.run_benchmark_comparators(
        dataset,
        assignments,
        "truth",
        "cluster",
        tmp_path,
        "human",
        "general",
        methods=("celltypepilot",),
    )
    assert len(resumed_predictions) == 4
    assert resumed_status["resumed_from_checkpoint"].fillna(False).all()


def test_running_checkpoint_is_machine_readable_before_comparator_returns(tmp_path, monkeypatch):
    dataset, assignments = _tiny_benchmark()

    def inspect_then_fail(paths, cluster_key, species, tissue):
        statuses = list((tmp_path / "checkpoints").glob("*.status.json"))
        assert statuses
        assert json.loads(statuses[0].read_text(encoding="utf-8"))["status"] == "running"
        raise RuntimeError("intentional interruption")

    monkeypatch.setattr(benchmark_runner, "run_celltypepilot_fold", inspect_then_fail)
    predictions, status = benchmark_runner.run_benchmark_comparators(
        dataset,
        assignments,
        "truth",
        "cluster",
        tmp_path,
        "human",
        "general",
        methods=("celltypepilot",),
        continue_on_unavailable=True,
    )
    assert predictions.empty
    assert set(status["status"]) == {"failed_or_unavailable"}
    assert "intentional interruption" in " ".join(status["error"].astype(str))


def test_native_comparator_methods_use_shared_fold_runner(tmp_path, monkeypatch):
    dataset, assignments = _tiny_benchmark()
    calls = []

    def fake_native(
        train,
        query,
        cluster_key,
        backend,
        run_dir,
        *,
        species,
        tissue,
        entry_overrides=None,
    ):
        test = ad.read_h5ad(query)
        calls.append(backend)
        return (
            pd.DataFrame(
                {
                    "cell_id": test.obs_names.astype(str),
                    "predicted_label": "A",
                    "confidence": 0.5,
                }
            ),
            {
                "implementation": "shared-native-test",
                "version": "1",
                "reference_policy": "fold_train_only",
                "confidence_semantics": "test",
            },
        )

    monkeypatch.setattr("celltypepilot.native_backends.run_fold_native_backend", fake_native)
    predictions, status = benchmark_runner.run_benchmark_comparators(
        dataset,
        assignments,
        "truth",
        "cluster",
        tmp_path,
        "human",
        "general",
        methods=("scanvi",),
    )

    assert calls == ["scanvi", "scanvi"]
    assert len(predictions) == 4
    assert set(status["implementation"]) == {"shared-native-test"}
