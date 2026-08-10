"""GPU batch isolation + distributed checkpoint worker contracts."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from celltypepilot import benchmark_runner
from celltypepilot.benchmark import build_holdout_assignments

REPO = Path(__file__).resolve().parents[1]
SCRATCH = REPO / "scratch" / "pytest_temp" / "gpu_batch"


@pytest.fixture(autouse=True)
def _restore_temp_env():
    """configure_benchmark_runtime rewires TMP; restore so other tests keep working."""
    old_tmp = os.environ.get("TMP")
    old_temp = os.environ.get("TEMP")
    old_numba = os.environ.get("NUMBA_CACHE_DIR")
    old_tempdir = tempfile.tempdir
    yield
    if old_tmp is None:
        os.environ.pop("TMP", None)
    else:
        os.environ["TMP"] = old_tmp
    if old_temp is None:
        os.environ.pop("TEMP", None)
    else:
        os.environ["TEMP"] = old_temp
    if old_numba is None:
        os.environ.pop("NUMBA_CACHE_DIR", None)
    else:
        os.environ["NUMBA_CACHE_DIR"] = old_numba
    tempfile.tempdir = old_tempdir


def _tiny():
    obs = pd.DataFrame(
        {
            "truth": ["A", "A", "B", "B"],
            "study": ["s"] * 4,
            "donor": ["d1", "d1", "d2", "d2"],
            "cluster": ["0", "0", "1", "1"],
        },
        index=["c1", "c2", "c3", "c4"],
    )
    data = ad.AnnData(X=np.ones((4, 2)), obs=obs, var=pd.DataFrame(index=["G1", "G2"]))
    assignments = build_holdout_assignments(obs, "study", "donor", "donor")
    return data, assignments


def test_worker_fold_filter_writes_only_assigned_checkpoints(monkeypatch):
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    out = SCRATCH / "worker_a"
    out.mkdir(parents=True)
    dataset, assignments = _tiny()
    folds = list(assignments["fold_id"].astype(str).unique())
    target = folds[0]

    def fake_run(paths, *args, **kwargs):
        test = ad.read_h5ad(paths["test"])
        return (
            pd.DataFrame(
                {
                    "cell_id": test.obs_names.astype(str),
                    "predicted_label": "A",
                    "confidence": 0.9,
                }
            ),
            {
                "implementation": "test",
                "version": "1",
                "reference_policy": "fold_train_only",
                "confidence_semantics": "test",
            },
        )

    monkeypatch.setattr(benchmark_runner, "run_celltypepilot_fold", fake_run)
    predictions, status = benchmark_runner.run_benchmark_comparators(
        dataset,
        assignments,
        "truth",
        "cluster",
        out,
        "human",
        "lung",
        methods=("celltypepilot",),
        fold_ids=(target,),
        write_aggregate_tables=False,
        worker_id="gpu-node-a",
        batch_id="gpu_popv_retrain_v1",
    )
    assert len(predictions) == 2
    assert set(status["fold_id"].astype(str)) == {target}
    # No global OOF rewrite in distributed mode
    assert not (out / "out_of_fold_predictions.csv").exists()
    ck = list((out / "checkpoints").glob("*.status.json"))
    assert len(ck) == 1
    payload = json.loads(ck[0].read_text(encoding="utf-8"))
    assert payload.get("worker_id") == "gpu-node-a"
    assert payload.get("batch_id") == "gpu_popv_retrain_v1"


def test_second_worker_does_not_clobber_sibling_checkpoints(monkeypatch):
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    out = SCRATCH / "shared"
    out.mkdir(parents=True)
    dataset, assignments = _tiny()
    folds = list(assignments["fold_id"].astype(str).unique())

    def fake_run(paths, *args, **kwargs):
        test = ad.read_h5ad(paths["test"])
        return (
            pd.DataFrame(
                {
                    "cell_id": test.obs_names.astype(str),
                    "predicted_label": "B",
                    "confidence": 0.5,
                }
            ),
            {
                "implementation": "test",
                "version": "1",
                "reference_policy": "fold_train_only",
                "confidence_semantics": "test",
            },
        )

    monkeypatch.setattr(benchmark_runner, "run_celltypepilot_fold", fake_run)
    benchmark_runner.run_benchmark_comparators(
        dataset,
        assignments,
        "truth",
        "cluster",
        out,
        "human",
        "lung",
        methods=("celltypepilot",),
        fold_ids=(folds[0],),
        write_aggregate_tables=False,
        worker_id="node-a",
        batch_id="gpu_popv_retrain_v1",
    )
    first_files = {p.name: p.read_bytes() for p in (out / "checkpoints").iterdir()}
    benchmark_runner.run_benchmark_comparators(
        dataset,
        assignments,
        "truth",
        "cluster",
        out,
        "human",
        "lung",
        methods=("celltypepilot",),
        fold_ids=(folds[1],),
        write_aggregate_tables=False,
        worker_id="node-b",
        batch_id="gpu_popv_retrain_v1",
    )
    # First worker's files unchanged
    for name, content in first_files.items():
        assert (out / "checkpoints" / name).read_bytes() == content
    assert len(list((out / "checkpoints").glob("*.status.json"))) == 2


def test_gpu_batch_manifest_isolation_policy():
    manifest = json.loads(
        (
            REPO
            / "benchmarks"
            / "public_v1"
            / "batches"
            / "gpu_popv_retrain_v1"
            / "batch_manifest.json"
        ).read_text(encoding="utf-8")
    )
    assert manifest["device_track"] == "gpu"
    assert "travaglini_lung_smartseq2_2020" in str(manifest["isolation"]["must_not_write_into"])
    assert len(manifest["required_folds"]) == 3
    assert manifest["methods"] == ["popv"]


def test_aggregate_script_read_only_merge():
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    batch = SCRATCH / "batch"
    run_ck = batch / "run" / "checkpoints"
    run_ck.mkdir(parents=True)
    # Minimal batch manifest
    (batch / "batch_manifest.json").write_text(
        json.dumps(
            {
                "batch_id": "gpu_popv_retrain_v1",
                "required_folds": ["donor=travaglini_2020::1", "donor=travaglini_2020::2"],
            }
        ),
        encoding="utf-8",
    )
    from celltypepilot.benchmark_runner import _safe_fold_name

    for fold in ["donor=travaglini_2020::1", "donor=travaglini_2020::2"]:
        stem = f"{_safe_fold_name(fold)}__popv"
        (run_ck / f"{stem}.status.json").write_text(
            json.dumps(
                {
                    "method": "popv",
                    "fold_id": fold,
                    "status": "completed",
                    "worker_id": "n1",
                    "batch_id": "gpu_popv_retrain_v1",
                }
            ),
            encoding="utf-8",
        )
        pd.DataFrame(
            {
                "cell_id": [f"{fold}-c1", f"{fold}-c2"],
                "predicted_label": ["epithelial", "endothelial"],
                "confidence": [0.8, 0.7],
                "method": ["popv", "popv"],
                "fold_id": [fold, fold],
            }
        ).to_csv(run_ck / f"{stem}.csv", index=False)

    import importlib.util
    import sys

    spec = importlib.util.spec_from_file_location(
        "aggregate_gpu_batch_checkpoints",
        REPO / "scripts" / "aggregate_gpu_batch_checkpoints.py",
    )
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    old = sys.argv
    try:
        sys.argv = [
            "aggregate_gpu_batch_checkpoints.py",
            "--batch-root",
            str(batch),
            "--require-complete",
        ]
        rc = mod.main()
    finally:
        sys.argv = old
    assert rc == 0
    assert (batch / "aggregate" / "out_of_fold_predictions.csv").is_file()
    pred = pd.read_csv(batch / "aggregate" / "out_of_fold_predictions.csv")
    assert len(pred) == 4
    report = json.loads((batch / "aggregate" / "aggregate_report.json").read_text(encoding="utf-8"))
    assert report["complete"] is True
    assert report["device_track"] == "gpu"
