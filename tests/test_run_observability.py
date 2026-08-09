"""Read-only run observability tests (never write live fold workspaces)."""

from __future__ import annotations

import json
import shutil
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from celltypepilot.run_observability import (
    ObservabilityError,
    build_observability_snapshot,
    estimate_fold_eta,
    load_checkpoint_records,
)

REPO = Path(__file__).resolve().parents[1]
LIVE_RUN = REPO / "benchmarks" / "public_v1" / "runs" / "travaglini_lung_smartseq2_2020"
SCRATCH = REPO / "scratch" / "pytest_temp" / "run_observability"


def _write_status(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def test_snapshot_reads_live_run_without_mutation():
    if not (LIVE_RUN / "checkpoints").is_dir():
        pytest.skip("live smartseq2 run not present")
    before = {
        p.name: p.stat().st_mtime_ns
        for p in (LIVE_RUN / "checkpoints").glob("*.status.json")
    }
    snapshot = build_observability_snapshot(LIVE_RUN, include_host=True)
    after = {
        p.name: p.stat().st_mtime_ns
        for p in (LIVE_RUN / "checkpoints").glob("*.status.json")
    }
    assert before == after, "observability must not touch checkpoint mtimes"
    assert snapshot["read_only"] is True
    assert snapshot["prediction_mutation_allowed"] is False
    assert snapshot["checkpoints"]["n_status_files"] >= 1
    assert "fold_eta" in snapshot
    assert "stale" in snapshot
    assert "products" in snapshot


def test_eta_from_completed_and_running():
    now = datetime.now(timezone.utc)
    records = [
        {
            "method": "celltypepilot",
            "fold_id": "f1",
            "status": "completed",
            "duration_seconds": 100.0,
            "elapsed_seconds": 100.0,
        },
        {
            "method": "celltypepilot",
            "fold_id": "f2",
            "status": "running",
            "duration_seconds": None,
            "elapsed_seconds": 40.0,
        },
    ]
    eta = estimate_fold_eta(records)
    assert eta["n_completed"] == 1
    assert eta["n_running"] == 1
    assert eta["mean_completed_duration_seconds"] == 100.0
    assert eta["estimated_remaining_seconds"] == pytest.approx(60.0)


def test_failure_reason_and_product_hash():
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    run = SCRATCH / "run_a"
    ck = run / "checkpoints"
    ck.mkdir(parents=True)
    _write_status(
        ck / "donor_x__celltypist.status.json",
        {
            "method": "celltypist",
            "fold_id": "donor=x",
            "status": "failed_or_unavailable",
            "started_at_utc": "2026-08-09T00:00:00+00:00",
            "failed_at_utc": "2026-08-09T00:01:00+00:00",
            "error": "dependency_unavailable",
        },
    )
    pred = ck / "donor_x__celltypepilot.csv"
    pred.write_text("cell_id,predicted_label,confidence\na,T,0.9\n", encoding="utf-8")
    _write_status(
        ck / "donor_x__celltypepilot.status.json",
        {
            "method": "celltypepilot",
            "fold_id": "donor=x",
            "status": "completed",
            "started_at_utc": "2026-08-09T00:00:00+00:00",
            "completed_at_utc": "2026-08-09T00:00:30+00:00",
            "provenance": {"implementation": "test"},
        },
    )
    (run / "artifact_status.json").write_text(
        json.dumps(
            {
                "review_state": "applied_overrides_artifacts_stale",
                "derived_artifacts_stale": True,
                "stale_reason": "test stale",
                "stale_artifacts": ["evidence_table.csv"],
            }
        ),
        encoding="utf-8",
    )
    snapshot = build_observability_snapshot(run, include_host=False)
    assert any(f["failure_reason"] == "dependency_unavailable" for f in snapshot["failures"])
    completed = [
        r for r in snapshot["checkpoints"]["records"] if r["method"] == "celltypepilot"
    ][0]
    assert completed["prediction_csv_present"] is True
    assert completed["prediction_sha256"]
    assert snapshot["stale"]["derived_artifacts_stale"] is True
    products = {p["name"]: p for p in snapshot["products"]}
    assert products["artifact_status.json"]["present"] is True
    assert products["artifact_status.json"]["sha256"]


def test_api_observability_read_only_and_mutation_blocked():
    from celltypepilot import web_inspector

    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    run = SCRATCH / "run_api"
    ck = run / "checkpoints"
    ck.mkdir(parents=True)
    _write_status(
        ck / "f1__celltypepilot.status.json",
        {
            "method": "celltypepilot",
            "fold_id": "f1",
            "status": "running",
            "started_at_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    web_inspector._output_dir = run
    web_inspector._run_dir = run
    web_inspector._overrides = {}
    web_inspector._adata_cache = None
    web_inspector._evidence_cache = None

    client = web_inspector.app.test_client()
    resp = client.get("/api/observability")
    assert resp.status_code == 200
    body = resp.get_json()
    assert body["ok"] is True
    assert body["observability"]["read_only"] is True
    assert body["observability"]["prediction_mutation_allowed"] is False
    assert body["observability"]["checkpoints"]["n_status_files"] == 1

    blocked = client.post(
        "/api/observability/predictions",
        json={"cell_id": "x", "predicted_label": "T"},
    )
    assert blocked.status_code == 405
    assert blocked.get_json()["prediction_mutation_allowed"] is False

    # Override into a pure benchmark run must be blocked (protect fold trees).
    ov = client.post(
        "/api/override",
        json={"cluster": "0", "new_type": "T cell", "reason": "should fail"},
    )
    assert ov.status_code == 403

    apply = client.post("/api/overrides/apply")
    assert apply.status_code == 403

    # No new files under checkpoints from the blocked write attempts.
    assert list(ck.glob("*"))  # still only status
    names = {p.name for p in ck.iterdir()}
    assert names == {"f1__celltypepilot.status.json"}


def test_missing_path_raises():
    with pytest.raises(ObservabilityError):
        build_observability_snapshot(REPO / "does_not_exist_observability")
