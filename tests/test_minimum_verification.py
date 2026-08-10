from __future__ import annotations

import json

import pandas as pd
import pytest
from scripts.build_minimum_verification import (
    checkpoint_runtimes,
    completed_fold_check,
)

from celltypepilot.benchmark import BenchmarkValidationError


def test_completed_fold_check_requires_exact_method_fold_sets() -> None:
    assignments = pd.DataFrame({"cell_id": ["a", "b"], "fold_id": ["donor=1", "donor=2"]})
    status = pd.DataFrame(
        [
            {"method": method, "fold_id": fold, "status": "completed"}
            for method in ("celltypepilot", "popv")
            for fold in ("donor=1", "donor=2")
        ]
    )

    completed_fold_check(status, assignments, ("celltypepilot", "popv"))
    with pytest.raises(BenchmarkValidationError, match="exact completed-fold"):
        completed_fold_check(
            status[~((status["method"] == "popv") & (status["fold_id"] == "donor=2"))],
            assignments,
            ("celltypepilot", "popv"),
        )


def test_checkpoint_runtimes_retains_only_completed(tmp_path) -> None:
    completed = {
        "method": "popv",
        "fold_id": "donor=1",
        "status": "completed",
        "started_at_utc": "2026-08-09T00:00:00+00:00",
        "completed_at_utc": "2026-08-09T00:01:30+00:00",
        "provenance": {
            "version": "0.6.1",
            "reference_policy": "fold_train_only",
            "confidence_semantics": "expert_agreement_not_probability",
        },
    }
    running = {
        "method": "popv",
        "fold_id": "donor=2",
        "status": "running",
        "started_at_utc": "2026-08-09T00:02:00+00:00",
    }
    (tmp_path / "a.status.json").write_text(json.dumps(completed), encoding="utf-8")
    (tmp_path / "b.status.json").write_text(json.dumps(running), encoding="utf-8")

    result = checkpoint_runtimes(tmp_path)

    assert result[["method", "fold_id"]].to_dict(orient="records") == [
        {"method": "popv", "fold_id": "donor=1"}
    ]
    assert result.loc[0, "wall_seconds"] == 90.0
    assert result.loc[0, "reference_policy"] == "fold_train_only"
