from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from scripts.build_three_way_product_acceptance import build_report


def _write_run(root: Path, methods: list[str]) -> None:
    root.mkdir(parents=True)
    (root / "benchmark_manifest.json").write_text(
        json.dumps(
            {
                "assignments_sha256": "assignments",
                "execution": {
                    "input_sha256": "input",
                    "cluster_map_sha256": "clusters",
                    "label_map_sha256": "labels",
                },
            }
        ),
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {
                "method": method,
                "status": "evaluated",
                "n_cells": 10,
                "accuracy": 0.8,
                "macro_f1": 0.7,
                "balanced_accuracy": 0.7,
                "coverage": 1.0 if method != "celltypepilot" else 0.2,
                "abstain_rate": 0.0 if method != "celltypepilot" else 0.8,
                "selective_accuracy": 0.8,
            }
            for method in methods
        ]
    ).to_csv(root / "benchmark_results.csv", index=False)
    pd.DataFrame(
        [
            {"method": method, "fold_id": f"donor={fold}", "status": "completed"}
            for method in methods
            for fold in (1, 2)
        ]
    ).to_csv(root / "comparator_status.csv", index=False)


def test_missing_agent_arms_rejects_product_claim(tmp_path: Path):
    traditional = tmp_path / "traditional"
    plugin = tmp_path / "plugin"
    output = tmp_path / "acceptance"
    _write_run(traditional, ["celltypist", "singler", "popv"])
    _write_run(plugin, ["celltypepilot"])

    report = build_report(traditional, plugin, output)

    assert report["decision"] == "rejected_not_claim_ready"
    assert report["comparison_contract"]["input_sha256"] is True
    assert report["arms"]["standalone_agent"]["status"] == "not_run"
    assert (
        report["arms"]["agent_plus_celltypepilot"]["status"]
        == "backend_completed_agent_product_not_run"
    )
    assert "standalone_agent_tasks_completed" in report["blocking_findings"]
    assert (output / "three_way_acceptance_report.json").is_file()
    assert (output / "blinded_agent_review_template.csv").is_file()
