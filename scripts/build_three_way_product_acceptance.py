"""Build a fail-closed traditional vs bare-Agent vs Agent+CellTypePilot report.

The script never substitutes deterministic backend metrics for Agent-product
outcomes. Missing Agent runs or blinded human review are explicit blockers and
produce a rejected/not-claim-ready decision.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SCHEMA = "celltypepilot.three-way-product-acceptance.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_run(root: Path) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    manifest = json.loads((root / "benchmark_manifest.json").read_text(encoding="utf-8"))
    results = pd.read_csv(root / "benchmark_results.csv")
    status = pd.read_csv(root / "comparator_status.csv")
    return manifest, results, status


def _metric_rows(frame: pd.DataFrame, methods: set[str]) -> list[dict]:
    columns = [
        "method",
        "status",
        "n_cells",
        "accuracy",
        "macro_f1",
        "balanced_accuracy",
        "coverage",
        "abstain_rate",
        "selective_accuracy",
    ]
    selected = frame[frame["method"].isin(methods)][columns].copy()
    return selected.where(pd.notna(selected), None).to_dict(orient="records")


def _agent_observations(path: Path | None, arm: str) -> dict:
    if path is None or not path.is_file():
        return {
            "status": "not_run",
            "path": str(path) if path else None,
            "n_tasks": 0,
            "n_expert_accepted": 0,
            "correction_minutes_available": False,
        }
    frame = pd.read_csv(path, dtype=str).fillna("")
    required = {
        "arm",
        "task_id",
        "completion_status",
        "provenance_complete",
        "correction_minutes",
        "expert_accept",
    }
    missing = required - set(frame)
    if missing:
        return {"status": "invalid", "path": str(path), "missing_columns": sorted(missing)}
    scoped = frame[frame["arm"] == arm]
    if scoped.empty:
        return {"status": "not_run", "path": str(path), "n_tasks": 0}
    complete = scoped["completion_status"].str.casefold().eq("completed")
    accepted = scoped["expert_accept"].str.casefold().isin({"true", "1", "yes"})
    times = pd.to_numeric(scoped["correction_minutes"], errors="coerce")
    return {
        "status": "completed" if complete.all() else "incomplete",
        "path": str(path),
        "n_tasks": int(len(scoped)),
        "n_completed": int(complete.sum()),
        "n_expert_accepted": int(accepted.sum()),
        "expert_accept_rate": float(accepted.mean()),
        "correction_minutes_available": bool(times.notna().all()),
        "median_correction_minutes": float(times.median()) if times.notna().all() else None,
        "all_provenance_complete": bool(
            scoped["provenance_complete"].str.casefold().isin({"true", "1", "yes"}).all()
        ),
    }


def build_report(
    traditional_root: Path,
    plugin_root: Path,
    output_dir: Path,
    bare_agent_csv: Path | None = None,
    plugin_agent_csv: Path | None = None,
) -> dict:
    old_manifest, old_results, old_status = _load_run(traditional_root)
    new_manifest, new_results, new_status = _load_run(plugin_root)
    traditional_methods = {"celltypist", "singler", "popv"}
    traditional_complete = all(
        set(old_status.loc[old_status["method"] == method, "status"]) == {"completed"}
        for method in traditional_methods
    )
    plugin_complete = set(new_status["status"]) == {"completed"}
    same_contract = {
        "input_sha256": old_manifest.get("execution", {}).get("input_sha256")
        == new_manifest.get("execution", {}).get("input_sha256"),
        "assignments_sha256": old_manifest.get("assignments_sha256")
        == new_manifest.get("assignments_sha256"),
        "cluster_map_sha256": old_manifest.get("execution", {}).get("cluster_map_sha256")
        == new_manifest.get("execution", {}).get("cluster_map_sha256"),
        "label_map_sha256": old_manifest.get("execution", {}).get("label_map_sha256")
        == new_manifest.get("execution", {}).get("label_map_sha256"),
    }
    bare = _agent_observations(bare_agent_csv, "standalone_agent")
    plugin_agent = _agent_observations(plugin_agent_csv, "agent_plus_celltypepilot")
    plugin_metrics = _metric_rows(new_results, {"celltypepilot"})
    plugin_metric = plugin_metrics[0]

    gates = {
        "same_locked_data_contract": all(same_contract.values()),
        "traditional_software_completed": traditional_complete,
        "fixed_celltypepilot_backend_completed": plugin_complete,
        "standalone_agent_tasks_completed": bare.get("status") == "completed",
        "agent_plus_celltypepilot_tasks_completed": plugin_agent.get("status") == "completed",
        "blinded_expert_acceptance_available": (
            bare.get("status") == "completed" and plugin_agent.get("status") == "completed"
        ),
        "correction_time_available_for_both_agent_arms": bool(
            bare.get("correction_minutes_available")
            and plugin_agent.get("correction_minutes_available")
        ),
        "celltypepilot_not_single_lineage_only": False,
    }
    blockers = [name for name, passed in gates.items() if not passed]
    report = {
        "schema_version": SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "decision": "rejected_not_claim_ready" if blockers else "accepted_within_scope",
        "claim_ready": not blockers,
        "comparison_contract": same_contract,
        "arms": {
            "traditional_software": {
                "status": "completed" if traditional_complete else "incomplete",
                "methods": sorted(traditional_methods),
                "metrics": _metric_rows(old_results, traditional_methods),
            },
            "standalone_agent": bare,
            "agent_plus_celltypepilot": {
                "status": (
                    "completed"
                    if plugin_agent.get("status") == "completed"
                    else "backend_completed_agent_product_not_run"
                ),
                "backend_metrics": plugin_metrics,
                "agent_observations": plugin_agent,
                "observed_prediction_scope": {
                    "accepted_cell_fraction": plugin_metric.get("coverage"),
                    "accepted_identity": "Endothelial cell only",
                    "warning": "Backend recovery is single-lineage and does not establish product utility.",
                },
            },
        },
        "gates": gates,
        "blocking_findings": blockers,
        "product_conclusion": (
            "The gene-identity defect is repaired, but current CellTypePilot does not yet "
            "solve general lung annotation: coverage remains low and accepted predictions "
            "collapse to one broad lineage. No superiority over traditional software or "
            "standalone Agents may be claimed without paired Agent runs and blinded review."
        ),
        "claim_boundary": (
            "Traditional backend accuracy, CellTypePilot backend accuracy, Agent usability, "
            "expert acceptance, and correction time are separate endpoints. Missing arms are "
            "never imputed."
        ),
        "inputs": {
            "traditional_run": str(traditional_root),
            "traditional_results_sha256": _sha256(traditional_root / "benchmark_results.csv"),
            "plugin_run": str(plugin_root),
            "plugin_results_sha256": _sha256(plugin_root / "benchmark_results.csv"),
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "three_way_acceptance_report.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    template = pd.DataFrame(
        columns=[
            "arm",
            "task_id",
            "completion_status",
            "proposed_label",
            "uncertainty_language_valid",
            "evidence_citations_n",
            "provenance_complete",
            "correction_minutes",
            "expert_accept",
            "expert_notes",
        ]
    )
    template.to_csv(output_dir / "blinded_agent_review_template.csv", index=False)
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--traditional-run", type=Path, required=True)
    parser.add_argument("--plugin-run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bare-agent-csv", type=Path)
    parser.add_argument("--plugin-agent-csv", type=Path)
    args = parser.parse_args()
    report = build_report(
        args.traditional_run,
        args.plugin_run,
        args.output,
        args.bare_agent_csv,
        args.plugin_agent_csv,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
