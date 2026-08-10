"""Agent lifecycle discrimination and host-acceptance contracts."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from celltypepilot.agent_lifecycle import (
    AgentLifecycleError,
    assert_agent_can_discriminate,
    classify_checkpoint_status,
    classify_release_status,
    scan_checkpoint_dir,
)
from celltypepilot.host_acceptance import (
    _normalize_cli_help,
    discover_cli_commands,
    discover_host_surfaces,
    discover_mcp_tools,
    exercise_lifecycle_discrimination,
    exercise_timeout_cancel_resume_semantics,
    run_host_acceptance,
    write_lifecycle_fixtures,
)

REPO = Path(__file__).resolve().parents[1]
SCRATCH = REPO / "scratch" / "pytest_temp" / "agent_lifecycle"


def test_cli_help_normalization_preserves_styled_options():
    styled = "\x1b[1;36m-\x1b[0m\x1b[1;36m-input\x1b[0m"
    assert _normalize_cli_help(styled) == "--input"


def test_checkpoint_state_matrix():
    cases = [
        ({"status": "running"}, "running"),
        ({"status": "completed"}, "completed"),
        (
            {
                "status": "failed_or_unavailable",
                "error": "dependency_unavailable: not_installed",
            },
            "unavailable",
        ),
        (
            {
                "status": "failed_or_unavailable",
                "error": "exceeded 900 seconds; terminated by the locked execution limit",
            },
            "failed",
        ),
        (
            {
                "status": "failed_or_unavailable",
                "error": "cancelled by host KeyboardInterrupt",
            },
            "cancelled",
        ),
        (
            {
                "status": "completed",
                "resumed_from_checkpoint": True,
                "previous_status": "running",
            },
            "resumed",
        ),
    ]
    for payload, expected in cases:
        got = classify_checkpoint_status(payload)["agent_state"]
        assert got == expected, (payload, got, expected)


def test_claim_ready_vs_incomplete():
    assert classify_release_status({"status": "claim_ready"})["agent_state"] == "claim_ready"
    assert (
        classify_release_status({"status": "incomplete_not_claim_ready"})["agent_state"]
        == "incomplete_not_claim_ready"
    )


def test_discrimination_matrix_distinct():
    root = SCRATCH / "disc"
    if root.exists():
        import shutil

        shutil.rmtree(root)
    result = exercise_lifecycle_discrimination(root)
    assert result["passed"]
    assert_agent_can_discriminate(result["states_observed"].values())
    core = set(result["core_distinct"])
    assert core == {"running", "completed", "failed", "unavailable", "claim_ready"}


def test_timeout_cancel_resume():
    root = SCRATCH / "tcr"
    if root.exists():
        import shutil

        shutil.rmtree(root)
    result = exercise_timeout_cancel_resume_semantics(root)
    assert result["passed"]


def test_scan_checkpoint_dir_rollup():
    root = SCRATCH / "scan"
    if root.exists():
        import shutil

        shutil.rmtree(root)
    write_lifecycle_fixtures(root)
    scan = scan_checkpoint_dir(root / "checkpoints")
    assert scan["n_status_files"] >= 5
    assert "running" in scan["counts_by_agent_state"]
    assert scan["rollup_agent_state"] in {
        "running",
        "incomplete_not_claim_ready",
        "completed",
        "failed",
        "unavailable",
    }


def test_missing_states_fail_closed():
    with pytest.raises(AgentLifecycleError):
        assert_agent_can_discriminate(["completed", "running"])


def test_host_surface_discovery():
    surfaces = discover_host_surfaces(REPO)
    assert surfaces["claude_plugin_manifest"]
    assert surfaces["codex_plugin_manifest"]
    assert surfaces["skill_md"]
    assert surfaces["mcp_json"]


def test_cli_and_mcp_discovery():
    cli = discover_cli_commands()
    assert cli["required_commands"]["doctor"]
    assert cli["required_commands"]["inspect"]
    assert cli["required_commands"]["benchmark-run"]
    mcp = discover_mcp_tools()
    assert mcp["required_tools_present"]["tool_inspect_h5ad"]
    assert mcp["required_tools_present"]["tool_read_artifact_status"]


def test_host_acceptance_harness_passes():
    report = run_host_acceptance(skip_worktree=True)
    assert report["schema_version"].startswith("celltypepilot.host-acceptance")
    assert report["overall_status"] == "passed", json.dumps(report["checks"], indent=2)
    disc = report["agent_discrimination"]["states_observed"]
    assert disc["running"] == "running"
    assert disc["completed"] == "completed"
    assert disc["failed"] == "failed"
    assert disc["unavailable"] == "unavailable"
    assert disc["claim_ready"] == "claim_ready"


def test_mcp_lifecycle_tool():
    from celltypepilot.mcp_server import tool_agent_lifecycle_status

    root = SCRATCH / "mcp_life"
    if root.exists():
        import shutil

        shutil.rmtree(root)
    write_lifecycle_fixtures(root)
    payload = tool_agent_lifecycle_status(str(root))
    assert payload["prediction_mutation_allowed"] is False
    assert "checkpoint_scan" in payload
    assert payload["checkpoint_scan"]["n_status_files"] >= 1
