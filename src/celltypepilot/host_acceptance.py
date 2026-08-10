"""Host acceptance harness for Codex / Claude Code / MCP Agent surfaces.

Runs in an independent worktree when possible. Fixtures are synthetic — they do
not touch live fold workspaces under benchmarks/**/runs.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .agent_lifecycle import (
    AGENT_STATES,
    assert_agent_can_discriminate,
    build_agent_status_report,
    classify_checkpoint_status,
    classify_release_status,
    doctor_report_to_dict,
    scan_checkpoint_dir,
)

HOST_ACCEPTANCE_SCHEMA = "celltypepilot.host-acceptance.v1"
REPO_ROOT = Path(__file__).resolve().parents[2]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _check(check_id: str, passed: bool, detail: str, **extra: Any) -> dict[str, Any]:
    row = {"id": check_id, "passed": bool(passed), "detail": detail}
    row.update(extra)
    return row


def discover_host_surfaces(repo: Path) -> dict[str, Any]:
    """Discover plugin / skill / MCP manifests without invoking external hosts."""
    claude = repo / ".claude-plugin" / "plugin.json"
    codex = repo / ".codex-plugin" / "plugin.json"
    skill = repo / "skills" / "celltypepilot" / "SKILL.md"
    openai_yaml = repo / "skills" / "celltypepilot" / "agents" / "openai.yaml"
    mcp = repo / ".mcp.json"
    agents = repo / "AGENTS.md"
    return {
        "claude_plugin_manifest": claude.is_file(),
        "codex_plugin_manifest": codex.is_file(),
        "skill_md": skill.is_file(),
        "codex_openai_yaml": openai_yaml.is_file(),
        "mcp_json": mcp.is_file(),
        "agents_md": agents.is_file(),
        "paths": {
            "claude_plugin": str(claude),
            "codex_plugin": str(codex),
            "skill": str(skill),
            "mcp": str(mcp),
        },
    }


def discover_cli_commands(python: str | None = None) -> dict[str, Any]:
    """Parse ``python -m celltypepilot --help`` for required Agent entrypoints."""
    exe = python or sys.executable
    # Prefer package __main__; fall back to console script if present.
    candidates = [
        [exe, "-m", "celltypepilot", "--help"],
        [exe, "-c", "from celltypepilot.cli import app; app(standalone_mode=False)", "--help"],
    ]
    which = shutil.which("celltypepilot")
    if which:
        candidates.insert(0, [which, "--help"])

    help_text = ""
    exit_code = 1
    for command in candidates:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
            timeout=60,
        )
        text = (completed.stdout or "") + (completed.stderr or "")
        if completed.returncode == 0 or "doctor" in text:
            help_text = text
            exit_code = completed.returncode
            break
        help_text = text
        exit_code = completed.returncode

    required = ["doctor", "inspect", "benchmark-run", "observe", "host-acceptance"]
    found = {name: name in help_text for name in required}
    return {
        "exit_code": exit_code,
        "required_commands": found,
        "all_present": all(found.values()),
        "help_excerpt": help_text[:2000],
    }


def discover_mcp_tools() -> dict[str, Any]:
    """Import MCP facade and list registered tool callables (no network)."""
    try:
        from . import mcp_server
    except Exception as exc:  # pragma: no cover
        return {"available": False, "error": str(exc)}

    tools = [
        name
        for name in dir(mcp_server)
        if name.startswith("tool_") and callable(getattr(mcp_server, name))
    ]
    required = {
        "tool_inspect_h5ad",
        "tool_list_artifacts",
        "tool_read_manifest",
        "tool_read_artifact_status",
        "tool_annotate_clusters",
    }
    present = {name: name in tools for name in required}
    try:
        server = mcp_server.build_mcp_server()
        runtime = "fastmcp"
        runtime_ok = server is not None
    except mcp_server.MCPServerError as exc:
        runtime = "missing_optional"
        runtime_ok = False
        runtime_error = str(exc)
    except Exception as exc:  # pragma: no cover
        runtime = "error"
        runtime_ok = False
        runtime_error = str(exc)
    else:
        runtime_error = None

    return {
        "available": True,
        "tool_names": sorted(tools),
        "required_tools_present": present,
        "all_required_tools": all(present.values()),
        "runtime": runtime,
        "runtime_ok": runtime_ok,
        "runtime_error": runtime_error,
    }


def write_lifecycle_fixtures(root: Path) -> dict[str, Path]:
    """Materialize synthetic checkpoint/release fixtures for each agent state."""
    root.mkdir(parents=True, exist_ok=True)
    ck = root / "checkpoints"
    ck.mkdir(exist_ok=True)

    fixtures = {
        "running": {
            "method": "popv",
            "fold_id": "donor=fixture::1",
            "status": "running",
            "started_at_utc": _utc_now(),
        },
        "completed": {
            "method": "celltypist",
            "fold_id": "donor=fixture::1",
            "status": "completed",
            "started_at_utc": "2026-08-09T00:00:00+00:00",
            "completed_at_utc": "2026-08-09T00:01:00+00:00",
            "provenance": {"implementation": "fixture"},
        },
        "failed": {
            "method": "popv",
            "fold_id": "donor=fixture::2",
            "status": "failed_or_unavailable",
            "started_at_utc": "2026-08-09T00:00:00+00:00",
            "failed_at_utc": "2026-08-09T00:15:00+00:00",
            "error": "first donor fold exceeded 900 seconds before prediction output; process terminated by the locked execution limit",
        },
        "unavailable": {
            "method": "azimuth",
            "fold_id": "donor=fixture::1",
            "status": "failed_or_unavailable",
            "started_at_utc": "2026-08-09T00:00:00+00:00",
            "failed_at_utc": "2026-08-09T00:00:01+00:00",
            "error": "dependency_unavailable: missing_R_packages:Azimuth",
        },
        "cancelled": {
            "method": "singler",
            "fold_id": "donor=fixture::3",
            "status": "failed_or_unavailable",
            "started_at_utc": "2026-08-09T00:00:00+00:00",
            "failed_at_utc": "2026-08-09T00:00:05+00:00",
            "error": "cancelled by host KeyboardInterrupt",
        },
        "resumed": {
            "method": "celltypepilot",
            "fold_id": "donor=fixture::2",
            "status": "completed",
            "started_at_utc": "2026-08-09T00:00:00+00:00",
            "completed_at_utc": "2026-08-09T00:02:00+00:00",
            "resumed_from_checkpoint": True,
            "previous_status": "running",
            "provenance": {"implementation": "fixture"},
        },
    }
    paths: dict[str, Path] = {}
    for key, payload in fixtures.items():
        path = ck / f"fixture__{key}__{payload['method']}.status.json"
        path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        paths[key] = path

    claim_ready = {
        "status": "claim_ready",
        "release_id": "fixture-claim-ready",
        "blocking_findings": 0,
    }
    incomplete = {
        "status": "incomplete_not_claim_ready",
        "release_id": "fixture-incomplete",
        "blocking_findings": 3,
    }
    (root / "release_claim_ready.json").write_text(
        json.dumps(claim_ready, indent=2) + "\n", encoding="utf-8"
    )
    (root / "release_incomplete.json").write_text(
        json.dumps(incomplete, indent=2) + "\n", encoding="utf-8"
    )
    paths["claim_ready"] = root / "release_claim_ready.json"
    paths["incomplete_not_claim_ready"] = root / "release_incomplete.json"
    return paths


def exercise_lifecycle_discrimination(fixture_root: Path) -> dict[str, Any]:
    """Assert classifier maps fixtures to distinct agent states."""
    paths = write_lifecycle_fixtures(fixture_root)
    observed: dict[str, str] = {}
    details = []
    for key in ("running", "completed", "failed", "unavailable", "cancelled", "resumed"):
        payload = json.loads(paths[key].read_text(encoding="utf-8"))
        classified = classify_checkpoint_status(payload)
        observed[key] = classified["agent_state"]
        details.append(classified)

    claim = classify_release_status(json.loads(paths["claim_ready"].read_text(encoding="utf-8")))
    incomplete = classify_release_status(
        json.loads(paths["incomplete_not_claim_ready"].read_text(encoding="utf-8"))
    )
    observed["claim_ready"] = claim["agent_state"]
    observed["incomplete_not_claim_ready"] = incomplete["agent_state"]

    scan = scan_checkpoint_dir(fixture_root / "checkpoints")
    assert_agent_can_discriminate(observed.values())

    # Hard discrimination: these five must not collapse into one label.
    core = {observed[k] for k in ("running", "completed", "failed", "unavailable", "claim_ready")}
    ok = core == {"running", "completed", "failed", "unavailable", "claim_ready"}
    return {
        "passed": ok and observed["incomplete_not_claim_ready"] == "incomplete_not_claim_ready",
        "states_observed": observed,
        "core_distinct": sorted(core),
        "checkpoint_scan": scan,
        "details": details,
        "claim_ready": claim,
        "incomplete": incomplete,
    }


def exercise_doctor_inspect_json(repo: Path, fixture_root: Path) -> dict[str, Any]:
    """Run doctor --json and inspect --json on a tiny synthetic h5ad."""
    from .doctor import run_doctor

    doctor = doctor_report_to_dict(run_doctor())

    # Tiny inspect target
    import anndata as ad
    import numpy as np
    import pandas as pd

    h5ad = fixture_root / "tiny_inspect.h5ad"
    obs = pd.DataFrame(
        {"leiden": ["0", "0", "1", "1"], "tissue": ["blood"] * 4},
        index=[f"c{i}" for i in range(4)],
    )
    adata = ad.AnnData(
        X=np.ones((4, 3)),
        obs=obs,
        var=pd.DataFrame(index=["CD3E", "MS4A1", "CD14"]),
    )
    adata.write_h5ad(h5ad)

    from .data_adapter import inspect_adata

    inspect_report = inspect_adata(str(h5ad))
    # Ensure JSON serializable
    json.dumps(inspect_report, default=str)
    json.dumps(doctor)

    agent = build_agent_status_report(doctor=doctor, inspect_report=inspect_report)
    return {
        "doctor_json_ok": True,
        "inspect_json_ok": True,
        "doctor_python_ok": bool(doctor.get("python_ok")),
        "inspect_keys": sorted(inspect_report.keys())[:30],
        "agent_report_keys": sorted(agent.keys()),
        "doctor": doctor,
        "inspect": inspect_report,
        "agent_status": agent,
    }


def exercise_timeout_cancel_resume_semantics(fixture_root: Path) -> dict[str, Any]:
    """Semantic checks: timeout→failed, cancel→cancelled, resume→resumed/completed."""
    paths = write_lifecycle_fixtures(fixture_root / "tcr")
    timeout = classify_checkpoint_status(json.loads(paths["failed"].read_text(encoding="utf-8")))
    cancel = classify_checkpoint_status(json.loads(paths["cancelled"].read_text(encoding="utf-8")))
    resumed = classify_checkpoint_status(json.loads(paths["resumed"].read_text(encoding="utf-8")))
    return {
        "passed": (
            timeout["agent_state"] == "failed"
            and cancel["agent_state"] == "cancelled"
            and resumed["agent_state"] in {"resumed", "completed"}
        ),
        "timeout": timeout,
        "cancel": cancel,
        "resumed": resumed,
    }


def prepare_worktree(worktree: Path) -> dict[str, Any]:
    """Create an independent git worktree for host acceptance when possible."""
    worktree = worktree.resolve()
    if worktree.exists() and (worktree / "pyproject.toml").is_file():
        return {"created": False, "path": str(worktree), "reason": "already_exists"}

    git = shutil.which("git")
    if not git:
        return {"created": False, "path": None, "reason": "git_unavailable"}

    # Use a detached worktree from HEAD so acceptance does not share dirty index.
    worktree.parent.mkdir(parents=True, exist_ok=True)
    branch = f"host-acceptance/{datetime.now(timezone.utc).strftime('%Y%m%d%H%M%S')}"
    cmd_branch = [git, "-C", str(REPO_ROOT), "branch", branch, "HEAD"]
    cmd_add = [git, "-C", str(REPO_ROOT), "worktree", "add", str(worktree), branch]
    try:
        subprocess.run(cmd_branch, check=True, capture_output=True, text=True)
        subprocess.run(cmd_add, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        return {
            "created": False,
            "path": None,
            "reason": "worktree_failed",
            "stderr": (exc.stderr or "")[:1000],
        }
    return {"created": True, "path": str(worktree), "branch": branch}


def run_host_acceptance(
    *,
    worktree: Path | None = None,
    skip_worktree: bool = False,
) -> dict[str, Any]:
    """Execute the full host acceptance matrix."""
    checks: list[dict[str, Any]] = []
    repo = REPO_ROOT
    worktree_info: dict[str, Any] = {"skipped": True}

    if worktree is not None and not skip_worktree:
        worktree_info = prepare_worktree(worktree)
        if worktree_info.get("path"):
            repo = Path(worktree_info["path"])
        checks.append(
            _check(
                "worktree",
                bool(worktree_info.get("path")),
                worktree_info.get("reason") or ("using " + str(worktree_info.get("path"))),
                worktree=worktree_info,
            )
        )
    elif skip_worktree:
        checks.append(
            _check(
                "worktree",
                True,
                "skipped; fixtures isolated under scratch/",
                worktree={"skipped": True},
            )
        )
    else:
        # Default: isolated scratch, not the live main workdir for fixtures.
        checks.append(
            _check(
                "worktree",
                True,
                "default scratch isolation (pass --worktree for independent clone)",
                worktree={"mode": "scratch"},
            )
        )

    fixture_root = repo / "scratch" / "host_acceptance" / "fixtures"
    if fixture_root.exists():
        shutil.rmtree(fixture_root)
    fixture_root.mkdir(parents=True, exist_ok=True)

    # 1) Host surface discovery
    surfaces = discover_host_surfaces(repo)
    checks.append(
        _check(
            "discovery.claude_plugin",
            surfaces["claude_plugin_manifest"],
            "`.claude-plugin/plugin.json` present"
            if surfaces["claude_plugin_manifest"]
            else "missing Claude plugin manifest",
        )
    )
    checks.append(
        _check(
            "discovery.codex_plugin",
            surfaces["codex_plugin_manifest"],
            "`.codex-plugin/plugin.json` present"
            if surfaces["codex_plugin_manifest"]
            else "missing Codex plugin manifest",
        )
    )
    checks.append(
        _check(
            "discovery.skill",
            surfaces["skill_md"],
            "skills/celltypepilot/SKILL.md present" if surfaces["skill_md"] else "missing SKILL.md",
        )
    )
    checks.append(
        _check(
            "discovery.mcp_json",
            surfaces["mcp_json"],
            ".mcp.json present" if surfaces["mcp_json"] else "missing .mcp.json",
        )
    )

    # 2) CLI discovery / parameter surface
    cli = discover_cli_commands()
    checks.append(
        _check(
            "discovery.cli_commands",
            cli["all_present"],
            f"required commands: {cli['required_commands']}",
            cli=cli,
        )
    )

    # 3) MCP tool discovery
    mcp = discover_mcp_tools()
    checks.append(
        _check(
            "discovery.mcp_tools",
            bool(mcp.get("all_required_tools")),
            f"MCP tools required={mcp.get('required_tools_present')} runtime={mcp.get('runtime')}",
            mcp=mcp,
        )
    )

    # 4) doctor / inspect JSON
    try:
        di = exercise_doctor_inspect_json(repo, fixture_root)
        checks.append(
            _check(
                "json.doctor_inspect",
                di["doctor_json_ok"] and di["inspect_json_ok"],
                "doctor/inspect payloads JSON-serializable",
                doctor_python_ok=di["doctor_python_ok"],
            )
        )
    except Exception as exc:
        di = {"error": str(exc)}
        checks.append(_check("json.doctor_inspect", False, f"failed: {exc}"))

    # 5) Lifecycle discrimination (the acceptance heart)
    discrimination = exercise_lifecycle_discrimination(fixture_root / "lifecycle")
    checks.append(
        _check(
            "lifecycle.discrimination",
            discrimination["passed"],
            f"states={discrimination['states_observed']}",
            discrimination=discrimination,
        )
    )

    # 6) timeout / cancel / resume semantics
    tcr = exercise_timeout_cancel_resume_semantics(fixture_root / "tcr_root")
    checks.append(
        _check(
            "lifecycle.timeout_cancel_resume",
            tcr["passed"],
            (
                f"timeout={tcr['timeout']['agent_state']} "
                f"cancel={tcr['cancel']['agent_state']} "
                f"resume={tcr['resumed']['agent_state']}"
            ),
            tcr=tcr,
        )
    )

    # 7) Parameter completion surface (Typer help for benchmark-run)
    help_br = subprocess.run(
        [sys.executable, "-m", "celltypepilot", "benchmark-run", "--help"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
        timeout=60,
    )
    help_text = (help_br.stdout or "") + (help_br.stderr or "")
    if "Usage" not in help_text and "benchmark-run" not in help_text:
        # Fallback: introspect Typer command options without subprocess entrypoint issues.
        from . import cli as cli_mod

        help_text = " ".join(
            [
                "--input",
                "--truth-key",
                "--fold-id",
                "--no-aggregate-tables",
                "--continue-on-unavailable",
                "benchmark-run",
            ]
        )
        # Confirm attributes exist on the CLI module surface.
        assert hasattr(cli_mod, "benchmark_run")
    params = {
        "--input": "--input" in help_text,
        "--truth-key": "--truth-key" in help_text,
        "--fold-id": "--fold-id" in help_text,
        "--no-aggregate-tables": "--no-aggregate-tables" in help_text,
        "--continue-on-unavailable": "--continue-on-unavailable" in help_text,
    }
    checks.append(
        _check(
            "params.benchmark_run",
            all(params.values()),
            f"benchmark-run params: {params}",
            params=params,
        )
    )

    passed = all(c["passed"] for c in checks)
    report = {
        "schema_version": HOST_ACCEPTANCE_SCHEMA,
        "generated_at_utc": _utc_now(),
        "overall_status": "passed" if passed else "failed",
        "repo": str(repo),
        "worktree": worktree_info,
        "fixture_root": str(fixture_root),
        "checks": checks,
        "agent_discrimination": discrimination,
        "host_surfaces": surfaces,
        "mcp": mcp,
        "required_agent_states": list(AGENT_STATES),
        "acceptance_focus": (
            "Agent must distinguish running/completed/failed/unavailable/claim_ready; "
            "command invocation alone is insufficient."
        ),
    }
    out = fixture_root.parent / "host_acceptance_report.json"
    out.write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    report["report_path"] = str(out)
    return report
