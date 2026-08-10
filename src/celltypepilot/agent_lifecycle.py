"""Agent-facing lifecycle status vocabulary and classifiers.

Hosts (Codex, Claude Code, MCP) must not treat “command returned” as success.
They must discriminate durable states:

- ``running`` — work in progress; do not claim results
- ``completed`` — unit finished; still not a public claim
- ``failed`` — executable error after start
- ``unavailable`` — dependency/adapter/reference missing; negative result
- ``claim_ready`` — release gate green for public robustness claims
- ``incomplete_not_claim_ready`` — protocol artifacts incomplete
- ``cancelled`` — cooperative stop before completion
- ``resumed`` — completed via checkpoint resume (informational overlay)

This module is pure classification over JSON/files. It never mutates fold workspaces.
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

LIFECYCLE_SCHEMA = "celltypepilot.agent-lifecycle.v1"

# Canonical states agents must understand.
AGENT_STATES = (
    "running",
    "completed",
    "failed",
    "unavailable",
    "claim_ready",
    "incomplete_not_claim_ready",
    "cancelled",
    "resumed",
    "unknown",
)

# Checkpoint / comparator raw statuses written by benchmark_runner.
CHECKPOINT_RAW_STATUSES = (
    "running",
    "completed",
    "failed_or_unavailable",
    "failed",
    "error",
    "unreadable_checkpoint",
)

CLAIM_READY = "claim_ready"
INCOMPLETE_NOT_CLAIM_READY = "incomplete_not_claim_ready"

UNAVAILABLE_MARKERS = (
    "dependency_unavailable",
    "not_installed",
    "missing_r_packages",
    "missing_package",
    "adapter_unavailable",
    "azimuth:dependency_unavailable",
    "no executable adapter",
    "ready_not_run",
    "blocked_overlap_audit",
    "not_provided",
)

CANCEL_MARKERS = (
    "cancelled",
    "canceled",
    "keyboardinterrupt",
    "timeout",
    "timed out",
    "exceeded",
    "terminated by the locked execution limit",
)


class AgentLifecycleError(ValueError):
    """Invalid lifecycle payload."""


def _lower(value: Any) -> str:
    return str(value or "").strip().lower()


def classify_checkpoint_status(payload: dict[str, Any]) -> dict[str, Any]:
    """Map one checkpoint/*.status.json payload to an agent lifecycle state."""
    raw = _lower(payload.get("status"))
    error = str(payload.get("error") or payload.get("detail") or "")
    error_l = error.lower()
    resumed = bool(payload.get("resumed_from_checkpoint")) or _lower(
        payload.get("previous_status")
    ) in {"running", "failed_or_unavailable"}

    if raw == "running":
        agent_state = "running"
    elif raw == "completed":
        agent_state = (
            "resumed" if resumed and payload.get("resumed_from_checkpoint") else "completed"
        )
    elif raw in {"failed", "error", "unreadable_checkpoint"}:
        agent_state = "failed"
    elif raw == "failed_or_unavailable":
        if any(marker in error_l for marker in CANCEL_MARKERS):
            agent_state = "cancelled" if "cancel" in error_l else "failed"
            if any(m in error_l for m in ("timeout", "timed out", "exceeded", "execution limit")):
                agent_state = "failed"
        elif any(marker in error_l for marker in UNAVAILABLE_MARKERS):
            agent_state = "unavailable"
        else:
            agent_state = "failed"
    elif raw in {"cancelled", "canceled"}:
        agent_state = "cancelled"
    else:
        agent_state = "unknown"

    # Timeout failures are failed (not unavailable): the adapter was present but did not finish.
    if (
        any(m in error_l for m in ("timeout", "timed out", "execution limit"))
        and agent_state == "unavailable"
    ):
        agent_state = "failed"

    return {
        "schema_version": LIFECYCLE_SCHEMA,
        "agent_state": agent_state,
        "raw_status": payload.get("status"),
        "method": payload.get("method"),
        "fold_id": payload.get("fold_id"),
        "error": error or None,
        "resumed_from_checkpoint": bool(payload.get("resumed_from_checkpoint")),
        "claim_language": _claim_language(agent_state),
        "agent_must_not": _must_not(agent_state),
    }


def _claim_language(agent_state: str) -> str:
    if agent_state == "claim_ready":
        return "May support predeclared public robustness claims for this release only."
    if agent_state == "completed":
        return "Unit finished; not a multi-cohort claim-ready release."
    if agent_state == "running":
        return "In progress; do not report results as final."
    if agent_state == "unavailable":
        return "Negative result: comparator/dependency missing; do not impute predictions."
    if agent_state == "failed":
        return "Execution failed; retain as negative/incomplete, do not claim success."
    if agent_state == "cancelled":
        return "Stopped before completion; not success."
    if agent_state == "incomplete_not_claim_ready":
        return "Release gate blocked; no public robustness claim."
    if agent_state == "resumed":
        return "Completed after resume; still not automatically claim-ready."
    return "Unknown state; fail closed — do not claim."


def _must_not(agent_state: str) -> list[str]:
    base = ["invent_predictions", "upgrade_to_claim_ready_without_release_gate"]
    if agent_state == "running":
        return base + ["report_final_metrics", "mark_fold_done"]
    if agent_state == "unavailable":
        return base + ["treat_as_completed", "fill_missing_predictions"]
    if agent_state == "failed":
        return base + ["retry_silently_without_status", "hide_error"]
    if agent_state == "incomplete_not_claim_ready":
        return base + ["publish_as_claim_ready"]
    if agent_state == "cancelled":
        return base + ["report_as_completed"]
    return base


def classify_release_status(payload: dict[str, Any]) -> dict[str, Any]:
    """Map benchmark-release readiness to agent states."""
    raw = _lower(payload.get("status") or payload.get("release_status"))
    if raw == CLAIM_READY or raw == "claim-ready":
        agent_state = CLAIM_READY
    elif raw in {INCOMPLETE_NOT_CLAIM_READY, "incomplete", "not_claim_ready"}:
        agent_state = INCOMPLETE_NOT_CLAIM_READY
    elif raw in {"running", "in_progress"}:
        agent_state = "running"
    else:
        agent_state = "unknown"
    return {
        "schema_version": LIFECYCLE_SCHEMA,
        "agent_state": agent_state,
        "raw_status": payload.get("status") or payload.get("release_status"),
        "claim_language": _claim_language(agent_state),
        "agent_must_not": _must_not(agent_state),
        "blocking_findings": payload.get("blocking_findings") or payload.get("n_blocking_findings"),
    }


def classify_doctor_capability(capabilities: dict[str, Any]) -> dict[str, Any]:
    """Map doctor capability strings to available/unavailable for agent routing."""
    classified = {}
    for name, value in capabilities.items():
        text = _lower(value)
        if text in {"available", "ok", "full", "yes", "true", "installed"}:
            state = "completed"  # capability ready (not a benchmark unit)
            agent = "available"
        elif text in {"missing", "unavailable", "not_installed", "false", "no"}:
            state = "unavailable"
            agent = "unavailable"
        elif text in {"degraded", "partial", "optional"}:
            state = "incomplete_not_claim_ready"
            agent = "degraded"
        else:
            state = "unknown"
            agent = "unknown"
        classified[name] = {
            "raw": value,
            "agent_capability": agent,
            "lifecycle_analogue": state,
        }
    return {
        "schema_version": LIFECYCLE_SCHEMA,
        "capabilities": classified,
        "note": (
            "Doctor capabilities are environment routing signals, not annotation accuracy "
            "or claim-ready evidence."
        ),
    }


def scan_checkpoint_dir(checkpoint_dir: str | Path) -> dict[str, Any]:
    """Classify all *.status.json files for an Agent host."""
    root = Path(checkpoint_dir)
    records = []
    counts: dict[str, int] = {state: 0 for state in AGENT_STATES}
    if root.is_dir():
        for path in sorted(root.glob("*.status.json")):
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                payload = {"status": "unreadable_checkpoint", "error": "unreadable"}
            if not isinstance(payload, dict):
                payload = {"status": "unreadable_checkpoint", "error": "not_object"}
            classified = classify_checkpoint_status(payload)
            classified["status_path"] = str(path)
            records.append(classified)
            counts[classified["agent_state"]] = counts.get(classified["agent_state"], 0) + 1

    # Aggregate agent view
    if counts.get("running", 0):
        rollup = "running"
    elif counts.get("failed", 0) and not counts.get("completed", 0):
        rollup = "failed"
    elif counts.get("unavailable", 0) and not counts.get("completed", 0):
        rollup = "unavailable"
    elif counts.get("completed", 0) or counts.get("resumed", 0):
        if any(counts.get(s, 0) for s in ("failed", "unavailable", "cancelled", "running")):
            rollup = "incomplete_not_claim_ready"
        else:
            rollup = "completed"
    else:
        rollup = "unknown"

    return {
        "schema_version": LIFECYCLE_SCHEMA,
        "checkpoint_dir": str(root),
        "n_status_files": len(records),
        "counts_by_agent_state": counts,
        "rollup_agent_state": rollup,
        "records": records,
        "claim_language": _claim_language(rollup),
        "agent_must_not": _must_not(rollup),
        "discrimination_required": list(AGENT_STATES),
    }


def build_agent_status_report(
    *,
    doctor: dict[str, Any] | None = None,
    inspect_report: dict[str, Any] | None = None,
    checkpoint_dir: str | Path | None = None,
    release_manifest: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compose a single Agent-facing status document for host acceptance."""
    sections: dict[str, Any] = {
        "schema_version": LIFECYCLE_SCHEMA,
        "purpose": (
            "Discriminate running/completed/failed/unavailable/claim_ready. "
            "Presence of a CLI exit code is not sufficient."
        ),
        "states": list(AGENT_STATES),
    }
    if doctor is not None:
        caps = doctor.get("capabilities") or {}
        sections["doctor"] = {
            "python_ok": doctor.get("python_ok"),
            "capability_classification": classify_doctor_capability(caps),
            "mcp_status": doctor.get("mcp_status"),
        }
    if inspect_report is not None:
        sections["inspect"] = {
            "species": inspect_report.get("species"),
            "annotation_supported": inspect_report.get("annotation_supported")
            or inspect_report.get("species_supported"),
            "cluster_keys": inspect_report.get("cluster_keys") or inspect_report.get("clusters"),
            "note": "Inspection is routing metadata, not claim-ready evidence.",
        }
    if checkpoint_dir is not None:
        sections["checkpoints"] = scan_checkpoint_dir(checkpoint_dir)
    if release_manifest is not None:
        readiness = release_manifest.get("readiness") or release_manifest
        if "status" not in readiness and release_manifest.get("status"):
            readiness = release_manifest
        sections["release"] = classify_release_status(readiness)
    return sections


def doctor_report_to_dict(report: Any) -> dict[str, Any]:
    """Serialize DoctorReport dataclass-like object to JSON dict."""
    deps = []
    for item in getattr(report, "dependencies", []) or []:
        deps.append(
            {
                "name": item.name,
                "installed": item.installed,
                "version": item.version,
                "required": item.required,
                "note": item.note,
            }
        )
    optional = []
    for item in getattr(report, "optional_deps", []) or []:
        optional.append(
            {
                "name": item.name,
                "installed": item.installed,
                "version": item.version,
                "required": item.required,
                "note": item.note,
            }
        )
    return {
        "schema_version": "celltypepilot.doctor.v1",
        "python_version": report.python_version,
        "python_ok": report.python_ok,
        "dependencies": deps,
        "optional_deps": optional,
        "capabilities": dict(report.capabilities or {}),
        "warnings": list(report.warnings or []),
        "mcp_status": dict(report.mcp_status or {}),
        "lifecycle": classify_doctor_capability(dict(report.capabilities or {})),
    }


def assert_agent_can_discriminate(states: Iterable[str]) -> None:
    """Fail closed if a host test matrix is missing required states."""
    required = {"running", "completed", "failed", "unavailable", "claim_ready"}
    have = set(states)
    missing = sorted(required - have)
    if missing:
        raise AgentLifecycleError(
            f"Host acceptance matrix missing required agent states: {missing}"
        )
