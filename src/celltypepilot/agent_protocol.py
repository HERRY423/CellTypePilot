"""Stable decision protocol shared by the four default Agent-host tools."""

from __future__ import annotations

from typing import Any

AGENT_DECISION_SCHEMA = "celltypepilot.agent-decision.v1"

GLOBAL_FORBIDDEN_CLAIMS = [
    "do_not_describe_ranking_scores_as_calibrated_probabilities",
    "do_not_use_cell_state_to_rescue_identity",
    "do_not_promote_free_text_to_marker_evidence",
    "do_not_name_novel_cell_types_automatically",
    "do_not_claim_external_accuracy_without_independent_benchmark_evidence",
]


def agent_decision(
    *,
    operation: str,
    status: str,
    decision_scope: str,
    blockers: list[str] | None = None,
    warnings: list[str] | None = None,
    evidence_summary: dict[str, Any] | None = None,
    allowed_next_actions: list[str] | None = None,
    forbidden_claims: list[str] | None = None,
    artifact_paths: dict[str, str] | None = None,
    human_action_required: bool,
    claim_boundary: str,
    **compatibility_fields: Any,
) -> dict[str, Any]:
    """Return one machine-stable envelope without discarding legacy fields."""
    merged_forbidden = list(GLOBAL_FORBIDDEN_CLAIMS)
    for item in forbidden_claims or []:
        if item not in merged_forbidden:
            merged_forbidden.append(item)
    payload = {
        "schema_version": AGENT_DECISION_SCHEMA,
        "operation": operation,
        "status": status,
        "decision_scope": decision_scope,
        "blockers": list(blockers or []),
        "warnings": list(warnings or []),
        "evidence_summary": dict(evidence_summary or {}),
        "allowed_next_actions": list(allowed_next_actions or []),
        "forbidden_claims": merged_forbidden,
        "artifact_paths": dict(artifact_paths or {}),
        "human_action_required": bool(human_action_required),
        "claim_boundary": claim_boundary,
    }
    for key, value in compatibility_fields.items():
        if key not in payload:
            payload[key] = value
    return payload


def validate_agent_decision(payload: dict[str, Any]) -> None:
    """Fail closed when an MCP result drifts from the decision protocol."""
    required = {
        "schema_version",
        "operation",
        "status",
        "decision_scope",
        "blockers",
        "warnings",
        "evidence_summary",
        "allowed_next_actions",
        "forbidden_claims",
        "artifact_paths",
        "human_action_required",
        "claim_boundary",
    }
    missing = required - set(payload)
    if missing:
        raise ValueError(f"Agent decision payload missing fields: {sorted(missing)}")
    if payload["schema_version"] != AGENT_DECISION_SCHEMA:
        raise ValueError(f"schema_version must be {AGENT_DECISION_SCHEMA}")
    for field in ("blockers", "warnings", "allowed_next_actions", "forbidden_claims"):
        if not isinstance(payload[field], list):
            raise ValueError(f"{field} must be a list")
    if not isinstance(payload["artifact_paths"], dict):
        raise ValueError("artifact_paths must be an object")
    if not isinstance(payload["human_action_required"], bool):
        raise ValueError("human_action_required must be boolean")
