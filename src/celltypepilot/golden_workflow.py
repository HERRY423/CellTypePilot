"""Four-step Agent-host workflow for governed cell annotation.

The public plugin surface intentionally exposes only four stateful operations:
prepare, annotate, review, and finalize. Lower-level scoring and governance
functions remain implementation details or advanced maintainer interfaces.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .agent_protocol import agent_decision

PLAN_SCHEMA = "celltypepilot.agent-annotation-plan.v1"
STATE_SCHEMA = "celltypepilot.agent-annotation-state.v1"
REVIEW_SCHEMA = "celltypepilot.agent-review-queue.v1"
PLAN_FILE = "annotation_plan.json"
STATE_FILE = "annotation_workflow_state.json"
REVIEW_FILE = "review_queue.json"


class GoldenWorkflowError(ValueError):
    """Raised when an Agent attempts to skip or mutate a locked workflow step."""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _canonical_hash(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode(
            "utf-8"
        )
    ).hexdigest()


def _write_json(path: Path, payload: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


def _load_plan(plan_path: str | Path) -> dict:
    path = Path(plan_path)
    if not path.is_file():
        raise GoldenWorkflowError(f"Annotation plan not found: {path}")
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != PLAN_SCHEMA:
        raise GoldenWorkflowError(f"Plan schema must be {PLAN_SCHEMA}")
    expected = plan.get("plan_sha256")
    unsigned = {key: value for key, value in plan.items() if key != "plan_sha256"}
    if not expected or _canonical_hash(unsigned) != expected:
        raise GoldenWorkflowError("Annotation plan hash mismatch; rerun prepare_annotation")
    if plan.get("status") != "ready":
        raise GoldenWorkflowError("Annotation plan is blocked and cannot be executed")
    return plan


def prepare_annotation(
    input_path: str,
    output_dir: str,
    *,
    cluster_key: str | None = None,
    species: str | None = None,
    tissue: str | None = None,
    embedding_key: str | None = None,
    packs: list[str] | None = None,
    marker_evidence_policy: str = "database",
    context_file_path: str | None = None,
    custom_markers_path: str | None = None,
) -> dict[str, Any]:
    """Step 1: inspect data and write an immutable, executable annotation plan."""
    from .data_adapter import compute_data_hash, inspect_adata
    from .mcp_server import tool_evidence_coverage
    from .pack_manager import list_installed_packs

    source = Path(input_path).resolve()
    root = Path(output_dir).resolve()
    if not source.is_file():
        raise GoldenWorkflowError(f"Input h5ad not found: {source}")
    inspection = inspect_adata(source, cluster_key, embedding_key)
    detected_species = species or inspection.get("species") or inspection.get("detected_species")
    detected_tissue = tissue or inspection.get("tissue") or "general"
    cluster_candidates = inspection.get("cluster_keys") or inspection.get("cluster_columns") or []
    embedding_candidates = inspection.get("embedding_keys") or []
    selected_cluster = cluster_key or (
        cluster_candidates[0] if len(cluster_candidates) == 1 else None
    )
    selected_embedding = embedding_key or (
        embedding_candidates[0] if len(embedding_candidates) == 1 else None
    )

    blockers: list[str] = []
    warnings: list[str] = []
    if detected_species not in {"human", "mouse"}:
        blockers.append(
            f"species {detected_species!r} is not supported by the runtime annotation atlas"
        )
    if not selected_cluster:
        blockers.append("cluster_key is required and could not be selected unambiguously")
    elif cluster_candidates and selected_cluster not in cluster_candidates:
        blockers.append(f"cluster_key {selected_cluster!r} was not reported by inspection")
    if context_file_path and not Path(context_file_path).is_file():
        blockers.append(f"context file not found: {context_file_path}")
    if custom_markers_path and not Path(custom_markers_path).is_file():
        blockers.append(f"custom marker file not found: {custom_markers_path}")

    installed = {item["name"]: item for item in list_installed_packs()}
    for name in packs or []:
        if name not in installed:
            blockers.append(f"evidence pack {name!r} is not installed")
        elif installed[name].get("trust") == "hypothesis":
            warnings.append(
                f"pack {name!r} is hypothesis-trust and cannot silently produce an accepted identity"
            )

    coverage = None
    if not blockers and detected_species and detected_tissue:
        coverage = tool_evidence_coverage(
            str(source),
            str(detected_species),
            str(detected_tissue),
            packs=list(packs or []),
            evidence_policy=marker_evidence_policy,
        )
        if coverage.get("gene_identity", {}).get("marker_overlap_after", 0) == 0:
            blockers.append("no runtime marker genes are addressable after identity normalization")

    plan = {
        "schema_version": PLAN_SCHEMA,
        "created_at": _utc_now(),
        "status": "blocked" if blockers else "ready",
        "input_path": str(source),
        "input_sha256": compute_data_hash(source),
        "output_dir": str(root),
        "cluster_key": selected_cluster,
        "species": detected_species,
        "tissue": detected_tissue,
        "embedding_key": selected_embedding,
        "packs": list(packs or []),
        "marker_evidence_policy": marker_evidence_policy,
        "context_file_path": str(Path(context_file_path).resolve()) if context_file_path else None,
        "custom_markers_path": (
            str(Path(custom_markers_path).resolve()) if custom_markers_path else None
        ),
        "blockers": blockers,
        "warnings": warnings,
        "inspection": inspection,
        "evidence_coverage": coverage,
        "claim_boundary": (
            "Preparation establishes executable scope and marker addressability; it does not "
            "make an annotation or biological accuracy claim."
        ),
    }
    plan["plan_sha256"] = _canonical_hash(plan)
    path = _write_json(root / PLAN_FILE, plan)
    _write_json(
        root / STATE_FILE,
        {
            "schema_version": STATE_SCHEMA,
            "phase": "prepared" if not blockers else "blocked",
            "updated_at": _utc_now(),
            "plan_sha256": plan["plan_sha256"],
            "next_action": "annotate_from_plan" if not blockers else "resolve_blockers",
        },
    )
    allowed = ["annotate_from_plan"] if not blockers else ["resolve_prepare_blockers"]
    return agent_decision(
        operation="prepare_annotation",
        status=plan["status"],
        decision_scope="input_and_annotation_plan_readiness",
        blockers=blockers,
        warnings=warnings,
        evidence_summary={
            "species": detected_species,
            "tissue": detected_tissue,
            "cluster_key": selected_cluster,
            "marker_overlap_after": (
                coverage.get("gene_identity", {}).get("marker_overlap_after")
                if coverage
                else None
            ),
        },
        allowed_next_actions=allowed,
        forbidden_claims=["do_not_treat_preparation_as_an_annotation_result"],
        artifact_paths={"annotation_plan": str(path), "workflow_state": str(root / STATE_FILE)},
        human_action_required=bool(blockers),
        claim_boundary=plan["claim_boundary"],
        plan_path=str(path),
        plan=plan,
    )


def annotate_from_plan(plan_path: str) -> dict[str, Any]:
    """Step 2: execute only the parameters locked by prepare_annotation."""
    from .data_adapter import compute_data_hash
    from .orchestrator import run_annotation_pipeline

    plan = _load_plan(plan_path)
    if compute_data_hash(plan["input_path"]) != plan["input_sha256"]:
        raise GoldenWorkflowError("Input h5ad changed after preparation")
    result = run_annotation_pipeline(
        input_path=plan["input_path"],
        cluster_key=plan["cluster_key"],
        output_dir=plan["output_dir"],
        species=plan["species"],
        tissue=plan["tissue"],
        embedding_key=plan["embedding_key"],
        marker_evidence_policy=plan["marker_evidence_policy"],
        context_file_path=plan.get("context_file_path"),
        custom_markers_path=plan.get("custom_markers_path"),
        packs=plan.get("packs"),
    )
    state = {
        "schema_version": STATE_SCHEMA,
        "phase": "annotated_draft",
        "updated_at": _utc_now(),
        "plan_sha256": plan["plan_sha256"],
        "next_action": "review_uncertain_clusters",
        "critic_summary": result.get("critic_summary"),
        "validation_scope": result.get("validation_scope"),
    }
    _write_json(Path(plan["output_dir"]) / STATE_FILE, state)
    paths = {key: str(value) for key, value in result.get("paths", {}).items()}
    critic_summary = result.get("critic_summary") or {}
    return agent_decision(
        operation="annotate_from_plan",
        status="annotated_draft",
        decision_scope="draft_identity_candidates_with_fail_closed_abstention",
        evidence_summary={
            "critic_summary": critic_summary,
            "validation_scope": result.get("validation_scope"),
        },
        allowed_next_actions=["review_uncertain_clusters"],
        forbidden_claims=[
            "do_not_treat_draft_annotations_as_human_approved",
            "do_not_replace_unknown_with_the_top_candidate_automatically",
        ],
        artifact_paths={**paths, "workflow_state": str(Path(plan["output_dir"]) / STATE_FILE)},
        human_action_required=True,
        claim_boundary=(
            "Annotation produces reviewable draft identities and explicit abstentions; human "
            "review remains required and the run is not external accuracy evidence."
        ),
        output_dir=plan["output_dir"],
        critic_summary=critic_summary,
        validation_scope=result.get("validation_scope"),
        paths=paths,
        next_action="review_uncertain_clusters",
    )


def review_uncertain_clusters(output_dir: str) -> dict[str, Any]:
    """Step 3: produce one bounded queue; never invent or apply a label."""
    import pandas as pd

    root = Path(output_dir).resolve()
    from .agent_evidence import load_agent_evidence_indexes

    contrast_index, gap_index = load_agent_evidence_indexes(root)
    evidence_path = root / "evidence_table.csv"
    if not evidence_path.is_file():
        raise GoldenWorkflowError(f"Evidence table not found: {evidence_path}")
    frame = pd.read_csv(evidence_path, dtype=str).fillna("")
    rows = []
    for record in frame.to_dict(orient="records"):
        decision = str(record.get("decision", "")).strip().casefold()
        flags = str(record.get("critic_flags", ""))
        confidence = str(record.get("critic_confidence", "")).strip().casefold()
        needs_review = (
            decision not in {"accepted", "supported"}
            or confidence in {"low", "needs_review", ""}
            or flags not in {"", "PASS"}
        )
        if not needs_review:
            continue
        cluster = str(record.get("cluster", ""))
        gap = gap_index.get(cluster, {})
        contrast = contrast_index.get(cluster, {})
        row_actions = list(gap.get("allowed_next_actions", []))
        if decision in {"accepted", "supported"}:
            row_actions.extend(["accept_existing_after_human_review", "human_override_with_reason"])
        else:
            row_actions.extend(["keep_unknown", "human_override_with_reason"])
        rows.append(
            {
                "cluster": cluster,
                "current_decision": record.get("decision", ""),
                "cell_type": record.get("cell_type", ""),
                "candidate_cell_type": record.get("candidate_cell_type", ""),
                "candidate_cl_id": record.get("candidate_cl_id", ""),
                "critic_confidence": record.get("critic_confidence", ""),
                "critic_flags": flags,
                "abstain_reason": record.get("abstain_reason", ""),
                "contrastive_evidence": contrast,
                "actionable_evidence_gaps": gap.get("gaps", []),
                "allowed_actions": list(dict.fromkeys(row_actions)),
                "forbidden_actions": gap.get(
                    "forbidden_actions",
                    ["do_not_replace_unknown_with_top_candidate_automatically"],
                ),
            }
        )
    payload = {
        "schema_version": REVIEW_SCHEMA,
        "created_at": _utc_now(),
        "status": "human_review_required" if rows else "no_flagged_clusters",
        "output_dir": str(root),
        "n_clusters_for_review": len(rows),
        "clusters": rows,
        "mutation_policy": "read_only_until_confirmed_human_finalize",
        "claim_boundary": (
            "The queue prioritizes review. It does not authorize the Agent to choose a label."
        ),
    }
    path = _write_json(root / REVIEW_FILE, payload)
    _write_json(
        root / STATE_FILE,
        {
            "schema_version": STATE_SCHEMA,
            "phase": "review_ready",
            "updated_at": _utc_now(),
            "review_queue": str(path),
            "next_action": "finalize_reviewed_annotations",
        },
    )
    artifacts = {
        "review_queue": str(path),
        "workflow_state": str(root / STATE_FILE),
        "evidence_table": str(evidence_path),
    }
    for name in ("evidence_gaps.json", "contrastive_evidence.csv"):
        artifact = root / name
        if artifact.is_file():
            artifacts[name.removesuffix(".json").removesuffix(".csv")] = str(artifact)
    return agent_decision(
        operation="review_uncertain_clusters",
        status=payload["status"],
        decision_scope="bounded_human_review_queue",
        evidence_summary={
            "n_clusters_for_review": len(rows),
            "n_unknown_clusters_with_actionable_gaps": len(gap_index),
        },
        allowed_next_actions=(
            ["inspect_actionable_evidence_gaps", "inspect_contrastive_evidence", "collect_human_decisions"]
            if rows
            else ["finalize_reviewed_annotations"]
        ),
        forbidden_claims=[
            "do_not_infer_or_apply_a_replacement_label_during_review",
            "do_not_interpret_top_two_score_margin_as_probability",
        ],
        artifact_paths=artifacts,
        human_action_required=bool(rows),
        claim_boundary=payload["claim_boundary"],
        review_queue_path=str(path),
        output_dir=payload["output_dir"],
        n_clusters_for_review=payload["n_clusters_for_review"],
        clusters=payload["clusters"],
        mutation_policy=payload["mutation_policy"],
        review_queue=payload,
    )


def _reconcile_human_overrides(root: Path, overrides: list[dict[str, Any]]) -> None:
    """Make the evidence table honest about human decisions before re-signing."""
    import pandas as pd

    path = root / "evidence_table.csv"
    frame = pd.read_csv(path, dtype=str).fillna("")
    for override in overrides:
        cluster = str(override["cluster"])
        mask = frame["cluster"].astype(str) == cluster
        if not mask.any():
            raise GoldenWorkflowError(f"Override cluster {cluster!r} absent from evidence table")
        frame.loc[mask, "pre_review_cell_type"] = frame.loc[mask, "cell_type"].astype(str)
        frame.loc[mask, "cell_type"] = str(override["new_type"])
        frame.loc[mask, "decision"] = "human_override"
        frame.loc[mask, "human_override_reason"] = str(override.get("reason", ""))
        frame.loc[mask, "evidence_alignment_status"] = "human_override_not_rescored"
        frame.loc[mask, "critic_flags"] = frame.loc[mask, "critic_flags"].map(
            lambda value: "; ".join(item for item in [str(value).strip(), "HUMAN_OVERRIDE"] if item)
        )
    frame.to_csv(path, index=False)


def finalize_reviewed_annotations(
    output_dir: str,
    *,
    human_overrides: list[dict[str, Any]] | None = None,
    confirm_human_review: bool = False,
    signer: str = "",
) -> dict[str, Any]:
    """Step 4: land explicit human decisions, reconcile artifacts, and re-sign."""
    from .mcp_server import tool_review_add_override, tool_review_apply_overrides
    from .review_resign import resign_review_outputs

    if not confirm_human_review:
        raise GoldenWorkflowError("confirm_human_review=True is required to finalize")
    if not str(signer).strip():
        raise GoldenWorkflowError("A human signer identity is required")
    root = Path(output_dir).resolve()
    overrides = list(human_overrides or [])
    for override in overrides:
        required = {"cluster", "new_type"}
        if required - set(override):
            raise GoldenWorkflowError("Every override requires cluster and new_type")
        staged = tool_review_add_override(
            str(root),
            str(override["cluster"]),
            str(override["new_type"]),
            str(override.get("reason", "")),
            confirm_human_review=True,
        )
        if staged.get("status") != "saved":
            raise GoldenWorkflowError(str(staged.get("error", "Override staging failed")))
    if overrides:
        applied = tool_review_apply_overrides(str(root), confirm_human_review=True)
        if applied.get("status") != "applied":
            raise GoldenWorkflowError(str(applied.get("error", "Override apply failed")))
        _reconcile_human_overrides(root, overrides)
    signed = resign_review_outputs(root, signer=signer, regenerate=True)
    state = {
        "schema_version": STATE_SCHEMA,
        "phase": "finalized_human_reviewed",
        "updated_at": _utc_now(),
        "signer": signer,
        "n_human_overrides": len(overrides),
        "signature_sha256": signed["signature"]["signature_sha256"],
        "next_action": "complete",
        "claim_boundary": (
            "Human finalization approves a review artifact; it does not convert marker scores "
            "into calibrated probabilities or establish external biological validity."
        ),
    }
    _write_json(root / STATE_FILE, state)
    return agent_decision(
        operation="finalize_reviewed_annotations",
        status=state["phase"],
        decision_scope="human_reviewed_artifact_finalization",
        evidence_summary={
            "n_human_overrides": len(overrides),
            "signer": signer,
            "signature_sha256": state["signature_sha256"],
        },
        allowed_next_actions=["export_human_reviewed_artifacts", "complete"],
        forbidden_claims=[
            "do_not_treat_human_review_as_external_benchmark_validation",
            "do_not_describe_human_overrides_as_rescored_marker_evidence",
        ],
        artifact_paths={
            "review_signature": str(root / "review_signature.json"),
            "artifact_status": str(root / "artifact_status.json"),
            "workflow_state": str(root / STATE_FILE),
        },
        human_action_required=False,
        claim_boundary=state["claim_boundary"],
        workflow_state=state,
        **signed,
    )
