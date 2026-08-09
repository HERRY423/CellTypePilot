"""Native MCP facade for CellTypePilot.

This server exposes deterministic CellTypePilot operations to Agent hosts.
It does not plan experiments, choose biological conclusions, or bypass the
same fail-closed gates used by the CLI.

Design constraints (see AGENTS.md):
- ``Plugin, not Agent``: bounded, task-scoped operations only; no autonomous
  planning, no self-directed biological claims.
- ``Governed context``: review actions only record/land a human's explicit
  annotation decision. They never invent, infer, or rescue a cell type.
- ``Review auditability``: manual review edits append to
  ``annotation_audit_log.jsonl`` and mark derived artifacts stale via
  ``artifact_status.json``.
- ``Fail closed``: unsupported species, missing artifacts, and invalid
  overrides return bounded error payloads; they never widen claim language.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__

# Artifact names that become stale after a human override is applied to the
# .h5ad file. Mirrors web_inspector.STALE_AFTER_OVERRIDE_APPLY so the MCP
# facade and the Web Inspector share one stale-by-default contract.
STALE_AFTER_OVERRIDE_APPLY = [
    "evidence_table.csv",
    "report_draft.html",
    "methodology_draft.txt",
    "manifest.json",
    "figures/",
]

AUDIT_LOG_FILENAME = "annotation_audit_log.jsonl"
ARTIFACT_STATUS_FILENAME = "artifact_status.json"
OVERRIDES_FILENAME = "annotation_overrides.json"  # pending (unapplied) overrides
OUTPUT_ANNOTATED = "data.annotated.h5ad"

# Resource names exposed by the MCP resource layer.
RESOURCE_NAMES = ("manifest", "evidence", "novelty", "artifact_status", "audit_log")


class MCPServerError(RuntimeError):
    """Raised when the optional MCP runtime is unavailable."""


def _jsonable(value: Any) -> Any:
    """Recursively normalize numpy/pandas/Path values to JSON-safe structures."""
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "to_dict"):
        try:
            return _jsonable(value.to_dict(orient="records"))
        except TypeError:
            return _jsonable(value.to_dict())
    if hasattr(value, "tolist"):
        return _jsonable(value.tolist())
    if hasattr(value, "item"):
        try:
            return _jsonable(value.item())
        except (TypeError, ValueError):
            pass
    return value


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _output_root(output_dir: str) -> Path:
    return Path(output_dir)


# ──────────────────────────────────────────────────────────────────────────
# Manifest / evidence / novelty / audit readers (deterministic, fail-closed)
# ──────────────────────────────────────────────────────────────────────────


def tool_list_artifacts(output_dir: str) -> dict[str, Any]:
    """List reviewable artifacts present in a CellTypePilot output directory.

    This is discovery-only and never changes claim language. It reports what
    durable files are available so an Agent can route a follow-up read.
    """
    root = _output_root(output_dir)
    available = {
        name: (root / name).is_file() if name not in {"figures/"} else True
        for name in (
            "data.annotated.h5ad",
            "evidence_table.csv",
            "contrastive_evidence.csv",
            "evidence_gaps.json",
            "novelty_results.csv",
            "manifest.json",
            AUDIT_LOG_FILENAME,
            ARTIFACT_STATUS_FILENAME,
            OVERRIDES_FILENAME,
            "figures/",
        )
    }
    return {
        "schema_version": "celltypepilot.artifacts.v1",
        "output_dir": str(root),
        "exists": root.is_dir(),
        "available": available,
        "note": "Artifact presence is discovery metadata, not a validity claim.",
    }


def tool_read_manifest(output_dir: str) -> dict[str, Any]:
    """Read a CellTypePilot output manifest and related audit state if present."""
    from .provenance import load_manifest

    root = _output_root(output_dir)
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        return {
            "schema_version": "celltypepilot.manifest.v1",
            "status": "not_found",
            "message": "manifest.json was not found; rerun annotation with current CellTypePilot.",
        }
    payload: dict[str, Any] = {"manifest": load_manifest(manifest_path)}
    for name in (ARTIFACT_STATUS_FILENAME, OVERRIDES_FILENAME):
        path = root / name
        if path.is_file():
            payload[name.removesuffix(".json")] = json.loads(path.read_text(encoding="utf-8"))
    audit = _load_audit_tail(root)
    if audit:
        payload["audit_log_tail"] = audit
    return {
        "schema_version": "celltypepilot.manifest.v1",
        "status": "available",
        "celltypepilot_version": __version__,
        **payload,
    }


def tool_read_evidence(output_dir: str, cluster: str | None = None) -> dict[str, Any]:
    """Read per-cluster annotation evidence from evidence_table.csv.

    Args:
        output_dir: CellTypePilot output directory.
        cluster: Optional cluster id to filter a single row.

    Provides the uncertainty-language contract alongside raw rows so Agents
    cannot widen ``evidence_score`` / ``combined_score`` into probabilities.
    """
    import pandas as pd

    root = _output_root(output_dir)
    path = root / "evidence_table.csv"
    if not path.is_file():
        return {
            "schema_version": "celltypepilot.evidence.v1",
            "status": "not_found",
            "message": "evidence_table.csv was not found; rerun annotation with current CellTypePilot.",
        }
    frame = pd.read_csv(path)
    # Normalize cluster ids to strings so the schema is type-stable regardless
    # of how pandas inferred the CSV column type on disk.
    if "cluster" in frame.columns:
        frame["cluster"] = frame["cluster"].astype(str)
    if cluster is not None:
        frame = frame[frame["cluster"] == str(cluster)]
        if frame.empty:
            return {
                "schema_version": "celltypepilot.evidence.v1",
                "status": "cluster_not_found",
                "cluster": str(cluster),
                "message": f"Cluster '{cluster}' not found in evidence_table.csv.",
            }
    from .agent_evidence import load_agent_evidence_indexes

    contrast_index, gap_index = load_agent_evidence_indexes(root)
    if cluster is not None:
        contrast_index = (
            {str(cluster): contrast_index[str(cluster)]}
            if str(cluster) in contrast_index
            else {}
        )
        gap_index = (
            {str(cluster): gap_index[str(cluster)]} if str(cluster) in gap_index else {}
        )
    return {
        "schema_version": "celltypepilot.evidence.v1",
        "status": "available",
        "rows": _jsonable(frame),
        "cluster_filter": cluster,
        "n_rows": int(len(frame)),
        "contrastive_evidence": list(contrast_index.values()),
        "actionable_evidence_gaps": list(gap_index.values()),
        "claim_boundary": (
            "evidence_score / combined_score are evidence-ranking signals, "
            "not calibrated probabilities; critic_confidence is a rule-based "
            "review category, not a probability."
        ),
    }


def tool_read_novelty_results(output_dir: str) -> dict[str, Any]:
    """Read novelty/OOD review results from a CellTypePilot output directory.

    Novelty/OOD output is a review-priority signal for human sign-off, not
    validated new cell-type discovery. This reader never renames identity or
    assigns new ontology terms.
    """
    import pandas as pd

    root = _output_root(output_dir)
    path = root / "novelty_results.csv"
    if not path.is_file():
        return {
            "schema_version": "celltypepilot.novelty-ood.v1",
            "status": "not_found",
            "message": "novelty_results.csv was not found; rerun annotation with current CellTypePilot.",
        }
    frame = pd.read_csv(path)
    return {
        "schema_version": "celltypepilot.novelty-ood.v1",
        "status": "available",
        "rows": _jsonable(frame),
        "decision_counts": (
            frame.get("novelty_decision", pd.Series(dtype=str)).value_counts().to_dict()
        ),
        "claim_boundary": (
            "Novelty/OOD output is a review-priority signal, not validated new "
            "cell-type discovery. Atlas-gap / OOD candidates require artifact/QC "
            "review, external evidence, and human sign-off."
        ),
    }


def tool_read_artifact_status(output_dir: str) -> dict[str, Any]:
    """Return whether derived artifacts are current or stale after review edits."""
    root = _output_root(output_dir)
    status = _load_artifact_status(root)
    return {
        "schema_version": "celltypepilot.artifact-status.v1",
        "output_dir": str(root),
        "artifact_status": status,
    }


def tool_qc_diagnostics(
    input_path: str,
    cluster_key: str | None = None,
    doublet_table_path: str | None = None,
    ambient_table_path: str | None = None,
) -> dict[str, Any]:
    """Assemble QC diagnostic axes without changing identity labels.

    Missing metadata → not_assessed (never clean). Doublet/ambient external tables
    are optional diagnostic inputs only.
    """
    from .data_adapter import load_h5ad
    from .qc_diagnostics import (
        QCDiagnosticError,
        assemble_qc_diagnostics,
        load_external_tool_table,
    )

    try:
        adata = load_h5ad(input_path)
        doublet = (
            load_external_tool_table(doublet_table_path, axis="doublet")
            if doublet_table_path
            else None
        )
        ambient = (
            load_external_tool_table(ambient_table_path, axis="ambient_rna")
            if ambient_table_path
            else None
        )
        report = assemble_qc_diagnostics(
            adata,
            cluster_key=cluster_key,
            doublet_table=doublet,
            ambient_table=ambient,
        )
    except (QCDiagnosticError, OSError, FileNotFoundError, ValueError) as exc:
        return {
            "schema_version": "celltypepilot.qc-diagnostics.v1",
            "status": "error",
            "error": str(exc),
            "can_rescue_identity": False,
        }
    return _jsonable(report)


def tool_agent_lifecycle_status(output_dir: str) -> dict[str, Any]:
    """Classify run/output lifecycle for Agent hosts (not a biological claim).

    Returns agent_state discrimination over checkpoints and artifact staleness:
    running / completed / failed / unavailable / claim_ready / incomplete_not_claim_ready.
    """
    from .agent_lifecycle import build_agent_status_report, scan_checkpoint_dir

    root = _output_root(output_dir)
    checkpoint_dir = root / "checkpoints"
    release_path = root / "release_manifest.json"
    release = None
    if release_path.is_file():
        try:
            release = json.loads(release_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            release = {"status": "unknown", "error": "release_manifest_unreadable"}
    report = build_agent_status_report(
        checkpoint_dir=checkpoint_dir if checkpoint_dir.is_dir() else None,
        release_manifest=release,
    )
    report["output_dir"] = str(root)
    report["artifact_status"] = _load_artifact_status(root)
    if checkpoint_dir.is_dir():
        report["checkpoint_scan"] = scan_checkpoint_dir(checkpoint_dir)
    report["prediction_mutation_allowed"] = False
    report["note"] = (
        "Lifecycle status is for Agent routing. completed ≠ claim_ready. "
        "Do not invent predictions for unavailable/failed folds."
    )
    return _jsonable(report)


def tool_read_audit_log(
    output_dir: str, limit: int = 20, event_type: str | None = None
) -> dict[str, Any]:
    """Read the tail of the append-only Web Review audit log."""
    root = _output_root(output_dir)
    rows = _load_audit_tail(root, limit=limit)
    if event_type:
        rows = [row for row in rows if row.get("event_type") == event_type]
    return {
        "schema_version": "celltypepilot.audit-log.v1",
        "output_dir": str(root),
        "events": rows,
        "n_events": len(rows),
        "note": "Audit log is append-only; events are immutable review records.",
    }


# ──────────────────────────────────────────────────────────────────────────
# Review actions (record / land a *human* annotation decision only)
#
# These operations never generate a biological conclusion. They record an
# explicit human decision, and applying overrides writes that decision back
# to the .h5ad file (with a timestamped backup) and marks derived artifacts
# stale. An Agent must not use these to invent or infer cell types.
# ──────────────────────────────────────────────────────────────────────────


def tool_review_list_overrides(output_dir: str) -> dict[str, Any]:
    """List pending (unapplied) human annotation overrides for an output dir."""
    root = _output_root(output_dir)
    overrides = _load_overrides(root)
    return {
        "schema_version": "celltypepilot.review.v1",
        "output_dir": str(root),
        "overrides": _jsonable(overrides),
        "count": len(overrides),
        "note": (
            "Overrides are pending human decisions; they are not written to the "
            ".h5ad file until tool_review_apply_overrides is called."
        ),
    }


def tool_review_add_override(
    output_dir: str,
    cluster: str,
    new_type: str,
    reason: str = "",
    confirm_human_review: bool = False,
) -> dict[str, Any]:
    """Record one human annotation override (does NOT modify the .h5ad).

    Args:
        output_dir: CellTypePilot output directory.
        cluster: Cluster id to override.
        new_type: The human-approved cell-type label.
        reason: Free-text human justification (provenance only, never evidence).
        confirm_human_review: MUST be True. This flag documents that the label
            came from a human decision, not from Agent inference.

    This is strictly a record-keeping action. It fails closed unless
    ``confirm_human_review`` is True, so Agents cannot silently stage their own
    annotation choices as if they were human review.
    """
    root = _output_root(output_dir)
    if not confirm_human_review:
        return {
            "schema_version": "celltypepilot.review.v1",
            "status": "error",
            "error": (
                "confirm_human_review must be True: this action records a human "
                "annotation decision, not an Agent-inferred conclusion."
            ),
        }
    if not str(cluster).strip() or not str(new_type).strip():
        return {
            "schema_version": "celltypepilot.review.v1",
            "status": "error",
            "error": "Both 'cluster' and 'new_type' must be non-empty.",
        }
    if not (root / OUTPUT_ANNOTATED).is_file():
        return {
            "schema_version": "celltypepilot.review.v1",
            "status": "error",
            "error": f"{OUTPUT_ANNOTATED} not found in output dir {root}.",
        }

    overrides = _load_overrides(root)
    overrides[str(cluster)] = {
        "new_type": str(new_type),
        "reason": str(reason),
        "timestamp": _utc_now(),
    }
    _save_overrides(root, overrides)
    _append_audit_event(
        root,
        "override_saved",
        {"cluster": str(cluster), "new_type": str(new_type), "reason": str(reason)},
    )
    return {
        "schema_version": "celltypepilot.review.v1",
        "status": "saved",
        "cluster": str(cluster),
        "new_type": str(new_type),
        "total_overrides": len(overrides),
        "note": "Override staged for review; call tool_review_apply_overrides to land it.",
    }


def tool_review_delete_override(output_dir: str, cluster: str) -> dict[str, Any]:
    """Remove one pending (unapplied) human override."""
    root = _output_root(output_dir)
    overrides = _load_overrides(root)
    if str(cluster) not in overrides:
        return {
            "schema_version": "celltypepilot.review.v1",
            "status": "not_found",
            "cluster": str(cluster),
            "error": "No pending override for that cluster.",
        }
    removed = overrides.pop(str(cluster))
    _save_overrides(root, overrides)
    _append_audit_event(root, "override_deleted", {"cluster": str(cluster), "override": removed})
    return {
        "schema_version": "celltypepilot.review.v1",
        "status": "deleted",
        "cluster": str(cluster),
        "total_overrides": len(overrides),
    }


def tool_review_clear_overrides(output_dir: str) -> dict[str, Any]:
    """Clear all pending overrides (does NOT revert already-applied .h5ad edits)."""
    root = _output_root(output_dir)
    cleared = len(_load_overrides(root))
    _save_overrides(root, {})
    _append_audit_event(root, "overrides_cleared", {"cleared": cleared})
    return {
        "schema_version": "celltypepilot.review.v1",
        "status": "cleared",
        "cleared": cleared,
        "note": "Only pending overrides were cleared; applied .h5ad edits are not reverted.",
    }


def tool_review_apply_overrides(
    output_dir: str, confirm_human_review: bool = False
) -> dict[str, Any]:
    """Apply all pending human overrides to data.annotated.h5ad.

    Creates a timestamped backup before modifying, appends an audit event, and
    marks derived artifacts (evidence, figures, report, methodology, manifest)
    as stale so downstream automated use cannot treat them as current.

    Args:
        output_dir: CellTypePilot output directory.
        confirm_human_review: MUST be True. This lands human-approved labels
            only; it never lets an Agent write its own conclusions.
    """
    root = _output_root(output_dir)
    if not confirm_human_review:
        return {
            "schema_version": "celltypepilot.review.v1",
            "status": "error",
            "error": (
                "confirm_human_review must be True: applying overrides writes "
                "human-approved labels into the data file."
            ),
        }
    overrides = _load_overrides(root)
    if not overrides:
        return {
            "schema_version": "celltypepilot.review.v1",
            "status": "error",
            "error": "No pending overrides to apply.",
        }
    h5ad_path = root / OUTPUT_ANNOTATED
    if not h5ad_path.is_file():
        return {
            "schema_version": "celltypepilot.review.v1",
            "status": "error",
            "error": f"{OUTPUT_ANNOTATED} not found in output dir {root}.",
        }

    from .orchestrator import apply_overrides_to_h5ad

    try:
        result = apply_overrides_to_h5ad(h5ad_path, overrides)
    except Exception as exc:  # bounded: surface apply failure, no claim widening
        return {
            "schema_version": "celltypepilot.review.v1",
            "status": "error",
            "error": f"Override apply failed: {exc}",
        }

    status = _mark_artifacts_stale_after_apply(root, result)
    _append_audit_event(root, "overrides_applied", {"result": result, "artifact_status": status})
    _save_overrides(root, {})  # applied overrides are no longer pending

    return {
        "schema_version": "celltypepilot.review.v1",
        "status": "applied",
        "result": _jsonable(result),
        "artifact_status": status,
        "note": (
            "Derived artifacts are now stale and must be regenerated before "
            "publication or downstream automated use."
        ),
    }


# ──────────────────────────────────────────────────────────────────────────
# Underlying file helpers (shared by tools; deterministic and fail-closed)
# ──────────────────────────────────────────────────────────────────────────


def _overrides_path(root: Path) -> Path:
    return root / OVERRIDES_FILENAME


def _audit_log_path(root: Path) -> Path:
    return root / AUDIT_LOG_FILENAME


def _artifact_status_path(root: Path) -> Path:
    return root / ARTIFACT_STATUS_FILENAME


def _load_overrides(root: Path) -> dict[str, Any]:
    path = _overrides_path(root)
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (json.JSONDecodeError, ValueError):
        return {}


def _save_overrides(root: Path, overrides: dict[str, Any]) -> None:
    path = _overrides_path(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(overrides, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _append_audit_event(
    root: Path, event_type: str, payload: dict[str, Any] | None = None
) -> dict[str, Any]:
    event = {
        "schema_version": "celltypepilot.web-audit.v1",
        "timestamp": _utc_now(),
        "event_type": event_type,
        "payload": payload or {},
    }
    path = _audit_log_path(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(event, ensure_ascii=False) + "\n")
    return event


def _load_audit_tail(root: Path, limit: int = 20) -> list[dict[str, Any]]:
    path = _audit_log_path(root)
    if not path.is_file():
        return []
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines()[-limit:]:
        if not line.strip():
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            rows.append({"event_type": "unparseable_audit_line", "raw": line})
    return rows


def _default_artifact_status() -> dict[str, Any]:
    return {
        "schema_version": "celltypepilot.artifact-status.v1",
        "updated_at": _utc_now(),
        "review_state": "current",
        "stale_artifacts": [],
        "current_artifacts": [OUTPUT_ANNOTATED],
        "message": "No applied review overrides have marked derived artifacts stale.",
    }


def _load_artifact_status(root: Path) -> dict[str, Any]:
    path = _artifact_status_path(root)
    if not path.is_file():
        return _default_artifact_status()
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        status = _default_artifact_status()
        status["review_state"] = "status_unreadable"
        status["message"] = f"{ARTIFACT_STATUS_FILENAME} is not valid JSON"
        return status


def _mark_artifacts_stale_after_apply(root: Path, result: dict[str, Any]) -> dict[str, Any]:
    status = {
        "schema_version": "celltypepilot.artifact-status.v1",
        "review_state": "applied_overrides_artifacts_stale",
        "stale_artifacts": STALE_AFTER_OVERRIDE_APPLY,
        "current_artifacts": [OUTPUT_ANNOTATED],
        "last_apply_result": _jsonable(result),
        "message": (
            "Human overrides were written to data.annotated.h5ad. Derived evidence, "
            "figures, report, methodology, and manifest should be regenerated before "
            "publication or downstream automated use."
        ),
    }
    path = _artifact_status_path(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return status


# ──────────────────────────────────────────────────────────────────────────
# Existing deterministic tools (inspect, markers, governance, critic, annotate)
# ──────────────────────────────────────────────────────────────────────────


def tool_inspect_h5ad(
    input_path: str,
    cluster_key: str | None = None,
    embedding_key: str | None = None,
) -> Any:
    """Inspect an h5ad file and report support boundaries."""
    from .data_adapter import inspect_adata

    return _jsonable(inspect_adata(input_path, cluster_key, embedding_key))


def tool_list_marker_scope(
    species: str = "human",
    tissue: str | None = None,
    packs: list[str] | None = None,
) -> Any:
    """List supported tissues or marker definitions for one tissue."""
    from .data_adapter import get_all_markers_for_tissue, load_marker_atlas

    atlas = load_marker_atlas(species)
    warnings: list[str] = []
    if packs:
        from .pack_manager import merge_marker_atlas, resolve_extension_packs

        records, pack_warnings = resolve_extension_packs(packs, species)
        atlas, merge_warnings = merge_marker_atlas(atlas, records, species)
        warnings.extend(pack_warnings)
        warnings.extend(merge_warnings)
    if tissue is None:
        return {
            "species": species,
            "available_tissues": sorted(atlas.get("tissues", {})),
            "warnings": warnings,
        }
    markers = get_all_markers_for_tissue(atlas, tissue)
    return {
        "species": species,
        "tissue": tissue,
        "markers": _jsonable(markers),
        "warnings": warnings,
    }


def tool_atlas_governance(include_packs: bool = True) -> Any:
    """Return the offline atlas governance report."""
    from .atlas_governance import build_atlas_governance_report

    return _jsonable(build_atlas_governance_report(include_packs=include_packs))


def tool_uncertainty_language(
    calibration_policy_path: str | None = None,
    uses_reference: bool = False,
) -> Any:
    """Return CellTypePilot's score/confidence/probability claim contract."""
    from .uncertainty import build_uncertainty_language_manifest

    policy = None
    if calibration_policy_path:
        policy = json.loads(Path(calibration_policy_path).read_text(encoding="utf-8"))
    return _jsonable(
        build_uncertainty_language_manifest(
            calibration_policy=policy,
            uses_reference=uses_reference,
        )
    )


def tool_evidence_coverage(
    input_path: str,
    species: str,
    tissue: str,
    packs: list[str] | None = None,
    evidence_policy: str = "database",
) -> dict[str, Any]:
    """Preflight gene/marker reachability without making annotations."""
    from .data_adapter import load_h5ad, load_marker_atlas
    from .identity_contract import (
        apply_gene_identity_contract,
        collect_pack_identity_contract,
        compose_marker_definitions,
        restore_original_gene_identifiers,
    )

    atlas = load_marker_atlas(species)
    records: list[dict] = []
    warnings: list[str] = []
    if packs:
        from .pack_manager import merge_marker_atlas, resolve_extension_packs

        records, pack_warnings = resolve_extension_packs(packs, species)
        atlas, merge_warnings = merge_marker_atlas(atlas, records, species)
        warnings.extend(pack_warnings)
        warnings.extend(merge_warnings)
    contract = collect_pack_identity_contract(records)
    runtime_markers, scope = compose_marker_definitions(
        atlas, tissue, evidence_policy=evidence_policy, pack_contract=contract
    )
    markers, candidate_scope = compose_marker_definitions(
        atlas,
        tissue,
        evidence_policy=evidence_policy,
        pack_contract=contract,
        include_unverified_candidates=True,
    )
    marker_universe = {
        gene
        for info in markers.values()
        for gene in (
            list(info.get("positive_markers", [])) + list(info.get("negative_markers", []))
        )
    }
    adata = load_h5ad(input_path)
    audit = apply_gene_identity_contract(adata, marker_universe)
    active = set(map(str, adata.var_names))
    per_type = []
    for cell_type, info in markers.items():
        expected = list(info.get("positive_markers", []))
        reachable = [gene for gene in expected if gene in active]
        runtime_expected = list(runtime_markers.get(cell_type, {}).get("positive_markers", []))
        runtime_reachable = [gene for gene in runtime_expected if gene in active]
        per_type.append(
            {
                "cell_type": cell_type,
                "cl_id": info.get("cl_id", ""),
                "evidence_tissue": info.get("evidence_tissue", ""),
                "n_positive_markers": len(expected),
                "n_reachable_positive_markers": len(reachable),
                "reachable_fraction": len(reachable) / len(expected) if expected else 0.0,
                "unreachable_positive_markers": [gene for gene in expected if gene not in active],
                "n_runtime_eligible_positive_markers": len(runtime_expected),
                "n_reachable_runtime_eligible_positive_markers": len(runtime_reachable),
                "runtime_reachable_fraction": (
                    len(runtime_reachable) / len(runtime_expected) if runtime_expected else 0.0
                ),
            }
        )
    restore_original_gene_identifiers(adata)
    return _jsonable(
        {
            "schema_version": "celltypepilot.evidence-coverage.v1",
            "input_path": str(Path(input_path)),
            "species": species,
            "tissue": tissue,
            "evidence_policy": evidence_policy,
            "gene_identity": audit,
            "identity_scope": scope,
            "candidate_inventory_scope": candidate_scope,
            "cell_types": per_type,
            "warnings": warnings,
            "claim_boundary": (
                "Candidate inventory coverage measures addressability, not annotation accuracy. "
                "Only runtime-eligible markers satisfying evidence_policy may enter scoring."
            ),
        }
    )


def tool_evidence_trace(
    cell_type: str,
    species: str = "human",
    tissue: str = "general",
    packs: list[str] | None = None,
    cl_id: str = "",
) -> dict[str, Any]:
    """Trace one requested identity to canonical label and marker-edge sources."""
    from .data_adapter import load_marker_atlas
    from .identity_contract import (
        build_identity_resolver,
        collect_pack_identity_contract,
        compose_marker_definitions,
        resolve_identity_label,
    )

    atlas = load_marker_atlas(species)
    records: list[dict] = []
    warnings: list[str] = []
    if packs:
        from .pack_manager import merge_marker_atlas, resolve_extension_packs

        records, pack_warnings = resolve_extension_packs(packs, species)
        atlas, merge_warnings = merge_marker_atlas(atlas, records, species)
        warnings.extend(pack_warnings)
        warnings.extend(merge_warnings)
    contract = collect_pack_identity_contract(records)
    markers, scope = compose_marker_definitions(atlas, tissue, pack_contract=contract)
    resolver = build_identity_resolver(atlas, scope["active_tissues"], contract)
    resolution = resolve_identity_label(cell_type, resolver, cl_id)
    info = markers.get(resolution["canonical_label"])
    return _jsonable(
        {
            "schema_version": "celltypepilot.evidence-trace.v1",
            "resolution": resolution,
            "identity_scope": scope,
            "marker_definition": info,
            "status": "resolved" if info else "unresolved_no_runtime_evidence",
            "warnings": warnings,
            "claim_boundary": (
                "A trace reports provenance and scope. Aggregate or co-occurrence sources are "
                "not primary-source validation."
            ),
        }
    )


def tool_resolve_evidence_packs(
    species: str,
    tissue: str,
    disease: str | None = None,
) -> dict[str, Any]:
    """List compatible packs without installing or trusting them."""
    from .pack_manager import list_installed_packs

    candidates = []
    for pack in list_installed_packs():
        if species not in pack.get("species", []):
            continue
        if tissue not in pack.get("tissues", []):
            continue
        diseases = [str(value).casefold() for value in pack.get("diseases", [])]
        if disease and diseases and str(disease).casefold() not in diseases:
            continue
        candidates.append(pack)
    return _jsonable(
        {
            "schema_version": "celltypepilot.pack-resolution.v1",
            "species": species,
            "tissue": tissue,
            "disease": disease,
            "candidates": candidates,
            "selection_policy": "compatible_only_no_auto_install_no_trust_upgrade",
        }
    )


def tool_evidence_gap_queue(output_dir: str) -> dict[str, Any]:
    """Read annotation gaps and return a deterministic curation queue."""
    import pandas as pd

    root = _output_root(output_dir)
    actionable_path = root / "evidence_gaps.json"
    if actionable_path.is_file():
        payload = json.loads(actionable_path.read_text(encoding="utf-8"))
        return _jsonable(
            {
                **payload,
                "path": str(actionable_path),
                "mutation_policy": "read_only_bounded_actions_no_label_selection",
            }
        )
    path = root / "evidence_table.csv"
    if not path.is_file():
        return {"status": "missing", "path": str(path), "gaps": []}
    frame = pd.read_csv(path, dtype=str).fillna("")
    gap_flags = {
        "NO_MARKERS",
        "UNKNOWN_ATLAS_LABEL",
        "LOW_DE_SUPPORT",
        "AGGREGATE_PROVENANCE_ONLY",
    }
    gaps = []
    for row in frame.to_dict(orient="records"):
        flags = set(str(row.get("critic_flags", "")).split("; "))
        active = sorted(flags & gap_flags)
        if not active:
            continue
        gaps.append(
            {
                "cluster": row.get("cluster", ""),
                "candidate_cell_type": row.get("candidate_cell_type") or row.get("cell_type", ""),
                "candidate_cl_id": row.get("candidate_cl_id", ""),
                "flags": active,
                "priority": 100 if "NO_MARKERS" in active else 75,
                "decision": row.get("decision", ""),
            }
        )
    gaps.sort(key=lambda item: (-item["priority"], item["cluster"]))
    return _jsonable(
        {
            "schema_version": "celltypepilot.evidence-gap-queue.v1",
            "status": "review_required" if gaps else "no_declared_gaps",
            "path": str(path),
            "gaps": gaps,
            "mutation_policy": "read_only_proposals_only",
        }
    )


def tool_benchmark_card(output_dir: str) -> dict[str, Any]:
    """Summarize benchmark endpoints while preserving claim boundaries."""
    import numpy as np
    import pandas as pd

    root = _output_root(output_dir)
    manifest_path = root / "benchmark_manifest.json"
    manifest = (
        json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    )

    def records(name: str) -> list[dict]:
        path = root / name
        return (
            pd.read_csv(path).replace({np.nan: None}).to_dict(orient="records")
            if path.is_file()
            else []
        )

    cell_results = records("benchmark_results.csv")
    cluster_results = records("cluster_track_results.csv")
    diagnostics = records("cluster_track_diagnostics.csv")
    release_manifest = root / "release_manifest.json"
    return _jsonable(
        {
            "schema_version": "celltypepilot.benchmark-card.v1",
            "manifest": manifest,
            "cell_track": cell_results,
            "cluster_track": cluster_results,
            "cluster_diagnostics": diagnostics,
            "claim_ready": release_manifest.is_file(),
            "claim_boundary": (
                "Cell and cluster tracks are separate endpoints. Presence of result tables "
                "does not establish claim readiness without a complete release manifest."
            ),
        }
    )


def tool_critic_review_h5ad(
    input_path: str,
    cluster_key: str,
    focus_cluster: str,
    species: str | None = None,
    tissue: str | None = None,
) -> Any:
    """Run deterministic critic review for one cluster in a local h5ad file."""
    from .data_adapter import load_h5ad
    from .orchestrator import critic_review

    adata = load_h5ad(input_path)
    return _jsonable(
        critic_review(
            adata,
            cluster_key,
            focus_cluster,
            species=species,
            tissue=tissue,
        )
    )


def tool_annotate_clusters(
    input_path: str,
    cluster_key: str,
    output_dir: str,
    species: str | None = None,
    tissue: str | None = None,
    embedding_key: str | None = None,
    no_figures: bool = False,
) -> Any:
    """Run the bounded annotation pipeline and write reviewable artifacts."""
    from .orchestrator import run_annotation_pipeline

    result = run_annotation_pipeline(
        input_path=input_path,
        cluster_key=cluster_key,
        output_dir=output_dir,
        species=species,
        tissue=tissue,
        embedding_key=embedding_key,
        no_figures=no_figures,
    )
    return _jsonable(
        {
            "species": result["species"],
            "tissue": result["tissue"],
            "critic_summary": result["critic_summary"],
            "validation_scope": result.get("validation_scope"),
            "novelty_decision_counts": result.get("manifest", {})
            .get("parameters", {})
            .get("novelty_ood", {})
            .get("decision_counts", {}),
            "paths": result["paths"],
            "manifest": result["manifest"],
        }
    )


# ──────────────────────────────────────────────────────────────────────────
# Resource schema (read-only exposure of reviewable artifacts)
# ──────────────────────────────────────────────────────────────────────────


def read_output_resource(output_dir: str, artifact: str) -> str:
    """Return one reviewable artifact as a JSON string (shared core logic).

    Used by the MCP resource layer. Kept framework-free so it is unit-testable
    without the optional MCP runtime installed. Only read-only artifacts are
    exposed as resources; mutation is restricted to explicit tool calls.
    """
    if artifact not in RESOURCE_NAMES:
        return json.dumps(
            {
                "schema_version": "celltypepilot.resource.v1",
                "status": "error",
                "error": f"Unknown artifact '{artifact}'. Available: {list(RESOURCE_NAMES)}",
            },
            ensure_ascii=False,
        )

    dispatch = {
        "manifest": tool_read_manifest,
        "evidence": tool_read_evidence,
        "novelty": tool_read_novelty_results,
        "artifact_status": tool_read_artifact_status,
        "audit_log": tool_read_audit_log,
    }
    return json.dumps(dispatch[artifact](output_dir), ensure_ascii=False, indent=2)


def build_mcp_server():
    """Create the FastMCP server, or raise an actionable dependency error.

    Exposes read-only resources for reviewable artifacts and deterministic
    tools. Review tools are gated on ``confirm_human_review`` so an Agent can
    only land labels that a human explicitly approved, never its own inference.
    """
    try:
        from fastmcp import FastMCP
    except ImportError as exc:
        raise MCPServerError(
            "CellTypePilot MCP support requires the optional MCP runtime. "
            "Install with: pip install -e .[mcp]"
        ) from exc

    mcp = FastMCP("CellTypePilot")

    # Read-only resources for reviewable artifacts.
    mcp.resource("celltypepilot://output/{output_dir}/{artifact}")(read_output_resource)

    # The default product surface is intentionally four stateful operations.
    # Maintainers can opt into legacy/diagnostic primitives with
    # CELLTYPEPILOT_MCP_SURFACE=advanced; ordinary Agents should not assemble
    # biological workflows from dozens of low-level calls.
    from .golden_workflow import (
        annotate_from_plan,
        finalize_reviewed_annotations,
        prepare_annotation,
        review_uncertain_clusters,
    )

    mcp.tool()(prepare_annotation)
    mcp.tool()(annotate_from_plan)
    mcp.tool()(review_uncertain_clusters)
    mcp.tool()(finalize_reviewed_annotations)

    if os.getenv("CELLTYPEPILOT_MCP_SURFACE", "golden").strip().casefold() == "advanced":
        # Deterministic maintainer / compatibility tools.
        mcp.tool()(tool_inspect_h5ad)
        mcp.tool()(tool_list_marker_scope)
        mcp.tool()(tool_atlas_governance)
        mcp.tool()(tool_uncertainty_language)
        mcp.tool()(tool_list_artifacts)
        mcp.tool()(tool_read_manifest)
        mcp.tool()(tool_read_evidence)
        mcp.tool()(tool_read_novelty_results)
        mcp.tool()(tool_read_artifact_status)
        mcp.tool()(tool_agent_lifecycle_status)
        mcp.tool()(tool_qc_diagnostics)
        mcp.tool()(tool_read_audit_log)
        mcp.tool()(tool_evidence_coverage)
        mcp.tool()(tool_evidence_trace)
        mcp.tool()(tool_resolve_evidence_packs)
        mcp.tool()(tool_evidence_gap_queue)
        mcp.tool()(tool_benchmark_card)
        mcp.tool()(tool_critic_review_h5ad)
        mcp.tool()(tool_annotate_clusters)
        mcp.tool()(tool_review_list_overrides)
        mcp.tool()(tool_review_add_override)
        mcp.tool()(tool_review_delete_override)
        mcp.tool()(tool_review_clear_overrides)
        mcp.tool()(tool_review_apply_overrides)

    return mcp


def main() -> None:
    """Run the CellTypePilot MCP server over stdio."""
    server = build_mcp_server()
    server.run()


if __name__ == "__main__":
    main()
