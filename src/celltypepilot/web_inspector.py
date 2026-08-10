"""CellTypePilot — Web Inspector: lightweight interactive review panel (Review Cockpit).

A Flask app serving an interactive Review Cockpit for reviewing, auditing,
and correcting cell-type annotations.
"""

from __future__ import annotations

import importlib.resources
import json
from datetime import datetime, timezone
from pathlib import Path

from flask import Flask, jsonify, make_response, render_template, request

from . import __version__

# Frontend assets live in package resources (templates/ + static/),
# so the Flask app resolves them regardless of the working directory.
_pkg_root = importlib.resources.files("celltypepilot")

app = Flask(
    __name__,
    template_folder=str(_pkg_root / "templates"),
    static_folder=str(_pkg_root / "static"),
)

# Global state — set by `celltypepilot inspect-web`
_output_dir: Path | None = None
# Optional dedicated benchmark-run root for read-only observability
# (defaults to _output_dir; never written by observability endpoints).
_run_dir: Path | None = None
_adata_cache = None
_evidence_cache = None

# Server-side override store (persists across page reloads)
_overrides: dict = {}  # {cluster_id: {new_type, reason, timestamp}}

AUDIT_LOG_FILENAME = "annotation_audit_log.jsonl"
ARTIFACT_STATUS_FILENAME = "artifact_status.json"
CHECKLIST_FILENAME = "review_checklist.json"
CLUSTER_REVIEWS_FILENAME = "cluster_reviews.json"
SIGNOFF_FILENAME = "review_signoff.json"

STALE_AFTER_OVERRIDE_APPLY = [
    "evidence_table.csv",
    "report_draft.html",
    "methodology_draft.txt",
    "manifest.json",
    "figures/",
]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _load_data():
    """Lazy-load annotation data from output directory."""
    global _adata_cache, _evidence_cache
    if _adata_cache is not None:
        return _adata_cache, _evidence_cache

    import pandas as pd
    import scanpy as sc

    adata_path = _output_dir / "data.annotated.h5ad"
    evidence_path = _output_dir / "evidence_table.csv"

    if not adata_path.exists():
        raise FileNotFoundError(f"No annotated data found at {adata_path}")

    _adata_cache = sc.read_h5ad(adata_path)

    _evidence_cache = pd.read_csv(evidence_path) if evidence_path.exists() else pd.DataFrame()

    # Load any existing overrides from disk
    _load_overrides_from_disk()

    return _adata_cache, _evidence_cache


def _overrides_path() -> Path:
    """Path to the server-side overrides file."""
    return _output_dir / "annotation_overrides.json"


def _audit_log_path() -> Path:
    return _output_dir / AUDIT_LOG_FILENAME


def _artifact_status_path() -> Path:
    return _output_dir / ARTIFACT_STATUS_FILENAME


def _checklist_path() -> Path:
    return _output_dir / CHECKLIST_FILENAME


def _cluster_reviews_path() -> Path:
    return _output_dir / CLUSTER_REVIEWS_FILENAME


def _signoff_path() -> Path:
    return _output_dir / SIGNOFF_FILENAME


def _append_audit_event(event_type: str, payload: dict | None = None) -> dict:
    """Append one immutable Web Review event to the audit log."""
    event = {
        "schema_version": "celltypepilot.web-audit.v1",
        "timestamp": _utc_now(),
        "event_type": event_type,
        "payload": payload or {},
    }
    path = _audit_log_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(event, ensure_ascii=False) + "\n")
    return event


def _audit_tail(limit: int = 20) -> list[dict]:
    path = _audit_log_path()
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


def _default_artifact_status() -> dict:
    return {
        "schema_version": "celltypepilot.artifact-status.v1",
        "updated_at": _utc_now(),
        "review_state": "current",
        "stale_artifacts": [],
        "current_artifacts": ["data.annotated.h5ad"],
        "message": "No applied Web Review overrides have marked derived artifacts stale.",
    }


def _load_artifact_status() -> dict:
    path = _artifact_status_path()
    if not path.is_file():
        return _default_artifact_status()
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        status = _default_artifact_status()
        status["review_state"] = "status_unreadable"
        status["message"] = f"{ARTIFACT_STATUS_FILENAME} is not valid JSON"
        return status


def _save_artifact_status(status: dict) -> None:
    status["updated_at"] = _utc_now()
    _artifact_status_path().write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _default_checklist() -> dict:
    return {
        "schema_version": "celltypepilot.review-checklist.v1",
        "updated_at": _utc_now(),
        "manual_items": {
            "marker_literature_alignment": {
                "label": "Lineage marker expression aligned with reference literature",
                "completed": False,
            },
            "sample_context_verified": {
                "label": "Sample background, species, and tissue context verified",
                "completed": False,
            },
            "doublet_and_neg_conflict_checked": {
                "label": "Potential doublets and negative marker conflicts inspected",
                "completed": False,
            },
            "pi_protocol_reviewed": {
                "label": "Final review protocol ready for PI sign-off",
                "completed": False,
            },
        },
    }


def _load_checklist() -> dict:
    path = _checklist_path()
    if not path.is_file():
        return _default_checklist()
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return _default_checklist()


def _save_checklist(data: dict) -> None:
    data["updated_at"] = _utc_now()
    _checklist_path().write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _load_cluster_reviews() -> dict:
    path = _cluster_reviews_path()
    if not path.is_file():
        return {"clusters": {}}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {"clusters": {}}


def _save_cluster_reviews(data: dict) -> None:
    _cluster_reviews_path().write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _load_signoff() -> dict:
    path = _signoff_path()
    if not path.is_file():
        return {
            "signed_off": False,
            "decision": "UNREVIEWED",
            "reviewer_name": "",
            "reviewer_role": "",
            "notes": "",
            "timestamp": None,
        }
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {
            "signed_off": False,
            "decision": "UNREVIEWED",
            "reviewer_name": "",
            "reviewer_role": "",
            "notes": "",
            "timestamp": None,
        }


def _save_signoff(data: dict) -> None:
    data["timestamp"] = _utc_now()
    _signoff_path().write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _compute_checklist_status() -> dict:
    checklist = _load_checklist()
    _, evidence = _load_data()

    flagged_clusters = []
    low_conf_clusters = []
    if not evidence.empty:
        for _, r in evidence.iterrows():
            cluster_str = str(r.get("cluster", ""))
            flags = str(r.get("critic_flags", "PASS"))
            conf = str(r.get("critic_confidence", "needs_review"))
            if flags != "PASS":
                flagged_clusters.append(cluster_str)
            if conf in ["low", "needs_review"]:
                low_conf_clusters.append(cluster_str)

    reviews = _load_cluster_reviews().get("clusters", {})

    def _are_reviewed(cids: list[str]) -> bool:
        return all(
            reviews.get(cluster_id, {}).get("status", "unreviewed") != "unreviewed"
            for cluster_id in cids
        )

    auto_items = {
        "critic_flags_reviewed": {
            "label": f"Critic flags reviewed ({len(flagged_clusters)} cluster(s) flagged)",
            "completed": len(flagged_clusters) == 0 or _are_reviewed(flagged_clusters),
            "flagged_clusters": flagged_clusters,
        },
        "low_confidence_verified": {
            "label": f"Low/Needs Review confidence clusters verified ({len(low_conf_clusters)} cluster(s))",
            "completed": len(low_conf_clusters) == 0 or _are_reviewed(low_conf_clusters),
            "low_conf_clusters": low_conf_clusters,
        },
        "pending_overrides_checked": {
            "label": f"Pending overrides reviewed ({len(_overrides)} pending)",
            "completed": True,
            "pending_count": len(_overrides),
        },
    }

    manual_items = checklist.get("manual_items", {})

    total_checks = len(auto_items) + len(manual_items)
    passed_checks = sum(1 for v in auto_items.values() if v.get("completed")) + sum(
        1 for v in manual_items.values() if v.get("completed")
    )

    readiness_pct = int((passed_checks / total_checks) * 100) if total_checks > 0 else 0

    return {
        "readiness_pct": readiness_pct,
        "passed_checks": passed_checks,
        "total_checks": total_checks,
        "automated_items": auto_items,
        "manual_items": manual_items,
        "is_ready_for_signoff": readiness_pct == 100,
    }


def _calculate_override_diff() -> dict:
    adata, evidence = _load_data()
    obs = adata.obs
    total_cells = len(obs)

    diff_rows = []
    modified_clusters_count = 0
    affected_cells_count = 0

    if not evidence.empty:
        for _, row in evidence.iterrows():
            cluster_id = str(row.get("cluster", ""))
            orig_type = str(row.get("cell_type", "Unknown"))
            n_cells = int(row.get("n_cells", 0))

            pending_override = _overrides.get(cluster_id)
            if pending_override:
                new_type = pending_override.get("new_type", "")
                reason = pending_override.get("reason", "")
                is_changed = new_type != orig_type
                if is_changed:
                    modified_clusters_count += 1
                    affected_cells_count += n_cells
                diff_rows.append(
                    {
                        "cluster": cluster_id,
                        "original_type": orig_type,
                        "pending_type": new_type,
                        "reason": reason,
                        "n_cells": n_cells,
                        "status": "pending_change" if is_changed else "pending_no_change",
                    }
                )
            else:
                diff_rows.append(
                    {
                        "cluster": cluster_id,
                        "original_type": orig_type,
                        "pending_type": orig_type,
                        "reason": "-",
                        "n_cells": n_cells,
                        "status": "unchanged",
                    }
                )

    affected_pct = round((affected_cells_count / total_cells * 100), 2) if total_cells > 0 else 0.0

    return {
        "total_clusters": len(diff_rows),
        "modified_clusters": modified_clusters_count,
        "affected_cells": affected_cells_count,
        "affected_cells_pct": affected_pct,
        "total_cells": total_cells,
        "diff_rows": diff_rows,
    }


def _generate_review_packet() -> dict:
    adata, evidence = _load_data()
    obs = adata.obs

    checklist_status = _compute_checklist_status()
    diff_summary = _calculate_override_diff()
    signoff = _load_signoff()
    audit_tail = _audit_tail(50)
    reviews = _load_cluster_reviews()
    artifact_status = _load_artifact_status()

    clusters_detail = []
    if not evidence.empty:
        for _, r in evidence.iterrows():
            cid = str(r.get("cluster", ""))
            clusters_detail.append(
                {
                    "cluster": cid,
                    "cell_type": r.get("cell_type", "Unknown"),
                    "cl_id": r.get("cl_id", ""),
                    "combined_score": float(r.get("combined_score", 0)),
                    "critic_confidence": r.get("critic_confidence", "unknown"),
                    "critic_flags": r.get("critic_flags", "PASS"),
                    "n_cells": int(r.get("n_cells", 0)),
                    "pending_override": _overrides.get(cid),
                    "review_state": reviews.get("clusters", {}).get(cid, {}),
                }
            )

    return {
        "schema_version": "celltypepilot.review-packet.v1",
        "generated_at": _utc_now(),
        "version": __version__,
        "dataset_summary": {
            "total_cells": len(obs),
            "total_clusters": obs["ctp_cl_id"].nunique()
            if "ctp_cl_id" in obs.columns
            else len(clusters_detail),
            "artifact_status": artifact_status,
        },
        "signoff_certificate": signoff,
        "readiness_checklist": checklist_status,
        "override_diff": diff_summary,
        "clusters": clusters_detail,
        "pending_overrides": _overrides,
        "audit_log_tail": audit_tail,
    }


def _mark_artifacts_stale_after_apply(result: dict) -> dict:
    status = {
        "schema_version": "celltypepilot.artifact-status.v1",
        "review_state": "applied_overrides_artifacts_stale",
        "stale_artifacts": STALE_AFTER_OVERRIDE_APPLY,
        "current_artifacts": ["data.annotated.h5ad"],
        "last_apply_result": result,
        "message": (
            "Manual overrides were written to data.annotated.h5ad. Derived evidence, "
            "figures, report, methodology, and manifest should be regenerated before "
            "publication or downstream automated use."
        ),
    }
    _save_artifact_status(status)
    return status


def _load_overrides_from_disk():
    """Load existing overrides from disk if present."""
    global _overrides
    opath = _overrides_path()
    if opath.exists():
        try:
            _overrides = json.loads(opath.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, ValueError):
            _overrides = {}


def _save_overrides_to_disk():
    """Persist current overrides to disk."""
    opath = _overrides_path()
    opath.write_text(json.dumps(_overrides, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _apply_overrides_to_h5ad() -> dict:
    """Apply all overrides to the .h5ad file and back up the original.

    Delegates to :func:`celltypepilot.orchestrator.apply_overrides_to_h5ad`
    so the CLI and the Web Inspector share a single implementation.

    Returns:
        Summary dict with counts and details.
    """
    from .orchestrator import apply_overrides_to_h5ad

    adata_path = _output_dir / "data.annotated.h5ad"
    result = apply_overrides_to_h5ad(adata_path, _overrides)

    # Clear cache so next load gets fresh data
    global _adata_cache
    _adata_cache = None

    return result


# ──────────────────────────────────────────────
# Routes
# ──────────────────────────────────────────────


def _observability_root() -> Path | None:
    """Directory used for read-only checkpoint / product observability."""
    if _run_dir is not None:
        return Path(_run_dir)
    if _output_dir is not None:
        return Path(_output_dir)
    return None


def _annotation_writes_blocked() -> bool:
    """Block override/apply when the output dir is a benchmark fold workspace.

    Prevents the dashboard from writing audit/override files into an active
    ``checkpoints/`` run tree that has no annotated review products.
    """
    if _output_dir is None:
        return True
    root = Path(_output_dir)
    has_checkpoints = (root / "checkpoints").is_dir()
    has_annotated = (root / "data.annotated.h5ad").is_file()
    return has_checkpoints and not has_annotated


def _build_observability_snapshot() -> dict:
    from .run_observability import ObservabilityError, build_observability_snapshot

    root = _observability_root()
    if root is None:
        return {
            "schema_version": "celltypepilot.run-observability.v1",
            "error": "No output/run directory configured",
            "read_only": True,
            "prediction_mutation_allowed": False,
        }
    try:
        return build_observability_snapshot(root)
    except ObservabilityError as exc:
        return {
            "schema_version": "celltypepilot.run-observability.v1",
            "error": str(exc),
            "run_root": str(root),
            "read_only": True,
            "prediction_mutation_allowed": False,
        }


@app.route("/")
def dashboard():
    """Main dashboard (Review Cockpit)."""
    stats = {
        "total_clusters": 0,
        "high": 0,
        "medium": 0,
        "low": 0,
        "needs_review": 0,
        "flagged": 0,
    }
    annotations = []
    evidence_dict = {}
    annotation_mode = True
    load_error = None

    try:
        adata, evidence = _load_data()
        obs = adata.obs
        stats = {
            "total_clusters": obs["ctp_cl_id"].nunique() if "ctp_cl_id" in obs.columns else 0,
            "high": int((obs["ctp_confidence"] == "high").sum())
            if "ctp_confidence" in obs.columns
            else 0,
            "medium": int((obs["ctp_confidence"] == "medium").sum())
            if "ctp_confidence" in obs.columns
            else 0,
            "low": int((obs["ctp_confidence"] == "low").sum())
            if "ctp_confidence" in obs.columns
            else 0,
            "needs_review": int((obs["ctp_confidence"] == "needs_review").sum())
            if "ctp_confidence" in obs.columns
            else 0,
            "flagged": 0,
        }
        cluster_reviews = _load_cluster_reviews().get("clusters", {})
        if not evidence.empty:
            for _, row in evidence.iterrows():
                cluster = str(row.get("cluster", ""))
                flags = row.get("critic_flags", "PASS")
                if flags != "PASS":
                    stats["flagged"] += 1

                review_info = cluster_reviews.get(cluster, {"status": "unreviewed", "notes": []})

                annotations.append(
                    {
                        "cluster": cluster,
                        "cell_type": row.get("cell_type", "Unknown"),
                        "cl_id": row.get("cl_id", ""),
                        "score": float(row.get("combined_score", 0)),
                        "confidence": row.get("critic_confidence", "unknown"),
                        "flags": flags,
                        "novelty_decision": row.get("novelty_decision", "not_assessed"),
                        "novelty_score": float(row.get("novelty_score", 0.0) or 0.0),
                        "n_cells": int(row.get("n_cells", 0)),
                        "review_status": review_info.get("status", "unreviewed"),
                        "note_count": len(review_info.get("notes", [])),
                    }
                )

                evidence_dict[str(cluster)] = {
                    "cell_type": row.get("cell_type", ""),
                    "combined_score": float(row.get("combined_score", 0)),
                    "pct_overlap": float(row.get("pct_overlap", 0)),
                    "top_markers": [],
                    "novelty_decision": row.get("novelty_decision", "not_assessed"),
                    "novelty_score": float(row.get("novelty_score", 0.0) or 0.0),
                    "top_unmapped_markers": row.get("top_unmapped_markers", ""),
                    "alternative_explanations": row.get("alternative_explanations", ""),
                    "recommended_next_actions": row.get("recommended_next_actions", ""),
                    "critic_notes": row.get("critic_notes", ""),
                }
    except FileNotFoundError as exc:
        # Benchmark-run directories may only have checkpoints — observability still works.
        annotation_mode = False
        load_error = str(exc)
        cluster_reviews = {}

    checklist_status = (
        _compute_checklist_status()
        if annotation_mode
        else {
            "readiness_pct": 0,
            "is_ready_for_signoff": False,
            "automated_items": {},
            "manual_items": {},
        }
    )
    diff_summary = (
        _calculate_override_diff()
        if annotation_mode
        else {
            "modified_clusters": 0,
            "rows": [],
        }
    )
    signoff = (
        _load_signoff()
        if annotation_mode
        else {
            "signed_off": False,
            "decision": "NOT_SIGNED",
        }
    )
    observability = _build_observability_snapshot()

    return render_template(
        "dashboard.html",
        version=__version__,
        stats=stats,
        annotations=annotations,
        evidence_json=json.dumps(evidence_dict),
        artifact_status=_load_artifact_status()
        if annotation_mode
        else {
            "review_state": "not_applicable",
            "message": load_error or "Annotation products not present; observability-only mode.",
            "stale_artifacts": [],
        },
        audit_tail=_audit_tail() if annotation_mode else [],
        checklist_status=checklist_status,
        diff_summary=diff_summary,
        signoff=signoff,
        cluster_reviews=cluster_reviews if annotation_mode else {},
        observability=observability,
        annotation_mode=annotation_mode,
        load_error=load_error,
    )


@app.route("/api/observability", methods=["GET"])
def api_observability():
    """Read-only run observability snapshot (checkpoints, ETA, host, hashes, stale).

    Never mutates predictions or fold workspaces. Manual overrides remain on the
    append-only audit log + apply-overrides path only.
    """
    snapshot = _build_observability_snapshot()
    # Defense in depth: refuse any write-shaped keys clients might mis-use.
    snapshot["prediction_mutation_allowed"] = False
    snapshot["read_only"] = True
    return jsonify({"ok": True, "observability": snapshot})


@app.route("/api/observability/checkpoints", methods=["GET"])
def api_observability_checkpoints():
    """Read-only checkpoint/*.status.json summary."""
    snapshot = _build_observability_snapshot()
    return jsonify(
        {
            "ok": True,
            "read_only": True,
            "prediction_mutation_allowed": False,
            "checkpoints": snapshot.get("checkpoints"),
            "fold_eta": snapshot.get("fold_eta"),
            "failures": snapshot.get("failures"),
        }
    )


@app.route("/api/observability/mutate", methods=["POST", "PUT", "PATCH", "DELETE"])
@app.route("/api/observability/predictions", methods=["POST", "PUT", "PATCH", "DELETE"])
def api_observability_mutation_blocked():
    """Explicitly reject prediction mutations via the observability surface."""
    return (
        jsonify(
            {
                "ok": False,
                "error": (
                    "Observability is read-only. Manual overrides must use "
                    "POST /api/override and POST /api/overrides/apply, which "
                    "append to annotation_audit_log.jsonl and mark artifacts stale."
                ),
                "prediction_mutation_allowed": False,
            }
        ),
        405,
    )


@app.route("/api/evidence/<cluster_id>")
def api_evidence(cluster_id):
    """Get detailed evidence for a cluster."""
    _, evidence = _load_data()
    if evidence.empty:
        return jsonify({})
    row = evidence[evidence["cluster"].astype(str) == str(cluster_id)]
    if row.empty:
        return jsonify({})
    return jsonify(row.iloc[0].to_dict())


@app.route("/api/clusters/<cluster_id>/review-panel", methods=["GET"])
def api_cluster_review_panel(cluster_id):
    """Identity × State × Novelty panel with markers, neighbors, strata, literature."""
    from .review_panel import build_cluster_review_panel, load_optional_csv

    try:
        adata, evidence = _load_data()
    except FileNotFoundError as exc:
        return jsonify({"ok": False, "error": str(exc)}), 404

    state = load_optional_csv(_output_dir / "state_results.csv")
    novelty = load_optional_csv(_output_dir / "novelty_results.csv")
    # Prefer ctp_cl_id if written; fall back to common cluster keys.
    cluster_key = "ctp_cl_id"
    if adata is not None and cluster_key not in adata.obs.columns:
        for key in ("leiden", "louvain", "cluster", "ctp_cluster"):
            if key in adata.obs.columns:
                cluster_key = key
                break

    panel = build_cluster_review_panel(
        cluster=str(cluster_id),
        evidence=evidence,
        state_results=state,
        novelty_results=novelty,
        adata=adata,
        cluster_key=cluster_key,
        pending_override=_overrides.get(str(cluster_id)),
        audit_events=[
            ev
            for ev in _audit_tail(100)
            if str(ev.get("payload", {}).get("cluster", "")) == str(cluster_id)
        ],
    )
    panel["ok"] = True
    panel["artifact_status"] = _load_artifact_status()
    return jsonify(panel)


@app.route("/api/review/resign", methods=["POST"])
def api_review_resign():
    """Regenerate derived artifacts and re-sign after human review (append-only audit)."""
    if _annotation_writes_blocked():
        return (
            jsonify(
                {
                    "ok": False,
                    "error": "Resign blocked on benchmark-run directory.",
                    "prediction_mutation_allowed": False,
                }
            ),
            403,
        )
    data = request.get_json(force=True, silent=True) or {}
    signer = str(data.get("signer") or "web_reviewer").strip()
    regenerate = bool(data.get("regenerate", True))
    try:
        from .review_resign import ReviewResignError, resign_review_outputs

        result = resign_review_outputs(_output_dir, signer=signer, regenerate=regenerate)
    except ReviewResignError as exc:
        return jsonify({"ok": False, "error": str(exc)}), 400
    except Exception as exc:  # pragma: no cover
        return jsonify({"ok": False, "error": f"Resign failed: {exc}"}), 500
    return jsonify({"ok": True, **result})


@app.route("/api/stats")
def api_stats():
    """Get summary statistics."""
    adata, evidence = _load_data()
    obs = adata.obs
    stats = {}
    if "ctp_confidence" in obs.columns:
        stats["confidence_distribution"] = obs["ctp_confidence"].value_counts().to_dict()
    if "ctp_cell_type" in obs.columns:
        stats["cell_type_counts"] = obs["ctp_cell_type"].value_counts().to_dict()
    stats["total_cells"] = len(obs)
    stats["total_clusters"] = obs["ctp_cl_id"].nunique() if "ctp_cl_id" in obs.columns else 0
    return jsonify(stats)


# ──────────────────────────────────────────────
# Checklist & Protocol API Routes
# ──────────────────────────────────────────────


@app.route("/api/checklist", methods=["GET"])
def api_get_checklist():
    """Get overall checklist readiness status."""
    return jsonify({"ok": True, "checklist": _compute_checklist_status()})


@app.route("/api/checklist", methods=["POST"])
def api_update_checklist():
    """Update manual checklist items.

    Expected JSON body:
        {"item_key": "marker_literature_alignment", "completed": true}
    """
    if _annotation_writes_blocked():
        return (
            jsonify(
                {
                    "ok": False,
                    "error": "Write blocked on benchmark-run directory (read-only observability).",
                    "prediction_mutation_allowed": False,
                }
            ),
            403,
        )
    data = request.get_json(force=True, silent=True) or {}
    item_key = data.get("item_key", "").strip()
    completed = bool(data.get("completed", False))

    if not item_key:
        return jsonify({"ok": False, "error": "Missing 'item_key'"}), 400

    checklist = _load_checklist()
    manual = checklist.setdefault("manual_items", {})
    if item_key not in manual:
        manual[item_key] = {"label": item_key.replace("_", " ").title(), "completed": completed}
    else:
        manual[item_key]["completed"] = completed

    _save_checklist(checklist)
    _append_audit_event("checklist_updated", {"item_key": item_key, "completed": completed})
    return jsonify({"ok": True, "checklist": _compute_checklist_status()})


# ──────────────────────────────────────────────
# Cluster-level History & Notes API Routes
# ──────────────────────────────────────────────


@app.route("/api/clusters/<cluster_id>/history", methods=["GET"])
def api_cluster_history(cluster_id):
    """Get complete lifecycle history and review notes for a cluster."""
    cluster_id = str(cluster_id)
    _, evidence = _load_data()

    base_row = {}
    if not evidence.empty:
        rows = evidence[evidence["cluster"].astype(str) == cluster_id]
        if not rows.empty:
            base_row = rows.iloc[0].to_dict()

    reviews = _load_cluster_reviews().get("clusters", {}).get(cluster_id, {})
    all_audit_events = _audit_tail(200)

    # Filter audit events for this cluster
    cluster_events = [
        ev for ev in all_audit_events if str(ev.get("payload", {}).get("cluster", "")) == cluster_id
    ]

    return jsonify(
        {
            "ok": True,
            "cluster": cluster_id,
            "baseline": base_row,
            "pending_override": _overrides.get(cluster_id),
            "review_status": reviews.get("status", "unreviewed"),
            "notes": reviews.get("notes", []),
            "audit_history": cluster_events,
        }
    )


@app.route("/api/clusters/<cluster_id>/note", methods=["POST"])
def api_add_cluster_note(cluster_id):
    """Add a review note/comment for a cluster.

    Expected JSON body:
        {"author": "Dr. Smith", "text": "Confirmed CD4 expression with qPCR panel"}
    """
    if _annotation_writes_blocked():
        return (
            jsonify(
                {
                    "ok": False,
                    "error": "Write blocked on benchmark-run directory (read-only observability).",
                    "prediction_mutation_allowed": False,
                }
            ),
            403,
        )
    cluster_id = str(cluster_id)
    data = request.get_json(force=True, silent=True) or {}
    text = data.get("text", "").strip()
    author = data.get("author", "Reviewer").strip()

    if not text:
        return jsonify({"ok": False, "error": "Missing note text"}), 400

    reviews_data = _load_cluster_reviews()
    clusters = reviews_data.setdefault("clusters", {})
    c_info = clusters.setdefault(cluster_id, {"status": "unreviewed", "notes": []})

    note_entry = {
        "id": len(c_info.get("notes", [])) + 1,
        "author": author,
        "text": text,
        "timestamp": _utc_now(),
    }
    c_info.setdefault("notes", []).append(note_entry)

    _save_cluster_reviews(reviews_data)
    _append_audit_event(
        "cluster_note_added",
        {"cluster": cluster_id, "author": author, "text": text},
    )

    return jsonify({"ok": True, "cluster": cluster_id, "notes": c_info["notes"]})


@app.route("/api/clusters/<cluster_id>/status", methods=["POST"])
def api_update_cluster_status(cluster_id):
    """Update cluster review status ('unreviewed' | 'reviewed' | 'flagged').

    Expected JSON body:
        {"status": "reviewed"}
    """
    if _annotation_writes_blocked():
        return (
            jsonify(
                {
                    "ok": False,
                    "error": "Write blocked on benchmark-run directory (read-only observability).",
                    "prediction_mutation_allowed": False,
                }
            ),
            403,
        )
    cluster_id = str(cluster_id)
    data = request.get_json(force=True, silent=True) or {}
    status = data.get("status", "").strip().lower()

    if status not in ["unreviewed", "reviewed", "flagged"]:
        return jsonify({"ok": False, "error": "Invalid status value"}), 400

    reviews_data = _load_cluster_reviews()
    clusters = reviews_data.setdefault("clusters", {})
    c_info = clusters.setdefault(cluster_id, {"status": "unreviewed", "notes": []})
    c_info["status"] = status

    _save_cluster_reviews(reviews_data)
    _append_audit_event("cluster_status_updated", {"cluster": cluster_id, "status": status})

    return jsonify({"ok": True, "cluster": cluster_id, "status": status})


# ──────────────────────────────────────────────
# Override Diff API Routes
# ──────────────────────────────────────────────


@app.route("/api/overrides/diff", methods=["GET"])
def api_override_diff():
    """Get differential impact view of pending overrides."""
    return jsonify({"ok": True, "diff": _calculate_override_diff()})


# ──────────────────────────────────────────────
# Sign-off & Review Packet API Routes
# ──────────────────────────────────────────────


@app.route("/api/signoff", methods=["GET"])
def api_get_signoff():
    """Get formal sign-off certificate status."""
    return jsonify({"ok": True, "signoff": _load_signoff()})


@app.route("/api/signoff", methods=["POST"])
def api_submit_signoff():
    """Submit formal sign-off certificate.

    Expected JSON body:
        {
            "reviewer_name": "Dr. Jane Doe",
            "reviewer_role": "Lead Bioinformatician",
            "decision": "APPROVED",  // "APPROVED" | "APPROVED_WITH_RESERVATIONS" | "REVISIONS_REQUIRED"
            "notes": "Cell types verified with lineage markers.",
            "force": false
        }
    """
    if _annotation_writes_blocked():
        return (
            jsonify(
                {
                    "ok": False,
                    "error": "Write blocked on benchmark-run directory (read-only observability).",
                    "prediction_mutation_allowed": False,
                }
            ),
            403,
        )
    data = request.get_json(force=True, silent=True) or {}
    name = data.get("reviewer_name", "").strip()
    role = data.get("reviewer_role", "").strip()
    decision = data.get("decision", "").strip().upper()
    notes = data.get("notes", "").strip()
    force = bool(data.get("force", False))

    if not name:
        return jsonify({"ok": False, "error": "Missing 'reviewer_name'"}), 400
    if decision not in ["APPROVED", "APPROVED_WITH_RESERVATIONS", "REVISIONS_REQUIRED"]:
        return jsonify({"ok": False, "error": "Invalid decision value"}), 400

    checklist_status = _compute_checklist_status()
    if not checklist_status["is_ready_for_signoff"] and not force:
        return (
            jsonify(
                {
                    "ok": False,
                    "error": (
                        f"Checklist incomplete ({checklist_status['readiness_pct']}% complete). "
                        "Complete all checks or pass 'force': true to override."
                    ),
                    "checklist_status": checklist_status,
                }
            ),
            400,
        )

    signoff = {
        "signed_off": True,
        "decision": decision,
        "reviewer_name": name,
        "reviewer_role": role,
        "notes": notes,
        "timestamp": _utc_now(),
        "checklist_snapshot": checklist_status,
        "diff_snapshot": _calculate_override_diff(),
    }

    _save_signoff(signoff)
    _append_audit_event("review_signed_off", signoff)

    return jsonify({"ok": True, "signoff": signoff})


@app.route("/api/export/review-packet", methods=["GET"])
def api_export_review_packet():
    """Export complete, self-contained JSON review packet bundle."""
    packet = _generate_review_packet()
    packet_json = json.dumps(packet, indent=2, ensure_ascii=False)

    response = make_response(packet_json)
    response.headers["Content-Type"] = "application/json; charset=utf-8"
    response.headers["Content-Disposition"] = "attachment; filename=review_packet.json"
    return response


# ──────────────────────────────────────────────
# Override API Routes
# ──────────────────────────────────────────────


@app.route("/api/overrides", methods=["GET"])
def api_get_overrides():
    """Get all current overrides."""
    return jsonify({"ok": True, "overrides": _overrides, "count": len(_overrides)})


@app.route("/api/override", methods=["POST"])
def api_add_override():
    """Add or update a single override."""
    global _overrides
    if _annotation_writes_blocked():
        return (
            jsonify(
                {
                    "ok": False,
                    "error": (
                        "Override writes blocked: this directory looks like a benchmark "
                        "run (checkpoints/ without data.annotated.h5ad). Observability is "
                        "read-only here; open an annotation output for review overrides."
                    ),
                    "prediction_mutation_allowed": False,
                }
            ),
            403,
        )
    data = request.get_json(force=True, silent=True)
    if not data:
        return jsonify({"ok": False, "error": "Invalid JSON body"}), 400

    cluster = str(data.get("cluster", "")).strip()
    new_type = data.get("new_type", "").strip()
    reason = data.get("reason", "").strip()

    if not cluster:
        return jsonify({"ok": False, "error": "Missing 'cluster' field"}), 400
    if not new_type:
        return jsonify({"ok": False, "error": "Missing 'new_type' field"}), 400

    _overrides[cluster] = {
        "new_type": new_type,
        "reason": reason,
        "timestamp": _utc_now(),
    }

    # Persist to disk immediately
    _save_overrides_to_disk()
    _append_audit_event(
        "override_saved",
        {"cluster": cluster, "new_type": new_type, "reason": reason},
    )

    return jsonify(
        {
            "ok": True,
            "cluster": cluster,
            "new_type": new_type,
            "total_overrides": len(_overrides),
        }
    )


@app.route("/api/override/<cluster_id>", methods=["DELETE"])
def api_delete_override(cluster_id):
    """Remove a single override."""
    global _overrides
    cluster_id = str(cluster_id)
    if cluster_id in _overrides:
        removed = _overrides[cluster_id]
        del _overrides[cluster_id]
        _save_overrides_to_disk()
        _append_audit_event("override_deleted", {"cluster": cluster_id, "override": removed})
        return jsonify({"ok": True, "removed": cluster_id})
    return jsonify({"ok": False, "error": "Override not found"}), 404


@app.route("/api/overrides/apply", methods=["POST"])
def api_apply_overrides():
    """Apply all overrides to the .h5ad file (audit-logged; marks artifacts stale)."""
    global _overrides
    if _annotation_writes_blocked():
        return (
            jsonify(
                {
                    "ok": False,
                    "error": (
                        "Apply blocked on benchmark-run directories. Observability never "
                        "mutates fold predictions. Use an annotation output directory."
                    ),
                    "prediction_mutation_allowed": False,
                }
            ),
            403,
        )
    if not _overrides:
        return jsonify({"ok": False, "error": "No overrides to apply"}), 400

    try:
        result = _apply_overrides_to_h5ad()
        status = _mark_artifacts_stale_after_apply(result)
        _append_audit_event("overrides_applied", {"result": result, "artifact_status": status})
        _overrides = {}
        _save_overrides_to_disk()
        return jsonify({"ok": True, "result": result, "artifact_status": status})
    except FileNotFoundError as e:
        return jsonify({"ok": False, "error": str(e)}), 404
    except Exception as e:
        return jsonify({"ok": False, "error": f"Apply failed: {e}"}), 500


@app.route("/api/overrides/clear", methods=["POST"])
def api_clear_overrides():
    """Clear all overrides (does NOT revert .h5ad changes)."""
    global _overrides
    count = len(_overrides)
    _overrides = {}
    _save_overrides_to_disk()
    _append_audit_event("overrides_cleared", {"cleared": count})
    return jsonify({"ok": True, "cleared": count})


@app.route("/api/audit", methods=["GET"])
def api_audit_log():
    """Return the tail of the append-only Web Review audit log."""
    limit = request.args.get("limit", default=20, type=int)
    limit = max(1, min(limit, 200))
    return jsonify({"ok": True, "events": _audit_tail(limit), "path": str(_audit_log_path())})


@app.route("/api/artifact-status", methods=["GET"])
def api_artifact_status():
    """Return whether derived artifacts are current or stale after review edits."""
    return jsonify({"ok": True, "artifact_status": _load_artifact_status()})


@app.route("/api/novelty/verify/<cluster_id>", methods=["GET"])
def api_novelty_verify(cluster_id: str):
    """Run 5-gate audit protocol on a focus cluster for Web Cockpit."""
    try:
        from .data_adapter import detect_species, load_marker_atlas
        from .novelty_verification import verify_novelty_candidate

        adata, evidence_df = _load_data()
        focus_row = evidence_df[evidence_df["cluster"].astype(str) == str(cluster_id)]
        row_dict = focus_row.iloc[0].to_dict() if not focus_row.empty else {}
        atlas = load_marker_atlas(detect_species(adata))
        tissue = str(adata.obs["tissue"].iloc[0]) if "tissue" in adata.obs else "general"

        packet = verify_novelty_candidate(adata, "cluster", cluster_id, row_dict, atlas, tissue)
        return jsonify({"ok": True, "packet": packet})
    except Exception as exc:
        return jsonify({"ok": False, "error": str(exc)}), 400


@app.route("/api/novelty/adjudicate", methods=["POST"])
def api_novelty_adjudicate():
    """Post human expert adjudication verdict."""
    if _annotation_writes_blocked():
        return (
            jsonify(
                {
                    "ok": False,
                    "error": "Write blocked on benchmark-run directory (read-only observability).",
                    "prediction_mutation_allowed": False,
                }
            ),
            403,
        )
    data = request.json or {}
    cluster = data.get("cluster")
    verdict = data.get("verdict")
    reviewer = data.get("reviewer", "web_inspector")
    notes = data.get("notes")
    pmid = data.get("pmid")

    if not cluster or not verdict:
        return jsonify({"ok": False, "error": "cluster and verdict are required"}), 400

    try:
        from .novelty_verification import log_novelty_adjudication

        entry = log_novelty_adjudication(
            _output_dir, cluster, verdict, reviewer, notes=notes, pmid=pmid
        )
        _append_audit_event("novelty_adjudicated", entry)
        return jsonify({"ok": True, "entry": entry})
    except Exception as exc:
        return jsonify({"ok": False, "error": str(exc)}), 400


def run_inspector(
    output_dir: str | Path,
    host: str = "127.0.0.1",
    port: int = 8765,
    run_dir: str | Path | None = None,
):
    """Start the Review Cockpit.

    Parameters
    ----------
    output_dir:
        Annotation review directory (evidence / h5ad / audit log).
    run_dir:
        Optional benchmark-run directory containing ``checkpoints/``. When
        omitted, observability reads from ``output_dir``. Observability never
        writes fold files or prediction tables.
    """
    global _output_dir, _run_dir, _adata_cache, _evidence_cache, _overrides
    _output_dir = Path(output_dir).resolve()
    _run_dir = Path(run_dir).resolve() if run_dir is not None else _output_dir
    _adata_cache = None
    _evidence_cache = None
    _overrides = {}
    if _output_dir.is_dir():
        _load_overrides_from_disk()

    print(f"  Output dir: {_output_dir}")
    print(f"  Observability run dir (read-only): {_run_dir}")
    print(f"  URL: http://{host}:{port}")
    print("  Observability never mutates predictions; overrides use audit log.")
    print("  Press Ctrl+C to stop\n")

    app.run(host=host, port=port, debug=False)
