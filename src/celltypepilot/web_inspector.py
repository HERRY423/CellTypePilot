"""CellTypePilot — Web Inspector: lightweight interactive review panel.

A minimal Flask app that serves an interactive dashboard for reviewing
and correcting cell-type annotations. Designed for local use — runs
on localhost, no authentication, no deployment complexity.

Override workflow:
1. User corrects a cell-type label in the browser
2. POST /api/override — stores override in memory + server-side JSON
3. POST /api/overrides/apply — writes overrides back to .h5ad file
   and regenerates figures/reports
"""

from __future__ import annotations

import importlib.resources
import json
from datetime import datetime, timezone
from pathlib import Path

from flask import Flask, jsonify, render_template, request

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
_adata_cache = None
_evidence_cache = None

# Server-side override store (persists across page reloads)
_overrides: dict = {}  # {cluster_id: {new_type, reason, timestamp}}

AUDIT_LOG_FILENAME = "annotation_audit_log.jsonl"
ARTIFACT_STATUS_FILENAME = "artifact_status.json"
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


@app.route("/")
def dashboard():
    """Main dashboard."""
    adata, evidence = _load_data()

    # Build stats
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

    # Build annotation rows
    annotations = []
    evidence_dict = {}

    if not evidence.empty:
        for _, row in evidence.iterrows():
            cluster = row.get("cluster", "")
            flags = row.get("critic_flags", "PASS")
            if flags != "PASS":
                stats["flagged"] += 1

            annotations.append(
                {
                    "cluster": cluster,
                    "cell_type": row.get("cell_type", "Unknown"),
                    "cl_id": row.get("cl_id", ""),
                    "score": float(row.get("combined_score", 0)),
                    "confidence": row.get("critic_confidence", "unknown"),
                    "flags": flags,
                    "n_cells": int(row.get("n_cells", 0)),
                }
            )

            evidence_dict[str(cluster)] = {
                "cell_type": row.get("cell_type", ""),
                "combined_score": float(row.get("combined_score", 0)),
                "pct_overlap": float(row.get("pct_overlap", 0)),
                "top_markers": [],
                "critic_notes": row.get("critic_notes", ""),
            }

    return render_template(
        "dashboard.html",
        version=__version__,
        stats=stats,
        annotations=annotations,
        evidence_json=json.dumps(evidence_dict),
        artifact_status=_load_artifact_status(),
        audit_tail=_audit_tail(),
    )


@app.route("/api/evidence/<cluster_id>")
def api_evidence(cluster_id):
    """Get detailed evidence for a cluster."""
    _, evidence = _load_data()
    if evidence.empty:
        return jsonify({})
    row = evidence[evidence["cluster"] == cluster_id]
    if row.empty:
        return jsonify({})
    return jsonify(row.iloc[0].to_dict())


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
# Override API Routes
# ──────────────────────────────────────────────


@app.route("/api/overrides", methods=["GET"])
def api_get_overrides():
    """Get all current overrides."""
    return jsonify({"ok": True, "overrides": _overrides, "count": len(_overrides)})


@app.route("/api/override", methods=["POST"])
def api_add_override():
    """Add or update a single override.

    Expected JSON body:
        {"cluster": "5", "new_type": "CD4 naive T cell", "reason": "Manual review"}
    """
    global _overrides
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
    if cluster_id in _overrides:
        removed = _overrides[cluster_id]
        del _overrides[cluster_id]
        _save_overrides_to_disk()
        _append_audit_event("override_deleted", {"cluster": cluster_id, "override": removed})
        return jsonify({"ok": True, "removed": cluster_id})
    return jsonify({"ok": False, "error": "Override not found"}), 404


@app.route("/api/overrides/apply", methods=["POST"])
def api_apply_overrides():
    """Apply all overrides to the .h5ad file.

    Creates a timestamped backup of the original file before modifying.
    Returns a summary of applied/skipped overrides.
    """
    global _overrides
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


def run_inspector(output_dir: str | Path, host: str = "127.0.0.1", port: int = 8765):
    """Launch the web inspector.

    Args:
        output_dir: Path to CellTypePilot output directory
        host: Host to bind to (default: localhost)
        port: Port to listen on (default: 8765)
    """
    global _output_dir
    _output_dir = Path(output_dir)

    if not (_output_dir / "data.annotated.h5ad").exists():
        raise FileNotFoundError(f"No annotated data found in {_output_dir}")

    print("CellTypePilot Web Inspector")
    print(f"  Output dir: {_output_dir}")
    print(f"  URL: http://{host}:{port}")
    print("  Press Ctrl+C to stop\n")

    app.run(host=host, port=port, debug=False)
