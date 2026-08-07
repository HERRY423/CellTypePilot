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

import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional

from flask import Flask, render_template_string, jsonify, request

from . import __version__

app = Flask(__name__)

# Global state — set by `celltypepilot inspect-web`
_output_dir: Optional[Path] = None
_adata_cache = None
_evidence_cache = None

# Server-side override store (persists across page reloads)
_overrides: dict = {}  # {cluster_id: {new_type, reason, timestamp}}


def _load_data():
    """Lazy-load annotation data from output directory."""
    global _adata_cache, _evidence_cache
    if _adata_cache is not None:
        return _adata_cache, _evidence_cache

    import scanpy as sc
    import pandas as pd

    adata_path = _output_dir / "data.annotated.h5ad"
    evidence_path = _output_dir / "evidence_table.csv"

    if not adata_path.exists():
        raise FileNotFoundError(f"No annotated data found at {adata_path}")

    _adata_cache = sc.read_h5ad(adata_path)

    if evidence_path.exists():
        _evidence_cache = pd.read_csv(evidence_path)
    else:
        _evidence_cache = pd.DataFrame()

    # Load any existing overrides from disk
    _load_overrides_from_disk()

    return _adata_cache, _evidence_cache


def _overrides_path() -> Path:
    """Path to the server-side overrides file."""
    return _output_dir / "annotation_overrides.json"


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
    opath.write_text(json.dumps(_overrides, indent=2), encoding="utf-8")


def _apply_overrides_to_h5ad() -> dict:
    """Apply all overrides to the .h5ad file and back up the original.

    Returns:
        Summary dict with counts and details.
    """
    import scanpy as sc

    adata_path = _output_dir / "data.annotated.h5ad"
    if not adata_path.exists():
        raise FileNotFoundError(f"No annotated data at {adata_path}")

    # Backup original
    backup_path = _output_dir / f"data.annotated.backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}.h5ad"
    shutil.copy2(adata_path, backup_path)

    # Reload fresh
    adata = sc.read_h5ad(adata_path)
    obs = adata.obs

    applied = 0
    skipped = 0
    details = []

    for cluster_id, override in _overrides.items():
        new_type = override.get("new_type", "")
        reason = override.get("reason", "")
        if not new_type:
            skipped += 1
            continue

        # Find cells in this cluster
        cluster_col = "ctp_cl_id" if "ctp_cl_id" in obs.columns else None
        if cluster_col is None:
            # Try numeric cluster column
            for col in obs.columns:
                if "cluster" in col.lower() or "cl_id" in col.lower():
                    cluster_col = col
                    break

        if cluster_col is None:
            skipped += 1
            details.append({"cluster": cluster_id, "status": "error", "reason": "No cluster column found"})
            continue

        mask = obs[cluster_col].astype(str) == str(cluster_id)
        n_cells = mask.sum()

        if n_cells == 0:
            skipped += 1
            details.append({"cluster": cluster_id, "status": "skipped", "reason": "No cells found"})
            continue

        old_type = obs.loc[mask, "ctp_cell_type"].iloc[0] if "ctp_cell_type" in obs.columns else "Unknown"

        # Apply override
        if "ctp_cell_type" in obs.columns:
            adata.obs.loc[mask, "ctp_cell_type"] = new_type
        if "ctp_override_reason" not in obs.columns:
            adata.obs["ctp_override_reason"] = ""
        adata.obs.loc[mask, "ctp_override_reason"] = reason
        if "ctp_overridden" not in obs.columns:
            adata.obs["ctp_overridden"] = False
        adata.obs.loc[mask, "ctp_overridden"] = True

        applied += 1
        details.append({
            "cluster": cluster_id,
            "old_type": old_type,
            "new_type": new_type,
            "n_cells": int(n_cells),
            "reason": reason,
            "status": "applied",
        })

    # Save updated h5ad
    adata.write(adata_path)

    # Clear cache so next load gets fresh data
    global _adata_cache
    _adata_cache = None

    return {
        "applied": applied,
        "skipped": skipped,
        "total": len(_overrides),
        "backup": str(backup_path),
        "details": details,
    }


# ──────────────────────────────────────────────
# HTML Template
# ──────────────────────────────────────────────

DASHBOARD_HTML = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CellTypePilot — Web Inspector</title>
    <style>
        :root {
            --bg: #fafafa;
            --card: #ffffff;
            --border: #e0e0e0;
            --text: #2c3e50;
            --muted: #7f8c8d;
            --accent: #3498db;
            --success: #27ae60;
            --warning: #f39c12;
            --danger: #e74c3c;
            --info: #9b59b6;
        }
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg);
            color: var(--text);
            line-height: 1.6;
        }
        .container { max-width: 1400px; margin: 0 auto; padding: 20px; }
        header {
            background: var(--card);
            border-bottom: 1px solid var(--border);
            padding: 16px 24px;
            margin-bottom: 24px;
        }
        header h1 { font-size: 1.5rem; font-weight: 600; }
        header .subtitle { color: var(--muted); font-size: 0.9rem; }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 16px;
            margin-bottom: 24px;
        }
        .stat-card {
            background: var(--card);
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 16px;
        }
        .stat-card .label { color: var(--muted); font-size: 0.8rem; text-transform: uppercase; }
        .stat-card .value { font-size: 2rem; font-weight: 700; margin-top: 4px; }
        .stat-card .value.high { color: var(--success); }
        .stat-card .value.medium { color: var(--warning); }
        .stat-card .value.low { color: var(--danger); }
        .stat-card .value.needs_review { color: var(--info); }
        .panel {
            background: var(--card);
            border: 1px solid var(--border);
            border-radius: 8px;
            margin-bottom: 24px;
        }
        .panel-header {
            padding: 12px 16px;
            border-bottom: 1px solid var(--border);
            font-weight: 600;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        .panel-body { padding: 16px; overflow-x: auto; }
        table { width: 100%; border-collapse: collapse; font-size: 0.9rem; }
        th, td { padding: 10px 12px; text-align: left; border-bottom: 1px solid var(--border); }
        th { font-weight: 600; color: var(--muted); font-size: 0.8rem; text-transform: uppercase; }
        tr:hover { background: #f8f9fa; }
        .badge {
            display: inline-block;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 0.75rem;
            font-weight: 600;
        }
        .badge-pass { background: #d4edda; color: #155724; }
        .badge-low_evidence { background: #fff3cd; color: #856404; }
        .badge-partial_evidence { background: #fff3cd; color: #856404; }
        .badge-neg_marker_conflict { background: #f8d7da; color: #721c24; }
        .badge-possible_doublet { background: #d1ecf1; color: #0c5460; }
        .badge-high { background: #d4edda; color: #155724; }
        .badge-medium { background: #fff3cd; color: #856404; }
        .badge-low { background: #f8d7da; color: #721c24; }
        .badge-needs_review { background: #d1ecf1; color: #0c5460; }
        .btn {
            display: inline-block;
            padding: 6px 12px;
            border: 1px solid var(--border);
            border-radius: 4px;
            background: var(--card);
            cursor: pointer;
            font-size: 0.85rem;
        }
        .btn:hover { background: #f0f0f0; }
        .btn-primary { background: var(--accent); color: white; border-color: var(--accent); }
        .btn-primary:hover { background: #2980b9; }
        .modal-overlay {
            display: none;
            position: fixed;
            top: 0; left: 0; right: 0; bottom: 0;
            background: rgba(0,0,0,0.5);
            z-index: 1000;
        }
        .modal {
            background: var(--card);
            border-radius: 8px;
            max-width: 600px;
            margin: 50px auto;
            padding: 24px;
        }
        .modal h3 { margin-bottom: 16px; }
        .form-group { margin-bottom: 12px; }
        .form-group label { display: block; font-weight: 600; margin-bottom: 4px; }
        .form-group input, .form-group select {
            width: 100%;
            padding: 8px;
            border: 1px solid var(--border);
            border-radius: 4px;
        }
        .evidence-list { list-style: none; padding: 0; }
        .evidence-list li { padding: 4px 0; font-size: 0.9rem; }
        .evidence-list .gene { font-weight: 600; color: var(--accent); }
        .filter-bar {
            display: flex;
            gap: 12px;
            margin-bottom: 16px;
            flex-wrap: wrap;
        }
        .filter-bar select, .filter-bar input {
            padding: 6px 10px;
            border: 1px solid var(--border);
            border-radius: 4px;
        }
        .save-bar {
            display: none;
            position: fixed;
            bottom: 0;
            left: 0;
            right: 0;
            background: var(--card);
            border-top: 2px solid var(--accent);
            padding: 12px 24px;
            z-index: 999;
            align-items: center;
            justify-content: space-between;
        }
        .save-bar .info { font-weight: 600; }
        .save-bar .actions { display: flex; gap: 8px; }
    </style>
</head>
<body>
    <header>
        <h1>CellTypePilot Web Inspector</h1>
        <div class="subtitle">v{{ version }} — Interactive annotation review panel</div>
    </header>

    <div class="container">
        <!-- Stats -->
        <div class="stats-grid">
            <div class="stat-card">
                <div class="label">Total Clusters</div>
                <div class="value">{{ stats.total_clusters }}</div>
            </div>
            <div class="stat-card">
                <div class="label">High Confidence</div>
                <div class="value high">{{ stats.high }}</div>
            </div>
            <div class="stat-card">
                <div class="label">Medium Confidence</div>
                <div class="value medium">{{ stats.medium }}</div>
            </div>
            <div class="stat-card">
                <div class="label">Low / Needs Review</div>
                <div class="value low">{{ stats.low + stats.needs_review }}</div>
            </div>
            <div class="stat-card">
                <div class="label">Critic Flags</div>
                <div class="value" style="color: var(--info);">{{ stats.flagged }}</div>
            </div>
        </div>

        <!-- Annotation Table -->
        <div class="panel">
            <div class="panel-header">
                <span>Annotation Results</span>
                <button class="btn" onclick="exportOverrides()">Export JSON</button>
                <button class="btn btn-primary" onclick="applyOverrides()" style="margin-left:8px">Apply to .h5ad</button>
            </div>
            <div class="panel-body">
                <div class="filter-bar">
                    <select id="confidenceFilter" onchange="filterTable()">
                        <option value="">All confidence levels</option>
                        <option value="high">High</option>
                        <option value="medium">Medium</option>
                        <option value="low">Low</option>
                        <option value="needs_review">Needs Review</option>
                    </select>
                    <select id="flagFilter" onchange="filterTable()">
                        <option value="">All critic flags</option>
                        <option value="PASS">PASS</option>
                        <option value="LOW_EVIDENCE">LOW_EVIDENCE</option>
                        <option value="PARTIAL_EVIDENCE">PARTIAL_EVIDENCE</option>
                        <option value="NEG_MARKER_CONFLICT">NEG_MARKER_CONFLICT</option>
                        <option value="POSSIBLE_DOUBLET">POSSIBLE_DOUBLET</option>
                    </select>
                    <input type="text" id="searchBox" placeholder="Search cell type..." oninput="filterTable()">
                </div>
                <table id="annotationTable">
                    <thead>
                        <tr>
                            <th>Cluster</th>
                            <th>Cell Type</th>
                            <th>CL ID</th>
                            <th>Score</th>
                            <th>Confidence</th>
                            <th>Critic Flags</th>
                            <th>Cells</th>
                            <th>Actions</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for row in annotations %}
                        <tr data-confidence="{{ row.confidence }}" data-flag="{{ row.flags }}" data-celltype="{{ row.cell_type|lower }}">
                            <td>{{ row.cluster }}</td>
                            <td><strong>{{ row.cell_type }}</strong></td>
                            <td>{{ row.cl_id }}</td>
                            <td>{{ "%.3f"|format(row.score) }}</td>
                            <td><span class="badge badge-{{ row.confidence }}">{{ row.confidence }}</span></td>
                            <td><span class="badge badge-{{ row.flags|lower }}">{{ row.flags }}</span></td>
                            <td>{{ row.n_cells }}</td>
                            <td>
                                <button class="btn" onclick="showEvidence({{ row.cluster }})">Evidence</button>
                                <button class="btn" onclick="showOverride({{ row.cluster }}, '{{ row.cell_type }}')">Override</button>
                            </td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </div>
        </div>
    </div>

    <!-- Save Bar -->
    <div class="save-bar" id="saveBar">
        <div class="info"><span id="overrideCount">0</span> override(s) pending</div>
        <div class="actions">
            <button class="btn" onclick="exportOverrides()">Export JSON</button>
            <button class="btn btn-primary" onclick="applyOverrides()">Apply to .h5ad</button>
            <button class="btn" onclick="hideSaveBar()" style="color:var(--muted)">Dismiss</button>
        </div>
    </div>

    <!-- Evidence Modal -->
    <div class="modal-overlay" id="evidenceModal">
        <div class="modal">
            <h3>Cluster <span id="evidenceCluster"></span> — Evidence</h3>
            <div id="evidenceContent"></div>
            <br>
            <button class="btn" onclick="closeModal('evidenceModal')">Close</button>
        </div>
    </div>

    <!-- Override Modal -->
    <div class="modal-overlay" id="overrideModal">
        <div class="modal">
            <h3>Override Annotation for Cluster <span id="overrideCluster"></span></h3>
            <div class="form-group">
                <label>Current: <span id="currentAnnotation"></span></label>
            </div>
            <div class="form-group">
                <label>New Cell Type</label>
                <input type="text" id="newCellType" placeholder="e.g., CD4 naive T cell">
            </div>
            <div class="form-group">
                <label>Reason</label>
                <input type="text" id="overrideReason" placeholder="e.g., Manual review based on marker expression">
            </div>
            <br>
            <button class="btn btn-primary" onclick="saveOverride()">Save Override</button>
            <button class="btn" onclick="closeModal('overrideModal')">Cancel</button>
        </div>
    </div>

    <script>
        const evidenceData = {{ evidence_json|safe }};
        const overrides = {};

        function filterTable() {
            const confFilter = document.getElementById('confidenceFilter').value;
            const flagFilter = document.getElementById('flagFilter').value;
            const search = document.getElementById('searchBox').value.toLowerCase();
            const rows = document.querySelectorAll('#annotationTable tbody tr');

            rows.forEach(row => {
                const conf = row.dataset.confidence;
                const flag = row.dataset.flag;
                const ct = row.dataset.celltype;
                const show = (!confFilter || conf === confFilter)
                          && (!flagFilter || flag === flagFilter)
                          && (!search || ct.includes(search));
                row.style.display = show ? '' : 'none';
            });
        }

        function showEvidence(cluster) {
            document.getElementById('evidenceCluster').textContent = cluster;
            const ev = evidenceData[cluster] || {};
            let html = '<p><strong>Cell Type:</strong> ' + (ev.cell_type || 'N/A') + '</p>';
            html += '<p><strong>Score:</strong> ' + (ev.combined_score || 0).toFixed(3) + '</p>';
            html += '<p><strong>Marker Overlap:</strong> ' + ((ev.pct_overlap || 0) * 100).toFixed(0) + '%</p>';
            html += '<h4 style="margin-top:12px">Top Markers</h4><ul class="evidence-list">';
            (ev.top_markers || []).forEach(m => {
                html += '<li><span class="gene">' + m.gene + '</span> (pct=' + (m.pct * 100).toFixed(0) + '%, FC=' + m.fc.toFixed(2) + ')</li>';
            });
            html += '</ul>';
            if (ev.critic_notes) {
                html += '<p style="margin-top:12px"><strong>Critic Notes:</strong> ' + ev.critic_notes + '</p>';
            }
            document.getElementById('evidenceContent').innerHTML = html;
            document.getElementById('evidenceModal').style.display = 'block';
        }

        function showOverride(cluster, currentType) {
            document.getElementById('overrideCluster').textContent = cluster;
            document.getElementById('currentAnnotation').textContent = currentType;
            document.getElementById('newCellType').value = overrides[cluster]?.new_type || '';
            document.getElementById('overrideReason').value = overrides[cluster]?.reason || '';
            document.getElementById('overrideModal').style.display = 'block';
        }

        function saveOverride() {
            const cluster = document.getElementById('overrideCluster').textContent;
            const newType = document.getElementById('newCellType').value;
            const reason = document.getElementById('overrideReason').value;
            if (!newType) { alert('Please enter a new cell type'); return; }

            // Send to server (not just browser memory)
            fetch('/api/override', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({cluster: cluster, new_type: newType, reason: reason})
            })
            .then(r => r.json())
            .then(data => {
                if (data.ok) {
                    overrides[cluster] = { new_type: newType, reason: reason };
                    closeModal('overrideModal');
                    showSaveBar();
                } else {
                    alert('Error: ' + (data.error || 'Unknown error'));
                }
            })
            .catch(err => alert('Network error: ' + err));
        }

        function showSaveBar() {
            const bar = document.getElementById('saveBar');
            const count = Object.keys(overrides).length;
            bar.style.display = 'flex';
            document.getElementById('overrideCount').textContent = count;
        }

        function hideSaveBar() {
            document.getElementById('saveBar').style.display = 'none';
        }

        function applyOverrides() {
            if (Object.keys(overrides).length === 0) {
                alert('No overrides to apply.');
                return;
            }
            if (!confirm('Apply ' + Object.keys(overrides).length + ' override(s) to .h5ad?\nThis will modify the data file and create a backup.')) return;

            fetch('/api/overrides/apply', { method: 'POST' })
            .then(r => r.json())
            .then(data => {
                if (data.ok) {
                    alert('Applied ' + data.result.applied + ' override(s).\nBackup: ' + data.result.backup);
                    location.reload();  // Reload to show updated data
                } else {
                    alert('Error: ' + (data.error || 'Apply failed'));
                }
            })
            .catch(err => alert('Network error: ' + err));
        }

        function exportOverrides() {
            if (Object.keys(overrides).length === 0) {
                alert('No overrides to export.');
                return;
            }
            // Download from server (authoritative copy)
            fetch('/api/overrides')
            .then(r => r.json())
            .then(data => {
                const blob = new Blob([JSON.stringify(data.overrides, null, 2)], { type: 'application/json' });
                const url = URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = 'annotation_overrides.json';
                a.click();
                URL.revokeObjectURL(url);
            });
        }

        // Close modals on overlay click
        document.querySelectorAll('.modal-overlay').forEach(el => {
            el.addEventListener('click', e => { if (e.target === el) el.style.display = 'none'; });
        });

        function closeModal(id) {
            document.getElementById(id).style.display = 'none';
        }

        // Load existing overrides from server on page load
        fetch('/api/overrides').then(r => r.json()).then(data => {
            Object.assign(overrides, data.overrides || {});
            if (Object.keys(overrides).length > 0) showSaveBar();
        });
    </script>
</body>
</html>
"""


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
        "high": int((obs["ctp_confidence"] == "high").sum()) if "ctp_confidence" in obs.columns else 0,
        "medium": int((obs["ctp_confidence"] == "medium").sum()) if "ctp_confidence" in obs.columns else 0,
        "low": int((obs["ctp_confidence"] == "low").sum()) if "ctp_confidence" in obs.columns else 0,
        "needs_review": int((obs["ctp_confidence"] == "needs_review").sum()) if "ctp_confidence" in obs.columns else 0,
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

            annotations.append({
                "cluster": cluster,
                "cell_type": row.get("cell_type", "Unknown"),
                "cl_id": row.get("cl_id", ""),
                "score": float(row.get("combined_score", 0)),
                "confidence": row.get("critic_confidence", "unknown"),
                "flags": flags,
                "n_cells": int(row.get("n_cells", 0)),
            })

            evidence_dict[str(cluster)] = {
                "cell_type": row.get("cell_type", ""),
                "combined_score": float(row.get("combined_score", 0)),
                "pct_overlap": float(row.get("pct_overlap", 0)),
                "top_markers": [],
                "critic_notes": row.get("critic_notes", ""),
            }

    return render_template_string(
        DASHBOARD_HTML,
        version=__version__,
        stats=stats,
        annotations=annotations,
        evidence_json=json.dumps(evidence_dict),
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
        "timestamp": datetime.now().isoformat(),
    }

    # Persist to disk immediately
    _save_overrides_to_disk()

    return jsonify({
        "ok": True,
        "cluster": cluster,
        "new_type": new_type,
        "total_overrides": len(_overrides),
    })


@app.route("/api/override/<cluster_id>", methods=["DELETE"])
def api_delete_override(cluster_id):
    """Remove a single override."""
    global _overrides
    if cluster_id in _overrides:
        del _overrides[cluster_id]
        _save_overrides_to_disk()
        return jsonify({"ok": True, "removed": cluster_id})
    return jsonify({"ok": False, "error": "Override not found"}), 404


@app.route("/api/overrides/apply", methods=["POST"])
def api_apply_overrides():
    """Apply all overrides to the .h5ad file.

    Creates a timestamped backup of the original file before modifying.
    Returns a summary of applied/skipped overrides.
    """
    if not _overrides:
        return jsonify({"ok": False, "error": "No overrides to apply"}), 400

    try:
        result = _apply_overrides_to_h5ad()
        return jsonify({"ok": True, "result": result})
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
    return jsonify({"ok": True, "cleared": count})


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

    print(f"CellTypePilot Web Inspector")
    print(f"  Output dir: {_output_dir}")
    print(f"  URL: http://{host}:{port}")
    print(f"  Press Ctrl+C to stop\n")

    app.run(host=host, port=port, debug=False)
