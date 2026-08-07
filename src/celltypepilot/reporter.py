"""Report generator — HTML report and evidence table output."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional
from datetime import datetime, timezone

import pandas as pd

from . import __version__, MKG_VERSION
from .constants import OUTPUT_EVIDENCE, OUTPUT_REPORT, OUTPUT_FIGURES_DIR


def save_evidence_table(
    critic_results: pd.DataFrame,
    output_dir: str | Path,
) -> Path:
    """Save the evidence table as CSV."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / OUTPUT_EVIDENCE
    critic_results.to_csv(path, index=False)
    return path


def generate_html_report(
    annotations: pd.DataFrame,
    critic_results: pd.DataFrame,
    critic_summary: dict,
    manifest: dict,
    figure_paths: list[str],
    output_dir: str | Path,
) -> Path:
    """Generate a comprehensive HTML report."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / OUTPUT_REPORT

    # Build HTML
    html_parts = [
        _html_header(),
        _html_overview(manifest, critic_summary),
        _html_annotation_table(critic_results),
        _html_figures(figure_paths, output_dir),
        _html_critic_details(critic_results),
        _html_footer(),
    ]

    html = "\n".join(html_parts)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)

    return report_path


def _html_header() -> str:
    return """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>CellTypePilot Annotation Report</title>
<style>
  :root {
    --bg: #ffffff; --fg: #1a1a2e; --accent: #0072B2;
    --success: #009E73; --warning: #E69F00; --danger: #D55E00;
    --muted: #6c757d; --border: #dee2e6; --card-bg: #f8f9fa;
  }
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    color: var(--fg); background: var(--bg); line-height: 1.6;
    max-width: 1200px; margin: 0 auto; padding: 2rem;
  }
  h1 { color: var(--accent); margin-bottom: 0.5rem; font-size: 1.8rem; }
  h2 { color: var(--accent); margin: 2rem 0 1rem; font-size: 1.3rem;
       border-bottom: 2px solid var(--accent); padding-bottom: 0.3rem; }
  h3 { margin: 1.5rem 0 0.5rem; font-size: 1.1rem; }
  .subtitle { color: var(--muted); font-size: 0.9rem; margin-bottom: 2rem; }
  table { width: 100%; border-collapse: collapse; margin: 1rem 0; font-size: 0.85rem; }
  th { background: var(--accent); color: white; padding: 0.6rem; text-align: left; }
  td { padding: 0.5rem; border-bottom: 1px solid var(--border); }
  tr:nth-child(even) { background: var(--card-bg); }
  .badge {
    display: inline-block; padding: 0.15rem 0.5rem; border-radius: 3px;
    font-size: 0.75rem; font-weight: 600; color: white;
  }
  .badge-high { background: var(--success); }
  .badge-medium { background: var(--warning); }
  .badge-low { background: var(--danger); }
  .badge-review { background: #CC0000; }
  .badge-pass { background: var(--success); }
  .badge-flag { background: var(--danger); }
  .card {
    background: var(--card-bg); border: 1px solid var(--border);
    border-radius: 6px; padding: 1rem; margin: 0.5rem 0;
  }
  .card-header { font-weight: 600; margin-bottom: 0.5rem; }
  .stats-grid {
    display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
    gap: 1rem; margin: 1rem 0;
  }
  .stat-card {
    background: var(--card-bg); border: 1px solid var(--border);
    border-radius: 6px; padding: 1rem; text-align: center;
  }
  .stat-value { font-size: 1.8rem; font-weight: 700; color: var(--accent); }
  .stat-label { font-size: 0.8rem; color: var(--muted); }
  .figure-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 1rem; }
  .figure-card { border: 1px solid var(--border); border-radius: 6px; overflow: hidden; }
  .figure-card img { width: 100%; display: block; }
  .figure-card .caption { padding: 0.5rem; font-size: 0.85rem; text-align: center; color: var(--muted); }
  .warning-box {
    background: #fff3cd; border: 1px solid #ffc107; border-radius: 6px;
    padding: 0.8rem; margin: 0.5rem 0; font-size: 0.85rem;
  }
  .evidence-text { font-size: 0.8rem; color: var(--muted); margin-top: 0.3rem; }
</style>
</head>
<body>
<h1>CellTypePilot Annotation Report</h1>"""


def _html_overview(manifest: dict, critic_summary: dict) -> str:
    params = manifest.get("parameters", {})
    total = critic_summary.get("total_clusters", 0)
    passed = critic_summary.get("pass", 0)
    flagged = critic_summary.get("flagged", 0)
    conf_dist = critic_summary.get("confidence_distribution", {})

    html = f"""
<p class="subtitle">
  Generated: {manifest.get('timestamp', 'N/A')} |
  CellTypePilot v{manifest.get('celltypepilot_version', '?')} |
  MKG: {manifest.get('mkg_version', '?')}
</p>

<h2>Overview</h2>
<div class="stats-grid">
  <div class="stat-card">
    <div class="stat-value">{manifest.get('input', {}).get('path', 'N/A').split('/')[-1].split(chr(92))[-1][:30]}</div>
    <div class="stat-label">Input File</div>
  </div>
  <div class="stat-card">
    <div class="stat-value">{total}</div>
    <div class="stat-label">Total Clusters</div>
  </div>
  <div class="stat-card">
    <div class="stat-value" style="color: var(--success)">{passed}</div>
    <div class="stat-label">Passed Critic</div>
  </div>
  <div class="stat-card">
    <div class="stat-value" style="color: var(--danger)">{flagged}</div>
    <div class="stat-label">Flagged for Review</div>
  </div>
</div>

<div class="stats-grid">
  <div class="stat-card">
    <div class="stat-value" style="color: var(--success)">{conf_dist.get('high', 0)}</div>
    <div class="stat-label">High Confidence</div>
  </div>
  <div class="stat-card">
    <div class="stat-value" style="color: var(--warning)">{conf_dist.get('medium', 0)}</div>
    <div class="stat-label">Medium Confidence</div>
  </div>
  <div class="stat-card">
    <div class="stat-value" style="color: var(--danger)">{conf_dist.get('low', 0)}</div>
    <div class="stat-label">Low Confidence</div>
  </div>
  <div class="stat-card">
    <div class="stat-value" style="color: #CC0000">{conf_dist.get('needs_review', 0)}</div>
    <div class="stat-label">Needs Review</div>
  </div>
</div>

<p><strong>Species:</strong> {params.get('species', 'N/A')} |
   <strong>Tissue:</strong> {params.get('tissue', 'N/A')} |
   <strong>Cluster key:</strong> {params.get('cluster_key', 'N/A')}</p>
"""
    return html


def _html_annotation_table(results: pd.DataFrame) -> str:
    rows_html = []
    for _, row in results.iterrows():
        conf = row.get("critic_confidence", "unknown")
        badge_class = f"badge-{conf.replace('_', '')}" if conf in ("high", "medium", "low", "needs_review") else ""
        flags = row.get("critic_flags", "PASS")
        flag_badge = "badge-pass" if flags == "PASS" else "badge-flag"

        rows_html.append(f"""<tr>
  <td><strong>{row.get('cluster', '')}</strong></td>
  <td>{row.get('cell_type', '')}</td>
  <td>{row.get('cl_id', '')}</td>
  <td>{row.get('combined_score', 0):.3f}</td>
  <td><span class="badge {badge_class}">{conf}</span></td>
  <td><span class="badge {flag_badge}">{flags[:40]}</span></td>
</tr>""")

    return f"""
<h2>Annotation Results</h2>
<table>
<thead>
<tr>
  <th>Cluster</th><th>Cell Type</th><th>CL ID</th>
  <th>Score</th><th>Confidence</th><th>Critic</th>
</tr>
</thead>
<tbody>
{''.join(rows_html)}
</tbody>
</table>
"""


def _html_figures(figure_paths: list[str], output_dir: Path) -> str:
    if not figure_paths:
        return ""

    figures_dir = output_dir / OUTPUT_FIGURES_DIR
    cards = []
    for fp in figure_paths:
        fname = Path(fp).name
        rel_path = f"figures/{fname}"
        title = fname.replace(".png", "").replace("_", " ").title()
        cards.append(f"""<div class="figure-card">
  <img src="{rel_path}" alt="{title}" loading="lazy">
  <div class="caption">{title}</div>
</div>""")

    return f"""
<h2>Figures</h2>
<div class="figure-grid">
{''.join(cards)}
</div>
"""


def _html_critic_details(results: pd.DataFrame) -> str:
    flagged = results[results.get("critic_flags", pd.Series(["PASS"])).apply(lambda x: x != "PASS")]
    if flagged.empty:
        return """
<h2>Critic Details</h2>
<div class="card"><p>All clusters passed the critic review. No flags raised.</p></div>
"""

    cards = []
    for _, row in flagged.iterrows():
        cards.append(f"""<div class="card">
  <div class="card-header">Cluster: {row.get('cluster', '?')} → {row.get('cell_type', '?')}</div>
  <p><strong>Flags:</strong> {row.get('critic_flags', '')}</p>
  <p class="evidence-text"><strong>Evidence:</strong> {row.get('critic_evidence', '')}</p>
  <p class="evidence-text"><strong>Notes:</strong> {row.get('critic_notes', '')}</p>
</div>""")

    return f"""
<h2>Critic Details — Flagged Clusters</h2>
<div class="warning-box">
  {len(flagged)} cluster(s) flagged for review. See details below.
</div>
{''.join(cards)}
"""


def _html_footer() -> str:
    return f"""
<hr style="margin: 2rem 0; border: none; border-top: 1px solid var(--border);">
<p style="text-align: center; color: var(--muted); font-size: 0.8rem;">
  Generated by CellTypePilot v{__version__} | MKG {MKG_VERSION} | {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}
</p>
</body>
</html>"""


def generate_methodology_text(
    manifest: dict,
    critic_summary: dict,
    annotations: pd.DataFrame,
) -> str:
    """Generate a draft methodology paragraph for papers."""
    params = manifest.get("parameters", {})
    total = critic_summary.get("total_clusters", 0)
    high = critic_summary.get("confidence_distribution", {}).get("high", 0)
    med = critic_summary.get("confidence_distribution", {}).get("medium", 0)
    flagged = critic_summary.get("flagged", 0)

    text = (
        f"Cell type annotation was performed using CellTypePilot (v{manifest.get('celltypepilot_version', '?')}), "
        f"an evidence-driven annotation pipeline with built-in critic review. "
        f"Marker gene evidence was sourced from the CellTypePilot Marker Knowledge Graph "
        f"(MKG {manifest.get('mkg_version', '?')}), a curated atlas integrating PanglaoDB, "
        f"CellMarker, and Cell Ontology resources. "
        f"For each of the {total} clusters identified by {params.get('cluster_key', 'clustering')} clustering, "
        f"marker gene overlap, expression specificity, fold-change magnitude, and negative marker "
        f"conflict were scored to generate candidate annotations with confidence levels. "
        f"An independent Annotation Critic module reviewed each assignment for evidence sufficiency, "
        f"negative marker conflicts, potential doublet signatures, and ontology consistency. "
        f"Of {total} clusters, {high} were assigned high confidence, {med} medium confidence, "
        f"and {flagged} were flagged for manual review. "
        f"Species: {params.get('species', 'N/A')}; tissue context: {params.get('tissue', 'N/A')}."
    )
    return text
