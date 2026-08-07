"""Report generator — HTML report and evidence table output.

HTML layout lives in Jinja2 templates (``celltypepilot/templates/report/``)
and styles in ``celltypepilot/static/report.css``, loaded via
``importlib.resources`` so the frontend code stays editable as real
HTML/CSS files.
"""

from __future__ import annotations

import importlib.resources
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from jinja2 import Environment, PackageLoader, select_autoescape

from . import MKG_VERSION, __version__
from .constants import OUTPUT_EVIDENCE, OUTPUT_FIGURES_DIR, OUTPUT_REPORT

# Jinja2 environment backed by package templates
_env = Environment(
    loader=PackageLoader("celltypepilot", "templates/report"),
    autoescape=select_autoescape(["html"]),
)


def _load_report_css() -> str:
    """Load the report stylesheet from package resources."""
    css_resource = importlib.resources.files("celltypepilot").joinpath("static/report.css")
    return css_resource.read_text(encoding="utf-8")


def _render(template_name: str, **context) -> str:
    return _env.get_template(template_name).render(**context)


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
    return _render("header.html", css=_load_report_css())


def _html_overview(manifest: dict, critic_summary: dict) -> str:
    params = manifest.get("parameters", {})
    input_path = manifest.get("input", {}).get("path", "N/A")
    input_filename = Path(input_path).name[:30]

    return _render(
        "overview.html",
        manifest=manifest,
        input_filename=input_filename,
        total=critic_summary.get("total_clusters", 0),
        passed=critic_summary.get("pass", 0),
        flagged=critic_summary.get("flagged", 0),
        conf_dist=critic_summary.get("confidence_distribution", {}),
        narrative=critic_summary.get("narrative", ""),
        params=params,
    )


def _annotation_rows(results: pd.DataFrame) -> list[dict]:
    """Convert critic results into template-friendly row dicts."""
    rows = []
    for _, row in results.iterrows():
        conf = row.get("critic_confidence", "unknown")
        badge_class = (
            f"badge-{conf.replace('_', '')}"
            if conf in ("high", "medium", "low", "needs_review")
            else ""
        )
        flags = row.get("critic_flags", "PASS")
        flag_badge = "badge-pass" if flags == "PASS" else "badge-flag"
        rows.append(
            {
                "cluster": row.get("cluster", ""),
                "cell_type": row.get("cell_type", ""),
                "cl_id": row.get("cl_id", ""),
                "score": float(row.get("combined_score", 0)),
                "confidence": conf,
                "badge_class": badge_class,
                "flags": flags,
                "flag_badge": flag_badge,
                "evidence_summary": row.get("evidence_summary", ""),
            }
        )
    return rows


def _html_annotation_table(results: pd.DataFrame) -> str:
    return _render("annotation_table.html", rows=_annotation_rows(results))


def _html_figures(figure_paths: list[str], output_dir: Path) -> str:
    if not figure_paths:
        return ""

    figures = []
    for fp in figure_paths:
        fname = Path(fp).name
        figures.append(
            {
                "rel_path": f"{OUTPUT_FIGURES_DIR}/{fname}",
                "title": fname.replace(".png", "").replace("_", " ").title(),
            }
        )
    return _render("figures.html", figures=figures)


def _html_critic_details(results: pd.DataFrame) -> str:
    flags_col = results.get("critic_flags", pd.Series(["PASS"] * len(results)))
    flagged = results[flags_col.apply(lambda x: x != "PASS")]
    if flagged.empty:
        return _render("critic_details.html", flagged_rows=[])

    flagged_rows = []
    for _, row in flagged.iterrows():
        flagged_rows.append(
            {
                "cluster": row.get("cluster", "?"),
                "cell_type": row.get("cell_type", "?"),
                "critic_flags": row.get("critic_flags", ""),
                "critic_evidence": row.get("critic_evidence", ""),
                "critic_notes": row.get("critic_notes", ""),
                "evidence_summary": row.get("evidence_summary", ""),
            }
        )
    return _render("critic_details.html", flagged_rows=flagged_rows)


def _html_footer() -> str:
    return _render(
        "footer.html",
        version=__version__,
        mkg_version=MKG_VERSION,
        generated_utc=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
    )


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
