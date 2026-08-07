"""Tests for the reporter module — HTML report and methodology text."""

import tempfile
from pathlib import Path

import pandas as pd

from celltypepilot.reporter import (
    _html_annotation_table,
    _html_critic_details,
    _html_figures,
    _html_footer,
    _html_header,
    _html_overview,
    generate_html_report,
    generate_methodology_text,
    save_evidence_table,
)


class TestSaveEvidenceTable:
    def test_creates_csv(self, sample_critic_results, tmp_output_dir):
        path = save_evidence_table(sample_critic_results, tmp_output_dir)
        assert path.exists()
        assert path.suffix == ".csv"
        assert path.stat().st_size > 0

    def test_csv_content_matches(self, sample_critic_results, tmp_output_dir):
        path = save_evidence_table(sample_critic_results, tmp_output_dir)
        loaded = pd.read_csv(path)
        assert len(loaded) == len(sample_critic_results)
        assert "cluster" in loaded.columns

    def test_creates_output_dir(self, sample_critic_results):
        with tempfile.TemporaryDirectory() as tmpdir:
            new_dir = Path(tmpdir) / "nested" / "output"
            path = save_evidence_table(sample_critic_results, new_dir)
            assert path.exists()


class TestHtmlComponents:
    def test_html_header(self):
        html = _html_header()
        assert "<!DOCTYPE html>" in html
        assert "CellTypePilot" in html
        assert "<html" in html

    def test_html_footer(self):
        html = _html_footer()
        assert "</html>" in html
        assert "CellTypePilot" in html

    def test_html_overview(self):
        manifest = {
            "timestamp": "2026-01-01T00:00:00",
            "celltypepilot_version": "0.1.0",
            "mkg_version": "mkg-2026.08",
            "input": {"path": "/data/test.h5ad"},
            "parameters": {"species": "human", "tissue": "blood", "cluster_key": "leiden"},
        }
        critic_summary = {
            "total_clusters": 10,
            "pass": 8,
            "flagged": 2,
            "confidence_distribution": {"high": 5, "medium": 3, "low": 1, "needs_review": 1},
        }
        html = _html_overview(manifest, critic_summary)
        assert "10" in html  # total clusters
        assert "8" in html  # passed
        assert "2" in html  # flagged
        assert "human" in html
        assert "blood" in html

    def test_html_annotation_table(self, sample_critic_results):
        html = _html_annotation_table(sample_critic_results)
        assert "<table>" in html
        assert "Cluster" in html
        assert "Cell Type" in html

    def test_html_figures_empty(self):
        html = _html_figures([], Path("/tmp"))
        assert html == ""

    def test_html_figures_with_paths(self):
        paths = ["/tmp/figures/umap_clusters.png", "/tmp/figures/dotplot.png"]
        html = _html_figures(paths, Path("/tmp"))
        assert "umap_clusters" in html
        assert "dotplot" in html

    def test_html_critic_details_all_pass(self, sample_critic_results):
        html = _html_critic_details(sample_critic_results)
        assert "All clusters passed" in html

    def test_html_critic_details_with_flags(self):
        results = pd.DataFrame(
            {
                "cluster": ["0", "1"],
                "cell_type": ["T cell", "B cell"],
                "cl_id": ["CL:001", ""],
                "combined_score": [0.8, 0.3],
                "critic_flags": ["PASS", "LOW_EVIDENCE"],
                "critic_evidence": ["Good", "Weak"],
                "critic_confidence": ["high", "needs_review"],
                "critic_notes": ["", "Review needed"],
            }
        )
        html = _html_critic_details(results)
        assert "Flagged" in html
        assert "LOW_EVIDENCE" in html


class TestGenerateHtmlReport:
    def test_full_report(self, sample_critic_results, tmp_output_dir):
        manifest = {
            "timestamp": "2026-01-01",
            "celltypepilot_version": "0.1.0",
            "mkg_version": "mkg-2026.08",
            "input": {"path": "test.h5ad"},
            "parameters": {"species": "human", "tissue": "blood", "cluster_key": "leiden"},
        }
        critic_summary = {
            "total_clusters": len(sample_critic_results),
            "pass": len(sample_critic_results),
            "flagged": 0,
            "confidence_distribution": {"high": len(sample_critic_results)},
        }
        path = generate_html_report(
            sample_critic_results,
            sample_critic_results,
            critic_summary,
            manifest,
            [],
            tmp_output_dir,
        )
        assert path.exists()
        content = path.read_text(encoding="utf-8")
        assert "<!DOCTYPE html>" in content
        assert "CellTypePilot" in content


class TestMethodologyText:
    def test_basic_methodology(self):
        manifest = {
            "celltypepilot_version": "0.1.0",
            "mkg_version": "mkg-2026.08",
            "parameters": {"cluster_key": "leiden", "species": "human", "tissue": "blood"},
        }
        critic_summary = {
            "total_clusters": 10,
            "confidence_distribution": {"high": 6, "medium": 3, "low": 1, "needs_review": 0},
            "flagged": 1,
        }
        annotations = pd.DataFrame({"cluster": range(10)})

        text = generate_methodology_text(manifest, critic_summary, annotations)
        assert "CellTypePilot" in text
        assert "0.1.0" in text
        assert "mkg-2026.08" in text
        assert "10" in text  # total clusters
        assert "human" in text
        assert "blood" in text
        assert len(text) > 100  # Should be a substantial paragraph
