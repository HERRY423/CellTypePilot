"""Tests for CLI commands using Typer's CliRunner."""

import json

from typer.testing import CliRunner

from celltypepilot.cli import app

runner = CliRunner()


class TestCLIBasic:
    def test_version(self):
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert "CellTypePilot" in result.output

    def test_help(self):
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        assert "annotate" in result.output
        assert "inspect" in result.output
        assert "doctor" in result.output

    def test_doctor(self):
        result = runner.invoke(app, ["doctor"])
        assert result.exit_code == 0

    def test_markers_list_tissues(self):
        result = runner.invoke(app, ["markers"])
        assert result.exit_code == 0
        assert "Available tissues" in result.output or "blood" in result.output.lower()

    def test_markers_specific_tissue(self):
        result = runner.invoke(app, ["markers", "--tissue", "blood"])
        assert result.exit_code == 0
        assert "T cell" in result.output or "cell types" in result.output.lower()

    def test_markers_json(self):
        result = runner.invoke(app, ["markers", "--tissue", "blood", "--json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert isinstance(data, dict)
        assert len(data) > 0

    def test_markers_unknown_tissue(self):
        result = runner.invoke(app, ["markers", "--tissue", "nonexistent_tissue"])
        assert result.exit_code == 0
        # Falls back to "general" tissue markers, so should still produce output
        assert len(result.output) > 20


class TestCLIInspect:
    def test_inspect_json(self, h5ad_path):
        result = runner.invoke(app, ["inspect", "-i", str(h5ad_path), "--json"])
        assert result.exit_code == 0
        # Output may have prefix text before JSON; extract the JSON block
        output = result.output
        json_start = output.index("{")
        # Use strict=False to tolerate control characters in Rich output
        data = json.loads(output[json_start:], strict=False)
        assert "species" in data
        assert "cluster_keys" in data

    def test_inspect_text(self, h5ad_path):
        result = runner.invoke(app, ["inspect", "-i", str(h5ad_path)])
        assert result.exit_code == 0
        assert len(result.output) > 50

    def test_inspect_missing_file(self):
        result = runner.invoke(app, ["inspect", "-i", "/nonexistent/path.h5ad"])
        assert result.exit_code != 0


class TestCLIAnnotate:
    def test_annotate_basic(self, h5ad_path, tmp_output_dir):
        result = runner.invoke(
            app,
            [
                "annotate",
                "-i",
                str(h5ad_path),
                "-k",
                "leiden",
                "-o",
                str(tmp_output_dir),
                "--species",
                "human",
                "--tissue",
                "blood",
                "--no-figures",
            ],
        )
        assert result.exit_code == 0
        assert "Done" in result.output

    def test_annotate_json_output(self, h5ad_path, tmp_output_dir):
        result = runner.invoke(
            app,
            [
                "annotate",
                "-i",
                str(h5ad_path),
                "-k",
                "leiden",
                "-o",
                str(tmp_output_dir),
                "--species",
                "human",
                "--tissue",
                "blood",
                "--no-figures",
                "--json",
            ],
        )
        assert result.exit_code == 0
        assert "annotations" in result.output

    def test_annotate_bad_cluster_key(self, h5ad_path, tmp_output_dir):
        result = runner.invoke(
            app,
            [
                "annotate",
                "-i",
                str(h5ad_path),
                "-k",
                "nonexistent_key",
                "-o",
                str(tmp_output_dir),
                "--species",
                "human",
                "--tissue",
                "blood",
                "--no-figures",
            ],
        )
        assert result.exit_code == 1


class TestCLILiterature:
    def test_literature_no_markers(self):
        result = runner.invoke(app, ["literature", "-c", "T cell"])
        assert result.exit_code == 0
        assert (
            "queries" in result.output.lower()
            or "query" in result.output.lower()
            or "search" in result.output.lower()
        )

    def test_literature_json(self):
        result = runner.invoke(app, ["literature", "-c", "T cell", "--json"])
        assert result.exit_code == 0
        # Extract JSON block from output (may have prefix text)
        output = result.output
        json_start = output.index("{")
        data = json.loads(output[json_start:])
        assert "queries" in data or "mcp_status" in data
