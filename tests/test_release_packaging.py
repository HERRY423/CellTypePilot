"""Release packaging contract tests."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).parents[1]
BUNDLE_SCRIPT = REPO_ROOT / "scripts/build_plugin_bundle.py"


def test_pyproject_keeps_extras_out_of_project_urls():
    pyproject = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    optional_block = pyproject.split("[project.optional-dependencies]", 1)[1].split(
        "[project.urls]", 1
    )[0]
    urls_block = pyproject.split("[project.urls]", 1)[1].split("[project.scripts]", 1)[0]
    assert "all = [" in optional_block
    assert "dev = [" in optional_block
    assert "all = [" not in urls_block
    assert "dev = [" not in urls_block


def test_plugin_bundle_contains_agent_surface_and_installable_backend(tmp_path):
    result = subprocess.run(
        [sys.executable, str(BUNDLE_SCRIPT), "--output-dir", str(tmp_path), "--tag", "v0.3.0"],
        check=True,
        capture_output=True,
        text=True,
    )
    bundle = Path(result.stdout.strip())
    assert bundle.is_file()

    prefix = "celltypepilot-plugin-0.3.0/"
    with zipfile.ZipFile(bundle) as archive:
        names = set(archive.namelist())
        required = {
            ".codex-plugin/plugin.json",
            ".claude-plugin/plugin.json",
            "skills/celltypepilot/SKILL.md",
            "src/celltypepilot/__init__.py",
            "src/celltypepilot/data/packs/premium/marker_atlas.json",
            "pyproject.toml",
            "BUNDLE_MANIFEST.json",
        }
        assert {prefix + path for path in required}.issubset(names)
        assert not any("/tests/" in name or "/benchmarks/" in name for name in names)
        assert not any(name.endswith((".h5ad", ".pyc")) for name in names)

        manifest = json.loads(archive.read(prefix + "BUNDLE_MANIFEST.json"))
        assert manifest["schema_version"] == "celltypepilot.plugin-bundle.v1"
        assert manifest["version"] == "0.3.0"
        assert manifest["distribution"] == "agent_plugin_bundle"
        for record in manifest["files"]:
            payload = archive.read(prefix + record["path"])
            assert len(payload) == record["size"]
            assert hashlib.sha256(payload).hexdigest() == record["sha256"]


def test_plugin_bundle_rejects_mismatched_release_tag(tmp_path):
    result = subprocess.run(
        [sys.executable, str(BUNDLE_SCRIPT), "--output-dir", str(tmp_path), "--tag", "v9.9.9"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "does not match project version v0.3.0" in result.stderr
