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
CHECKSUM_SCRIPT = REPO_ROOT / "scripts/build_release_checksums.py"


def _project_version() -> str:
    for line in (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8").splitlines():
        if line.startswith('version = "'):
            return line.split('"', 2)[1]
    raise AssertionError("project version not found")


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
    version = _project_version()
    result = subprocess.run(
        [
            sys.executable,
            str(BUNDLE_SCRIPT),
            "--output-dir",
            str(tmp_path),
            "--tag",
            f"v{version}",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    bundle = Path(result.stdout.strip())
    assert bundle.is_file()

    prefix = f"celltypepilot-plugin-{version}/"
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
        assert manifest["version"] == version
        assert manifest["distribution"] == "agent_plugin_bundle"
        for record in manifest["files"]:
            payload = archive.read(prefix + record["path"])
            assert len(payload) == record["size"]
            assert hashlib.sha256(payload).hexdigest() == record["sha256"]


def test_plugin_bundle_rejects_mismatched_release_tag(tmp_path):
    version = _project_version()
    result = subprocess.run(
        [sys.executable, str(BUNDLE_SCRIPT), "--output-dir", str(tmp_path), "--tag", "v9.9.9"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert f"does not match project version v{version}" in result.stderr


def test_release_checksums_cover_exact_versioned_artifacts(tmp_path):
    version = _project_version()
    dist_dir = tmp_path / "dist"
    plugin_dir = tmp_path / "release-assets"
    dist_dir.mkdir()
    plugin_dir.mkdir()
    artifacts = {
        dist_dir / f"celltypepilot-{version}-py3-none-any.whl": b"wheel payload",
        dist_dir / f"celltypepilot-{version}.tar.gz": b"sdist payload",
        plugin_dir / f"celltypepilot-plugin-{version}.zip": b"plugin payload",
    }
    for path, payload in artifacts.items():
        path.write_bytes(payload)

    output_path = plugin_dir / "SHA256SUMS"
    result = subprocess.run(
        [
            sys.executable,
            str(CHECKSUM_SCRIPT),
            "--dist-dir",
            str(dist_dir),
            "--plugin-dir",
            str(plugin_dir),
            "--output",
            str(output_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert Path(result.stdout.strip()) == output_path
    expected = [
        f"{hashlib.sha256(payload).hexdigest()}  {path.name}"
        for path, payload in sorted(artifacts.items())
    ]
    assert output_path.read_text(encoding="utf-8").splitlines() == expected
    assert all("dist" not in line and "release-assets" not in line for line in expected)


def test_release_checksums_fail_closed_when_artifact_is_missing(tmp_path):
    version = _project_version()
    dist_dir = tmp_path / "dist"
    plugin_dir = tmp_path / "release-assets"
    dist_dir.mkdir()
    plugin_dir.mkdir()
    (dist_dir / f"celltypepilot-{version}.tar.gz").write_bytes(b"sdist payload")
    (plugin_dir / f"celltypepilot-plugin-{version}.zip").write_bytes(b"plugin payload")

    result = subprocess.run(
        [
            sys.executable,
            str(CHECKSUM_SCRIPT),
            "--dist-dir",
            str(dist_dir),
            "--plugin-dir",
            str(plugin_dir),
            "--output",
            str(plugin_dir / "SHA256SUMS"),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "expected exactly one wheel" in result.stderr


def test_release_workflow_uses_tested_checksum_builder():
    workflow = (REPO_ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    assert "python scripts/build_release_checksums.py" in workflow
    assert "> ../release-assets/SHA256SUMS" not in workflow
    assert ") >> SHA256SUMS" not in workflow
