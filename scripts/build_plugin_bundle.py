"""Build a deterministic CellTypePilot Agent plugin release bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
VERSION_PATTERN = re.compile(r'^version = "([^"]+)"$', re.MULTILINE)
INIT_VERSION_PATTERN = re.compile(r'^__version__ = "([^"]+)"$', re.MULTILINE)
BUNDLE_COMPONENTS = (
    Path(".claude-plugin"),
    Path(".codex-plugin"),
    Path("skills"),
    Path("commands"),
    Path("hooks"),
    Path("rules"),
    Path("src/celltypepilot"),
    Path(".mcp.json"),
    Path("AGENTS.md"),
    Path("CHANGELOG.md"),
    Path("LICENSE"),
    Path("README.md"),
    Path("pyproject.toml"),
)
EXCLUDED_PARTS = {"__pycache__", ".pytest_cache", ".ruff_cache"}
EXCLUDED_SUFFIXES = {".pyc", ".pyo"}
ZIP_TIMESTAMP = (2020, 1, 1, 0, 0, 0)


def _read_version(root: Path = REPO_ROOT) -> str:
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    match = VERSION_PATTERN.search(pyproject)
    if not match:
        raise ValueError("project version not found in pyproject.toml")
    return match.group(1)


def _validate_versions(root: Path, version: str, tag: str | None) -> None:
    init_text = (root / "src/celltypepilot/__init__.py").read_text(encoding="utf-8")
    init_match = INIT_VERSION_PATTERN.search(init_text)
    if not init_match or init_match.group(1) != version:
        raise ValueError("src/celltypepilot/__init__.py version does not match pyproject.toml")

    for manifest_path in (Path(".codex-plugin/plugin.json"), Path(".claude-plugin/plugin.json")):
        manifest = json.loads((root / manifest_path).read_text(encoding="utf-8"))
        if manifest.get("version") != version:
            raise ValueError(f"{manifest_path} version does not match pyproject.toml")

    if tag is not None and tag != f"v{version}":
        raise ValueError(f"release tag {tag!r} does not match project version v{version}")
    if tag is not None:
        changelog = (root / "CHANGELOG.md").read_text(encoding="utf-8")
        dated_heading = re.compile(
            rf"^## \[{re.escape(version)}\] - \d{{4}}-\d{{2}}-\d{{2}}$", re.MULTILINE
        )
        if not dated_heading.search(changelog):
            raise ValueError(f"CHANGELOG.md needs a dated {version} release heading")


def _bundle_files(root: Path) -> list[Path]:
    files: set[Path] = set()
    for component in BUNDLE_COMPONENTS:
        source = root / component
        if not source.exists():
            raise FileNotFoundError(f"required bundle component is missing: {component}")
        candidates = [source] if source.is_file() else source.rglob("*")
        for candidate in candidates:
            if not candidate.is_file():
                continue
            relative = candidate.relative_to(root)
            if any(part in EXCLUDED_PARTS for part in relative.parts):
                continue
            if relative.suffix in EXCLUDED_SUFFIXES:
                continue
            if candidate.is_symlink():
                raise ValueError(f"symbolic links are not allowed in the plugin bundle: {relative}")
            files.add(relative)
    return sorted(files, key=lambda path: path.as_posix())


def _zip_info(name: str) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(name, date_time=ZIP_TIMESTAMP)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.external_attr = 0o644 << 16
    return info


def build_bundle(output_dir: Path, tag: str | None = None, root: Path = REPO_ROOT) -> Path:
    """Build the plugin bundle and return its path."""
    version = _read_version(root)
    _validate_versions(root, version, tag)
    payload_files = _bundle_files(root)
    archive_root = f"celltypepilot-plugin-{version}"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"celltypepilot-plugin-{version}.zip"

    manifest_files = []
    payloads: list[tuple[Path, bytes]] = []
    for relative in payload_files:
        data = (root / relative).read_bytes()
        payloads.append((relative, data))
        manifest_files.append(
            {
                "path": relative.as_posix(),
                "sha256": hashlib.sha256(data).hexdigest(),
                "size": len(data),
            }
        )

    manifest = {
        "schema_version": "celltypepilot.plugin-bundle.v1",
        "name": "celltypepilot",
        "version": version,
        "distribution": "agent_plugin_bundle",
        "backend_install": "python -m pip install .",
        "pypi_backend": f"celltypepilot=={version}",
        "files": manifest_files,
    }
    manifest_data = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()

    with zipfile.ZipFile(output_path, "w") as archive:
        for relative, data in payloads:
            archive.writestr(_zip_info(f"{archive_root}/{relative.as_posix()}"), data)
        archive.writestr(_zip_info(f"{archive_root}/BUNDLE_MANIFEST.json"), manifest_data)

    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("release-assets"))
    parser.add_argument("--tag", help="Release tag; must equal v<project.version>")
    args = parser.parse_args()
    print(build_bundle(args.output_dir, tag=args.tag))


if __name__ == "__main__":
    main()
