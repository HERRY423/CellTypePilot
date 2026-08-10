"""Build a deterministic SHA256SUMS file for CellTypePilot release artifacts."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
VERSION_PATTERN = re.compile(r'^version = "([^"]+)"$', re.MULTILINE)


def _read_version(root: Path = REPO_ROOT) -> str:
    pyproject = (root / "pyproject.toml").read_text(encoding="utf-8")
    match = VERSION_PATTERN.search(pyproject)
    if not match:
        raise ValueError("project version not found in pyproject.toml")
    return match.group(1)


def _require_single(directory: Path, pattern: str, label: str) -> Path:
    matches = sorted(path for path in directory.glob(pattern) if path.is_file())
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one {label} matching {directory / pattern}, found {len(matches)}"
        )
    return matches[0]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_checksums(
    dist_dir: Path,
    plugin_dir: Path,
    output_path: Path,
    *,
    root: Path = REPO_ROOT,
) -> Path:
    """Hash the one wheel, sdist, and plugin bundle for the current project version."""
    version = _read_version(root)
    artifacts = [
        _require_single(dist_dir, f"celltypepilot-{version}-*.whl", "wheel"),
        _require_single(dist_dir, f"celltypepilot-{version}.tar.gz", "source distribution"),
        _require_single(
            plugin_dir,
            f"celltypepilot-plugin-{version}.zip",
            "Agent plugin bundle",
        ),
    ]

    names = [artifact.name for artifact in artifacts]
    if len(names) != len(set(names)):
        raise ValueError("release artifact basenames must be unique")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{_sha256(artifact)}  {artifact.name}\n" for artifact in sorted(artifacts)]
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.writelines(lines)
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dist-dir", type=Path, default=Path("dist"))
    parser.add_argument("--plugin-dir", type=Path, default=Path("release-assets"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("release-assets/SHA256SUMS"),
    )
    args = parser.parse_args()
    print(build_checksums(args.dist_dir, args.plugin_dir, args.output))


if __name__ == "__main__":
    main()
