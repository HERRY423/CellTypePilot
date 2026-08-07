"""Provenance tracking — version, parameter, and data hash management."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from . import MKG_VERSION, __version__
from .constants import OUTPUT_MANIFEST


def create_manifest(
    input_path: str,
    data_hash: str,
    cluster_key: str,
    species: str,
    tissue: str,
    parameters: dict,
    output_dir: str | Path,
) -> dict:
    """Create a run manifest for reproducibility."""
    manifest = {
        "celltypepilot_version": __version__,
        "mkg_version": MKG_VERSION,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "input": {
            "path": str(input_path),
            "sha256": data_hash,
        },
        "parameters": {
            "cluster_key": cluster_key,
            "species": species,
            "tissue": tissue,
            **parameters,
        },
        "outputs": {},
    }
    return manifest


def update_manifest_outputs(manifest: dict, output_dir: str | Path) -> dict:
    """Update manifest with output file hashes."""
    output_dir = Path(output_dir)
    for fpath in output_dir.rglob("*"):
        if fpath.is_file() and fpath.name != OUTPUT_MANIFEST:
            rel_path = str(fpath.relative_to(output_dir))
            file_hash = _file_hash(fpath)
            manifest["outputs"][rel_path] = {
                "sha256": file_hash,
                "size_bytes": fpath.stat().st_size,
            }
    return manifest


def save_manifest(manifest: dict, output_dir: str | Path) -> Path:
    """Save manifest to JSON file."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / OUTPUT_MANIFEST
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    return manifest_path


def load_manifest(path: str | Path) -> dict:
    """Load a previously saved manifest."""
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _file_hash(path: Path) -> str:
    """Compute SHA-256 hash of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def format_manifest_summary(manifest: dict) -> str:
    """Format manifest as human-readable summary."""
    lines = [
        "=" * 50,
        "CellTypePilot — Run Manifest",
        "=" * 50,
        f"Version:     {manifest['celltypepilot_version']}",
        f"MKG Version: {manifest['mkg_version']}",
        f"Timestamp:   {manifest['timestamp']}",
        "",
        "Input:",
        f"  Path:      {manifest['input']['path']}",
        f"  SHA-256:   {manifest['input']['sha256'][:16]}...",
        "",
        "Parameters:",
    ]
    for k, v in manifest["parameters"].items():
        lines.append(f"  {k}: {v}")

    if manifest.get("outputs"):
        lines.append("")
        lines.append(f"Outputs ({len(manifest['outputs'])} files):")
        for name, info in manifest["outputs"].items():
            lines.append(f"  {name} ({info['size_bytes']:,} bytes)")

    lines.append("=" * 50)
    return "\n".join(lines)
