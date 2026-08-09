"""Extension Pack management — data-only domain packs (spec v1).

Extension packs let third-party labs extend CellTypePilot with domain
marker atlases and state atlases WITHOUT executable code. A pack is a
directory containing:

    pack.json           manifest (celltypepilot.pack.v1)
    marker_atlas.json   optional, bundled-atlas schema with provenance
    state_atlas.json    optional, celltypepilot.state-atlas.v1 schema

Design principles (mirrors the reference_registry fail-closed contracts):

- Data only: packs never contain Python; loading them cannot execute code.
- Fail closed: packs failing schema/provenance validation are rejected at
  install time unless the user explicitly opts into hypothesis trust.
- Trust tiers: provenance-verified packs merge at "atlas" trust (same
  rights as the built-in MKG); incomplete packs merge at "hypothesis"
  trust and can never silently produce an accepted identity.
- Auditable: installed packs record per-file sha256; runs record packs in
  the manifest and re-verify hashes at load time.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

from .atlas_conflict import detect_marker_conflicts
from .constants import FIRST_PARTY_PACKS_DIR

PACK_SCHEMA_VERSION = "celltypepilot.pack.v1"
PACKS_ENV_VAR = "CELLTYPEPILOT_PACKS_DIR"
INSTALLED_METADATA_FILE = "_installed.json"
ATLAS_FILE = "marker_atlas.json"
STATE_FILE = "state_atlas.json"
ONTOLOGY_MAP_FILE = "ontology_map.json"
REFERENCE_MANIFEST_FILE = "reference_manifest.json"
ALLOWED_PACK_FILES = {
    ATLAS_FILE,
    STATE_FILE,
    ONTOLOGY_MAP_FILE,
    REFERENCE_MANIFEST_FILE,
}
TRUST_ATLAS = "atlas"
TRUST_HYPOTHESIS = "hypothesis"
TRUST_LEVELS = (TRUST_ATLAS, TRUST_HYPOTHESIS)
PACK_LICENSE_TIERS = ("community", "academic", "commercial")
LICENSE_TIER_RANK = {"free": 0, "trial": 0, "academic": 1, "commercial": 2}
PACK_TIER_RANK = {"community": 0, "academic": 1, "commercial": 2}
SPECIES_VALUES = ("human", "mouse")


class PackError(ValueError):
    """Raised when an extension pack is malformed, unauthorized, or missing."""


def packs_dir() -> Path:
    """User pack installation directory (override via env var for tests)."""
    override = os.environ.get(PACKS_ENV_VAR)
    if override:
        return Path(override)
    return Path.home() / ".celltypepilot" / "packs"


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _is_git_url(source: str) -> bool:
    return source.startswith(("http://", "https://", "git@", "ssh://")) or source.endswith(
        ".git"
    )


def _current_license_rank() -> int:
    from .license_manager import load_license

    return LICENSE_TIER_RANK.get(load_license().tier.value, 0)


def _check_pack_license(manifest: dict, action: str) -> None:
    required = str(manifest.get("license_tier", "community"))
    if required not in PACK_TIER_RANK:
        raise PackError(f"Pack declares unknown license_tier {required!r}")
    if PACK_TIER_RANK[required] > _current_license_rank():
        raise PackError(
            f"Pack {manifest.get('name', '?')!r} requires a {required} license "
            f"to {action}. Upgrade at https://celltypepilot.io/license"
        )


def read_pack_manifest(pack_dir: str | Path) -> dict:
    """Load and minimally parse a pack.json manifest."""
    manifest_path = Path(pack_dir) / "pack.json"
    if not manifest_path.is_file():
        raise PackError(f"pack.json not found in {pack_dir}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise PackError(f"pack.json is not valid JSON: {exc}") from exc
    if not isinstance(manifest, dict):
        raise PackError("pack.json must be a JSON object")
    return manifest


def validate_pack_manifest(manifest: dict) -> list[str]:
    """Return schema problems for a pack manifest (empty list = valid)."""
    issues = []
    schema = manifest.get("schema_version")
    if schema != PACK_SCHEMA_VERSION:
        issues.append(f"schema_version must be {PACK_SCHEMA_VERSION!r}, got {schema!r}")
    name = manifest.get("name", "")
    if not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", str(name)):
        issues.append(f"name must match [a-z0-9][a-z0-9_-]*, got {name!r}")
    if not str(manifest.get("version", "")).strip():
        issues.append("version is required")
    species = manifest.get("species", ["human"])
    if not isinstance(species, list) or not species:
        issues.append("species must be a non-empty list")
    elif any(item not in SPECIES_VALUES for item in species):
        issues.append(f"species must be a subset of {list(SPECIES_VALUES)}")
    tissues = manifest.get("tissues", [])
    if not isinstance(tissues, list):
        issues.append("tissues must be a list of tissue keys")
    diseases = manifest.get("diseases", [])
    if diseases is not None and not isinstance(diseases, list):
        issues.append("diseases must be a list when present")
    pack_kind = manifest.get("pack_kind", "evidence")
    if pack_kind not in ("evidence", "reference", "mixed"):
        issues.append("pack_kind must be evidence, reference, or mixed")
    if manifest.get("license_tier", "community") not in PACK_LICENSE_TIERS:
        issues.append("license_tier must be one of community/academic/commercial")
    if manifest.get("trust", TRUST_ATLAS) not in TRUST_LEVELS:
        issues.append("trust must be atlas or hypothesis")
    files = manifest.get("files", [])
    if not isinstance(files, list) or not files:
        issues.append("files must be a non-empty list of pack data files")
    elif any(item not in ALLOWED_PACK_FILES for item in files):
        issues.append(f"files must be a subset of {sorted(ALLOWED_PACK_FILES)}")
    # Runtime gate declaration is documentation only; all packs still hit gates.
    gates = manifest.get("runtime_gates")
    if gates is not None and (not isinstance(gates, list) or not gates):
        issues.append("runtime_gates must be a non-empty list when present")
    return issues


def validate_pack(pack_dir: str | Path) -> list[str]:
    """Validate a pack directory: data-only + manifest + atlas provenance gates.

    Marker atlases must pass the same provenance validator as the built-in
    MKG; state atlases must pass the state provenance validator. Executable
    code is always rejected.
    """
    pack_dir = Path(pack_dir)
    try:
        from .pack_signing import scan_pack_for_code

        code_issues = scan_pack_for_code(pack_dir)
    except Exception:
        code_issues = []
    try:
        manifest = read_pack_manifest(pack_dir)
    except PackError as exc:
        return code_issues + [str(exc)]
    issues = code_issues + validate_pack_manifest(manifest)

    for filename in ALLOWED_PACK_FILES:
        path = pack_dir / filename
        if not path.exists():
            continue
        try:
            content = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            issues.append(f"{filename} is not valid JSON: {exc}")
            continue
        if filename == ATLAS_FILE:
            from .data_adapter import validate_atlas_provenance

            if not str(content.get("version", "")).strip():
                issues.append(f"{filename} missing top-level version")
            issues.extend(validate_atlas_provenance(content))
        elif filename == STATE_FILE:
            from .state_scorer import validate_state_atlas

            if not str(content.get("version", "")).strip():
                issues.append(f"{filename} missing top-level version")
            issues.extend(validate_state_atlas(content))
        elif filename == ONTOLOGY_MAP_FILE:
            if content.get("schema_version") != "celltypepilot.ontology-map.v1":
                issues.append(
                    f"{filename} schema_version must be celltypepilot.ontology-map.v1"
                )
            for field in ("aliases", "safe_parent_fallbacks"):
                if field in content and not isinstance(content[field], dict):
                    issues.append(f"{filename} {field} must be an object")
            if "include_tissues" in content and not isinstance(
                content["include_tissues"], list
            ):
                issues.append(f"{filename} include_tissues must be a list")
        elif filename == REFERENCE_MANIFEST_FILE:
            if not isinstance(content, dict) or not content.get("schema_version"):
                issues.append(f"{filename} requires a schema_version")

    declared = set(manifest.get("files", []))
    present = {name for name in ALLOWED_PACK_FILES if (pack_dir / name).exists()}
    if declared and not (declared & ALLOWED_PACK_FILES) & present:
        issues.append("no declared pack data file is present on disk")
    return issues


def _materialize_source(source: str) -> tuple[Path, Path | None]:
    """Resolve an install source to a local directory.

    Returns (pack_dir, temp_root). temp_root is set when the source was a
    git URL and the clone must be cleaned up by the caller.
    """
    if not _is_git_url(source):
        local = Path(source)
        if not local.is_dir():
            raise PackError(f"Pack source directory not found: {local}")
        return local, None

    if shutil.which("git") is None:
        raise PackError("git is required to install packs from URLs but was not found")
    temp_root = Path(tempfile.mkdtemp(prefix="celltypepilot-pack-"))
    clone_dir = temp_root / "repo"
    try:
        subprocess.run(
            ["git", "clone", "--depth", "1", source, str(clone_dir)],
            check=True,
            capture_output=True,
            text=True,
            timeout=300,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        detail = getattr(exc, "stderr", "") or str(exc)
        shutil.rmtree(temp_root, ignore_errors=True)
        raise PackError(f"git clone failed for {source}: {detail}") from exc
    return clone_dir, temp_root


def install_pack(
    source: str,
    trust: str = TRUST_ATLAS,
    force: bool = False,
) -> dict:
    """Install a pack from a local directory or git URL.

    Fails closed: a pack with validation issues is rejected unless
    ``trust="hypothesis"`` is explicitly requested, in which case the pack
    is recorded at hypothesis trust and can never act as atlas-level
    evidence.

    Returns a summary dict describing the installed pack.
    """
    if trust not in TRUST_LEVELS:
        raise PackError(f"trust must be one of {list(TRUST_LEVELS)}, got {trust!r}")

    pack_dir, temp_root = _materialize_source(source)
    try:
        manifest = read_pack_manifest(pack_dir)
        issues = validate_pack(pack_dir)
        if issues and trust == TRUST_ATLAS:
            raise PackError(
                "Pack failed provenance validation (fail closed). "
                "Fix the pack or install with --trust hypothesis: " + "; ".join(issues[:5])
            )
        _check_pack_license(manifest, "install")

        name = str(manifest["name"])
        if (FIRST_PARTY_PACKS_DIR / name).exists():
            raise PackError(f"Pack name {name!r} is reserved for a first-party pack")

        content_files = {
            filename: pack_dir / filename
            for filename in ALLOWED_PACK_FILES
            if (pack_dir / filename).exists()
        }
        if not content_files:
            raise PackError("Pack contains no marker_atlas.json or state_atlas.json")

        destination = packs_dir() / name
        if destination.exists():
            if not force:
                raise PackError(
                    f"Pack {name!r} is already installed; use --force to reinstall"
                )
            shutil.rmtree(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)

        def _ignore_git(directory, entries):
            return [".git"] if Path(directory) == pack_dir else []

        shutil.copytree(pack_dir, destination, ignore=_ignore_git)
        metadata = {
            "trust": trust,
            "installed_at": datetime.now(timezone.utc).isoformat(),
            "source": source,
            "validation_issues_at_install": issues,
            "content_sha256": {
                filename: _file_sha256(path) for filename, path in content_files.items()
            },
        }
        (destination / INSTALLED_METADATA_FILE).write_text(
            json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        return {
            "name": name,
            "version": str(manifest.get("version", "")),
            "trust": trust,
            "path": str(destination),
            "tissues": list(manifest.get("tissues", [])),
            "species": list(manifest.get("species", ["human"])),
            "validation_issues": issues,
        }
    finally:
        if temp_root is not None:
            shutil.rmtree(temp_root, ignore_errors=True)


def _installed_metadata(pack_dir: Path) -> dict:
    path = pack_dir / INSTALLED_METADATA_FILE
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def list_installed_packs(include_first_party: bool = True) -> list[dict]:
    """List first-party and user-installed packs (metadata only, no loading)."""
    entries = []
    if include_first_party and FIRST_PARTY_PACKS_DIR.is_dir():
        for child in sorted(FIRST_PARTY_PACKS_DIR.iterdir()):
            if not (child / "pack.json").is_file():
                continue
            try:
                manifest = read_pack_manifest(child)
            except PackError:
                continue
            entries.append(_describe_pack(child, manifest, origin="first_party"))
    root = packs_dir()
    if root.is_dir():
        for child in sorted(root.iterdir()):
            if not (child / "pack.json").is_file():
                continue
            try:
                manifest = read_pack_manifest(child)
            except PackError:
                continue
            entries.append(_describe_pack(child, manifest, origin="user"))
    return entries


def _describe_pack(pack_dir: Path, manifest: dict, origin: str) -> dict:
    metadata = _installed_metadata(pack_dir)
    sig_status = "unsigned"
    sig_path = pack_dir / "pack.sig.json"
    if sig_path.is_file():
        try:
            from .pack_signing import verify_pack_signature

            sig_status = verify_pack_signature(pack_dir).get("status", "unknown")
        except Exception:
            sig_status = "unreadable"
    return {
        "name": str(manifest.get("name", pack_dir.name)),
        "version": str(manifest.get("version", "")),
        "description": str(manifest.get("description", "")),
        "species": list(manifest.get("species", ["human"])),
        "tissues": list(manifest.get("tissues", [])),
        "diseases": list(manifest.get("diseases", []) or []),
        "pack_kind": str(manifest.get("pack_kind", "evidence")),
        "license": str(manifest.get("license") or manifest.get("license_spdx") or ""),
        "license_tier": str(manifest.get("license_tier", "community")),
        "ontology": manifest.get("ontology") or {},
        "provenance": manifest.get("provenance") or {},
        "runtime_gates": list(
            manifest.get("runtime_gates")
            or ["marker_evidence", "critic", "abstention", "conflict_detection"]
        ),
        "signature_status": sig_status,
        "trust": metadata.get("trust", manifest.get("trust", TRUST_ATLAS)),
        "origin": origin,
        "path": str(pack_dir),
        "code_policy": "data_only_no_executables",
    }


def remove_pack(name: str) -> Path:
    """Remove a user-installed pack. First-party packs cannot be removed."""
    if (FIRST_PARTY_PACKS_DIR / name / "pack.json").is_file():
        raise PackError(f"Pack {name!r} is first-party and cannot be removed")
    target = packs_dir() / name
    if not target.is_dir():
        raise PackError(f"Pack {name!r} is not installed")
    shutil.rmtree(target)
    return target


def _verify_installed_hashes(pack_dir: Path, metadata: dict) -> None:
    for filename, expected in (metadata.get("content_sha256") or {}).items():
        path = pack_dir / filename
        if not path.is_file():
            raise PackError(f"Installed pack file missing after install: {path}")
        if _file_sha256(path) != expected:
            raise PackError(
                f"Installed pack file was modified after install: {path}. "
                "Reinstall the pack to restore its recorded hash."
            )


def _load_pack_content(pack_dir: Path, manifest: dict, trust: str) -> dict:
    """Load and re-validate pack data files (fail closed at load time)."""
    record = {
        "name": str(manifest.get("name", pack_dir.name)),
        "version": str(manifest.get("version", "")),
        "trust": trust,
        "license_tier": str(manifest.get("license_tier", "community")),
        "files": {},
        "marker_atlas": None,
        "state_atlas": None,
        "ontology_map": None,
        "reference_manifest": None,
    }
    issues = validate_pack(pack_dir)
    if issues and trust == TRUST_ATLAS:
        raise PackError(
            f"Pack {record['name']!r} no longer passes validation: " + "; ".join(issues[:5])
        )
    for filename in ALLOWED_PACK_FILES:
        path = pack_dir / filename
        if not path.is_file():
            continue
        record["files"][filename] = _file_sha256(path)
        content = json.loads(path.read_text(encoding="utf-8"))
        if filename == ATLAS_FILE:
            record["marker_atlas"] = content
        elif filename == STATE_FILE:
            record["state_atlas"] = content
        elif filename == ONTOLOGY_MAP_FILE:
            record["ontology_map"] = content
        elif filename == REFERENCE_MANIFEST_FILE:
            record["reference_manifest"] = content
    if record["marker_atlas"] is None and record["state_atlas"] is None:
        raise PackError(f"Pack {record['name']!r} provides no atlas data")
    return record


def resolve_extension_packs(
    names: list[str], species: str
) -> tuple[list[dict], list[str]]:
    """Resolve requested pack names to loaded, license-checked records.

    Returns (records, warnings). Packs whose species scope does not include
    the run species are skipped with a warning; unknown names fail closed.
    """
    warnings: list[str] = []
    if not names:
        return [], warnings
    index = {entry["name"]: entry for entry in list_installed_packs()}
    records = []
    seen = set()
    for name in names:
        if name in seen:
            continue
        seen.add(name)
        entry = index.get(name)
        if entry is None:
            available = sorted(index) or ["(none installed)"]
            raise PackError(
                f"Extension pack {name!r} is not installed. Available: {', '.join(available)}"
            )
        if species not in entry["species"]:
            warnings.append(
                f"Pack {name!r} skipped: species scope {entry['species']} "
                f"does not include {species!r}"
            )
            continue
        pack_dir = Path(entry["path"])
        manifest = read_pack_manifest(pack_dir)
        _check_pack_license(manifest, "use")
        metadata = _installed_metadata(pack_dir)
        if metadata:
            _verify_installed_hashes(pack_dir, metadata)
        trust = metadata.get("trust", entry["trust"]) if metadata else entry["trust"]
        records.append(_load_pack_content(pack_dir, manifest, trust))
    return records, warnings


def _tag_hypothesis_cell_types(cell_types: dict) -> None:
    """Mark pack cell types as unreviewed context so the existing trust
    boundary applies: they can never silently produce an accepted identity."""
    for info in cell_types.values():
        positive = list(info.get("positive_markers", []))
        negative = list(info.get("negative_markers", []))
        info["context_origin"] = True
        info["context_review_status"] = "draft"
        info["context_positive_markers"] = positive
        info["context_negative_markers"] = negative
        # Force atlas_supporting to be empty so support from these markers is
        # treated as context-only support by the scorer and critic.
        info["atlas_positive_markers"] = []
        _tag_hypothesis_cell_types(info.get("subtypes", {}))


def merge_marker_atlas(
    base_atlas: dict,
    records: list[dict],
    species: str = "human",
) -> tuple[dict, list[str]]:
    """Merge extension-pack marker atlases into the built-in atlas.

    Deterministic precedence: built-in wins. Pack tissues absent from the
    base are added wholesale; for shared tissues only top-level cell types
    are merged, and a name collision leaves the built-in entry untouched
    (the pack entry is reported as shadowed).
    """
    from .data_adapter import convert_atlas_to_mouse

    merged = deepcopy(base_atlas)
    warnings: list[str] = []
    for record in records:
        pack_atlas = record.get("marker_atlas")
        if not pack_atlas:
            continue
        name = record["name"]
        if species == "mouse":
            pack_atlas = convert_atlas_to_mouse(pack_atlas)
        else:
            pack_atlas = deepcopy(pack_atlas)
        for tissue, tissue_data in pack_atlas.get("tissues", {}).items():
            if tissue not in merged.get("tissues", {}):
                added = deepcopy(tissue_data)
                if record["trust"] == TRUST_HYPOTHESIS:
                    _tag_hypothesis_cell_types(added.get("cell_types", {}))
                merged.setdefault("tissues", {})[tissue] = added
                continue
            existing = merged["tissues"][tissue].setdefault("cell_types", {})
            for ct_name, ct_info in tissue_data.get("cell_types", {}).items():
                if ct_name in existing:
                    warnings.append(
                        f"Pack {name!r}: cell type {ct_name!r} in tissue {tissue!r} "
                        "is shadowed by the built-in atlas"
                    )
                    continue
                entry = deepcopy(ct_info)
                if record["trust"] == TRUST_HYPOTHESIS:
                    _tag_hypothesis_cell_types({ct_name: entry})
                existing[ct_name] = entry

    conflicts = detect_marker_conflicts(merged)
    for c in conflicts:
        if c.severity == "high":
            warnings.append(f"High severity conflict introduced by pack: {c.conflict_type} on {c.gene} between {c.cell_type_a} and {c.cell_type_b}")

    return merged, warnings


def pack_conflict_report(pack_dir: str | Path, base_atlas: dict) -> dict:
    """Generate a conflict report for a pack against a base atlas."""
    pack_dir = Path(pack_dir)
    manifest = read_pack_manifest(pack_dir)
    record = _load_pack_content(pack_dir, manifest, TRUST_ATLAS)

    merged, _ = merge_marker_atlas(base_atlas, [record])
    conflicts = detect_marker_conflicts(merged)

    return {
        "pack_name": record["name"],
        "conflicts": [vars(c) for c in conflicts]
    }


def collect_pack_state_definitions_input(records: list[dict]) -> list[dict]:
    """Extract state atlases from pack records for load_state_definitions."""
    return [
        {
            "pack_name": record["name"],
            "pack_version": record["version"],
            "trust": record["trust"],
            "atlas": record["state_atlas"],
        }
        for record in records
        if record.get("state_atlas")
    ]


def pack_manifest_parameters(records: list[dict]) -> list[dict]:
    """Reproducibility-critical pack fields for the run manifest."""
    return [
        {
            "name": record["name"],
            "version": record["version"],
            "trust": record["trust"],
            "license_tier": record["license_tier"],
            "files": record["files"],
        }
        for record in records
    ]
