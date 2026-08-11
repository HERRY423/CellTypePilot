"""Immutable release freeze for CellTypePilot's governed decision layer."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

GOVERNANCE_FREEZE_SCHEMA = "celltypepilot.governance-freeze.v2"
GOVERNANCE_HASH_SEMANTICS = "utf8_text_lf_v1"


class GovernanceFreezeError(ValueError):
    """Raised when a frozen governance artifact is incomplete or changed."""


def _normalized_content(path: Path) -> bytes:
    """Return governed UTF-8 text with platform-independent LF endings."""
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise GovernanceFreezeError(
            f"Governed file must be UTF-8 text for {GOVERNANCE_HASH_SEMANTICS}: {path}"
        ) from exc
    return text.encode("utf-8")


def _sha256(path: Path) -> str:
    return hashlib.sha256(_normalized_content(path)).hexdigest()


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def governed_paths(root: Path | None = None) -> list[Path]:
    """Return release-governed code/data paths, excluding the freeze itself."""
    base = (root or _repo_root()).resolve()
    fixed = [
        "src/celltypepilot/calibration.py",
        "src/celltypepilot/calibration_split.py",
        "src/celltypepilot/candidate_backends.py",
        "src/celltypepilot/critic.py",
        "src/celltypepilot/hierarchical_selector.py",
        "src/celltypepilot/identity_contract.py",
        "src/celltypepilot/lineage_coverage.py",
        "src/celltypepilot/novelty_detector.py",
        "src/celltypepilot/pack_manager.py",
        "src/celltypepilot/state_scorer.py",
        "src/celltypepilot/validation_domains.py",
        "src/celltypepilot/data/marker_atlas.json",
        "src/celltypepilot/data/state_atlas.json",
        "src/celltypepilot/data/validation_domains.json",
    ]
    paths = [base / value for value in fixed]
    pack_root = base / "src/celltypepilot/data/packs"
    if pack_root.is_dir():
        paths.extend(
            child
            for child in pack_root.rglob("*")
            if child.is_file()
            and child.name
            in {
                "pack.json",
                "marker_atlas.json",
                "state_atlas.json",
                "ontology_map.json",
                "reference_manifest.json",
            }
        )
    return sorted(set(paths))


def build_governance_freeze(
    output_path: str | Path,
    *,
    root: str | Path | None = None,
    release_id: str,
) -> dict:
    """Hash the current decision layer into a deterministic release artifact."""
    base = Path(root).resolve() if root else _repo_root()
    records = []
    for path in governed_paths(base):
        if not path.is_file():
            raise GovernanceFreezeError(f"Governed file is missing: {path}")
        records.append(
            {
                "path": path.relative_to(base).as_posix(),
                "normalized_size": len(_normalized_content(path)),
                "sha256": _sha256(path),
            }
        )
    payload = {
        "schema_version": GOVERNANCE_FREEZE_SCHEMA,
        "hash_semantics": GOVERNANCE_HASH_SEMANTICS,
        "release_id": release_id,
        "files": records,
        "frozen_invariants": {
            "marker_scorer_role": "evidence_only_not_candidate_vote",
            "candidate_decision_backends": [
                "celltypist",
                "popv",
                "singler",
                "scanvi",
                "custom_reference",
            ],
            "llm_role": "hypothesis_only_by_default",
            "minimum_independent_backend_groups": 2,
            "selector_outputs": ["accepted_leaf", "accepted_ancestor", "abstain"],
            "calibration_effect": "downgrade_only",
            "identity_state_novelty": "independent_axes",
            "human_signoff": "required",
        },
        "claim_boundary": (
            "This freezes software governance and content hashes. It does not freeze or establish "
            "biological accuracy, selective-risk calibration, or domain validation."
        ),
    }
    canonical = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    payload["freeze_sha256"] = hashlib.sha256(canonical).hexdigest()
    output = Path(output_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return payload


def verify_governance_freeze(
    freeze_path: str | Path,
    *,
    root: str | Path | None = None,
) -> dict:
    """Fail closed if a governed release file differs from the frozen hash."""
    freeze_path = Path(freeze_path).resolve()
    payload = json.loads(freeze_path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != GOVERNANCE_FREEZE_SCHEMA:
        raise GovernanceFreezeError(f"Freeze schema must be {GOVERNANCE_FREEZE_SCHEMA}")
    if payload.get("hash_semantics") != GOVERNANCE_HASH_SEMANTICS:
        raise GovernanceFreezeError(f"Freeze hash semantics must be {GOVERNANCE_HASH_SEMANTICS}")
    unsigned = {key: value for key, value in payload.items() if key != "freeze_sha256"}
    canonical = json.dumps(
        unsigned, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    if hashlib.sha256(canonical).hexdigest() != payload.get("freeze_sha256"):
        raise GovernanceFreezeError("Freeze manifest hash mismatch")
    expected = {record["path"]: record for record in payload.get("files", [])}
    if root is not None:
        base = Path(root).resolve()
        current_paths = {path.relative_to(base).as_posix(): path for path in governed_paths(base)}
    else:
        source_root = _repo_root()
        if (source_root / "src/celltypepilot").is_dir():
            current_paths = {
                path.relative_to(source_root).as_posix(): path
                for path in governed_paths(source_root)
            }
        else:
            package_root = Path(__file__).resolve().parent
            prefix = "src/celltypepilot/"
            current_paths = {
                relative: package_root / relative.removeprefix(prefix)
                for relative in expected
                if relative.startswith(prefix)
            }
    if set(expected) != set(current_paths):
        missing = sorted(set(expected) - set(current_paths))
        added = sorted(set(current_paths) - set(expected))
        raise GovernanceFreezeError(
            f"Governed file set changed; missing={missing[:3]} added={added[:3]}"
        )
    changed = []
    for relative, record in expected.items():
        path = current_paths[relative]
        normalized = _normalized_content(path)
        if (
            len(normalized) != record["normalized_size"]
            or hashlib.sha256(normalized).hexdigest() != record["sha256"]
        ):
            changed.append(relative)
    if changed:
        raise GovernanceFreezeError(f"Governed files changed after freeze: {changed[:5]}")
    return {
        "schema_version": GOVERNANCE_FREEZE_SCHEMA,
        "hash_semantics": GOVERNANCE_HASH_SEMANTICS,
        "status": "verified",
        "release_id": payload["release_id"],
        "freeze_sha256": payload["freeze_sha256"],
        "n_files": len(expected),
    }
