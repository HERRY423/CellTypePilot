"""Atlas governance reporting for CellTypePilot.

The governance report is an offline, machine-readable asset summary. It is
designed for plugin hosts and reviewers that need to know what the atlas can
support, what remains draft-level, and which maintenance checks are pending.
"""

from __future__ import annotations

import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from .constants import ANNOTATION_SUPPORTED_SPECIES, ATLAS_PATH, FIRST_PARTY_PACKS_DIR
from .data_adapter import summarize_atlas_evidence, validate_atlas_provenance
from .atlas_conflict import detect_marker_conflicts
from .atlas_lifecycle import sunset_check
from .atlas_curation import build_curation_queue

GOVERNANCE_SCHEMA_VERSION = "celltypepilot.atlas-governance.v1"


def _iter_nodes(atlas: dict):
    def walk(cell_types: dict, tissue: str, parents: tuple[str, ...]):
        for name, info in cell_types.items():
            path = (*parents, name)
            yield tissue, path, info
            yield from walk(info.get("subtypes", {}), tissue, path)

    for tissue, tissue_info in atlas.get("tissues", {}).items():
        yield from walk(tissue_info.get("cell_types", {}), tissue, ())


def _atlas_counts(atlas: dict) -> dict:
    node_counts = Counter()
    relationship_counts = Counter()
    status_counts = Counter()
    species_counts = Counter()
    max_depth = 0
    missing_cl_id = 0
    deprecated_nodes = 0
    deprecated_edges = 0

    for tissue, path, info in _iter_nodes(atlas):
        node_counts[tissue] += 1
        max_depth = max(max_depth, len(path))
        if not info.get("cl_id"):
            missing_cl_id += 1
        if info.get("deprecated"):
            deprecated_nodes += 1
        for record in info.get("marker_evidence", []):
            relationship_counts[tissue] += 1
            if record.get("deprecated"):
                deprecated_edges += 1
            status_counts[record.get("verification_status", "missing")] += 1
            for species in record.get("species", []):
                species_counts[str(species)] += 1

    return {
        "tissues": sorted(atlas.get("tissues", {})),
        "n_tissues": len(atlas.get("tissues", {})),
        "n_cell_type_nodes": sum(node_counts.values()),
        "n_marker_relationships": sum(relationship_counts.values()),
        "max_cell_type_depth": max_depth,
        "missing_cl_id_nodes": missing_cl_id,
        "deprecated_nodes": deprecated_nodes,
        "deprecated_edges": deprecated_edges,
        "nodes_by_tissue": dict(sorted(node_counts.items())),
        "relationships_by_tissue": dict(sorted(relationship_counts.items())),
        "verification_status_counts": dict(sorted(status_counts.items())),
        "relationship_species_counts": dict(sorted(species_counts.items())),
    }


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _describe_atlas(path: Path, label: str, origin: str) -> dict:
    atlas = _load_json(path)
    provenance_issues = validate_atlas_provenance(atlas)
    evidence = summarize_atlas_evidence(atlas)
    counts = _atlas_counts(atlas)
    
    current_version = atlas.get("version", "1.0")
    sunsets = sunset_check(atlas, current_version)
    
    conflicts = detect_marker_conflicts(atlas)
    high_conflicts = [c for c in conflicts if c.severity == "high"]
    
    curation_q = build_curation_queue(atlas)
    high_priority_curation = len(curation_q[curation_q["priority"] >= 90]) if not curation_q.empty else 0
    
    health_score = 100
    health_score -= min(30, len(provenance_issues) * 2)
    health_score -= min(40, len(high_conflicts) * 5)
    health_score -= min(20, len(sunsets) * 2)
    health_score -= min(10, high_priority_curation // 10)
    health_score = max(0, health_score)

    return {
        "label": label,
        "origin": origin,
        "path": str(path),
        "version": atlas.get("version", ""),
        "schema_version": atlas.get("schema_version", ""),
        "health_score": health_score,
        "provenance_validation": "passed" if not provenance_issues else "failed",
        "provenance_issue_count": len(provenance_issues),
        "provenance_issue_examples": provenance_issues[:10],
        "evidence_summary": evidence,
        "counts": counts,
        "conflicts": {
            "total": len(conflicts),
            "high_severity": len(high_conflicts),
            "examples": [vars(c) for c in high_conflicts[:5]]
        },
        "lifecycle": {
            "sunset_issues": len(sunsets),
            "sunset_examples": sunsets[:5]
        },
        "curation_queue": {
            "total_pending": len(curation_q),
            "high_priority": high_priority_curation
        }
    }


def _first_party_pack_atlases() -> list[tuple[str, Path]]:
    if not FIRST_PARTY_PACKS_DIR.is_dir():
        return []
    entries = []
    for child in sorted(FIRST_PARTY_PACKS_DIR.iterdir()):
        atlas_path = child / "marker_atlas.json"
        if atlas_path.is_file():
            entries.append((f"pack:{child.name}", atlas_path))
    return entries


def _pack_inventory() -> list[dict]:
    try:
        from .pack_manager import list_installed_packs, validate_pack
    except Exception as exc:
        return [{"status": "unavailable", "error": str(exc)}]

    entries = []
    for item in list_installed_packs():
        issues = validate_pack(Path(item["path"]))
        entry = dict(item)
        entry["validation"] = "passed" if not issues else "failed"
        entry["validation_issue_count"] = len(issues)
        entry["validation_issue_examples"] = issues[:5]
        entries.append(entry)
    return entries


def _ontology_section() -> dict:
    try:
        from .ontology import (
            check_atlas_ontology,
            load_ontology,
            ontology_cache_status,
            summarize_findings,
        )
    except Exception as exc:
        return {"available": False, "error": str(exc)}

    status = ontology_cache_status()
    section = {"available": True, "cache": status, "checks": {}}
    if not status.get("cached"):
        section["claim_boundary"] = "Ontology validation not run because no cache is present."
        return section
    try:
        service = load_ontology()
        targets = [("builtin", ATLAS_PATH), *_first_party_pack_atlases()]
        for label, path in targets:
            atlas = _load_json(path)
            findings = check_atlas_ontology(service, atlas)
            section["checks"][label] = {
                "summary": summarize_findings(findings),
                "finding_examples": findings[:10],
            }
    except Exception as exc:
        section["error"] = str(exc)
    return section


def build_atlas_governance_report(include_packs: bool = True) -> dict:
    """Build an offline governance report for bundled and installed atlas assets."""
    atlas_assets = [_describe_atlas(ATLAS_PATH, "builtin", "builtin")]
    if include_packs:
        for label, path in _first_party_pack_atlases():
            atlas_assets.append(_describe_atlas(path, label, "first_party_pack"))

    aggregate_status = Counter()
    aggregate_relationships = 0
    for asset in atlas_assets:
        counts = asset["counts"]["verification_status_counts"]
        aggregate_status.update(counts)
        aggregate_relationships += asset["counts"]["n_marker_relationships"]

    return {
        "schema_version": GOVERNANCE_SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "supported_annotation_species": list(ANNOTATION_SUPPORTED_SPECIES),
        "claim_boundary": (
            "This report summarizes atlas governance state. It is not biological "
            "validation, batch robustness evidence, or primary-source review."
        ),
        "governance_health_score": float(sum(asset["health_score"] for asset in atlas_assets) / len(atlas_assets)) if atlas_assets else 100.0,
        "atlas_assets": atlas_assets,
        "aggregate": {
            "n_assets": len(atlas_assets),
            "n_marker_relationships": aggregate_relationships,
            "verification_status_counts": dict(sorted(aggregate_status.items())),
            "needs_edge_curation": int(
                aggregate_status.get("aggregate_source_only_not_edge_verified", 0)
            ),
        },
        "ontology": _ontology_section(),
        "extension_packs": _pack_inventory() if include_packs else [],
        "governance_actions": [
            "Run ontology update/check before releases.",
            "Publish curation queue and sweep reports with every atlas version.",
            "Keep unsupported species fail-closed until species-scoped atlas packs pass governance.",
            "Treat literature co-occurrence as supporting discussion, not marker specificity proof.",
        ],
    }


def write_atlas_governance_report(output_path: str | Path, include_packs: bool = True) -> Path:
    """Write the governance report to JSON and return its path."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    report = build_atlas_governance_report(include_packs=include_packs)
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path
