"""Add explicit, honest per-marker provenance records to bundled atlases."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ATLAS_PATHS = (
    ROOT / "src" / "celltypepilot" / "data" / "marker_atlas.json",
    ROOT / "src" / "celltypepilot" / "data" / "packs" / "premium" / "marker_atlas.json",
)

SOURCE_REGISTRY = {
    "cellmarker_2_0": {
        "source_id": "cellmarker_2_0",
        "name": "CellMarker 2.0",
        "pmid": "36300619",
        "doi": "10.1093/nar/gkac947",
        "url": "https://pubmed.ncbi.nlm.nih.gov/36300619/",
    },
    "panglaodb": {
        "source_id": "panglaodb",
        "name": "PanglaoDB",
        "pmid": "30951143",
        "doi": "10.1093/database/baz046",
        "url": "https://pubmed.ncbi.nlm.nih.gov/30951143/",
    },
    "cell_ontology": {
        "source_id": "cell_ontology",
        "name": "Cell Ontology 2016",
        "pmid": "27377652",
        "doi": "10.1186/s13326-016-0088-7",
        "url": "https://pubmed.ncbi.nlm.nih.gov/27377652/",
    },
}


def _state_for_tissue(tissue: str) -> str:
    return {
        "tumor_microenvironment": "tumor",
        "developing_brain": "developmental",
        "inflamed_tissue": "inflamed",
        "immune_activation": "activated",
    }.get(tissue, "baseline_or_unspecified")


def _walk_cell_types(cell_types: dict):
    for cell_type, info in cell_types.items():
        yield cell_type, info
        yield from _walk_cell_types(info.get("subtypes", {}))


def _relationship_record(
    gene: str,
    polarity: str,
    tissue: str,
    version: str,
) -> dict:
    return {
        "gene": gene,
        "polarity": polarity,
        "species": ["human", "mouse"],
        "tissue": tissue,
        "state": _state_for_tissue(tissue),
        "atlas_version": version,
        "sources": [SOURCE_REGISTRY["cellmarker_2_0"], SOURCE_REGISTRY["panglaodb"]],
        "evidence_scope": "database_level_source",
        "verification_status": "aggregate_source_only_not_edge_verified",
    }


def enrich(path: Path) -> None:
    atlas = json.loads(path.read_text(encoding="utf-8"))
    version = atlas["version"]
    atlas["schema_version"] = "celltypepilot.marker-atlas.v2"
    atlas["source_registry"] = SOURCE_REGISTRY
    atlas["provenance_policy"] = {
        "edge_fields_required": [
            "gene",
            "polarity",
            "species",
            "tissue",
            "state",
            "atlas_version",
            "sources",
            "evidence_scope",
            "verification_status",
        ],
        "aggregate_source_only_not_edge_verified": (
            "The cited database papers establish source provenance but do not prove that "
            "each relationship was independently traced to a marker-specific experiment."
        ),
        "database_record_verified": (
            "A curator verified the exact relationship against a stable database record."
        ),
        "primary_source_verified": (
            "A curator verified the exact relationship at a primary-source evidence locator."
        ),
    }
    for tissue, tissue_info in atlas.get("tissues", {}).items():
        for _cell_type, info in _walk_cell_types(tissue_info.get("cell_types", {})):
            existing = {
                (record.get("gene"), record.get("polarity")): record
                for record in info.get("marker_evidence", [])
            }
            relationships = []
            for polarity, key in (
                ("positive", "positive_markers"),
                ("negative", "negative_markers"),
            ):
                for gene in info.get(key, []):
                    prior = existing.get((gene, polarity))
                    relationships.append(
                        prior
                        if prior
                        and prior.get("verification_status")
                        in {"database_record_verified", "primary_source_verified"}
                        else _relationship_record(gene, polarity, tissue, version)
                    )
            info["marker_evidence"] = relationships
            if info.get("cl_id"):
                info["ontology_evidence"] = {
                    "cl_id": info["cl_id"],
                    "source": SOURCE_REGISTRY["cell_ontology"],
                    "verification_status": "identifier_format_recorded_not_live_resolved",
                }
    path.write_text(json.dumps(atlas, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    for atlas_path in ATLAS_PATHS:
        enrich(atlas_path)
