"""Migrate premium identities and normalize identity/state relation versions."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ATLAS_PATH = ROOT / "src" / "celltypepilot" / "data" / "premium_atlas.json"
STATE_ATLAS_PATH = ROOT / "src" / "celltypepilot" / "data" / "state_atlas.json"
VERSION = "premium-atlas-2026.08.1"

MAPPINGS = {
    "cancer_associated_fibroblast": ("fibroblast", "CL:0000057", "cancer_associated", "role"),
    "tumor_associated_macrophage_M1": (
        "macrophage",
        "CL:0000235",
        "tumor_associated_inflammatory",
        "activation",
    ),
    "tumor_associated_macrophage_M2": (
        "macrophage",
        "CL:0000235",
        "tumor_associated_immunosuppressive",
        "activation",
    ),
    "regulatory_T_cell_tumor": (
        "regulatory T cell",
        "CL:0000815",
        "tumor_associated",
        "context",
    ),
    "exhausted_T_cell": ("T cell", "CL:0000084", "exhausted", "activation"),
    "intermediate_progenitor": (
        "neural progenitor cell",
        "CL:0011020",
        "intermediate",
        "development",
    ),
    "migrating_neuroblast": (
        "neuroblast (sensu Vertebrata)",
        "CL:0000031",
        "migrating",
        "development",
    ),
    "microglia_developing": (
        "microglial cell",
        "CL:0000129",
        "developmental",
        "development",
    ),
    "inflammatory_macrophage": (
        "macrophage",
        "CL:0000235",
        "inflammatory",
        "activation",
    ),
    "plasma_cell_inflamed": ("plasma cell", "CL:0000786", "inflamed", "context"),
    "activated_fibroblast": ("fibroblast", "CL:0000057", "activated", "activation"),
    "activated_CD4_T_cell": (
        "CD4-positive, alpha-beta T cell",
        "CL:0000624",
        "activated",
        "activation",
    ),
    "activated_CD8_T_cell": (
        "CD8-positive, alpha-beta T cell",
        "CL:0000625",
        "activated",
        "activation",
    ),
    "germinal_center_B_cell": (
        "germinal center B cell",
        "CL:0000844",
        "germinal_center",
        "context",
    ),
    "NK_cell_activated": ("natural killer cell", "CL:0000623", "activated", "activation"),
}


def main() -> None:
    atlas = json.loads(ATLAS_PATH.read_text(encoding="utf-8"))
    seen = set()
    atlas["version"] = VERSION
    atlas["schema_version"] = "celltypepilot.marker-atlas.v2.1"
    for tissue in atlas.get("tissues", {}).values():
        for node_name, node in tissue.get("cell_types", {}).items():
            if node_name in MAPPINGS:
                base_label, cl_id, state_label, state_axis = MAPPINGS[node_name]
                seen.add(node_name)
                node["base_cell_type"] = base_label
                node["cl_id"] = cl_id
                node["state_label"] = state_label
                node["state_axis"] = state_axis
                node["display_label"] = node_name
                ontology = node.setdefault("ontology_evidence", {})
                ontology["cl_id"] = cl_id
                ontology["label"] = base_label
                ontology["verification_status"] = "identifier_format_recorded_not_live_resolved"
            for record in node.get("marker_evidence", []):
                record["atlas_version"] = VERSION
    missing = sorted(set(MAPPINGS) - seen)
    if missing:
        raise RuntimeError(f"Premium atlas nodes not found: {missing}")
    ATLAS_PATH.write_text(
        json.dumps(atlas, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    state_atlas = json.loads(STATE_ATLAS_PATH.read_text(encoding="utf-8"))
    for state_name, definition in state_atlas.get("states", {}).items():
        for record in definition.get("marker_evidence", []):
            record["species"] = list(definition.get("species", []))
            record["tissues"] = list(definition.get("tissues", []))
            record["state"] = state_name
            record["atlas_version"] = state_atlas["version"]
    STATE_ATLAS_PATH.write_text(
        json.dumps(state_atlas, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
