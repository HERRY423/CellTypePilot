"""One-off: batch-fix incorrect CL identifiers found by live ontology check.

Every fix below was flagged by `celltypepilot ontology check` as a label
mismatch whose current CL term is semantically unrelated to the atlas key
(e.g. Astrocyte resolving to "natural killer cell"). Replacement terms were
resolved against the live Cell Ontology (cl.obo, downloaded 2026-08).

Keys whose CL term is a genuine refinement (e.g. "Alveolar type 1 cell" vs
"pulmonary alveolar type 1 cell") are intentionally left as warnings.

Run: python scripts/fix_atlas_cl_ids_batch.py
"""

import json
from pathlib import Path

from celltypepilot.ontology import load_ontology

MARKER_ATLAS = Path("src/celltypepilot/data/marker_atlas.json")
PREMIUM_ATLAS = Path("src/celltypepilot/data/packs/premium/marker_atlas.json")

# (tissue, cell_path, new_cl_id, expected_label_fragment)
MARKER_FIXES = [
    ("blood", ["T cell", "CD4+ T cell", "Regulatory T cell"], "CL:0000815", "regulatory T cell"),
    ("blood", ["T cell", "CD4+ T cell", "Th17 cell"], "CL:0000899", "T-helper 17"),
    ("blood", ["Dendritic cell", "cDC1"], "CL:0002394", "CD141-positive"),
    ("blood", ["Dendritic cell", "cDC2"], "CL:0002399", "CD1c-positive"),
    ("blood", ["Monocyte", "Non-classical monocyte"], "CL:0000875", "non-classical monocyte"),
    ("lung", ["Basal cell"], "CL:0002633", "respiratory basal cell"),
    ("liver", ["Liver endothelial cell"], "CL:1000398", "hepatic sinusoid"),
    ("brain", ["Excitatory neuron"], "CL:0000679", "glutamatergic"),
    ("brain", ["Inhibitory neuron"], "CL:0000617", "GABAergic"),
    ("brain", ["Astrocyte"], "CL:0000127", "astrocyte"),
    ("brain", ["Oligodendrocyte"], "CL:0000128", "oligodendrocyte"),
    ("kidney", ["Collecting duct cell"], "CL:1000454", "collecting duct"),
    ("gut", ["Goblet cell (intestinal)"], "CL:0000160", "goblet cell"),
    ("gut", ["Intestinal stem cell"], "CL:0002250", "intestinal crypt stem cell"),
    ("skin", ["Dermal fibroblast"], "CL:0002551", "fibroblast of dermis"),
    ("heart", ["Cardiac endothelial cell"], "CL:0010006", "cardiac blood vessel"),
    ("heart", ["Cardiac smooth muscle cell"], "CL:0000359", "vascular associated smooth muscle"),
    ("general", ["Pericyte"], "CL:0000669", "pericyte"),
    ("general", ["Mast cell"], "CL:0000097", "mast cell"),
    ("muscle", ["Myofiber"], "CL:0008002", "skeletal muscle fiber"),
]

PREMIUM_FIXES = [
    ("tumor_microenvironment", ["myeloid_derived_suppressor_cell"], "CL:0000889", "myeloid suppressor"),
    ("developing_brain", ["radial_glia"], "CL:0000681", "radial glial"),
]


def apply_fixes(path: Path, fixes: list, service) -> int:
    atlas = json.loads(path.read_text(encoding="utf-8"))
    changed = 0
    for tissue, cell_path, new_id, fragment in fixes:
        node = atlas["tissues"][tissue]["cell_types"]
        for step in cell_path[:-1]:
            node = node[step]["subtypes"]
        info = node[cell_path[-1]]
        old_id = info.get("cl_id")
        term = service.resolve(new_id)
        assert term is not None and fragment.lower() in term.name.lower(), (
            f"replacement check failed for {new_id}: {term and term.name}"
        )
        assert old_id != new_id, f"{tissue}/{cell_path}: already {new_id}"
        info["cl_id"] = new_id
        ontology = info.get("ontology_evidence", {})
        if ontology.get("cl_id") == old_id:
            ontology["cl_id"] = new_id
        print(f"{tissue}/{'/'.join(cell_path)}: {old_id} -> {new_id} ({term.name})")
        changed += 1
    path.write_text(json.dumps(atlas, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return changed


def main() -> None:
    service = load_ontology()
    marker_changed = apply_fixes(MARKER_ATLAS, MARKER_FIXES, service)
    premium_changed = 0
    if PREMIUM_ATLAS.is_file():
        premium_changed = apply_fixes(PREMIUM_ATLAS, PREMIUM_FIXES, service)
    print(f"marker_atlas fixes: {marker_changed}; premium fixes: {premium_changed}")


if __name__ == "__main__":
    main()
