"""One-off: fix broken CL identifiers found by live ontology check.

Findings from `celltypepilot ontology check` against cl.obo (2026-08):

- liver/Cholangiocyte used CL:0000063, which the live Cell Ontology lists
  as "obsolete cell by histology". The correct current term for
  cholangiocyte is CL:1000488.
- kidney/Loop of Henle cell used CL:1000935, which does not exist in the
  current Cell Ontology. The correct general term is CL:1000909
  (kidney loop of Henle epithelial cell).

Run: python scripts/fix_atlas_cl_ids.py
"""

import json
from pathlib import Path

ATLAS = Path("src/celltypepilot/data/marker_atlas.json")

FIXES = [
    ("liver", "Cholangiocyte", "CL:0000063", "CL:1000488"),
    ("kidney", "Loop of Henle cell", "CL:1000935", "CL:1000909"),
]


def main() -> None:
    atlas = json.loads(ATLAS.read_text(encoding="utf-8"))
    for tissue, cell_type, old_id, new_id in FIXES:
        info = atlas["tissues"][tissue]["cell_types"][cell_type]
        assert info["cl_id"] == old_id, f"unexpected cl_id for {tissue}/{cell_type}"
        info["cl_id"] = new_id
        ontology = info.get("ontology_evidence", {})
        if ontology.get("cl_id") == old_id:
            ontology["cl_id"] = new_id
        print(f"{tissue}/{cell_type}: {old_id} -> {new_id}")
    ATLAS.write_text(json.dumps(atlas, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print("atlas written")


if __name__ == "__main__":
    main()
