"""Export unverified marker edges as a curator-owned CSV work queue."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ATLASES = (
    ROOT / "src" / "celltypepilot" / "data" / "marker_atlas.json",
    ROOT / "src" / "celltypepilot" / "data" / "packs" / "premium" / "marker_atlas.json",
)


def iter_nodes(cell_types: dict, parents: tuple[str, ...] = ()):
    for name, info in cell_types.items():
        path = (*parents, name)
        yield path, info
        yield from iter_nodes(info.get("subtypes", {}), path)


def export_queue(atlas_paths: tuple[Path, ...], output: Path) -> int:
    rows = []
    for atlas_path in atlas_paths:
        atlas = json.loads(atlas_path.read_text(encoding="utf-8"))
        for tissue, tissue_info in atlas.get("tissues", {}).items():
            for cell_path, info in iter_nodes(tissue_info.get("cell_types", {})):
                for record in info.get("marker_evidence", []):
                    if (
                        record.get("verification_status")
                        != "aggregate_source_only_not_edge_verified"
                    ):
                        continue
                    rows.append(
                        {
                            "atlas_file": atlas_path.name,
                            "atlas_version": atlas.get("version", ""),
                            "tissue": tissue,
                            "cell_type_path": " > ".join(cell_path),
                            "gene": record.get("gene", ""),
                            "polarity": record.get("polarity", ""),
                            "species": ";".join(record.get("species", [])),
                            "state": record.get("state", ""),
                            "target_status": "",
                            "source_type": "",
                            "source_name": "",
                            "pmid": "",
                            "doi": "",
                            "url": "",
                            "source_record_id": "",
                            "source_record_url": "",
                            "evidence_locator": "",
                            "curator": "",
                            "verified_at": "",
                            "curator_notes": "",
                        }
                    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["atlas_file"])
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "atlas_curation_queue.csv")
    parser.add_argument("--atlas", action="append", type=Path)
    args = parser.parse_args()
    paths = tuple(args.atlas) if args.atlas else DEFAULT_ATLASES
    count = export_queue(paths, args.output)
    print(f"exported {count} unverified marker edges to {args.output}")


if __name__ == "__main__":
    main()
