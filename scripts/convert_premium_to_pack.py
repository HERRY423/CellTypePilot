"""One-time migration: convert data/premium_atlas.json into a first-party
extension pack at data/packs/premium/ (pack spec v1).

The license gate moves from the atlas payload (license_required) into the
pack manifest (license_tier), which pack_manager enforces at install/resolve
time. The atlas payload keeps its provenance fields untouched so existing
validator contracts continue to pass.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "celltypepilot" / "data" / "premium_atlas.json"
DEST_DIR = ROOT / "src" / "celltypepilot" / "data" / "packs" / "premium"

sys.path.insert(0, str(ROOT / "src"))


def main() -> None:
    atlas = json.loads(SRC.read_text(encoding="utf-8"))
    license_tier = atlas.pop("license_required", "academic")

    DEST_DIR.mkdir(parents=True, exist_ok=True)

    manifest = {
        "schema_version": "celltypepilot.pack.v1",
        "name": "premium",
        "version": atlas["version"],
        "description": atlas.get(
            "description",
            "Premium identity atlas: TME, developing brain, inflamed tissue, immune activation",
        ),
        "species": ["human", "mouse"],
        "tissues": sorted(atlas.get("tissues", {})),
        "license_tier": license_tier,
        "files": ["marker_atlas.json"],
    }
    (DEST_DIR / "pack.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (DEST_DIR / "marker_atlas.json").write_text(
        json.dumps(atlas, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    # Verify pack passes the same gates the manager applies at install time.
    from celltypepilot.pack_manager import validate_pack

    issues = validate_pack(DEST_DIR)
    if issues:
        print("VALIDATION FAILED:")
        for issue in issues:
            print(" -", issue)
        raise SystemExit(1)
    print(f"Pack written to {DEST_DIR} and validated (0 issues).")


if __name__ == "__main__":
    main()
