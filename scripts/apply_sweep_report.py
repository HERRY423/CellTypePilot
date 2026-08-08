"""Apply a previously-generated sweep report to the atlas.

The CLI's `--apply` flag applies the in-memory sweep results from the same
run. This script lets us apply a saved sweep report (e.g. after review)
without re-running the PubMed queries.

Usage: python scripts/apply_sweep_report.py docs/curate/sweep_full.json
"""

import json
import sys
from pathlib import Path

from celltypepilot.atlas_curation import apply_sweep_results
from celltypepilot.constants import ATLAS_PATH


def main(report_path: str, new_version: str) -> None:
    report = json.loads(Path(report_path).read_text(encoding="utf-8"))
    atlas = json.loads(ATLAS_PATH.read_text(encoding="utf-8"))
    updated, applied = apply_sweep_results(atlas, report["results"], new_version)
    ATLAS_PATH.write_text(
        json.dumps(updated, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"Applied {applied} upgrades to {ATLAS_PATH} (version {new_version})")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <sweep_report.json> <new_version>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
