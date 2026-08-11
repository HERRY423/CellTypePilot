"""Build the deterministic CellTypePilot governance freeze."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from celltypepilot.governance_freeze import build_governance_freeze  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-id", required=True)
    parser.add_argument(
        "--output",
        default="src/celltypepilot/data/governance_freeze_v1.json",
    )
    args = parser.parse_args()
    output = Path(args.output)
    if not output.is_absolute():
        output = REPO_ROOT / output
    build_governance_freeze(output, root=REPO_ROOT, release_id=args.release_id)
    print(output)


if __name__ == "__main__":
    main()
