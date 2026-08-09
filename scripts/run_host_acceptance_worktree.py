#!/usr/bin/env python3
"""Create an independent worktree and run host acceptance there.

Usage:
  python scripts/run_host_acceptance_worktree.py
  python scripts/run_host_acceptance_worktree.py --worktree ../ctp-host-acceptance
  python scripts/run_host_acceptance_worktree.py --skip-worktree
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from celltypepilot.host_acceptance import run_host_acceptance  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--worktree",
        type=Path,
        default=REPO.parent / "celltypepilot-host-acceptance",
        help="Independent worktree path (default: sibling directory)",
    )
    parser.add_argument(
        "--skip-worktree",
        action="store_true",
        help="Run fixtures under scratch/ in the current clone",
    )
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    report = run_host_acceptance(
        worktree=None if args.skip_worktree else args.worktree,
        skip_worktree=args.skip_worktree,
    )
    if args.json:
        print(json.dumps(report, indent=2, default=str))
    else:
        print(f"overall_status={report['overall_status']}")
        print(f"report_path={report.get('report_path')}")
        for check in report["checks"]:
            flag = "PASS" if check["passed"] else "FAIL"
            print(f"  [{flag}] {check['id']}: {check['detail']}")
        disc = (report.get("agent_discrimination") or {}).get("states_observed") or {}
        print(f"discrimination={disc}")
    return 0 if report["overall_status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
