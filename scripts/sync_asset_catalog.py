"""Verify and materialize the immutable asset catalog without touching fold runs.

Never writes under benchmarks/**/runs/**. Use this to:
  - validate catalog + storage policy schemas;
  - verify SHA-256 of local file: sources and object caches;
  - materialize label maps / Dockerfiles into content-addressed objects/.

Remote CDN upload is out of scope for this script (publish pipeline separate).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Allow running without install when repo root is on PYTHONPATH.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from celltypepilot.asset_catalog import (  # noqa: E402
    AssetCatalogError,
    load_asset_catalog,
    load_storage_policy,
    materialize_source_to_object_cache,
    summarize_catalog,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog",
        default=str(REPO_ROOT / "benchmarks" / "assets" / "catalog.json"),
        help="Path to immutable asset catalog JSON",
    )
    parser.add_argument(
        "--policy",
        default=str(REPO_ROOT / "benchmarks" / "assets" / "storage_policy.json"),
        help="Path to object-store / CDN storage policy JSON",
    )
    parser.add_argument(
        "--materialize",
        action="store_true",
        help="Copy file: sources into content-addressed objects/ under the catalog root",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="With --materialize, report actions without writing",
    )
    parser.add_argument(
        "--kind",
        action="append",
        default=[],
        help="Limit to asset kind(s); repeatable",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print machine-readable summary",
    )
    args = parser.parse_args()

    try:
        catalog, catalog_path = load_asset_catalog(args.catalog)
        policy, _policy_path = load_storage_policy(args.policy)
    except (OSError, json.JSONDecodeError, AssetCatalogError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    catalog_root = catalog_path.parent
    assets = catalog["assets"]
    if args.kind:
        allowed = set(args.kind)
        assets = [asset for asset in assets if asset["kind"] in allowed]
        catalog = {**catalog, "assets": assets}

    materialize_rows = []
    if args.materialize:
        for asset in assets:
            try:
                row = materialize_source_to_object_cache(
                    asset,
                    catalog_root,
                    policy,
                    dry_run=args.dry_run,
                )
            except AssetCatalogError as exc:
                print(f"ERROR: {exc}", file=sys.stderr)
                return 2
            materialize_rows.append(row)

    summary = summarize_catalog(
        catalog,
        policy=policy,
        catalog_root=catalog_root,
        verify_local=True,
    )
    summary["materialize"] = materialize_rows
    summary["catalog_path"] = str(catalog_path)
    summary["catalog_root"] = str(catalog_root)
    summary["fold_workspace_protection"] = "benchmarks/**/runs/** writes denied"

    failures = [
        row
        for row in summary.get("local_verification", [])
        if row["status"]
        in {
            "sha256_mismatch",
            "byte_size_mismatch",
            "missing_local",
        }
        and row["availability"] == "available"
    ]
    materialize_failures = [
        row
        for row in materialize_rows
        if row["status"]
        in {
            "source_missing",
            "source_sha256_mismatch",
            "source_byte_size_mismatch",
            "target_exists_with_mismatch",
            "write_verify_failed",
        }
    ]

    if args.json:
        print(json.dumps(summary, indent=2))
    else:
        print(f"catalog: {catalog.get('catalog_id')} ({catalog_path})")
        print(f"assets: {summary['n_assets']}")
        print("availability:", json.dumps(summary["by_availability"], sort_keys=True))
        print("by_kind:", json.dumps(summary["by_kind"], sort_keys=True))
        print(
            f"local verified: {summary.get('local_verified_count', 0)}; "
            f"local hard failures: {summary.get('local_failure_count', 0)}"
        )
        if materialize_rows:
            counts: dict[str, int] = {}
            for row in materialize_rows:
                counts[row["status"]] = counts.get(row["status"], 0) + 1
            print("materialize:", json.dumps(counts, sort_keys=True))
        if failures:
            print("available-asset verification failures:")
            for row in failures:
                print(f"  - {row['asset_id']}: {row['status']}")
        if materialize_failures:
            print("materialize failures:")
            for row in materialize_failures:
                print(f"  - {row['asset_id']}: {row['status']}")

    if failures or materialize_failures:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
