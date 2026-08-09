"""Download or verify immutable public benchmark assets from the locked registry."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import urllib.request
from pathlib import Path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify(path: Path, cohort: dict[str, object]) -> None:
    expected_bytes = int(cohort["expected_bytes"])
    if path.stat().st_size != expected_bytes:
        raise ValueError(
            f"{cohort['cohort_id']}: bytes={path.stat().st_size}, expected={expected_bytes}"
        )
    observed_hash = file_sha256(path)
    if observed_hash != str(cohort["expected_sha256"]):
        raise ValueError(
            f"{cohort['cohort_id']}: SHA-256={observed_hash}, "
            f"expected={cohort['expected_sha256']}"
        )


def download(target: Path, cohort: dict[str, object]) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f"{target.name}.partial.{os.getpid()}")
    request = urllib.request.Request(
        str(cohort["dataset_url"]),
        headers={"User-Agent": "CellTypePilot-public-benchmark/1.0"},
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response, temporary.open(
            "xb"
        ) as output:
            digest = hashlib.sha256()
            n_bytes = 0
            while chunk := response.read(8 * 1024 * 1024):
                output.write(chunk)
                digest.update(chunk)
                n_bytes += len(chunk)
        if n_bytes != int(cohort["expected_bytes"]):
            raise ValueError(
                f"{cohort['cohort_id']}: downloaded {n_bytes} bytes, "
                f"expected {cohort['expected_bytes']}"
            )
        if digest.hexdigest() != str(cohort["expected_sha256"]):
            raise ValueError(f"{cohort['cohort_id']}: downloaded SHA-256 mismatch")
        if target.exists():
            raise FileExistsError(
                f"{target} appeared during download; verified partial retained at {temporary}"
            )
        temporary.replace(target)
    except Exception:
        if temporary.exists():
            print(f"Partial download retained for audit: {temporary}")
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", required=True)
    parser.add_argument("--cohort", action="append", default=[])
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()

    registry_path = Path(args.registry).resolve()
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    selected = set(args.cohort)
    results = []
    for cohort in registry["cohorts"]:
        cohort_id = str(cohort["cohort_id"])
        if selected and cohort_id not in selected:
            continue
        target = (registry_path.parent / cohort["local_path"]).resolve()
        if target.exists():
            verify(target, cohort)
            status = "verified_existing"
        elif args.verify_only:
            status = "missing"
        else:
            download(target, cohort)
            verify(target, cohort)
            status = "downloaded_and_verified"
        results.append({"cohort_id": cohort_id, "status": status, "path": str(target)})
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
