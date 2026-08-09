"""Rescore raw checkpoint labels after a documented label-map amendment."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--label-map", required=True)
    parser.add_argument("--amendment-id", required=True)
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    map_path = Path(args.label_map).resolve()
    label_map = pd.read_csv(map_path, dtype=str)
    mapping = {
        (str(row.method), str(row.raw_label)): str(row.canonical_label)
        for row in label_map.itertuples(index=False)
    }
    changed_rows = 0
    changed_files = 0
    for path in sorted((run_dir / "checkpoints").glob("*.csv")):
        frame = pd.read_csv(path, dtype={"cell_id": str})
        if "raw_predicted_label" not in frame or "method" not in frame:
            continue
        keys = list(
            zip(
                frame["method"].astype(str),
                frame["raw_predicted_label"].astype(str),
                strict=True,
            )
        )
        missing = sorted(set(keys) - set(mapping))
        if missing:
            raise ValueError(f"{path.name}: missing mappings {missing[:5]}")
        remapped = [mapping[key] for key in keys]
        differences = frame["predicted_label"].astype(str).ne(remapped)
        if differences.any():
            changed_rows += int(differences.sum())
            changed_files += 1
            frame["predicted_label"] = remapped
            temporary = path.with_name(path.name + ".tmp")
            frame.to_csv(temporary, index=False, lineterminator="\n")
            os.replace(temporary, path)

    audit = {
        "schema_version": "celltypepilot.checkpoint-remap.v1",
        "amendment_id": args.amendment_id,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "label_map": str(map_path),
        "label_map_sha256": sha256(map_path),
        "changed_files": changed_files,
        "changed_rows": changed_rows,
        "source": "raw_predicted_label",
        "comparator_reexecution": False,
    }
    audit_path = run_dir / "checkpoint_remap_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
