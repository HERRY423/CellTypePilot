"""Read-only aggregation of GPU batch atomic checkpoints.

Does not re-run models. Writes only under batch_root/aggregate/.
Never writes into the CPU three-fold run tree.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))

DEFAULT_BATCH = REPO / "benchmarks" / "public_v1" / "batches" / "gpu_popv_retrain_v1"
CPU_RUN_FORBIDDEN = REPO / "benchmarks" / "public_v1" / "runs" / "travaglini_lung_smartseq2_2020"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-root", type=Path, default=DEFAULT_BATCH)
    parser.add_argument(
        "--require-complete",
        action="store_true",
        help="Exit non-zero unless every required fold has a completed popV checkpoint",
    )
    args = parser.parse_args()

    batch_root = args.batch_root.resolve()
    manifest = json.loads((batch_root / "batch_manifest.json").read_text(encoding="utf-8"))
    required = [str(x) for x in manifest["required_folds"]]
    method = "popv"
    checkpoint_dir = batch_root / "run" / "checkpoints"
    aggregate_dir = batch_root / "aggregate"
    if CPU_RUN_FORBIDDEN.resolve() in aggregate_dir.parents or aggregate_dir == CPU_RUN_FORBIDDEN.resolve():
        raise SystemExit("Refusing to aggregate into CPU run tree")

    aggregate_dir.mkdir(parents=True, exist_ok=True)

    if not checkpoint_dir.is_dir():
        raise SystemExit(f"No checkpoints directory: {checkpoint_dir}")

    from celltypepilot.benchmark_runner import _safe_fold_name

    records = []
    frames = []
    missing = []
    for fold_id in required:
        stem = f"{_safe_fold_name(fold_id)}__{method}"
        status_path = checkpoint_dir / f"{stem}.status.json"
        pred_path = checkpoint_dir / f"{stem}.csv"
        if not status_path.is_file():
            missing.append({"fold_id": fold_id, "reason": "status_missing"})
            continue
        status = json.loads(status_path.read_text(encoding="utf-8"))
        row = {
            "method": method,
            "fold_id": fold_id,
            "status": status.get("status"),
            "status_path": str(status_path),
            "status_sha256": sha256_file(status_path),
            "prediction_path": str(pred_path) if pred_path.is_file() else None,
            "prediction_sha256": sha256_file(pred_path) if pred_path.is_file() else None,
            "worker_id": status.get("worker_id")
            or (status.get("provenance") or {}).get("worker_id"),
            "batch_id": status.get("batch_id")
            or (status.get("provenance") or {}).get("batch_id"),
            "error": status.get("error"),
            "started_at_utc": status.get("started_at_utc"),
            "completed_at_utc": status.get("completed_at_utc"),
        }
        records.append(row)
        if status.get("status") == "completed" and pred_path.is_file():
            frame = pd.read_csv(pred_path, dtype={"cell_id": str})
            if "method" not in frame.columns:
                frame["method"] = method
            if "fold_id" not in frame.columns:
                frame["fold_id"] = fold_id
            frames.append(frame)
        elif status.get("status") != "completed":
            missing.append(
                {
                    "fold_id": fold_id,
                    "reason": f"status={status.get('status')}",
                    "error": status.get("error"),
                }
            )
        else:
            missing.append({"fold_id": fold_id, "reason": "prediction_csv_missing"})

    status_frame = pd.DataFrame(records)
    status_frame.to_csv(aggregate_dir / "comparator_status.csv", index=False, lineterminator="\n")

    if frames:
        predictions = pd.concat(frames, ignore_index=True)
        predictions.to_csv(
            aggregate_dir / "out_of_fold_predictions.csv",
            index=False,
            lineterminator="\n",
        )
    else:
        predictions = pd.DataFrame()

    report = {
        "schema_version": "celltypepilot.gpu-batch-aggregate.v1",
        "batch_id": manifest["batch_id"],
        "aggregated_at_utc": datetime.now(timezone.utc).isoformat(),
        "checkpoint_dir": str(checkpoint_dir),
        "aggregate_dir": str(aggregate_dir),
        "required_folds": required,
        "n_completed_prediction_frames": len(frames),
        "n_prediction_rows": int(len(predictions)),
        "records": records,
        "missing_or_incomplete": missing,
        "complete": len(missing) == 0 and len(frames) == len(required),
        "read_only_inputs": True,
        "cpu_run_isolation": str(CPU_RUN_FORBIDDEN),
        "device_track": "gpu",
        "claim": (
            "GPU batch aggregate is a separate execution track from the CPU three-fold run. "
            "Do not pool metrics across tracks without an explicit multi-track analysis plan."
        ),
    }
    (aggregate_dir / "aggregate_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: report[k] for k in (
        "batch_id",
        "complete",
        "n_completed_prediction_frames",
        "n_prediction_rows",
        "missing_or_incomplete",
    )}, indent=2))

    if args.require_complete and not report["complete"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
