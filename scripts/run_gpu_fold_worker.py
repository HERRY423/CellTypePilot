"""Run one or more donor folds for the frozen GPU popV batch.

Workers write only atomic checkpoints under the GPU batch run tree.
They never write into the CPU three-fold run directory.
"""

from __future__ import annotations

import argparse
import json
import os
import socket
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))

CPU_RUN_FORBIDDEN = REPO / "benchmarks" / "public_v1" / "runs" / "travaglini_lung_smartseq2_2020"
DEFAULT_BATCH = REPO / "benchmarks" / "public_v1" / "batches" / "gpu_popv_retrain_v1"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-root", type=Path, default=DEFAULT_BATCH)
    parser.add_argument(
        "--fold-id",
        action="append",
        default=[],
        help="Fold id to execute (repeatable). Use --all-folds for the three required donors.",
    )
    parser.add_argument("--all-folds", action="store_true")
    parser.add_argument("--worker-id", default=None, help="Defaults to hostname")
    parser.add_argument(
        "--write-aggregate",
        action="store_true",
        help="Also rewrite global OOF tables (single-node only; off by default for multi-node)",
    )
    parser.add_argument(
        "--continue-on-unavailable",
        action="store_true",
        help="Record failures without aborting remaining assigned folds",
    )
    args = parser.parse_args()

    batch_root = args.batch_root.resolve()
    manifest_path = batch_root / "batch_manifest.json"
    if not manifest_path.is_file():
        raise SystemExit(f"Missing batch manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    batch_id = str(manifest["batch_id"])
    cohort = manifest["cohort"]
    required = [str(x) for x in manifest["required_folds"]]

    if args.all_folds:
        fold_ids = required
    elif args.fold_id:
        fold_ids = [str(x) for x in args.fold_id]
        unknown = sorted(set(fold_ids) - set(required))
        if unknown:
            raise SystemExit(f"fold_id not in required_folds for this batch: {unknown}")
    else:
        raise SystemExit("Pass --fold-id ... and/or --all-folds")

    output_root = (batch_root / "run").resolve()
    # Hard isolation from CPU tree
    if output_root == CPU_RUN_FORBIDDEN.resolve() or CPU_RUN_FORBIDDEN.resolve() in output_root.parents:
        raise SystemExit("Refusing to use CPU three-fold run tree as GPU batch output")
    for forbidden in manifest.get("isolation", {}).get("must_not_write_into", []):
        forbidden_path = (REPO / forbidden).resolve() if not Path(forbidden).is_absolute() else Path(forbidden)
        if output_root == forbidden_path or forbidden_path in output_root.parents:
            raise SystemExit(f"Output root collides with forbidden path {forbidden}")

    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "checkpoints").mkdir(parents=True, exist_ok=True)

    # Resolve cohort assets relative to repo root.
    def resolve_asset(rel: str) -> Path:
        path = Path(rel)
        if not path.is_absolute():
            path = REPO / path
        return path.resolve()

    input_h5ad = resolve_asset(cohort["input_h5ad"])
    cluster_map = resolve_asset(cohort["cluster_map"])
    label_map = resolve_asset(cohort["label_map"])
    comparator = resolve_asset(manifest["comparator_config"])
    # Prefer copying holdout assignments into batch run without rewriting source.
    assignments_src = resolve_asset(cohort["holdout_assignments"])

    worker_id = args.worker_id or socket.gethostname()

    # Import heavy deps only after path guards.
    import pandas as pd

    from celltypepilot.benchmark import (
        apply_truth_label_map,
        build_holdout_assignments,
        save_benchmark_plan,
    )
    from celltypepilot.benchmark_runner import CommandComparator, run_benchmark_comparators
    from celltypepilot.data_adapter import load_h5ad

    adata = load_h5ad(str(input_h5ad))
    study_key = "ctp_study_id"
    if study_key not in adata.obs:
        adata.obs[study_key] = str(cohort["constant_study_id"])
    cluster_key = str(cohort["cluster_key"])
    if cluster_key not in adata.obs:
        cluster_frame = pd.read_csv(cluster_map, dtype=str)
        cluster_values = cluster_frame.set_index("cell_id")["cluster"]
        cluster_values.index = cluster_values.index.astype(str)
        expected = pd.Index(adata.obs_names.astype(str))
        adata.obs[cluster_key] = cluster_values.reindex(expected).to_numpy()

    # Rebuild assignments from live obs so fold ids match the locked plan; then
    # optionally verify against the frozen assignments file.
    assignments = build_holdout_assignments(
        adata.obs, study_key, str(cohort["donor_key"]), str(cohort["strategy"])
    )
    if assignments_src.is_file():
        locked = pd.read_csv(assignments_src, dtype=str)
        live_folds = set(assignments["fold_id"].astype(str))
        locked_folds = set(locked["fold_id"].astype(str))
        if not set(required).issubset(live_folds) or not set(required).issubset(locked_folds):
            raise SystemExit(
                "Required folds missing from live or locked assignments; "
                f"live={sorted(live_folds)} locked={sorted(locked_folds)}"
            )

    save_benchmark_plan(
        assignments,
        output_root,
        study_key,
        str(cohort["donor_key"]),
        str(cohort["strategy"]),
    )
    # Copy-forward isolation note without touching CPU artifacts.
    note = {
        "batch_id": batch_id,
        "worker_id": worker_id,
        "fold_ids": fold_ids,
        "write_aggregate_tables": bool(args.write_aggregate),
        "cpu_run_isolation": str(CPU_RUN_FORBIDDEN),
        "device_track": "gpu",
    }
    (output_root / "worker_invocation.json").write_text(
        json.dumps(note, indent=2) + "\n", encoding="utf-8"
    )

    map_frame = pd.read_csv(label_map, dtype=str)
    specs = (CommandComparator.from_json(comparator),)
    os.environ.setdefault("CTP_POPV_GPU_IMAGE", "celltypepilot-popv-gpu:0.6.1-cu124")

    predictions, status = run_benchmark_comparators(
        adata,
        assignments,
        str(cohort["truth_key"]),
        cluster_key,
        output_root,
        str(cohort["species"]),
        str(cohort["tissue"]),
        methods=("popv",),
        command_specs=specs,
        label_map=map_frame,
        continue_on_unavailable=args.continue_on_unavailable,
        fold_ids=tuple(fold_ids),
        write_aggregate_tables=bool(args.write_aggregate),
        worker_id=worker_id,
        batch_id=batch_id,
    )

    print(
        json.dumps(
            {
                "batch_id": batch_id,
                "worker_id": worker_id,
                "fold_ids": fold_ids,
                "n_prediction_rows": int(len(predictions)),
                "status_rows": int(len(status)),
                "checkpoint_dir": str(output_root / "checkpoints"),
                "write_aggregate_tables": bool(args.write_aggregate),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
