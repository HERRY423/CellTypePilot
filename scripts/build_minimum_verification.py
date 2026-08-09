"""Build a single-cohort technical verification bundle from locked OOF predictions."""

from __future__ import annotations

import argparse
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from celltypepilot.benchmark import (
    BenchmarkValidationError,
    apply_truth_label_map,
    evaluate_holdout_predictions,
    validate_out_of_fold_predictions,
)
from celltypepilot.benchmark_release import (
    _markdown_table,
    _normalise_metadata,
    _read_obs,
    _sha256,
)
from celltypepilot.benchmark_runner import apply_locked_label_map
from celltypepilot.robustness import (
    batch_sensitivity,
    evaluate_by_independent_unit,
    merge_prediction_metadata,
    paired_method_comparisons,
    qc_stratified_performance,
    sample_enrichment_diagnostics,
    summarize_independent_units,
)

METHODS = ("celltypepilot", "celltypist", "singler", "popv")


def atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    os.close(descriptor)
    try:
        frame.to_csv(temporary_name, index=False)
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def resolve(base: Path, value: str) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def completed_fold_check(
    status: pd.DataFrame, assignments: pd.DataFrame, methods: tuple[str, ...]
) -> None:
    expected = set(assignments["fold_id"].astype(str))
    for method in methods:
        rows = status[
            (status["method"].astype(str) == method)
            & (status["status"].astype(str) == "completed")
        ]
        completed = set(rows["fold_id"].astype(str))
        if completed != expected:
            raise BenchmarkValidationError(
                f"{method} does not have exact completed-fold provenance: "
                f"missing={sorted(expected - completed)}, extra={sorted(completed - expected)}"
            )


def checkpoint_runtimes(checkpoint_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for path in sorted(checkpoint_dir.glob("*.status.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("status") != "completed":
            continue
        started = datetime.fromisoformat(str(payload["started_at_utc"]))
        completed = datetime.fromisoformat(str(payload["completed_at_utc"]))
        provenance = payload.get("provenance", {})
        rows.append(
            {
                "method": payload.get("method"),
                "fold_id": payload.get("fold_id"),
                "started_at_utc": started.isoformat(),
                "completed_at_utc": completed.isoformat(),
                "wall_seconds": (completed - started).total_seconds(),
                "version": provenance.get("version"),
                "reference_policy": provenance.get("reference_policy"),
                "confidence_semantics": provenance.get("confidence_semantics"),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", required=True)
    parser.add_argument("--cohort-id", required=True)
    parser.add_argument("--azimuth-audit", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-boot", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    registry_path = Path(args.registry).resolve()
    base = registry_path.parent
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    matches = [row for row in registry["cohorts"] if row["cohort_id"] == args.cohort_id]
    if len(matches) != 1:
        raise BenchmarkValidationError(f"Expected one cohort named {args.cohort_id!r}")
    cohort = matches[0]
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)

    data_path = resolve(base, cohort["local_path"])
    label_map_path = resolve(base, cohort["label_map_path"])
    assignment_path = resolve(base, cohort["assignments_path"])
    prediction_path = resolve(base, cohort["predictions_path"])
    status_path = resolve(base, cohort["comparator_status_path"])
    azimuth_path = Path(args.azimuth_audit).resolve()
    if data_path.stat().st_size != int(cohort["expected_bytes"]):
        raise BenchmarkValidationError("Dataset byte size differs from the locked registry")
    if _sha256(data_path) != cohort["expected_sha256"]:
        raise BenchmarkValidationError("Dataset SHA-256 differs from the locked registry")
    if _sha256(label_map_path) != cohort["label_map_sha256"]:
        raise BenchmarkValidationError("Label-map SHA-256 differs from the locked registry")
    if not azimuth_path.exists():
        raise BenchmarkValidationError("Azimuth reference audit is required")

    assignments = pd.read_csv(assignment_path, dtype={"cell_id": str, "fold_id": str})
    predictions = pd.read_csv(prediction_path, dtype={"cell_id": str, "fold_id": str})
    predictions = predictions[predictions["method"].astype(str).isin(METHODS)].copy()
    status = pd.read_csv(status_path, dtype={"fold_id": str})
    completed_fold_check(status, assignments, METHODS)
    validate_out_of_fold_predictions(assignments, predictions)

    obs = _read_obs(data_path)
    metadata, truth = _normalise_metadata(obs, cohort)
    label_map = pd.read_csv(label_map_path, dtype=str)
    truth = apply_truth_label_map(truth, label_map)
    predictions = apply_locked_label_map(predictions, label_map)
    aggregate, _ = evaluate_holdout_predictions(
        truth,
        assignments,
        predictions,
        expected_methods=METHODS,
        bootstrap_ci=False,
    )
    merged = merge_prediction_metadata(truth, metadata, predictions)
    merged["__donor__"] = (
        str(cohort.get("donor_namespace", args.cohort_id))
        + "::"
        + merged["__donor__"].astype(str)
    )

    donor = evaluate_by_independent_unit(
        merged, study_key="__study__", donor_key="__donor__"
    )
    summary = summarize_independent_units(
        donor, n_boot=args.n_boot, seed=args.seed
    )
    comparisons = paired_method_comparisons(
        donor, reference_method="celltypepilot", seed=args.seed
    )
    batch_levels, batch_summary = batch_sensitivity(
        merged,
        study_key="__study__",
        donor_key="__donor__",
        axes={
            "platform": "__platform__" if "__platform__" in merged else None,
            "plate_batch": "__batch__" if "__batch__" in merged else None,
            "condition": "__condition__" if "__condition__" in merged else None,
        },
        metric="macro_f1",
    )
    qc = qc_stratified_performance(
        merged,
        study_key="__study__",
        donor_key="__donor__",
        diagnostics=cohort.get("diagnostics", {}),
    )
    enrichment_frames = []
    for method, frame in merged.groupby("method", sort=True):
        enrichment = sample_enrichment_diagnostics(
            frame,
            cluster_key="predicted_label",
            sample_key="__sample__" if "__sample__" in frame else None,
        )
        enrichment.insert(0, "method", method)
        enrichment_frames.append(enrichment)
    enrichment = pd.concat(enrichment_frames, ignore_index=True)
    runtimes = checkpoint_runtimes(status_path.parent / "checkpoints")
    azimuth = json.loads(azimuth_path.read_text(encoding="utf-8"))

    negative_rows: list[dict[str, Any]] = [
        {
            "category": "scope",
            "finding": "single_study_single_platform_technical_verification",
            "severity": "claim_boundary",
            "detail": "No multi-study or cross-platform robustness claim is supported.",
        },
        {
            "category": "azimuth",
            "finding": azimuth.get("primary_track_status", "audit_status_unknown"),
            "severity": "retained_negative_result",
            "detail": azimuth.get("primary_track_reason", ""),
        },
    ]
    celltypepilot_aggregate = aggregate[
        aggregate["method"].astype(str) == "celltypepilot"
    ]
    if (
        not celltypepilot_aggregate.empty
        and float(celltypepilot_aggregate.iloc[0]["coverage"]) == 0.0
    ):
        negative_rows.append(
            {
                "category": "celltypepilot",
                "finding": "complete_abstention_on_all_out_of_fold_cells",
                "severity": "retained_negative_result",
                "detail": "Coverage was 0/9409; this result is retained and was not replaced by candidate labels.",
            }
        )
    for row in batch_summary[batch_summary["status"] != "estimated_descriptive"].to_dict(
        orient="records"
    ):
        negative_rows.append(
            {
                "category": f"batch_sensitivity:{row.get('axis')}",
                "finding": row.get("status"),
                "severity": "diagnostic_unassessed",
                "detail": row.get("method", ""),
            }
        )
    for row in qc[qc["status"].astype(str).str.startswith("not_assessed")].to_dict(
        orient="records"
    ):
        negative_rows.append(
            {
                "category": row.get("diagnostic"),
                "finding": row.get("status"),
                "severity": "diagnostic_unassessed",
                "detail": "Missing diagnostics are not evidence of absence.",
            }
        )
    for row in enrichment[enrichment["flag"].astype(str) == "SAMPLE_ENRICHED"].to_dict(
        orient="records"
    ):
        negative_rows.append(
            {
                "category": "sample_enrichment",
                "finding": "predicted_label_concentrated_in_one_sample",
                "severity": "diagnostic_flag",
                "detail": (
                    f"{row.get('method')}:{row.get('cluster')} -> "
                    f"{row.get('dominant_sample')} "
                    f"({float(row.get('dominant_sample_fraction')):.3f})"
                ),
            }
        )
    negative_rows.append(
        {
            "category": "method_comparison",
            "finding": "minimal_inference_only_three_paired_donors",
            "severity": "underpowered",
            "detail": "Exact sign-flip tests and donor bootstrap intervals are reported, but n=3 donors cannot support a superiority claim.",
        }
    )
    for row in comparisons[comparisons["status"] != "estimated"].to_dict(
        orient="records"
    ):
        negative_rows.append(
            {
                "category": "method_comparison",
                "finding": row.get("status"),
                "severity": "underpowered",
                "detail": f"{row.get('method_a')} vs {row.get('method_b')} {row.get('metric')}",
            }
        )
    negative = pd.DataFrame(negative_rows)

    artifacts = {
        "aggregate_metrics.csv": aggregate,
        "donor_metrics.csv": donor,
        "donor_weighted_summary.csv": summary,
        "paired_method_comparisons.csv": comparisons,
        "batch_sensitivity_levels.csv": batch_levels,
        "batch_sensitivity_summary.csv": batch_summary,
        "qc_stratified_performance.csv": qc,
        "sample_enrichment.csv": enrichment,
        "runtime_by_fold.csv": runtimes,
        "negative_results.csv": negative,
    }
    for name, frame in artifacts.items():
        write_csv(frame, output / name)

    report = "\n".join(
        [
            f"# Minimum technical verification: {args.cohort_id}",
            "",
            "Status: **minimum_verification_complete**",
            "",
            "This is a real single-study, single-platform, three-donor technical verification. "
            "It is not the public multi-cohort benchmark release.",
            "",
            "## Donor-weighted results",
            "",
            _markdown_table(summary),
            "",
            "## Paired donor-level comparisons",
            "",
            _markdown_table(comparisons),
            "",
            "## Batch sensitivity",
            "",
            _markdown_table(batch_summary),
            "",
            "## QC-stratified performance",
            "",
            _markdown_table(qc),
            "",
            "The locked low-quality rule was `nGene <= 500`. The dataset minimum was 501, "
            "so no retained cell entered the flagged stratum; this is evidence of upstream "
            "filtering in this asset, not evidence that low-quality robustness was tested.",
            "",
            "## Sample-enrichment flags",
            "",
            _markdown_table(enrichment[enrichment["flag"].astype(str) != "PASS"]),
            "",
            "## QC and retained negative findings",
            "",
            _markdown_table(negative),
            "",
            "## Claim boundary",
            "",
            "Cell-level metrics are descriptive. Donors are the independent units. Three donors "
            "provide only minimal inferential support; one study and one platform cannot establish "
            "study, platform, tumor, inflamed-tissue, or general batch robustness. Missing doublet "
            "or ambient-RNA metadata is unassessed, not negative. Confidence values from each tool "
            "retain their method-specific, non-probabilistic semantics.",
            "",
        ]
    )
    report_path = output / "verification_report.md"
    atomic_write_text(report_path, report)

    copied_azimuth = output / "azimuth_reference_audit.json"
    atomic_write_text(copied_azimuth, json.dumps(azimuth, indent=2) + "\n")
    manifest_artifacts = [*(output / name for name in artifacts), report_path, copied_azimuth]
    manifest = {
        "schema_version": "celltypepilot.minimum-verification.v1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "cohort_id": args.cohort_id,
        "dataset_sha256": _sha256(data_path),
        "label_map_sha256": _sha256(label_map_path),
        "methods": list(METHODS),
        "n_donor_folds": int(assignments["fold_id"].nunique()),
        "statistical_unit": "donor",
        "random_seed": args.seed,
        "bootstrap_replicates": args.n_boot,
        "validation_scope": {
            "run_role": "single_cohort_minimum_technical_verification",
            "public_multi_cohort_release": False,
            "cell_weighted_inference_prohibited": True,
            "negative_results_retained": True,
        },
        "artifacts": {
            path.name: {"path": str(path), "sha256": _sha256(path)}
            for path in manifest_artifacts
        },
    }
    manifest_path = output / "verification_manifest.json"
    atomic_write_text(manifest_path, json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"status": "minimum_verification_complete", **manifest}, indent=2))


if __name__ == "__main__":
    main()
