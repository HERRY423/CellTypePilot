"""Build an auditable multi-cohort public benchmark release.

The release builder consumes immutable cohort metadata plus locked out-of-fold
predictions.  It does not download data, tune thresholds, or train comparators;
those actions remain separate so that outcome access cannot silently change the
analysis plan.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

from .benchmark import (
    BenchmarkValidationError,
    apply_truth_label_map,
    evaluate_holdout_predictions,
    validate_out_of_fold_predictions,
)
from .benchmark_runner import apply_locked_label_map
from .robustness import (
    batch_sensitivity,
    evaluate_by_independent_unit,
    merge_prediction_metadata,
    paired_method_comparisons,
    qc_stratified_performance,
    sample_enrichment_diagnostics,
    summarize_independent_units,
)

REGISTRY_SCHEMA = "celltypepilot.public-benchmark-registry.v1"
RELEASE_SCHEMA = "celltypepilot.public-benchmark-release.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve(base: Path, value: str | None) -> Path | None:
    if not value:
        return None
    path = Path(value)
    return path if path.is_absolute() else (base / path).resolve()


def load_benchmark_registry(path: str | Path) -> tuple[dict[str, Any], Path]:
    registry_path = Path(path).resolve()
    payload = json.loads(registry_path.read_text(encoding="utf-8"))
    validate_benchmark_registry(payload)
    return payload, registry_path


def validate_benchmark_registry(payload: dict[str, Any]) -> None:
    """Validate scientific-design and provenance fields before reading outcomes."""
    if payload.get("schema_version") != REGISTRY_SCHEMA:
        raise BenchmarkValidationError(f"Registry schema must be {REGISTRY_SCHEMA!r}")
    required_top = {"release_id", "analysis_plan", "required_methods", "cohorts"}
    missing = required_top - set(payload)
    if missing:
        raise BenchmarkValidationError(f"Registry missing fields: {sorted(missing)}")
    plan = payload["analysis_plan"]
    required_plan = {
        "primary_metric",
        "independent_unit",
        "holdout_policy",
        "multiplicity_correction",
        "negative_result_policy",
    }
    missing_plan = required_plan - set(plan)
    if missing_plan:
        raise BenchmarkValidationError(
            f"Analysis plan missing fields: {sorted(missing_plan)}"
        )
    if plan["independent_unit"] != "donor":
        raise BenchmarkValidationError("Public release independent_unit must be 'donor'")
    if len(payload["cohorts"]) < 2:
        raise BenchmarkValidationError("A multi-cohort release requires at least two cohorts")
    if len(set(payload["required_methods"])) != len(payload["required_methods"]):
        raise BenchmarkValidationError("required_methods must be unique")

    ids = []
    for cohort in payload["cohorts"]:
        required = {
            "cohort_id",
            "title",
            "species",
            "tissue",
            "collection_url",
            "dataset_url",
            "dataset_version_id",
            "expected_cells",
            "expected_bytes",
            "expected_sha256",
            "citation_doi",
            "truth_provenance",
            "label_map_path",
            "label_map_sha256",
            "metadata",
        }
        missing_cohort = required - set(cohort)
        if missing_cohort:
            raise BenchmarkValidationError(
                f"Cohort {cohort.get('cohort_id', '<unknown>')} missing: {sorted(missing_cohort)}"
            )
        ids.append(str(cohort["cohort_id"]))
        expected_sha256 = str(cohort["expected_sha256"])
        for field, digest in (
            ("expected_sha256", expected_sha256),
            ("label_map_sha256", str(cohort["label_map_sha256"])),
        ):
            if len(digest) != 64 or any(
                character not in "0123456789abcdef" for character in digest
            ):
                raise BenchmarkValidationError(
                    f"Cohort {cohort['cohort_id']} {field} must be lowercase SHA-256"
                )
        metadata = cohort["metadata"]
        for key in ("truth_key", "donor_key", "platform_key"):
            if not metadata.get(key):
                raise BenchmarkValidationError(f"Cohort {cohort['cohort_id']} needs {key}")
        if not metadata.get("study_key") and not cohort.get("constant_study_id"):
            raise BenchmarkValidationError(
                f"Cohort {cohort['cohort_id']} needs study_key or constant_study_id"
            )
        if cohort.get("dataset_role", "benchmark") == "calibration":
            raise BenchmarkValidationError("Calibration cohorts cannot enter the test registry")
    if len(ids) != len(set(ids)):
        raise BenchmarkValidationError("cohort_id values must be unique")


def public_cohort_inventory(payload: dict[str, Any], registry_path: Path) -> pd.DataFrame:
    rows = []
    base = registry_path.parent
    for cohort in payload["cohorts"]:
        local_path = _resolve(base, cohort.get("local_path"))
        predictions = _resolve(base, cohort.get("predictions_path"))
        assignments = _resolve(base, cohort.get("assignments_path"))
        cluster_map = _resolve(base, cohort.get("cluster_map_path"))
        comparator_status = _resolve(base, cohort.get("comparator_status_path"))
        label_map = _resolve(base, cohort.get("label_map_path"))
        row = {
                "cohort_id": cohort["cohort_id"],
                "title": cohort["title"],
                "species": cohort["species"],
                "tissue": cohort["tissue"],
                "dataset_version_id": cohort["dataset_version_id"],
                "expected_cells": cohort.get("expected_cells"),
                "expected_bytes": cohort.get("expected_bytes"),
                "expected_sha256": cohort.get("expected_sha256"),
                "citation_doi": cohort["citation_doi"],
                "collection_url": cohort["collection_url"],
                "dataset_url": cohort["dataset_url"],
                "truth_provenance": cohort["truth_provenance"],
                "local_data_status": "present" if local_path and local_path.exists() else "missing",
                "predictions_status": (
                    "present" if predictions and predictions.exists() else "missing"
                ),
                "assignments_status": (
                    "present" if assignments and assignments.exists() else "missing"
                ),
                "cluster_map_status": (
                    "present" if cluster_map and cluster_map.exists() else "missing"
                ),
                "comparator_status": (
                    "present" if comparator_status and comparator_status.exists() else "missing"
                ),
                "label_map_status": (
                    "verified"
                    if label_map
                    and label_map.exists()
                    and _sha256(label_map) == cohort.get("label_map_sha256")
                    else "mismatch"
                    if label_map and label_map.exists()
                    else "missing"
                ),
                "expected_label_map_sha256": cohort.get("label_map_sha256"),
                "actual_label_map_sha256": (
                    _sha256(label_map) if label_map and label_map.exists() else None
                ),
                "assignments_sha256": (
                    _sha256(assignments) if assignments and assignments.exists() else None
                ),
                "cluster_map_sha256": (
                    _sha256(cluster_map) if cluster_map and cluster_map.exists() else None
                ),
                "predictions_sha256": (
                    _sha256(predictions) if predictions and predictions.exists() else None
                ),
                "comparator_status_sha256": (
                    _sha256(comparator_status)
                    if comparator_status and comparator_status.exists()
                    else None
                ),
        }
        if local_path and local_path.exists():
            obs = _read_obs(local_path)
            actual_sha256 = _sha256(local_path)
            metadata = cohort["metadata"]
            truth_key = metadata["truth_key"]
            donor_key = metadata["donor_key"]
            platform_key = metadata["platform_key"]
            condition_key = metadata.get("condition_key")
            required = [truth_key, donor_key, platform_key]
            missing_metadata = [key for key in required if key not in obs]
            row.update(
                {
                    "actual_sha256": actual_sha256,
                    "actual_cells": len(obs),
                    "actual_donors": obs[donor_key].nunique() if donor_key in obs else None,
                    "actual_platforms": obs[platform_key].nunique()
                    if platform_key in obs
                    else None,
                    "actual_conditions": obs[condition_key].nunique()
                    if condition_key and condition_key in obs
                    else None,
                    "actual_truth_labels": obs[truth_key].nunique() if truth_key in obs else None,
                    "metadata_status": (
                        "verified"
                        if not missing_metadata
                        and (not cohort.get("expected_cells") or len(obs) == cohort["expected_cells"])
                        and actual_sha256 == cohort["expected_sha256"]
                        else "mismatch"
                    ),
                    "missing_required_metadata": ",".join(missing_metadata),
                }
            )
        rows.append(row)
    return pd.DataFrame(rows)


def _read_obs(path: Path) -> pd.DataFrame:
    dataset = ad.read_h5ad(path, backed="r")
    try:
        return dataset.obs.copy()
    finally:
        if dataset.file is not None:
            dataset.file.close()


def _normalise_metadata(obs: pd.DataFrame, cohort: dict[str, Any]) -> tuple[pd.DataFrame, pd.Series]:
    spec = cohort["metadata"]
    truth_key = spec["truth_key"]
    donor_key = spec["donor_key"]
    if truth_key not in obs or donor_key not in obs:
        raise BenchmarkValidationError(
            f"Cohort {cohort['cohort_id']} lacks truth_key or donor_key in obs"
        )
    normal = obs.copy()
    normal.index = normal.index.astype(str)
    normal["__cohort__"] = str(cohort["cohort_id"])
    study_key = spec.get("study_key")
    if study_key and study_key in normal:
        normal["__study__"] = normal[study_key].astype(str)
    elif cohort.get("constant_study_id"):
        normal["__study__"] = str(cohort["constant_study_id"])
    else:
        raise BenchmarkValidationError(f"Cohort {cohort['cohort_id']} study metadata is absent")
    normal["__donor__"] = normal[donor_key].astype(str)
    for target, source_key in (
        ("__platform__", "platform_key"),
        ("__sample__", "sample_key"),
        ("__condition__", "condition_key"),
        ("__batch__", "batch_key"),
        ("__cluster__", "cluster_key"),
    ):
        source = spec.get(source_key)
        if source and source in normal:
            normal[target] = normal[source]
    truth = normal[truth_key].astype(str).copy()
    return normal, truth


def _prefix_cohort_cells(frame: pd.DataFrame, cohort_id: str) -> pd.DataFrame:
    output = frame.copy()
    output["cell_id"] = cohort_id + "::" + output["cell_id"].astype(str)
    return output


def _write_csv(frame: pd.DataFrame, path: Path, artifacts: dict[str, Path]) -> None:
    frame.to_csv(path, index=False)
    artifacts[path.stem] = path


def _software_versions() -> dict[str, str]:
    packages = ["celltypepilot", "anndata", "numpy", "pandas", "scipy", "scanpy"]
    versions = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not_installed_as_distribution"
    return versions


def _markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "No rows."
    display = frame.fillna("").astype(str)
    header = "| " + " | ".join(display.columns) + " |"
    rule = "| " + " | ".join("---" for _ in display.columns) + " |"
    rows = [
        "| "
        + " | ".join(value.replace("|", "\\|").replace("\n", " ") for value in row)
        + " |"
        for row in display.to_numpy()
    ]
    return "\n".join([header, rule, *rows])


def _render_report(
    payload: dict[str, Any],
    inventory: pd.DataFrame,
    summary: pd.DataFrame,
    comparisons: pd.DataFrame,
    negative: pd.DataFrame,
    readiness: dict[str, Any],
) -> str:
    lines = [
        f"# CellTypePilot public benchmark release: {payload['release_id']}",
        "",
        f"Release status: **{readiness['status']}**",
        "",
        "The independent statistical unit is the donor. Cell-level metrics are descriptive; "
        "confidence intervals and method comparisons use donor-level observations.",
        "",
        "## Locked analysis plan",
        "",
        f"- Primary metric: `{payload['analysis_plan']['primary_metric']}`",
        f"- Holdout policy: `{payload['analysis_plan']['holdout_policy']}`",
        f"- Multiplicity correction: `{payload['analysis_plan']['multiplicity_correction']}`",
        f"- Negative-result policy: {payload['analysis_plan']['negative_result_policy']}",
        "",
        "## Cohort materialization",
        "",
        _markdown_table(inventory),
        "",
        "## Donor-weighted summary",
        "",
        _markdown_table(summary) if not summary.empty else "No evaluable predictions.",
        "",
        "## Paired donor comparisons",
        "",
        _markdown_table(comparisons) if not comparisons.empty else "No estimable comparisons.",
        "",
        "## Negative and incomplete findings",
        "",
        _markdown_table(negative) if not negative.empty else "No negative rows recorded.",
        "",
        "## Claim boundary",
        "",
        "Results apply only to the immutable dataset versions, truth-label provenance, label map, "
        "comparators, and holdout plan recorded here. Expert-curated labels are reference labels, "
        "not infallible biological ground truth. Missing ambient-RNA or doublet metadata is an "
        "unassessed diagnostic, not evidence that those artifacts are absent.",
        "",
    ]
    return "\n".join(lines)


def build_public_benchmark_release(
    registry: str | Path,
    output_dir: str | Path,
    *,
    allow_incomplete: bool = False,
    n_boot: int = 2000,
    seed: int = 42,
) -> dict[str, Path]:
    """Evaluate all materialized cohorts and write a release with retained failures."""
    payload, registry_path = load_benchmark_registry(registry)
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    inventory = public_cohort_inventory(payload, registry_path)
    base = registry_path.parent
    artifacts: dict[str, Path] = {}
    negative_rows: list[dict[str, Any]] = []
    cohort_result_rows = []
    merged_frames = []
    enrichment_frames = []
    qc_frames = []

    for cohort in payload["cohorts"]:
        cohort_id = str(cohort["cohort_id"])
        data_path = _resolve(base, cohort.get("local_path"))
        prediction_path = _resolve(base, cohort.get("predictions_path"))
        assignment_path = _resolve(base, cohort.get("assignments_path"))
        cluster_map_path = _resolve(base, cohort.get("cluster_map_path"))
        status_path = _resolve(base, cohort.get("comparator_status_path"))
        label_map_path = _resolve(base, cohort.get("label_map_path"))
        prerequisites = (
            (prediction_path, "predictions", "out_of_fold_predictions_missing"),
            (assignment_path, "independence", "locked_holdout_assignments_missing"),
            (status_path, "comparator_provenance", "comparator_status_missing"),
            (label_map_path, "label_harmonization", "predeclared_label_map_missing"),
        )
        for prerequisite, category, finding in prerequisites:
            if not prerequisite or not prerequisite.exists():
                negative_rows.append(
                    {
                        "cohort_id": cohort_id,
                        "category": category,
                        "finding": finding,
                        "severity": "blocks_release",
                    }
                )
        if not data_path or not data_path.exists():
            negative_rows.append(
                {
                    "cohort_id": cohort_id,
                    "category": "data",
                    "finding": "public_dataset_not_materialized",
                    "severity": "blocks_release",
                }
            )
            continue
        if cohort.get("expected_bytes") and data_path.stat().st_size != int(cohort["expected_bytes"]):
            raise BenchmarkValidationError(
                f"Cohort {cohort_id} byte size differs from immutable registry"
            )
        if _sha256(data_path) != str(cohort["expected_sha256"]):
            raise BenchmarkValidationError(
                f"Cohort {cohort_id} SHA-256 differs from immutable registry"
            )
        obs = _read_obs(data_path)
        normal, truth = _normalise_metadata(obs, cohort)
        if "__cluster__" not in normal and (
            not cluster_map_path or not cluster_map_path.exists()
        ):
            negative_rows.append(
                {
                    "cohort_id": cohort_id,
                    "category": "clustering",
                    "finding": "truth_blind_cluster_map_missing",
                    "severity": "blocks_release",
                }
            )
        if "__cluster__" not in normal and cluster_map_path and cluster_map_path.exists():
            cluster_map = pd.read_csv(cluster_map_path, dtype=str)
            if set(cluster_map) < {"cell_id", "cluster"}:
                raise BenchmarkValidationError(
                    f"Cohort {cohort_id} cluster map needs cell_id and cluster"
                )
            if cluster_map["cell_id"].duplicated().any():
                raise BenchmarkValidationError(f"Cohort {cohort_id} cluster map has duplicates")
            mapped = cluster_map.set_index("cell_id")["cluster"]
            mapped.index = mapped.index.astype(str)
            missing_clusters = normal.index.difference(mapped.index)
            if len(missing_clusters):
                raise BenchmarkValidationError(
                    f"Cohort {cohort_id} cluster map misses {len(missing_clusters)} cells"
                )
            normal["__cluster__"] = mapped.reindex(normal.index).to_numpy()
        if status_path and status_path.exists():
            status = pd.read_csv(status_path)
            status_methods = set(status["method"].astype(str))
            for missing_method in set(payload["required_methods"]) - status_methods:
                negative_rows.append(
                    {
                        "cohort_id": cohort_id,
                        "category": "comparator_provenance",
                        "finding": f"{missing_method}:status_missing",
                        "severity": "blocks_release",
                    }
                )
            failed = status[status["status"] != "completed"]
            for row in failed.to_dict(orient="records"):
                negative_rows.append(
                    {
                        "cohort_id": cohort_id,
                        "category": "comparator",
                        "finding": f"{row.get('method', 'unknown')}:{row.get('status')}",
                        "severity": (
                            "blocks_release"
                            if row.get("method") in payload["required_methods"]
                            else "retained_negative_result"
                        ),
                        "detail": row.get("error", row.get("detail", "")),
                    }
                )
        if label_map_path and label_map_path.exists() and _sha256(label_map_path) != str(
            cohort["label_map_sha256"]
        ):
            raise BenchmarkValidationError(
                f"Cohort {cohort_id} label-map SHA-256 differs from frozen registry"
            )
        if not prediction_path or not prediction_path.exists():
            continue
        if not assignment_path or not assignment_path.exists():
            continue

        predictions = pd.read_csv(prediction_path, dtype={"cell_id": str, "fold_id": str})
        assignments = pd.read_csv(assignment_path, dtype={"cell_id": str, "fold_id": str})
        if status_path and status_path.exists():
            expected_folds = set(assignments["fold_id"].astype(str))
            status = pd.read_csv(status_path, dtype={"fold_id": str})
            for method in payload["required_methods"]:
                completed_folds = set(
                    status.loc[
                        (status["method"].astype(str) == method)
                        & (status["status"].astype(str) == "completed"),
                        "fold_id",
                    ].astype(str)
                )
                if completed_folds != expected_folds:
                    negative_rows.append(
                        {
                            "cohort_id": cohort_id,
                            "category": "comparator_provenance",
                            "finding": f"{method}:incomplete_fold_provenance",
                            "severity": "blocks_release",
                            "detail": (
                                f"missing={len(expected_folds - completed_folds)},"
                                f"extra={len(completed_folds - expected_folds)}"
                            ),
                        }
                    )
        label_map = pd.read_csv(label_map_path, dtype=str) if label_map_path else None
        truth = apply_truth_label_map(truth, label_map)
        predictions = apply_locked_label_map(predictions, label_map)
        validate_out_of_fold_predictions(assignments, predictions)
        aggregate, _ = evaluate_holdout_predictions(
            truth,
            assignments,
            predictions,
            expected_methods=tuple(payload["required_methods"]),
            bootstrap_ci=False,
        )
        aggregate.insert(0, "cohort_id", cohort_id)
        cohort_result_rows.append(aggregate)
        expected_cells = set(assignments["cell_id"].astype(str))
        for method in payload["required_methods"]:
            method_cells = set(
                predictions.loc[predictions["method"].astype(str) == method, "cell_id"].astype(str)
            )
            if method_cells != expected_cells:
                negative_rows.append(
                    {
                        "cohort_id": cohort_id,
                        "category": "comparator_coverage",
                        "finding": f"{method}:incomplete_or_missing_predictions",
                        "severity": "blocks_release",
                        "detail": (
                            f"missing={len(expected_cells - method_cells)},"
                            f"extra={len(method_cells - expected_cells)}"
                        ),
                    }
                )

        merged = merge_prediction_metadata(truth, normal, predictions)
        donor_namespace = str(cohort.get("donor_namespace", cohort_id))
        merged["cell_id"] = cohort_id + "::" + merged["cell_id"].astype(str)
        merged["__donor__"] = donor_namespace + "::" + merged["__donor__"].astype(str)
        merged_frames.append(merged)
        for method, method_frame in merged.groupby("method", sort=True):
            enrichment = sample_enrichment_diagnostics(
                method_frame,
                cluster_key="predicted_label",
                sample_key="__sample__" if "__sample__" in method_frame else None,
            )
            enrichment.insert(0, "method", str(method))
            enrichment.insert(0, "cohort_id", cohort_id)
            enrichment_frames.append(enrichment)

        qc = qc_stratified_performance(
            merged,
            study_key="__study__",
            donor_key="__donor__",
            diagnostics=cohort.get("diagnostics", {}),
        )
        qc.insert(0, "cohort_id", cohort_id)
        qc_frames.append(qc)

    cohort_results = (
        pd.concat(cohort_result_rows, ignore_index=True)
        if cohort_result_rows
        else pd.DataFrame()
    )
    merged_all = pd.concat(merged_frames, ignore_index=True) if merged_frames else pd.DataFrame()
    donor_metrics = (
        evaluate_by_independent_unit(
            merged_all,
            study_key="__study__",
            donor_key="__donor__",
        )
        if not merged_all.empty
        else pd.DataFrame()
    )
    summary = (
        summarize_independent_units(donor_metrics, n_boot=n_boot, seed=seed)
        if not donor_metrics.empty
        else pd.DataFrame()
    )
    comparisons = (
        paired_method_comparisons(
            donor_metrics,
            reference_method=str(payload["analysis_plan"].get("reference_method", "celltypepilot")),
            seed=seed,
        )
        if not donor_metrics.empty
        else pd.DataFrame()
    )
    if not merged_all.empty:
        batch_levels, batch_summary = batch_sensitivity(
            merged_all,
            study_key="__study__",
            donor_key="__donor__",
            axes={
                "platform": "__platform__" if "__platform__" in merged_all else None,
                "batch": "__batch__" if "__batch__" in merged_all else None,
                "condition": "__condition__" if "__condition__" in merged_all else None,
            },
            metric=str(payload["analysis_plan"]["primary_metric"]),
        )
    else:
        batch_levels, batch_summary = pd.DataFrame(), pd.DataFrame()

    enrichment = (
        pd.concat(enrichment_frames, ignore_index=True) if enrichment_frames else pd.DataFrame()
    )
    qc = pd.concat(qc_frames, ignore_index=True) if qc_frames else pd.DataFrame()
    if not enrichment.empty:
        for row in enrichment[enrichment["flag"].isin(["SAMPLE_ENRICHED", "NOT_ASSESSED"])].to_dict(
            orient="records"
        ):
            negative_rows.append(
                {
                    "cohort_id": row.get("cohort_id"),
                    "category": "sample_enrichment",
                    "finding": row.get("flag"),
                    "severity": "diagnostic_flag",
                    "detail": row.get("cluster"),
                }
            )
    if not qc.empty:
        for row in qc[qc["status"].astype(str).str.startswith("not_assessed")].to_dict(
            orient="records"
        ):
            negative_rows.append(
                {
                    "cohort_id": row.get("cohort_id"),
                    "category": row.get("diagnostic"),
                    "finding": row.get("status"),
                    "severity": "diagnostic_unassessed",
                }
            )
    if not batch_summary.empty:
        for row in batch_summary[batch_summary["status"] != "estimated_descriptive"].to_dict(
            orient="records"
        ):
            negative_rows.append(
                {
                    "cohort_id": "__release__",
                    "category": f"batch_sensitivity:{row.get('axis')}",
                    "finding": row.get("status"),
                    "severity": "diagnostic_unassessed",
                }
            )
    if not comparisons.empty:
        for row in comparisons[comparisons["status"] != "estimated"].to_dict(orient="records"):
            negative_rows.append(
                {
                    "cohort_id": "__release__",
                    "category": "method_comparison",
                    "finding": row.get("status"),
                    "severity": "underpowered",
                    "detail": f"{row.get('method_a')} vs {row.get('method_b')} {row.get('metric')}",
                }
            )
    negative = pd.DataFrame(negative_rows)

    _write_csv(inventory, output / "cohort_inventory.csv", artifacts)
    _write_csv(cohort_results, output / "benchmark_results_by_cohort.csv", artifacts)
    _write_csv(donor_metrics, output / "benchmark_results_by_donor.csv", artifacts)
    _write_csv(summary, output / "donor_weighted_summary.csv", artifacts)
    _write_csv(comparisons, output / "method_comparisons.csv", artifacts)
    _write_csv(batch_levels, output / "batch_sensitivity_levels.csv", artifacts)
    _write_csv(batch_summary, output / "batch_sensitivity_summary.csv", artifacts)
    _write_csv(enrichment, output / "sample_enrichment.csv", artifacts)
    _write_csv(qc, output / "qc_stratified_performance.csv", artifacts)
    _write_csv(negative, output / "negative_results.csv", artifacts)

    blocking = (
        negative[negative["severity"] == "blocks_release"] if not negative.empty else negative
    )
    readiness = {
        "status": "claim_ready" if blocking.empty else "incomplete_not_claim_ready",
        "blocking_findings": int(len(blocking)),
        "n_registered_cohorts": int(len(inventory)),
        "n_evaluated_cohorts": int(cohort_results["cohort_id"].nunique())
        if not cohort_results.empty
        else 0,
        "n_independent_donors": int(donor_metrics["donor_unit"].nunique())
        if not donor_metrics.empty
        else 0,
        "negative_results_retained": True,
    }
    report_path = output / "benchmark_report.md"
    report_path.write_text(
        _render_report(payload, inventory, summary, comparisons, negative, readiness),
        encoding="utf-8",
    )
    artifacts["benchmark_report"] = report_path

    manifest = {
        "schema_version": RELEASE_SCHEMA,
        "release_id": payload["release_id"],
        "registry_path": str(registry_path),
        "registry_sha256": _sha256(registry_path),
        "analysis_plan": payload["analysis_plan"],
        "required_methods": payload["required_methods"],
        "cohorts": json.loads(inventory.to_json(orient="records")),
        "readiness": readiness,
        "software": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "packages": _software_versions(),
        },
        "random_seed": seed,
        "bootstrap_replicates": n_boot,
        "artifacts": {
            name: {"path": str(path), "sha256": _sha256(path)}
            for name, path in artifacts.items()
        },
        "validation_scope": {
            "run_role": "public_multi_cohort_independent_benchmark_release",
            "statistical_unit": "donor",
            "cell_weighted_inference_prohibited": True,
            "truth_semantics": "expert_curated_reference_labels_not_infallible_ground_truth",
            "missing_diagnostic_semantics": "not_assessed_not_negative",
            "negative_result_policy": "retain_all_failures_unavailable_methods_and_null_results",
            "claim_boundary": (
                "Only claim performance for materialized immutable cohorts and completed "
                "comparators recorded in this manifest."
            ),
        },
    }
    manifest_path = output / "release_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    artifacts["release_manifest"] = manifest_path
    if not blocking.empty and not allow_incomplete:
        raise BenchmarkValidationError(
            f"Release has {len(blocking)} blocking findings and is not claim-ready; "
            f"auditable artifacts were written to {output}. Use allow_incomplete only to "
            "accept an explicitly incomplete release."
        )
    return artifacts
