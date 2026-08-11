"""Recoverable execution plan for the three depth-validation domains."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd

from .validation_domains import load_validation_domains

DOMAIN_PLAN_SCHEMA = "celltypepilot.domain-validation-plan.v1"
DOMAIN_RUN_SCHEMA = "celltypepilot.domain-validation-run.v1"
PRODUCT_METHOD = "celltypepilot"
FOLD_NATIVE_TRUTH_LABEL_METHODS = {"scanvi", "custom_reference"}


class DomainValidationPipelineError(ValueError):
    """Raised when a domain-validation plan cannot be built or executed safely."""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _canonical_hash(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode(
            "utf-8"
        )
    ).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_json(path: Path, payload: dict) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


def _domain_for_cohort(cohort: dict) -> str | None:
    cohort_id = str(cohort.get("cohort_id", "")).casefold()
    title = str(cohort.get("title", "")).casefold()
    tissue = str(cohort.get("tissue", "")).casefold()
    if tissue == "lung" or "lung" in cohort_id:
        return "lung"
    if tissue in {"gut", "colon", "intestine"} or any(
        token in cohort_id for token in ("ibd", "colon", "gut")
    ):
        return "gut_ibd"
    if any(token in f"{cohort_id} {title}" for token in ("tumor", "cancer", "crc", "carcinoma")):
        return "tumor_microenvironment"
    return None


def _resolve_registry_path(registry_root: Path, value: str | None) -> Path | None:
    if not value:
        return None
    path = Path(value)
    return (path if path.is_absolute() else registry_root / path).resolve()


def _inspect_local_cohort(cohort: dict, registry_root: Path, required_methods: list[str]) -> dict:
    data_path = _resolve_registry_path(registry_root, cohort.get("local_path"))
    label_map_path = _resolve_registry_path(registry_root, cohort.get("label_map_path"))
    cluster_map_path = _resolve_registry_path(registry_root, cohort.get("cluster_map_path"))
    metadata = cohort.get("metadata", {})
    blockers: list[str] = []
    donors = 0
    platform_families = 0
    donor_units: list[str] = []
    platform_values: list[str] = []
    if data_path is None or not data_path.is_file():
        blockers.append("LOCAL_H5AD_MISSING")
    else:
        backed = ad.read_h5ad(data_path, backed="r")
        try:
            donor_key = metadata.get("donor_key")
            platform_key = metadata.get("platform_key")
            cluster_key = metadata.get("cluster_key")
            truth_key = metadata.get("truth_key")
            if not donor_key or donor_key not in backed.obs:
                blockers.append("DONOR_KEY_MISSING")
            else:
                namespace = str(
                    cohort.get("donor_namespace") or cohort.get("constant_study_id") or "cohort"
                )
                donor_units = sorted(
                    {f"{namespace}::{value}" for value in backed.obs[donor_key].astype(str)}
                )
                donors = len(donor_units)
            if platform_key and platform_key in backed.obs:
                platform_values = sorted(set(backed.obs[platform_key].astype(str)))
                platform_families = len(platform_values)
            if not truth_key or truth_key not in backed.obs:
                blockers.append("TRUTH_KEY_MISSING")
            if (not cluster_key or cluster_key not in backed.obs) and (
                cluster_map_path is None or not cluster_map_path.is_file()
            ):
                blockers.append("CLUSTER_ASSIGNMENT_MISSING")
        finally:
            if backed.file is not None:
                backed.file.close()
    mapped_methods: set[str] = set()
    if label_map_path is None or not label_map_path.is_file():
        blockers.append("LOCKED_LABEL_MAP_MISSING")
    else:
        label_map = pd.read_csv(label_map_path, dtype=str)
        if not {"method", "raw_label", "canonical_label"}.issubset(label_map.columns):
            blockers.append("LOCKED_LABEL_MAP_INVALID")
        else:
            mapped_methods = set(label_map["method"].dropna().astype(str).str.casefold())
            missing = sorted(set(required_methods) - mapped_methods)
            if missing:
                blockers.append("LABEL_MAP_METHODS_MISSING:" + ",".join(missing))
    return {
        "data_path": str(data_path) if data_path else None,
        "label_map_path": str(label_map_path) if label_map_path else None,
        "cluster_map_path": str(cluster_map_path) if cluster_map_path else None,
        "donors": donors,
        "donor_units": donor_units,
        "platform_families": platform_families,
        "platform_values": platform_values,
        "mapped_methods": sorted(mapped_methods),
        "blockers": blockers,
        "execution_status": "ready" if not blockers else "blocked",
    }


def _lock_extended_label_map(source: Path, destination: Path) -> Path:
    """Predeclare mappings for fold-native methods whose label set is fold truth itself."""
    if not source.is_file():
        return source
    frame = pd.read_csv(source, dtype=str)
    required = {"method", "raw_label", "canonical_label"}
    if not required.issubset(frame.columns):
        return source
    truth = frame[frame["method"].astype(str) == "__truth__"].copy()
    additions = []
    present = set(frame["method"].astype(str).str.casefold())
    for method in sorted(FOLD_NATIVE_TRUTH_LABEL_METHODS - present):
        copied = truth.copy()
        copied["method"] = method
        additions.append(copied)
    locked = pd.concat([frame, *additions], ignore_index=True) if additions else frame
    destination.parent.mkdir(parents=True, exist_ok=True)
    locked.to_csv(destination, index=False, lineterminator="\n")
    return destination


def build_domain_validation_plan(
    registry_path: str | Path,
    output_dir: str | Path,
) -> dict:
    """Compile current cohort assets into an immutable, fail-closed domain plan."""
    registry_path = Path(registry_path).resolve()
    if not registry_path.is_file():
        raise DomainValidationPipelineError(f"Public registry not found: {registry_path}")
    try:
        registry = json.loads(registry_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise DomainValidationPipelineError(f"Invalid public registry JSON: {exc}") from exc
    domains = load_validation_domains()
    root = Path(output_dir).resolve()
    registry_root = registry_path.parent
    cohorts = []
    for cohort in registry.get("cohorts", []):
        domain_id = _domain_for_cohort(cohort)
        if domain_id is None:
            continue
        required = [
            PRODUCT_METHOD,
            *domains["domains"][domain_id]["required_candidate_backends"],
        ]
        required_packs = list(domains["domains"][domain_id]["atlas_contract"]["required_packs"])
        original_label_map = _resolve_registry_path(registry_root, cohort.get("label_map_path"))
        cohort_for_plan = dict(cohort)
        if original_label_map is not None:
            locked_label_map = _lock_extended_label_map(
                original_label_map,
                root / "locked_label_maps" / f"{cohort['cohort_id']}.csv",
            )
            cohort_for_plan["label_map_path"] = str(locked_label_map)
        local = _inspect_local_cohort(cohort_for_plan, registry_root, required)
        data_sha256 = cohort.get("expected_sha256")
        label_map_sha256 = None
        if local["label_map_path"] and Path(local["label_map_path"]).is_file():
            label_map_sha256 = _sha256_file(Path(local["label_map_path"]))
        cohorts.append(
            {
                "domain_id": domain_id,
                "cohort_id": cohort["cohort_id"],
                "title": cohort.get("title"),
                "species": cohort.get("species"),
                "tissue": cohort.get("tissue"),
                "constant_study_id": cohort.get("constant_study_id"),
                "donor_namespace": cohort.get("donor_namespace"),
                "metadata": cohort.get("metadata", {}),
                "required_methods": required,
                "required_packs": required_packs,
                "expected_data_sha256": data_sha256,
                "locked_label_map_sha256": label_map_sha256,
                **local,
            }
        )

    summaries = {}
    for domain_id, domain in domains["domains"].items():
        selected = [item for item in cohorts if item["domain_id"] == domain_id]
        studies = {str(item["constant_study_id"]) for item in selected if item["constant_study_id"]}
        donor_units = {unit for item in selected for unit in item.get("donor_units", [])}
        platform_values = {value for item in selected for value in item.get("platform_values", [])}
        donors = len(donor_units)
        platforms = len(platform_values)
        minimum = domain["minimum_evidence"]
        blockers = list(domain["current_blockers"])
        if len(selected) < minimum["independent_public_cohorts"]:
            blockers.append("INSUFFICIENT_PUBLIC_COHORTS")
        if len(studies) < minimum["independent_studies"]:
            blockers.append("INSUFFICIENT_INDEPENDENT_STUDIES")
        if donors < minimum["donors"]:
            blockers.append("INSUFFICIENT_DONORS")
        if platforms < minimum["platform_families"]:
            blockers.append("INSUFFICIENT_PLATFORM_FAMILIES")
        if any(item["execution_status"] != "ready" for item in selected):
            blockers.append("COHORT_ASSETS_OR_LABEL_MAPS_INCOMPLETE")
        summaries[domain_id] = {
            "atlas_contract": domain["atlas_contract"],
            "registered_cohorts": len(selected),
            "independent_studies": len(studies),
            "observed_donors": donors,
            "observed_platform_families": platforms,
            "ready_cohorts": sum(item["execution_status"] == "ready" for item in selected),
            "claim_ready": False,
            "blockers": list(dict.fromkeys(blockers)),
            "minimum_evidence": minimum,
        }
    plan = {
        "schema_version": DOMAIN_PLAN_SCHEMA,
        "created_at_utc": _utc_now(),
        "registry_path": str(registry_path),
        "registry_sha256": hashlib.sha256(registry_path.read_bytes()).hexdigest(),
        "output_dir": str(root),
        "domains": summaries,
        "cohorts": cohorts,
        "execution_policy": {
            "holdout": "leave_one_study_out_else_leave_one_donor_out",
            "truth_exposure": "test_truth_and_label_like_obs_removed_before_backend_execution",
            "evaluation_units": ["cell", "cluster"],
            "continue_on_unavailable": True,
            "checkpoint_contract": "atomic_per_method_fold_v1",
            "llm_role": "excluded_from_accuracy_comparator_family_hypothesis_only",
        },
        "claim_boundary": (
            "This plan inventories and executes evidence work. It cannot mark a domain claim-ready "
            "without the predeclared cohort, calibration, adjudication, and robustness artifacts."
        ),
    }
    plan["plan_sha256"] = _canonical_hash(plan)
    plan_path = _write_json(root / "domain_validation_plan.json", plan)
    cohort_status = pd.DataFrame(cohorts)
    cohort_status.to_csv(root / "domain_cohort_status.csv", index=False)
    return {"plan": plan, "plan_path": plan_path, "cohort_status": cohort_status}


def _load_plan(path: str | Path) -> dict:
    plan_path = Path(path).resolve()
    if not plan_path.is_file():
        raise DomainValidationPipelineError(f"Domain validation plan not found: {plan_path}")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("schema_version") != DOMAIN_PLAN_SCHEMA:
        raise DomainValidationPipelineError(f"Plan schema must be {DOMAIN_PLAN_SCHEMA}")
    unsigned = {key: value for key, value in plan.items() if key != "plan_sha256"}
    if _canonical_hash(unsigned) != plan.get("plan_sha256"):
        raise DomainValidationPipelineError("Domain validation plan hash mismatch")
    registry = Path(plan["registry_path"])
    if hashlib.sha256(registry.read_bytes()).hexdigest() != plan["registry_sha256"]:
        raise DomainValidationPipelineError("Public registry changed after planning")
    return plan


def execute_domain_validation_plan(
    plan_path: str | Path,
    *,
    domain_ids: set[str] | None = None,
    cohort_ids: set[str] | None = None,
) -> dict:
    """Execute all ready cohorts and retain every unavailable/failed method status."""
    from .benchmark import (
        apply_truth_label_map,
        build_cluster_level_track,
        build_holdout_assignments,
        evaluate_holdout_predictions,
        save_benchmark_plan,
    )
    from .benchmark_runner import run_benchmark_comparators

    plan = _load_plan(plan_path)
    root = Path(plan["output_dir"])
    run_rows = []
    for cohort in plan["cohorts"]:
        if domain_ids and cohort["domain_id"] not in domain_ids:
            continue
        if cohort_ids and cohort["cohort_id"] not in cohort_ids:
            continue
        if cohort["execution_status"] != "ready":
            run_rows.append(
                {
                    "domain_id": cohort["domain_id"],
                    "cohort_id": cohort["cohort_id"],
                    "status": "blocked",
                    "detail": ";".join(cohort["blockers"]),
                }
            )
            continue
        cohort_output = root / cohort["domain_id"] / cohort["cohort_id"]
        cohort_output.mkdir(parents=True, exist_ok=True)
        try:
            data_path = Path(cohort["data_path"])
            expected_data_hash = cohort.get("expected_data_sha256")
            if expected_data_hash:
                observed_data_hash = _sha256_file(data_path)
                if observed_data_hash != expected_data_hash:
                    raise DomainValidationPipelineError("Cohort H5AD hash mismatch")
            label_map_path = Path(cohort["label_map_path"])
            if _sha256_file(label_map_path) != cohort.get("locked_label_map_sha256"):
                raise DomainValidationPipelineError("Locked label map hash mismatch")
            adata = ad.read_h5ad(cohort["data_path"])
            metadata = cohort["metadata"]
            truth_key = metadata["truth_key"]
            donor_key = metadata["donor_key"]
            cluster_key = metadata["cluster_key"]
            study_key = metadata.get("study_key") or "__ctp_locked_study__"
            if study_key not in adata.obs:
                adata.obs[study_key] = str(cohort["constant_study_id"])
            if cluster_key not in adata.obs:
                cluster_map = pd.read_csv(cohort["cluster_map_path"], dtype=str)
                mapping = cluster_map.set_index("cell_id")["cluster"]
                mapping.index = mapping.index.astype(str)
                expected = pd.Index(adata.obs_names.astype(str))
                if len(expected.difference(mapping.index)) or len(
                    mapping.index.difference(expected)
                ):
                    raise DomainValidationPipelineError("Cluster map cell set mismatch")
                adata.obs[cluster_key] = mapping.reindex(expected).to_numpy()
            strategy = "study" if adata.obs[study_key].astype(str).nunique() >= 2 else "donor"
            assignments = build_holdout_assignments(adata.obs, study_key, donor_key, strategy)
            save_benchmark_plan(assignments, cohort_output, study_key, donor_key, strategy)
            label_map = pd.read_csv(cohort["label_map_path"], dtype=str)
            evaluation_truth = apply_truth_label_map(adata.obs[truth_key], label_map)
            predictions, statuses = run_benchmark_comparators(
                adata,
                assignments,
                truth_key,
                cluster_key,
                cohort_output,
                cohort["species"],
                cohort["tissue"],
                packs=tuple(cohort.get("required_packs", ())),
                methods=tuple(cohort["required_methods"]),
                label_map=label_map,
                continue_on_unavailable=True,
            )
            completed_pairs = {
                (str(row.method), str(row.fold_id))
                for row in statuses.loc[statuses["status"] == "completed"].itertuples()
            }
            if not predictions.empty:
                aggregate, per_fold = evaluate_holdout_predictions(
                    evaluation_truth, assignments, predictions
                )
                aggregate.to_csv(cohort_output / "benchmark_results.csv", index=False)
                per_fold.to_csv(cohort_output / "benchmark_results_by_fold.csv", index=False)
                cluster_truth, cluster_assignments, cluster_predictions, diagnostics = (
                    build_cluster_level_track(
                        evaluation_truth,
                        assignments,
                        predictions,
                        adata.obs[cluster_key],
                    )
                )
                cluster_results, cluster_folds = evaluate_holdout_predictions(
                    cluster_truth, cluster_assignments, cluster_predictions
                )
                cluster_results.to_csv(cohort_output / "cluster_track_results.csv", index=False)
                cluster_folds.to_csv(
                    cohort_output / "cluster_track_results_by_fold.csv", index=False
                )
                diagnostics.to_csv(cohort_output / "cluster_track_diagnostics.csv", index=False)
            expected_pairs = {
                (method, fold_id)
                for method in cohort["required_methods"]
                for fold_id in assignments["fold_id"].drop_duplicates().astype(str)
            }
            run_rows.append(
                {
                    "domain_id": cohort["domain_id"],
                    "cohort_id": cohort["cohort_id"],
                    "status": ("completed" if expected_pairs <= completed_pairs else "incomplete"),
                    "detail": (
                        "" if expected_pairs <= completed_pairs else "METHOD_FOLDS_INCOMPLETE"
                    ),
                }
            )
        except Exception as exc:
            run_rows.append(
                {
                    "domain_id": cohort["domain_id"],
                    "cohort_id": cohort["cohort_id"],
                    "status": "failed",
                    "detail": " ".join(str(exc).split())[:1000],
                }
            )
    status = pd.DataFrame(run_rows)
    status_path = root / "domain_validation_run_status.csv"
    status.to_csv(status_path, index=False)
    manifest = {
        "schema_version": DOMAIN_RUN_SCHEMA,
        "completed_at_utc": _utc_now(),
        "plan_sha256": plan["plan_sha256"],
        "status_counts": status["status"].value_counts().to_dict() if not status.empty else {},
        "claim_ready": False,
        "claim_boundary": (
            "Execution completion is not domain validation. Separate calibration, expert "
            "adjudication, minimum independent evidence, and release audit remain mandatory."
        ),
    }
    manifest_path = _write_json(root / "domain_validation_run_manifest.json", manifest)
    return {
        "status": status,
        "status_path": status_path,
        "manifest": manifest,
        "manifest_path": manifest_path,
    }
