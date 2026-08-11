"""Outcome-blind donor-role locking for independent calibration."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import anndata as ad
import pandas as pd

CALIBRATION_SPLIT_SCHEMA = "celltypepilot.calibration-split.v1"


class CalibrationSplitError(ValueError):
    """Raised when calibration and evaluation roles cannot be made disjoint."""


def _canonical_hash(payload: dict) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode(
            "utf-8"
        )
    ).hexdigest()


def build_calibration_split(
    registry_path: str | Path,
    output_dir: str | Path,
    *,
    calibration_fraction: float = 0.2,
    seed: str = "celltypepilot-calibration-v1",
) -> dict:
    """Lock donor-disjoint calibration/evaluation roles without reading labels."""
    if not 0 < calibration_fraction < 0.5:
        raise CalibrationSplitError("calibration_fraction must be greater than 0 and below 0.5")
    registry_input = Path(registry_path)
    registry_display = registry_input.as_posix()
    registry_path = registry_input.resolve()
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    registry_root = registry_path.parent
    rows: list[dict] = []
    all_units: set[str] = set()
    cohort_units: dict[str, list[str]] = {}
    cohort_specs: dict[str, dict] = {}
    for cohort in registry.get("cohorts", []):
        local = Path(cohort["local_path"])
        if not local.is_absolute():
            local = registry_root / local
        if not local.is_file():
            raise CalibrationSplitError(f"Cohort data missing: {local}")
        donor_key = cohort.get("metadata", {}).get("donor_key")
        if not donor_key:
            raise CalibrationSplitError(f"Cohort {cohort['cohort_id']} lacks donor_key")
        dataset = ad.read_h5ad(local, backed="r")
        try:
            if donor_key not in dataset.obs:
                raise CalibrationSplitError(
                    f"Cohort {cohort['cohort_id']} lacks donor column {donor_key}"
                )
            namespace = str(
                cohort.get("donor_namespace")
                or cohort.get("constant_study_id")
                or cohort["cohort_id"]
            )
            units = sorted(
                {f"{namespace}::{value}" for value in dataset.obs[donor_key].astype(str)}
            )
        finally:
            if dataset.file is not None:
                dataset.file.close()
        if len(units) < 3:
            raise CalibrationSplitError(
                f"Cohort {cohort['cohort_id']} needs at least three donors for donor-disjoint roles"
            )
        cohort_units[str(cohort["cohort_id"])] = units
        cohort_specs[str(cohort["cohort_id"])] = cohort
        all_units.update(units)

    # Shared donor namespaces (for example two platforms from one study) receive
    # one global role, preventing cross-platform role leakage.
    scores = {unit: hashlib.sha256(f"{seed}\n{unit}".encode()).hexdigest() for unit in all_units}
    calibration_units: set[str] = set()
    for units in cohort_units.values():
        n_calibration = max(1, int(round(len(units) * calibration_fraction)))
        n_calibration = min(n_calibration, len(units) - 2)
        calibration_units.update(
            sorted(units, key=lambda unit: (scores[unit], unit))[:n_calibration]
        )

    # If a shared unit was selected by either platform, it is calibration in all
    # cohorts. Preserve at least two evaluation donors per cohort.
    for cohort_id, units in cohort_units.items():
        selected = [unit for unit in units if unit in calibration_units]
        if len(units) - len(selected) < 2:
            for unit in sorted(selected, key=lambda value: (scores[value], value), reverse=True):
                calibration_units.remove(unit)
                selected.remove(unit)
                if len(units) - len(selected) >= 2:
                    break
        cohort = cohort_specs[cohort_id]
        for unit in units:
            rows.append(
                {
                    "cohort_id": cohort_id,
                    "dataset_version_id": cohort.get("dataset_version_id", ""),
                    "donor_unit": unit,
                    "role": "calibration" if unit in calibration_units else "evaluation",
                    "assignment_sha256": scores[unit],
                }
            )

    assignments = pd.DataFrame(rows).sort_values(["cohort_id", "donor_unit"])
    overlap = set(assignments.loc[assignments["role"] == "calibration", "donor_unit"]) & set(
        assignments.loc[assignments["role"] == "evaluation", "donor_unit"]
    )
    if overlap:
        raise CalibrationSplitError(f"Donor roles overlap: {sorted(overlap)[:3]}")
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    assignments_path = output / "donor_role_assignments.csv"
    assignments.to_csv(assignments_path, index=False, lineterminator="\n")
    manifest = {
        "schema_version": CALIBRATION_SPLIT_SCHEMA,
        "registry_path": registry_display,
        "registry_sha256": hashlib.sha256(registry_path.read_bytes()).hexdigest(),
        "calibration_fraction": calibration_fraction,
        "seed": seed,
        "independent_unit": "donor",
        "role_scope": "donor_independent_within_cohort_not_independent_cohort",
        "truth_access": "prohibited_not_read",
        "n_unique_donors": len(all_units),
        "n_calibration_donors": len(calibration_units),
        "n_evaluation_donors": len(all_units - calibration_units),
        "cohort_role_counts": assignments.groupby(["cohort_id", "role"])
        .size()
        .unstack(fill_value=0)
        .to_dict(orient="index"),
        "claim_boundary": (
            "This split supports donor-independent calibration. It is not an independent-cohort "
            "calibration claim and cannot satisfy domains that require a separate cohort."
        ),
    }
    manifest["split_sha256"] = _canonical_hash(manifest)
    manifest_path = output / "calibration_split_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {
        "assignments": assignments,
        "assignments_path": assignments_path,
        "manifest": manifest,
        "manifest_path": manifest_path,
    }


def fit_policy_from_locked_donors(
    registry_path: str | Path,
    assignments_path: str | Path,
    cohort_id: str,
    predictions_path: str | Path,
    label_map_path: str | Path,
    output_path: str | Path,
    *,
    max_selective_error: float = 0.25,
    min_coverage: float = 0.2,
) -> dict:
    """Fit a downgrade-only policy on donors preassigned to calibration."""
    from .benchmark import apply_truth_label_map
    from .benchmark_runner import apply_locked_label_map
    from .calibration import fit_abstention_policy, save_abstention_policy

    registry_path = Path(registry_path).resolve()
    assignments_path = Path(assignments_path).resolve()
    predictions_path = Path(predictions_path).resolve()
    label_map_path = Path(label_map_path).resolve()
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    matches = [item for item in registry.get("cohorts", []) if item["cohort_id"] == cohort_id]
    if len(matches) != 1:
        raise CalibrationSplitError(f"Expected one registry cohort named {cohort_id!r}")
    cohort = matches[0]
    roles = pd.read_csv(assignments_path, dtype=str)
    roles = roles[roles["cohort_id"] == cohort_id]
    calibration_units = set(roles.loc[roles["role"] == "calibration", "donor_unit"])
    evaluation_units = set(roles.loc[roles["role"] == "evaluation", "donor_unit"])
    if not calibration_units or not evaluation_units or calibration_units & evaluation_units:
        raise CalibrationSplitError("Locked calibration/evaluation donor roles are invalid")
    data_path = Path(cohort["local_path"])
    if not data_path.is_absolute():
        data_path = registry_path.parent / data_path
    dataset = ad.read_h5ad(data_path, backed="r")
    try:
        donor_key = cohort["metadata"]["donor_key"]
        truth_key = cohort["metadata"]["truth_key"]
        namespace = str(
            cohort.get("donor_namespace") or cohort.get("constant_study_id") or cohort_id
        )
        donor_units = dataset.obs[donor_key].astype(str).map(lambda donor: f"{namespace}::{donor}")
        calibration_cells = dataset.obs_names[donor_units.isin(calibration_units)].astype(str)
        truth = dataset.obs[truth_key].astype(str).copy()
        truth.index = dataset.obs_names.astype(str)
    finally:
        if dataset.file is not None:
            dataset.file.close()
    predictions = pd.read_csv(predictions_path, dtype={"cell_id": str})
    predictions = predictions[predictions["cell_id"].astype(str).isin(calibration_cells)].copy()
    if predictions.empty:
        raise CalibrationSplitError("No selector predictions belong to calibration donors")
    label_map = pd.read_csv(label_map_path, dtype=str)
    # v0.4 selectors emit canonical Atlas labels. Derive any absent product-map
    # row only through ontology equivalence to an already locked raw label; never
    # use truth frequencies, errors, or performance to choose a mapping.
    from .data_adapter import load_marker_atlas
    from .identity_contract import (
        build_identity_resolver,
        collect_pack_identity_contract,
        compose_marker_definitions,
        resolve_identity_label,
    )
    from .pack_manager import merge_marker_atlas, resolve_extension_packs
    from .validation_domains import assess_validation_domain

    domain = assess_validation_domain(cohort.get("tissue"))
    pack_names = domain.get("requirements", {}).get("atlas_contract", {}).get("required_packs", [])
    records, _ = resolve_extension_packs(pack_names, cohort.get("species", "human"))
    atlas = load_marker_atlas(cohort.get("species", "human"))
    atlas, _ = merge_marker_atlas(atlas, records, cohort.get("species", "human"))
    pack_contract = collect_pack_identity_contract(records)
    _, scope = compose_marker_definitions(
        atlas, cohort.get("tissue", "general"), pack_contract=pack_contract
    )
    resolver = build_identity_resolver(atlas, scope["active_tissues"], pack_contract)
    identity_to_evaluation: dict[str, set[str]] = {}
    for row in label_map.itertuples(index=False):
        resolution = resolve_identity_label(str(row.raw_label), resolver)
        if resolution["resolution"] in {"unresolved", "ambiguous_alias_unresolved"}:
            continue
        identity_to_evaluation.setdefault(resolution["canonical_label"], set()).add(
            str(row.canonical_label)
        )
    existing_product = set(
        label_map.loc[label_map["method"] == "celltypepilot", "raw_label"].astype(str)
    )
    additions = []
    for raw_label in sorted(set(predictions["predicted_label"].astype(str)) - existing_product):
        resolution = resolve_identity_label(raw_label, resolver)
        choices = identity_to_evaluation.get(resolution["canonical_label"], set())
        if len(choices) != 1:
            raise CalibrationSplitError(
                f"Cannot derive one ontology-equivalent evaluation mapping for {raw_label!r}"
            )
        additions.append(
            {
                "method": "celltypepilot",
                "raw_label": raw_label,
                "canonical_label": next(iter(choices)),
            }
        )
    extended_label_map = (
        pd.concat([label_map, pd.DataFrame(additions)], ignore_index=True)
        if additions
        else label_map
    )
    extended_path = Path(output_path).resolve().with_suffix(".label_map.csv")
    extended_path.parent.mkdir(parents=True, exist_ok=True)
    extended_label_map.to_csv(extended_path, index=False, lineterminator="\n")
    mapped_predictions = apply_locked_label_map(predictions, extended_label_map)
    mapped_truth = apply_truth_label_map(truth.loc[calibration_cells], label_map)
    policy = fit_abstention_policy(
        mapped_truth,
        mapped_predictions,
        "celltypepilot",
        max_selective_error=max_selective_error,
        min_coverage=min_coverage,
        dataset_role="calibration",
    )
    policy["independent_unit"] = "donor"
    policy["independence_scope"] = "donor_independent_within_cohort_not_independent_cohort"
    policy["cohort_id"] = cohort_id
    policy["calibration_donor_units"] = sorted(calibration_units)
    policy["n_calibration_donors"] = len(calibration_units)
    policy["n_calibration_cells"] = len(predictions)
    policy["source_sha256"] = {
        "registry": hashlib.sha256(registry_path.read_bytes()).hexdigest(),
        "assignments": hashlib.sha256(assignments_path.read_bytes()).hexdigest(),
        "predictions": hashlib.sha256(predictions_path.read_bytes()).hexdigest(),
        "label_map": hashlib.sha256(label_map_path.read_bytes()).hexdigest(),
        "extended_label_map": hashlib.sha256(extended_path.read_bytes()).hexdigest(),
    }
    policy["label_map_extension_policy"] = "ontology_equivalence_only_no_outcome_optimization"
    policy["extended_label_map_path"] = extended_path.name
    policy["claim_boundary"] = (
        "Threshold is fit on donor-disjoint calibration cells and is downgrade-only. "
        "It is not independent-cohort calibration or a guarantee on evaluation donors."
    )
    saved = save_abstention_policy(policy, output_path)
    return {"policy": policy, "policy_path": saved}
