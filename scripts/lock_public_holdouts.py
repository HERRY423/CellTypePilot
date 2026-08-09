"""Freeze donor-level test assignments for the public benchmark registry."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad

from celltypepilot.benchmark import build_holdout_assignments


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", required=True)
    parser.add_argument("--cohort", action="append", default=[])
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    registry_path = Path(args.registry).resolve()
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    selected = set(args.cohort)
    results: list[dict[str, object]] = []

    for cohort in registry["cohorts"]:
        cohort_id = str(cohort["cohort_id"])
        if selected and cohort_id not in selected:
            continue
        source = (registry_path.parent / cohort["local_path"]).resolve()
        output = (registry_path.parent / cohort["assignments_path"]).resolve()
        audit_path = output.with_name("holdout_plan.json")
        if output.exists() and not args.force:
            results.append({"cohort_id": cohort_id, "status": "already_present"})
            continue
        if source.stat().st_size != int(cohort["expected_bytes"]):
            raise ValueError(f"{cohort_id}: immutable byte-size check failed")

        dataset = ad.read_h5ad(source, backed="r")
        metadata = dataset.obs.copy()
        donor_key = str(cohort["metadata"]["donor_key"])
        declared_study_key = cohort["metadata"].get("study_key")
        study_key = str(declared_study_key or "__benchmark_study__")
        if declared_study_key is None:
            metadata[study_key] = str(cohort["constant_study_id"])
        assignments = build_holdout_assignments(
            metadata,
            study_key=study_key,
            donor_key=donor_key,
            strategy="donor",
        )

        output.parent.mkdir(parents=True, exist_ok=True)
        assignments.to_csv(output, index=False)
        audit = {
            "schema_version": "celltypepilot.public-holdout-plan.v1",
            "cohort_id": cohort_id,
            "registry_release_id": registry["release_id"],
            "input_path": str(source),
            "input_sha256": sha256(source),
            "strategy": "leave_one_donor_out",
            "study_key": study_key,
            "constant_study_id": cohort.get("constant_study_id"),
            "donor_key": donor_key,
            "n_cells": int(assignments["cell_id"].nunique()),
            "n_donors": int(assignments["held_out_donor"].nunique()),
            "n_folds": int(assignments["fold_id"].nunique()),
            "output_sha256": sha256(output),
        }
        audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
        results.append({"cohort_id": cohort_id, "status": "locked", **audit})

    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
