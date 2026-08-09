"""Create truth-blind, donor-local cluster maps for the public benchmark registry.

The immutable source H5AD files are never rewritten.  Clustering uses expression
only and is repeated independently inside each donor so that neither labels nor
other donors can influence the held-out donor's cluster boundaries.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def top_variance_mask(matrix: object, n_top_genes: int) -> np.ndarray:
    """Select expression-variance features without reading annotation labels."""
    if sparse.issparse(matrix):
        means = np.asarray(matrix.mean(axis=0)).ravel()
        squared_means = np.asarray(matrix.power(2).mean(axis=0)).ravel()
        variances = squared_means - means**2
    else:
        variances = np.var(np.asarray(matrix), axis=0)
    n_select = min(n_top_genes, len(variances))
    selected = np.argpartition(variances, -n_select)[-n_select:]
    mask = np.zeros(len(variances), dtype=bool)
    mask[selected] = True
    return mask


def cluster_donor(
    donor_data: ad.AnnData,
    *,
    seed: int,
) -> tuple[np.ndarray, dict[str, object]]:
    n_cells = donor_data.n_obs
    if n_cells < 50:
        return np.repeat("small_0", n_cells), {
            "status": "single_cluster_lt50_cells",
            "n_cells": n_cells,
            "n_clusters": 1,
        }
    working = donor_data.copy()
    n_components = min(15, working.n_obs - 1, working.n_vars - 1)
    if n_components < 2:
        return np.repeat("degenerate_0", n_cells), {
            "status": "single_cluster_insufficient_variable_genes",
            "n_cells": n_cells,
            "n_clusters": 1,
        }
    embedding = TruncatedSVD(
        n_components=n_components,
        algorithm="randomized",
        n_iter=3,
        random_state=seed,
    ).fit_transform(working.X)
    n_clusters = max(2, min(20, int(round(np.sqrt(n_cells / 100)))))
    labels = KMeans(
        n_clusters=n_clusters,
        random_state=seed,
        n_init=10,
    ).fit_predict(embedding).astype(str)
    return labels, {
        "status": "clustered_expression_only",
        "n_cells": n_cells,
        "n_clusters": int(pd.Series(labels).nunique()),
        "n_highly_variable_genes": int(working.n_vars),
        "n_pcs": n_components,
    }


def prepare_cohort(
    cohort: dict[str, object],
    base: Path,
    *,
    n_top_genes: int,
    seed: int,
    force: bool,
) -> dict[str, object]:
    source = (base / str(cohort["local_path"])).resolve()
    output = (base / str(cohort["cluster_map_path"])).resolve()
    if not source.exists():
        return {"cohort_id": cohort["cohort_id"], "status": "source_missing"}
    if output.exists() and not force:
        return {"cohort_id": cohort["cohort_id"], "status": "already_present"}
    if source.stat().st_size != int(cohort["expected_bytes"]):
        raise ValueError(f"{cohort['cohort_id']}: immutable byte-size check failed")

    dataset = ad.read_h5ad(source)
    donor_key = str(cohort["metadata"]["donor_key"])
    if donor_key not in dataset.obs:
        raise ValueError(f"{cohort['cohort_id']}: donor key {donor_key!r} missing")
    if dataset.obs[donor_key].isna().any():
        raise ValueError(f"{cohort['cohort_id']}: donor identifiers contain missing values")

    hvg_mask = top_variance_mask(dataset.X, n_top_genes)
    if hvg_mask.sum() < 50:
        raise ValueError(f"{cohort['cohort_id']}: fewer than 50 globally variable genes")

    assignments = pd.Series(index=dataset.obs_names.astype(str), dtype="object")
    donor_audit = []
    for donor in sorted(dataset.obs[donor_key].astype(str).unique()):
        mask = dataset.obs[donor_key].astype(str).to_numpy() == donor
        labels, audit = cluster_donor(
            dataset[mask, hvg_mask],
            seed=seed,
        )
        cell_ids = dataset.obs_names[mask].astype(str)
        assignments.loc[cell_ids] = [f"{donor}::{label}" for label in labels]
        donor_audit.append({"donor_id": donor, **audit})
    if assignments.isna().any():
        raise RuntimeError(f"{cohort['cohort_id']}: some cells were not clustered")

    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {"cell_id": assignments.index.astype(str), "cluster": assignments.to_numpy()}
    ).to_csv(output, index=False)
    audit_path = output.with_name("cluster_preparation.json")
    audit_payload = {
        "schema_version": "celltypepilot.truth-blind-clustering.v1",
        "cohort_id": cohort["cohort_id"],
        "input_path": str(source),
        "input_sha256": sha256(source),
        "output_path": str(output),
        "output_sha256": sha256(output),
        "truth_access": "prohibited_not_used",
        "metadata_used": [donor_key],
        "expression_used": "X",
        "feature_selection": "truth_blind_global_top_expression_variance",
        "cross_donor_clustering": False,
        "parameters": {
            "n_top_genes": n_top_genes,
            "seed": seed,
            "decomposition": "sklearn_TruncatedSVD_randomized_n_iter_3",
            "clustering": "sklearn_KMeans",
            "cluster_count_rule": "max(2,min(20,round(sqrt(n_cells/100))))",
        },
        "donors": donor_audit,
    }
    audit_path.write_text(json.dumps(audit_payload, indent=2), encoding="utf-8")
    return {
        "cohort_id": cohort["cohort_id"],
        "status": "prepared",
        "n_cells": dataset.n_obs,
        "n_donors": dataset.obs[donor_key].nunique(),
        "cluster_map": str(output),
        "audit": str(audit_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", required=True)
    parser.add_argument("--cohort", action="append", default=[])
    parser.add_argument("--n-top-genes", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    registry_path = Path(args.registry).resolve()
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    selected = set(args.cohort)
    results = []
    for cohort in registry["cohorts"]:
        if selected and cohort["cohort_id"] not in selected:
            continue
        results.append(
            prepare_cohort(
                cohort,
                registry_path.parent,
                n_top_genes=args.n_top_genes,
                seed=args.seed,
                force=args.force,
            )
        )
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
