"""Fold-train-only popV adapter for the CellTypePilot comparator contract."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import popv


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: run_popv.py TRAIN_H5AD TEST_H5AD OUTPUT_CSV")
    train_path, test_path, output_path = map(Path, sys.argv[1:])
    train = ad.read_h5ad(train_path)
    query = ad.read_h5ad(test_path)
    if "cell_type" not in train.obs:
        raise ValueError("fold-train H5AD lacks cell_type labels")
    if train.raw is None or query.raw is None:
        raise ValueError("popV requires integer raw counts in AnnData.raw")
    train = train.raw.to_adata()
    query = query.raw.to_adata()
    shared = train.var_names.intersection(query.var_names)
    if len(shared) < 50:
        raise ValueError("fewer than 50 shared genes")
    train = train[:, shared].copy()
    query = query[:, shared].copy()
    train.obs["__ctp_batch__"] = "reference"
    query.obs["__ctp_batch__"] = "query"

    processed = popv.preprocessing.Process_Query(
        query,
        train,
        query_batch_key="__ctp_batch__",
        ref_labels_key="cell_type",
        ref_batch_key="__ctp_batch__",
        cl_obo_folder=False,
        unknown_celltype_label="Unknown",
        save_path_trained_models=str(output_path.parent / "popv_models"),
        prediction_mode="retrain",
        hvg=min(4000, len(shared)),
    ).adata
    popv.annotation.annotate_data(
        processed,
        save_path=str(output_path.parent / "popv_output"),
    )
    query_rows = processed.obs["_dataset"].astype(str).eq("query")
    result = processed.obs.loc[query_rows]
    if "popv_prediction" not in result or "popv_prediction_score" not in result:
        raise RuntimeError("popV did not emit prediction and agreement score columns")
    prediction_keys = processed.uns.get("prediction_keys_seen", processed.uns.get("prediction_keys", []))
    denominator = max(1, len(prediction_keys))
    score = pd.to_numeric(result["popv_prediction_score"], errors="raise").to_numpy(float)
    if np.nanmax(score) > 1:
        score = score / denominator
    output = pd.DataFrame(
        {
            "cell_id": result.index.astype(str),
            "predicted_label": result["popv_prediction"].astype(str).to_numpy(),
            "confidence": np.clip(score, 0, 1),
        }
    )
    output.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()
