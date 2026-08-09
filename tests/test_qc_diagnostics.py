"""QC diagnostic contracts: composable, missing→not_assessed, no identity rescue."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from celltypepilot.qc_diagnostics import (
    FORBIDDEN_CLEAN_TOKENS,
    QCDiagnosticError,
    assemble_qc_diagnostics,
    assess_high_mito,
    assess_low_rna,
    assess_sample_enrichment,
    assess_tool_axis,
    load_external_tool_table,
    write_qc_diagnostics,
)
from celltypepilot.robustness import qc_stratified_performance, sample_enrichment_diagnostics

REPO = Path(__file__).resolve().parents[1]
SCRATCH = REPO / "scratch" / "pytest_temp" / "qc_diagnostics"


def _tiny_adata(*, with_qc: bool = False) -> ad.AnnData:
    obs = pd.DataFrame(
        {
            "leiden": ["0", "0", "1", "1"],
            "sample": ["s1", "s1", "s2", "s2"],
        },
        index=[f"c{i}" for i in range(4)],
    )
    if with_qc:
        obs["n_genes_by_counts"] = [50, 300, 400, 500]
        obs["pct_counts_mt"] = [30.0, 5.0, 8.0, np.nan]
        obs["doublet_score"] = [0.9, 0.1, 0.05, 0.02]
    return ad.AnnData(X=np.ones((4, 2)), obs=obs, var=pd.DataFrame(index=["G1", "G2"]))


def test_missing_metadata_is_not_assessed_never_clean():
    adata = _tiny_adata(with_qc=False)
    report = assemble_qc_diagnostics(adata, cluster_key="leiden")
    assert report["can_rescue_identity"] is False
    assert report["missing_metadata_policy"] == "not_assessed_never_clean"
    # low_rna / high_mito / doublet / ambient missing → not_assessed
    for axis in ("low_rna", "high_mito", "doublet", "ambient_rna"):
        status = report["axes"][axis]["status"]
        assert status.startswith("not_assessed"), (axis, status)
        assert "clean" not in status.lower()
        assert report["axes"][axis]["flag"] == "NOT_ASSESSED"
    # sample enrichment can assess with sample+cluster present
    assert report["axes"]["sample_enrichment"]["status"].startswith("assessed")
    # rollup must not claim clean
    assert "clean" not in report["rollup_status"].lower()
    assert "clean" not in report["rollup_flag"].lower()
    for token in FORBIDDEN_CLEAN_TOKENS:
        assert token not in report["rollup_status"].lower()


def test_external_doublet_and_ambient_are_diagnostic_only():
    root = SCRATCH
    if root.exists():
        import shutil

        shutil.rmtree(root)
    root.mkdir(parents=True, exist_ok=True)
    adata = _tiny_adata(with_qc=False)
    doublet_csv = root / "doublet.csv"
    ambient_csv = root / "ambient.csv"
    pd.DataFrame(
        {
            "cell_id": ["c0", "c1", "c2", "c3"],
            "doublet_score": [0.9, 0.1, 0.05, 0.02],
        }
    ).to_csv(doublet_csv, index=False)
    pd.DataFrame(
        {
            "cell_id": ["c0", "c1", "c2", "c3"],
            "ambient_score": [0.5, 0.01, 0.0, 0.1],
        }
    ).to_csv(ambient_csv, index=False)

    doublet = load_external_tool_table(
        doublet_csv, axis="doublet", score_column="doublet_score", threshold=0.25
    )
    ambient = load_external_tool_table(
        ambient_csv, axis="ambient_rna", score_column="ambient_score", threshold=0.2
    )
    identity_before = {"cell_type": ["T cell", "B cell"], "decision": ["accepted", "abstain"]}
    report = assemble_qc_diagnostics(
        adata,
        cluster_key="leiden",
        doublet_table=doublet,
        ambient_table=ambient,
        identity_snapshot=identity_before,
    )
    assert report["axes"]["doublet"]["n_cells_flagged"] >= 1
    assert report["axes"]["ambient_rna"]["n_cells_flagged"] >= 1
    assert report["can_rescue_identity"] is False
    assert report["identity_effect"] == "none"
    # Identity snapshot keys must remain untouched conceptually
    assert identity_before["decision"] == ["accepted", "abstain"]


def test_obs_qc_columns_assessed():
    adata = _tiny_adata(with_qc=True)
    low = assess_low_rna(adata.obs, threshold=200)
    assert low.status.startswith("assessed")
    assert low.n_cells_flagged >= 1
    mito = assess_high_mito(adata.obs, threshold=0.2)
    assert mito.status.startswith("assessed")
    # pct_counts_mt uses percent scale
    assert mito.n_cells_flagged >= 1
    dbl = assess_tool_axis(
        adata.obs, axis="doublet", obs_aliases=("doublet_score",), threshold=0.25
    )
    assert dbl.n_cells_flagged >= 1


def test_forbidden_clean_token_rejected():
    with pytest.raises(QCDiagnosticError):
        from celltypepilot.qc_diagnostics import _validate_no_clean_claim

        _validate_no_clean_claim("clean", "PASS")


def test_sample_enrichment_missing_is_not_assessed():
    table = sample_enrichment_diagnostics(
        pd.DataFrame({"leiden": [0, 1]}),
        cluster_key="leiden",
        sample_key=None,
    )
    assert str(table.iloc[0]["status"]).startswith("not_assessed")
    assert table.iloc[0]["flag"] == "NOT_ASSESSED"

    empty = assess_sample_enrichment(
        pd.DataFrame({"leiden": ["0", "1"]}),
        cluster_key="leiden",
        sample_key=None,
    )
    assert empty.status.startswith("not_assessed")


def test_qc_stratified_missing_never_clean():
    merged = pd.DataFrame(
        {
            "method": ["m"] * 2,
            "study": ["s"] * 2,
            "donor": ["d"] * 2,
            "__truth__": ["A", "B"],
            "predicted_label": ["A", "B"],
        }
    )
    result = qc_stratified_performance(
        merged,
        study_key="study",
        donor_key="donor",
        diagnostics={},  # nothing predeclared
    )
    assert not result.empty
    for status in result["status"].astype(str):
        assert status.startswith("not_assessed")
        assert "clean" not in status.lower()


def test_write_qc_artifacts():
    root = SCRATCH / "write"
    root.mkdir(parents=True, exist_ok=True)
    adata = _tiny_adata(with_qc=True)
    report = assemble_qc_diagnostics(adata, cluster_key="leiden")
    paths = write_qc_diagnostics(report, root)
    payload = json.loads(paths["qc_diagnostics_json"].read_text(encoding="utf-8"))
    assert payload["can_rescue_identity"] is False
    csv = pd.read_csv(paths["qc_diagnostics_csv"])
    assert set(csv["axis"]) >= {"low_rna", "high_mito", "doublet", "ambient_rna"}
    assert (csv["can_rescue_identity"] == False).all()  # noqa: E712
