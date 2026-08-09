"""Review panel + resign contracts."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from celltypepilot.review_panel import build_cluster_review_panel
from celltypepilot.review_resign import resign_review_outputs

REPO = Path(__file__).resolve().parents[1]
SCRATCH = REPO / "scratch" / "pytest_temp" / "review_panel"


def test_review_panel_three_axes():
    evidence = pd.DataFrame(
        [
            {
                "cluster": "0",
                "cell_type": "T cell",
                "candidate_cell_type": "T cell",
                "decision": "accepted",
                "cl_id": "CL:0000084",
                "combined_score": 0.8,
                "critic_confidence": "high",
                "critic_flags": "PASS",
                "pos_supporting_markers": "CD3E,CD3D",
                "neg_expressed_markers": "MS4A1",
                "pos_silent_markers": "CD2",
                "pos_missing_markers": "TRAC",
                "pct_overlap": 0.5,
                "novelty_decision": "known_supported",
                "novelty_score": 0.1,
                "marker_provenance_status": "literature_cooccurrence_supported",
                "marker_provenance_sources": "cellmarker_2_0",
            }
        ]
    )
    state = pd.DataFrame(
        [
            {
                "cluster": "0",
                "state_candidate": "activated",
                "state_decision": "supported",
                "state_score": 0.6,
                "state_confidence": "medium",
                "pos_supporting_markers": "IL2RA",
            }
        ]
    )
    novelty = pd.DataFrame(
        [
            {
                "cluster": "0",
                "novelty_decision": "known_supported",
                "novelty_score": 0.1,
                "top_unmapped_markers": "GENEX",
            }
        ]
    )
    obs = pd.DataFrame(
        {
            "ctp_cl_id": ["0", "0", "0"],
            "donor_id": ["d1", "d1", "d2"],
            "batch": ["b1", "b1", "b2"],
            "sample": ["s1", "s1", "s2"],
        },
        index=["c1", "c2", "c3"],
    )
    adata = ad.AnnData(X=np.ones((3, 2)), obs=obs, var=pd.DataFrame(index=["G1", "G2"]))
    panel = build_cluster_review_panel(
        cluster="0",
        evidence=evidence,
        state_results=state,
        novelty_results=novelty,
        adata=adata,
        cluster_key="ctp_cl_id",
    )
    assert panel["axes"]["identity"]["cell_type"] == "T cell"
    assert "CD3E" in panel["axes"]["identity"]["supporting_markers"]
    assert "MS4A1" in panel["axes"]["identity"]["opposing_markers"]
    assert panel["axes"]["state"]["state_candidate"] == "activated"
    assert panel["axes"]["novelty"]["novelty_decision"] == "known_supported"
    assert panel["donor_batch_strata"]["donors"]["status"] == "assessed"
    assert panel["edit_policy"]["append_only_audit"] is True
    assert panel["literature"]["status"] != ""


def test_resign_clears_stale():
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    out = SCRATCH / "run"
    out.mkdir(parents=True)
    # Minimal annotated products
    obs = pd.DataFrame({"ctp_cl_id": ["0"], "ctp_cell_type": ["T cell"]}, index=["c1"])
    ad.AnnData(X=np.ones((1, 2)), obs=obs, var=pd.DataFrame(index=["G1", "G2"])).write_h5ad(
        out / "data.annotated.h5ad"
    )
    pd.DataFrame(
        [
            {
                "cluster": "0",
                "cell_type": "T cell",
                "combined_score": 0.5,
                "critic_confidence": "medium",
                "critic_flags": "PASS",
            }
        ]
    ).to_csv(out / "evidence_table.csv", index=False)
    (out / "artifact_status.json").write_text(
        json.dumps(
            {
                "review_state": "applied_overrides_artifacts_stale",
                "stale_artifacts": ["report_draft.html", "manifest.json"],
                "message": "stale",
            }
        ),
        encoding="utf-8",
    )
    result = resign_review_outputs(out, signer="tester", regenerate=True)
    assert result["artifact_status"]["review_state"] == "current_after_resign"
    assert result["artifact_status"]["stale_artifacts"] == []
    assert (out / "review_signature.json").is_file()
    assert (out / "annotation_audit_log.jsonl").is_file()
    assert (out / "report_draft.html").is_file()
