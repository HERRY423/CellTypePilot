"""Unit and integration tests for Novelty Verification & Adjudication Protocol."""

from __future__ import annotations

import json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import pytest

from celltypepilot.novelty_verification import (
    verify_novelty_candidate,
    log_novelty_adjudication,
    verify_subclustering_homogeneity,
    _check_qc_and_batch_confounding,
    _check_doublet_signature,
    _evaluate_state_vs_identity,
    VERIFICATION_SCHEMA,
    ADJUDICATION_SCHEMA,
)


@pytest.fixture
def synthetic_adata():
    rng = np.random.default_rng(42)
    # 60 cells, 10 genes
    X = rng.normal(size=(60, 10))
    # Cluster 0: normal
    # Cluster 1: high mito
    # Cluster 2: batch confounded
    obs = pd.DataFrame({
        "cluster": ["0"] * 20 + ["1"] * 20 + ["2"] * 20,
        "n_genes_by_counts": [500] * 60,
        "pct_counts_mt": [2.0] * 20 + [25.0] * 20 + [3.0] * 20,  # cluster 1 high mito
        "batch": ["b1", "b2"] * 10 + ["b1", "b2"] * 10 + ["b1"] * 20,  # cluster 2 100% b1
    })
    var = pd.DataFrame(index=[f"g{i}" for i in range(10)])
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def sample_atlas():
    return {
        "version": "mkg-2026.08.1",
        "tissues": {
            "general": {
                "cell_types": {
                    "T cell": {
                        "positive_markers": ["g0", "g1"],
                        "negative_markers": ["g8"],
                    },
                    "B cell": {
                        "positive_markers": ["g8", "g9"],
                        "negative_markers": ["g0"],
                    },
                }
            }
        },
    }


def test_gate1_qc_and_batch_confounding(synthetic_adata):
    # Cluster 0: should pass QC & batch
    g1_cl0 = _check_qc_and_batch_confounding(synthetic_adata, "cluster", "0")
    assert g1_cl0["qc_passed"] is True
    assert g1_cl0["batch_passed"] is True

    # Cluster 1: should fail QC (high mito)
    g1_cl1 = _check_qc_and_batch_confounding(synthetic_adata, "cluster", "1")
    assert g1_cl1["qc_passed"] is False
    assert any("HIGH_MITO_PERCENT" in f for f in g1_cl1["flags"])

    # Cluster 2: should fail batch check (100% b1)
    g1_cl2 = _check_qc_and_batch_confounding(synthetic_adata, "cluster", "2")
    assert g1_cl2["batch_passed"] is False
    assert any("CONFOUNDED_BATCH" in f for f in g1_cl2["flags"])


def test_gate2_subclustering_homogeneity(synthetic_adata):
    res = verify_subclustering_homogeneity(synthetic_adata, "cluster", "0", resolution=0.5)
    assert "is_homogeneous" in res
    assert "n_subclusters" in res


def test_gate3_state_vs_identity():
    # Stress/cycling gene overlap
    res = _evaluate_state_vs_identity(["CDK1", "MKI67", "TOP2A", "CCNB1"])
    assert isinstance(res, dict)
    assert "is_cell_state_driven" in res


def test_verify_novelty_candidate_full_packet(synthetic_adata, sample_atlas):
    critic_row = {
        "cluster": "0",
        "cell_type": "T cell",
        "evidence_score": 0.30,
        "pct_overlap": 0.20,
        "top_unmapped_markers": "MKI67;CDK1",
    }
    packet = verify_novelty_candidate(
        synthetic_adata, "cluster", "0", critic_row, sample_atlas, "general"
    )
    assert packet["schema_version"] == VERIFICATION_SCHEMA
    assert "gates" in packet
    assert "suggested_classification" in packet
    assert packet["adjudication_status"] == "pending_human_review"


def test_log_novelty_adjudication(tmp_path):
    entry = log_novelty_adjudication(
        tmp_path,
        cluster="0",
        verdict="validated_novel_cell_type",
        reviewer="Dr. Smith",
        notes="Validated via flow cytometry and FISH",
        pmid="PMID:12345678",
    )
    assert entry["schema_version"] == ADJUDICATION_SCHEMA
    assert entry["cluster"] == "0"
    assert entry["verdict"] == "validated_novel_cell_type"

    # Check append-only log file
    log_file = tmp_path / "novelty_adjudication_log.jsonl"
    assert log_file.exists()
    content = log_file.read_text(encoding="utf-8")
    assert "validated_novel_cell_type" in content

    # Check artifact_status.json stale marking
    status_file = tmp_path / "artifact_status.json"
    assert status_file.exists()
    status_data = json.loads(status_file.read_text(encoding="utf-8"))
    assert status_data["novelty_review_status"] == "adjudicated"
    assert status_data["derived_artifacts_stale"] is True


def test_invalid_verdict_raises(tmp_path):
    with pytest.raises(ValueError, match="Invalid verdict"):
        log_novelty_adjudication(tmp_path, "0", "invalid_verdict_foo", "Dr. Smith")
