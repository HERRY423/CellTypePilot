"""Contracts for backend-neutral hierarchical selective identity decisions."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from celltypepilot.candidate_backends import (
    CandidateContractError,
    aggregate_cell_candidates,
    marker_scores_as_evidence,
    normalize_candidate_table,
)
from celltypepilot.hierarchical_selector import (
    SelectiveDecisionError,
    attach_marker_evidence,
    enforce_selective_decisions,
    select_hierarchical_identities,
    validate_policy,
)
from celltypepilot.identity_contract import build_identity_resolver
from celltypepilot.validation_domains import assess_validation_domain, load_validation_domains


@pytest.fixture
def blood_resolver(blood_atlas):
    return build_identity_resolver(blood_atlas, ["blood", "general"])


def _candidates(rows: list[dict], resolver: dict) -> pd.DataFrame:
    return normalize_candidate_table(pd.DataFrame(rows), resolver)


def test_two_independent_backends_accept_same_leaf(blood_resolver):
    candidates = _candidates(
        [
            {
                "cluster": "0",
                "backend": "celltypist",
                "cell_type": "CD4+ T cell",
                "score": 0.91,
            },
            {
                "cluster": "0",
                "backend": "popv",
                "cell_type": "CD4+ T cell",
                "score": 0.83,
            },
        ],
        blood_resolver,
    )

    row = select_hierarchical_identities(candidates, blood_resolver, ["0"]).iloc[0]

    assert row["selective_decision"] == "accepted_leaf"
    assert row["selected_cell_type"] == "CD4+ T cell"
    assert row["selected_cl_id"] == "CL:0000624"
    assert row["n_independent_backends"] == 2
    assert row["backend_agreement_fraction"] == 1.0
    assert pd.isna(row["calibrated_probability"])


def test_sibling_disagreement_collapses_only_to_governed_parent(blood_resolver):
    candidates = _candidates(
        [
            {
                "cluster": "0",
                "backend": "celltypist",
                "cell_type": "Naive CD4+ T cell",
            },
            {
                "cluster": "0",
                "backend": "singler",
                "cell_type": "Memory CD4+ T cell",
            },
        ],
        blood_resolver,
    )

    row = select_hierarchical_identities(candidates, blood_resolver, ["0"]).iloc[0]

    assert row["selective_decision"] == "accepted_ancestor"
    assert row["selected_cell_type"] == "CD4+ T cell"
    assert row["ancestor_hops"] == 1
    assert row["backend_agreement_fraction"] == 0.5


def test_cross_lineage_disagreement_abstains_with_candidate_set(blood_resolver):
    candidates = _candidates(
        [
            {"cluster": "0", "backend": "celltypist", "cell_type": "B cell"},
            {"cluster": "0", "backend": "scanvi", "cell_type": "T cell"},
        ],
        blood_resolver,
    )

    row = select_hierarchical_identities(candidates, blood_resolver, ["0"]).iloc[0]

    assert row["selective_decision"] == "abstain"
    assert row["selected_cell_type"] == "Unknown"
    assert row["selective_abstain_reason"] == "BACKEND_DISAGREEMENT_NO_SAFE_ANCESTOR"
    assert set(row["candidate_set"].split(";")) == {"B cell", "T cell"}


def test_marker_and_llm_are_candidates_but_not_default_votes(blood_resolver):
    marker_scores = pd.DataFrame(
        [
            {
                "cluster": "0",
                "cell_type": "B cell",
                "cl_id": "CL:0000236",
                "combined_score": 0.99,
                "rank": 1,
            }
        ]
    )
    marker = marker_scores_as_evidence(marker_scores, blood_resolver)
    llm = _candidates(
        [{"cluster": "0", "backend": "llm", "cell_type": "B cell", "score": 0.99}],
        blood_resolver,
    )
    candidates = pd.concat([marker, llm], ignore_index=True)

    row = select_hierarchical_identities(candidates, blood_resolver, ["0"]).iloc[0]

    assert set(candidates["decision_role"]) == {"evidence_only", "hypothesis_only"}
    assert row["selective_decision"] == "abstain"
    assert row["n_independent_backends"] == 0
    assert row["selective_candidate_cell_type"] == "B cell"


def test_same_family_adapters_cannot_manufacture_independent_consensus(blood_resolver):
    candidates = _candidates(
        [
            {"cluster": "0", "backend": "knn", "cell_type": "B cell"},
            {"cluster": "0", "backend": "correlation", "cell_type": "B cell"},
        ],
        blood_resolver,
    )

    row = select_hierarchical_identities(candidates, blood_resolver, ["0"]).iloc[0]

    assert row["n_independent_backends"] == 1
    assert row["selective_decision"] == "abstain"
    assert row["selective_abstain_reason"] == "INSUFFICIENT_INDEPENDENT_BACKENDS"


def test_cell_level_backend_output_aggregates_to_auditable_cluster_votes():
    frame = pd.DataFrame(
        [
            {"cell_id": "a", "backend": "popv", "predicted_label": "T cell"},
            {"cell_id": "b", "backend": "popv", "predicted_label": "T cell"},
            {"cell_id": "c", "backend": "popv", "predicted_label": "B cell"},
        ]
    )
    assignments = pd.Series(["0", "0", "0"], index=["a", "b", "c"])

    aggregated = aggregate_cell_candidates(frame, assignments).sort_values("rank")

    assert aggregated.iloc[0]["predicted_label"] == "T cell"
    assert aggregated.iloc[0]["score"] == pytest.approx(2 / 3)
    assert aggregated.iloc[0]["score_semantics"] == (
        "within_backend_cluster_top1_vote_fraction_not_probability"
    )
    assert aggregated.iloc[1]["rank"] == 2


def test_unknown_backend_and_fake_calibration_fail_closed(blood_resolver):
    with pytest.raises(CandidateContractError, match="Unsupported candidate backend"):
        _candidates(
            [{"cluster": "0", "backend": "mystery", "cell_type": "B cell"}],
            blood_resolver,
        )
    with pytest.raises(SelectiveDecisionError, match="cannot claim calibration"):
        validate_policy({"calibrated": True})


def test_selector_abstention_cannot_be_rescued_by_marker_critic(blood_resolver):
    candidates = _candidates(
        [{"cluster": "0", "backend": "celltypist", "cell_type": "B cell"}],
        blood_resolver,
    )
    decisions = select_hierarchical_identities(candidates, blood_resolver, ["0"])
    marker_scores = pd.DataFrame(
        [
            {
                "cluster": "0",
                "cell_type": "B cell",
                "cl_id": "CL:0000236",
                "combined_score": 0.99,
                "pct_overlap": 0.99,
                "neg_conflict": 0.0,
                "specificity": 0.99,
            }
        ]
    )
    critic_input = attach_marker_evidence(decisions, marker_scores, blood_resolver)
    critic_input["decision"] = "accepted"
    critic_input["candidate_cell_type"] = "B cell"
    critic_input["candidate_cl_id"] = "CL:0000236"
    critic_input["critic_flags"] = "PASS"
    critic_input["abstain_reason"] = ""

    enforced = enforce_selective_decisions(critic_input, decisions).iloc[0]

    assert enforced["cell_type"] == "Unknown"
    assert enforced["decision"] == "abstain"
    assert enforced["candidate_cell_type"] == "B cell"
    assert "INSUFFICIENT_INDEPENDENT_BACKENDS" in enforced["abstain_reason"]


def test_registry_defines_exactly_three_fail_closed_depth_domains():
    registry = load_validation_domains()

    assert set(registry["domains"]) == {"lung", "gut_ibd", "tumor_microenvironment"}
    assert all(domain["claim_ready"] is False for domain in registry["domains"].values())
    assert assess_validation_domain("lung")["domain_id"] == "lung"
    assert assess_validation_domain("gut")["domain_id"] == "gut_ibd"
    assert assess_validation_domain("brain")["status"] == "out_of_focus_exploratory"


@pytest.mark.slow
@pytest.mark.integration
def test_pipeline_writes_backend_and_hierarchical_decision_artifacts(tmp_path, monkeypatch):
    from celltypepilot.orchestrator import run_annotation_pipeline

    def fake_marker_scores(*_args, **_kwargs):
        rows = []
        for cluster, first, second in (
            ("t", ("T cell", "CL:0000084", 0.9), ("B cell", "CL:0000236", 0.1)),
            ("b", ("B cell", "CL:0000236", 0.9), ("T cell", "CL:0000084", 0.1)),
        ):
            for rank, (label, cl_id, score) in enumerate((first, second), start=1):
                rows.append(
                    {
                        "cluster": cluster,
                        "cell_type": label,
                        "cl_id": cl_id,
                        "rank": rank,
                        "combined_score": score,
                        "pct_overlap": score,
                        "mean_log2fc": score,
                        "pct_expressed": score,
                        "specificity": score,
                        "neg_conflict": 0.0,
                        "pos_supporting_markers": "CD3D;CD3E"
                        if label == "T cell"
                        else "CD19;MS4A1",
                        "context_only_support": False,
                        "context_review_status": "not_applicable",
                    }
                )
        return pd.DataFrame(rows)

    monkeypatch.setattr("celltypepilot.marker_scorer.compute_marker_scores", fake_marker_scores)

    genes = [
        "CD3D",
        "CD3E",
        "CD2",
        "TRAC",
        "IL7R",
        "CD19",
        "MS4A1",
        "CD79A",
        "CD79B",
        "PAX5",
        "BACKGROUND",
    ]
    rng = np.random.RandomState(7)
    matrix = rng.exponential(0.05, size=(60, len(genes)))
    matrix[:30, :5] += rng.exponential(3.0, size=(30, 5))
    matrix[30:, 5:10] += rng.exponential(3.0, size=(30, 5))
    adata = ad.AnnData(np.log1p(matrix).astype(np.float32))
    adata.var_names = genes
    adata.obs["cluster"] = pd.Categorical(["t"] * 30 + ["b"] * 30)
    input_path = tmp_path / "query.h5ad"
    adata.write(input_path)

    candidate_path = tmp_path / "candidates.csv"
    pd.DataFrame(
        [
            {"cluster": "t", "backend": "celltypist", "cell_type": "T cell"},
            {"cluster": "t", "backend": "popv", "cell_type": "T cell"},
            {"cluster": "b", "backend": "celltypist", "cell_type": "B cell"},
            {"cluster": "b", "backend": "popv", "cell_type": "B cell"},
        ]
    ).to_csv(candidate_path, index=False)

    result = run_annotation_pipeline(
        input_path,
        "cluster",
        tmp_path / "output",
        species="human",
        tissue="blood",
        no_figures=True,
        enable_states=False,
        candidate_artifact_paths=[candidate_path],
        progress=lambda step, total, message: print(step, total, message, flush=True),
    )

    assert set(result["critic_results"]["selective_decision"]) == {"accepted_leaf"}
    assert result["manifest"]["parameters"]["hierarchical_selector"]["marker_role"] == (
        "evidence_only_not_candidate_vote"
    )
    assert result["manifest"]["parameters"]["validation_domain"]["status"] == (
        "out_of_focus_exploratory"
    )
    assert result["paths"]["backend_candidates"].exists()
    assert result["paths"]["hierarchical_decisions"].exists()
    written = ad.read_h5ad(result["paths"]["annotated"])
    assert set(written.obs["ctp_selective_decision"]) == {"accepted_leaf"}
    assert set(written.obs["ctp_independent_backend_count"]) == {2}
