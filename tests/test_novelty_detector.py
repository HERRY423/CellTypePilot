import anndata as ad
import numpy as np
import pandas as pd

from celltypepilot.novelty_detector import (
    KNOWN_SUPPORTED,
    NOVELTY_SCORE_SEMANTICS,
    OOD_NOVEL_CANDIDATE,
    REVIEW_ARTIFACT_OR_MIXED,
    attach_novelty_results,
    build_novelty_manifest,
    score_novelty_candidates,
)


def _de(*genes):
    return pd.DataFrame(
        {
            "gene": list(genes),
            "logfoldchange": [1.4 for _ in genes],
            "pval_adj": [0.001 for _ in genes],
        }
    )


def test_novelty_detector_surfaces_ood_candidates_without_renaming_identity():
    adata = ad.AnnData(
        X=np.ones((6, 8)),
        obs=pd.DataFrame(
            {
                "cluster": ["0", "0", "1", "1", "2", "2"],
                "batch": ["b1", "b1", "b2", "b2", "b3", "b3"],
            },
            index=[f"c{i}" for i in range(6)],
        ),
        var=pd.DataFrame(index=["CD3E", "CD4", "MS4A1", "X1", "X2", "X3", "Y1", "Y2"]),
    )
    markers = {
        "T cell": {"positive_markers": ["CD3E", "CD4"], "negative_markers": ["MS4A1"]},
        "B cell": {"positive_markers": ["MS4A1"], "negative_markers": ["CD3E"]},
    }
    critic = pd.DataFrame(
        {
            "cluster": ["0", "1", "2"],
            "cell_type": ["T cell", "Unknown", "Unknown"],
            "candidate_cell_type": ["T cell", "T cell", "B cell"],
            "combined_score": [0.82, 0.18, 0.22],
            "evidence_score": [0.82, 0.18, 0.22],
            "pct_overlap": [0.60, 0.10, 0.10],
            "neg_conflict": [0.0, 0.0, 0.0],
            "decision": ["accepted", "abstain", "abstain"],
            "critic_flags": ["PASS", "LOW_EVIDENCE", "POSSIBLE_DOUBLET; LOW_EVIDENCE"],
        }
    )
    novelty = score_novelty_candidates(
        adata,
        "cluster",
        critic,
        markers,
        de_results={
            "0": _de("CD3E", "CD4"),
            "1": _de("X1", "X2", "X3"),
            "2": _de("Y1", "Y2", "X1"),
        },
    )

    decisions = novelty.set_index("cluster")["novelty_decision"].to_dict()
    assert decisions["0"] == KNOWN_SUPPORTED
    assert decisions["1"] == OOD_NOVEL_CANDIDATE
    assert decisions["2"] == REVIEW_ARTIFACT_OR_MIXED
    assert set(novelty["novelty_score_semantics"]) == {NOVELTY_SCORE_SEMANTICS}
    assert "X1;X2;X3" in novelty.set_index("cluster").loc["1", "top_unmapped_markers"]
    assert "batch_enriched" in novelty.set_index("cluster").loc["1", "alternative_explanations"]

    enriched = attach_novelty_results(critic, novelty)
    assert enriched["cell_type"].tolist() == critic["cell_type"].tolist()
    assert enriched.loc[1, "novelty_decision"] == OOD_NOVEL_CANDIDATE


def test_novelty_manifest_is_explicitly_review_only():
    novelty = pd.DataFrame(
        {
            "cluster": ["1"],
            "novelty_decision": [OOD_NOVEL_CANDIDATE],
        }
    )
    manifest = build_novelty_manifest(novelty)
    assert manifest["schema_version"] == "celltypepilot.novelty-ood.v1"
    assert manifest["identity_invariant"] is True
    assert manifest["decision_counts"][OOD_NOVEL_CANDIDATE] == 1
    assert "not a validated new cell-type discovery" in manifest["claim_boundary"]
