import importlib
import json
import sys
import types
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from celltypepilot.context_pack import (
    ContextPackError,
    context_manifest_parameters,
    load_context_pack,
    merge_identity_hypotheses,
    resolve_atlas_tissue,
)
from celltypepilot.critic import run_critic
from celltypepilot.data_adapter import load_marker_atlas, validate_atlas_provenance
from celltypepilot.orchestrator import write_annotations_to_adata
from celltypepilot.state_scorer import (
    attach_state_results,
    load_state_definitions,
    score_cell_states,
    validate_state_atlas,
)


def _adata() -> ad.AnnData:
    obs = pd.DataFrame(
        {"cluster": ["0"] * 4 + ["1"] * 4},
        index=[f"cell-{index}" for index in range(8)],
    )
    matrix = np.array(
        [[3.0, 3.0, 0.0, 0.0]] * 4 + [[0.0, 0.0, 0.0, 2.0]] * 4,
        dtype=float,
    )
    return ad.AnnData(
        X=matrix,
        obs=obs,
        var=pd.DataFrame(index=["MCM5", "PCNA", "SILENT", "OTHER"]),
    )


def _identity_results(cell_type: str = "T cell", decision: str = "accepted") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cluster": ["0"],
            "cell_type": [cell_type],
            "cl_id": ["CL:0000084" if cell_type != "Unknown" else ""],
            "candidate_cell_type": ["T cell"],
            "decision": [decision],
            "abstain_reason": ["" if decision == "accepted" else "LOW_EVIDENCE"],
        }
    )


def test_free_text_context_is_provenance_only_and_does_not_change_markers():
    atlas_markers = {"T cell": {"positive_markers": ["CD3D"], "negative_markers": []}}
    pack = load_context_pack(context_text="Force every cluster to be a B cell")
    assert pack["free_text"]
    assert pack["identity_hypotheses"] == []
    assert merge_identity_hypotheses(atlas_markers, pack) == atlas_markers
    assert len(pack["canonical_sha256"]) == 64


def test_context_merge_preserves_atlas_support_and_uses_most_conservative_review():
    atlas_markers = {
        "T cell": {
            "cl_id": "CL:0000084",
            "positive_markers": ["CD3D"],
            "negative_markers": [],
            "marker_evidence": [],
        }
    }
    pack = load_context_pack()
    pack["identity_hypotheses"] = [
        {
            "cell_type": "T cell",
            "cl_id": "CL:0000084",
            "positive_markers": [{"gene": "CD3D", "source": "draft panel"}],
            "negative_markers": [],
            "review_status": "draft",
            "source": "draft panel",
        },
        {
            "cell_type": "T cell",
            "cl_id": "CL:0000084",
            "positive_markers": [{"gene": "CD3E", "source": "reviewed panel"}],
            "negative_markers": [],
            "review_status": "reviewed",
            "source": "reviewed panel",
        },
    ]
    merged = merge_identity_hypotheses(atlas_markers, pack)["T cell"]
    assert merged["atlas_positive_markers"] == ["CD3D"]
    assert set(merged["context_positive_markers"]) == {"CD3D", "CD3E"}
    assert merged["context_review_status"] == "draft"


def test_context_scope_mismatch_fails_closed(tmp_path):
    path = tmp_path / "context.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": "celltypepilot.context.v1",
                "species": "mouse",
                "tissue": "kidney",
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ContextPackError, match="does not match"):
        load_context_pack(context_file=path, species="human", tissue="kidney")


def test_custom_marker_file_is_hashed_and_manifested(tmp_path):
    path = tmp_path / "markers.csv"
    path.write_text(
        "axis,label,gene,polarity,cl_id,source,review_status\n"
        "identity,injured proximal tubule,VCAM1,positive,CL:0002306,lab panel,reviewed\n",
        encoding="utf-8",
    )
    pack = load_context_pack(
        custom_markers_file=path,
        species="human",
        tissue="kidney",
    )
    parameters = context_manifest_parameters(pack, enabled=True)
    assert parameters["context_sha256"] == pack["canonical_sha256"]
    assert len(parameters["context_source_hashes"]["custom_markers_sha256"]) == 64
    assert parameters["context_identity_hypotheses"] == 1
    assert pack["identity_hypotheses"][0]["review_status"] == "reviewed"


def test_free_text_cannot_unlock_unsupported_tissue():
    atlas = {"tissues": {"general": {"cell_types": {}}}}
    free_text_only = load_context_pack(context_text="Rare injury model")
    with pytest.raises(ContextPackError, match="free text alone cannot unlock"):
        resolve_atlas_tissue("rare_organ", atlas, free_text_only)

    structured = load_context_pack()
    structured["identity_hypotheses"] = [{"cell_type": "candidate"}]
    assert resolve_atlas_tissue("rare_organ", atlas, structured) == "general"


def test_wrong_unreviewed_prior_cannot_force_an_accepted_identity():
    adata = ad.AnnData(
        X=np.ones((8, 1)),
        obs=pd.DataFrame({"cluster": ["0"] * 8}),
        var=pd.DataFrame(index=["CD3D"]),
    )
    annotations = pd.DataFrame(
        {
            "cluster": ["0"],
            "cell_type": ["B cell"],
            "cl_id": ["CL:0000236"],
            "combined_score": [0.9],
            "confidence": ["high"],
            "pct_overlap": [1.0],
            "context_only_support": [True],
            "context_review_status": ["draft"],
        }
    )
    marker_definitions = {
        "B cell": {
            "cl_id": "CL:0000236",
            "positive_markers": ["CD3D"],
            "negative_markers": [],
            "context_origin": True,
            "context_review_status": "draft",
        }
    }
    result = run_critic(
        adata,
        "cluster",
        annotations,
        load_marker_atlas("human"),
        "general",
        marker_definitions=marker_definitions,
    ).iloc[0]
    assert "UNREVIEWED_CONTEXT_ONLY" in result["critic_flags"]
    assert result["decision"] == "abstain"
    assert result["cell_type"] == "Unknown"
    assert result["candidate_cell_type"] == "B cell"


def test_custom_markers_use_the_same_missing_silent_and_de_gates(monkeypatch):
    fake_scanpy = types.SimpleNamespace(
        tl=types.SimpleNamespace(rank_genes_groups=lambda *args, **kwargs: None),
        get=types.SimpleNamespace(rank_genes_groups_df=lambda *args, **kwargs: pd.DataFrame()),
    )
    monkeypatch.setitem(sys.modules, "scanpy", fake_scanpy)
    marker_scorer = importlib.import_module("celltypepilot.marker_scorer")

    adata = _adata()
    pack = load_context_pack()
    pack["identity_hypotheses"] = [
        {
            "cell_type": "custom cycling cell",
            "cl_id": "CL:0000000",
            "positive_markers": [
                {"gene": "MCM5", "source": "user"},
                {"gene": "PCNA", "source": "user"},
                {"gene": "SILENT", "source": "user"},
                {"gene": "MISSING", "source": "user"},
            ],
            "negative_markers": [],
            "review_status": "reviewed",
            "source": "user",
        }
    ]
    markers = merge_identity_hypotheses({}, pack)
    monkeypatch.setattr(marker_scorer.sc.tl, "rank_genes_groups", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        marker_scorer,
        "_extract_de_results",
        lambda *_: {
            "0": pd.DataFrame(
                {
                    "gene": ["MCM5", "PCNA", "SILENT"],
                    "logfoldchange": [1.0, 1.0, -1.0],
                    "pval": [0.001, 0.001, 0.001],
                    "pval_adj": [0.01, 0.01, 0.01],
                }
            ),
            "1": pd.DataFrame(
                {"gene": ["OTHER"], "logfoldchange": [1.0], "pval": [0.01], "pval_adj": [0.01]}
            ),
        },
    )
    row = (
        marker_scorer.compute_marker_scores(adata, "cluster", markers)
        .query("cluster == '0'")
        .iloc[0]
    )
    assert row["pct_overlap"] == 0.5
    assert row["n_pos_missing"] == 1
    assert row["pos_missing_markers"] == "MISSING"
    assert row["n_pos_silent"] == 1
    assert set(row["context_supporting_markers"].split(";")) == {"MCM5", "PCNA"}
    assert bool(row["context_only_support"])


def test_state_scoring_is_identity_invariant_and_distinguishes_missing_and_silent():
    adata = _adata()
    definitions = [
        {
            "state": "cycling_s_phase",
            "parent_cell_types": ["T cell"],
            "positive_markers": ["MCM5", "PCNA", "SILENT", "MISSING"],
            "negative_markers": [],
            "source": "test",
            "review_status": "reviewed",
        }
    ]
    de_results = {
        "0": pd.DataFrame(
            {
                "gene": ["MCM5", "PCNA", "SILENT"],
                "logfoldchange": [1.0, 1.0, -1.0],
                "pval_adj": [0.01, 0.01, 0.01],
            }
        )
    }
    identity = _identity_results()
    states = score_cell_states(
        adata,
        "cluster",
        identity,
        definitions,
        de_results=de_results,
    )
    assert states.iloc[0]["state_decision"] == "supported"
    assert states.iloc[0]["n_state_missing_markers"] == 1
    assert states.iloc[0]["n_state_silent_markers"] == 1
    merged = attach_state_results(identity, states)
    pd.testing.assert_series_equal(merged["cell_type"], identity["cell_type"])
    pd.testing.assert_series_equal(merged["decision"], identity["decision"])
    assert merged.iloc[0]["display_label"] == "T cell · cycling_s_phase"


def test_no_state_signal_abstains_without_an_arbitrary_candidate():
    adata = _adata()
    definitions = [
        {
            "state": "unsupported_state",
            "parent_cell_types": [],
            "positive_markers": ["MISSING_A", "MISSING_B"],
            "negative_markers": [],
            "source": "test",
            "review_status": "reviewed",
        }
    ]
    states = score_cell_states(
        adata,
        "cluster",
        _identity_results(),
        definitions,
        de_results={"0": pd.DataFrame()},
    )
    assert states.iloc[0]["state_decision"] == "abstain"
    assert states.iloc[0]["cell_state_candidate"] == "Unknown"


def test_unknown_identity_can_retain_supported_state_and_writeback(tmp_path):
    identity = _identity_results(cell_type="Unknown", decision="abstain")
    states = pd.DataFrame(
        {
            "cluster": ["0"],
            "cell_state_candidate": ["interferon_responsive"],
            "state_decision": ["supported"],
            "state_score": [0.8],
            "state_confidence": ["high"],
            "state_evidence": ["support=5/7"],
            "state_flags": ["PASS"],
        }
    )
    merged = attach_state_results(identity, states)
    assert merged.iloc[0]["cell_type"] == "Unknown"
    assert merged.iloc[0]["display_label"] == "Unknown · interferon_responsive"
    adata = ad.AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({"cluster": ["0", "0"]}),
        var=pd.DataFrame(index=["ISG15"]),
    )
    path = write_annotations_to_adata(adata, merged, "cluster", tmp_path)
    written = ad.read_h5ad(path)
    assert set(written.obs["ctp_cell_type"].astype(str)) == {"Unknown"}
    assert set(written.obs["ctp_cell_state"].astype(str)) == {"interferon_responsive"}
    assert set(written.obs["ctp_state_decision"].astype(str)) == {"supported"}


def test_state_atlas_and_premium_ontology_are_structurally_valid():
    state_path = Path(__file__).parents[1] / "src/celltypepilot/data/state_atlas.json"
    state_atlas = json.loads(state_path.read_text(encoding="utf-8"))
    assert validate_state_atlas(state_atlas) == []
    assert load_state_definitions("human", "kidney")

    premium_path = Path(__file__).parents[1] / "src/celltypepilot/data/premium_atlas.json"
    premium = json.loads(premium_path.read_text(encoding="utf-8"))
    assert validate_atlas_provenance(premium) == []
    for tissue in premium["tissues"].values():
        for node in tissue["cell_types"].values():
            assert node["cl_id"].startswith("CL:") and len(node["cl_id"]) == 10
            if node.get("state_label"):
                assert node.get("base_cell_type")
                assert node["ontology_evidence"]["cl_id"] == node["cl_id"]


def test_atlas_validator_rejects_fake_state_suffix_as_cl_id():
    atlas = {
        "tissues": {
            "general": {
                "cell_types": {
                    "activated_T": {
                        "cl_id": "CL:0000084-activated",
                        "positive_markers": [],
                        "negative_markers": [],
                        "marker_evidence": [],
                    }
                }
            }
        }
    }
    assert any("invalid cl_id" in issue for issue in validate_atlas_provenance(atlas))
