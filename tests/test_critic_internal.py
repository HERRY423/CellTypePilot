"""Tests for critic internal functions — evidence, negative markers, doublet, ontology."""

import anndata as ad
import numpy as np
import pandas as pd

from celltypepilot.constants import (
    CONFIDENCE_HIGH,
    CONFIDENCE_LOW,
    CONFIDENCE_MEDIUM,
    CONFIDENCE_REVIEW,
)
from celltypepilot.critic import (
    _check_doublet_signal,
    _check_ensemble_agreement,
    _check_evidence_sufficiency,
    _check_negative_markers,
    _check_ontology_consistency,
    _recalibrate_confidence,
    format_evidence_summary,
    format_run_narrative,
    generate_critic_summary,
    run_critic,
)


def _make_expression_adata(high_genes: list[str], n_cells: int = 60) -> ad.AnnData:
    """One-cluster AnnData where ``high_genes`` are expressed in all cells."""
    X = np.zeros((n_cells, len(high_genes)), dtype=np.float32)
    X[:, :] = 3.0
    adata = ad.AnnData(X=X)
    adata.var_names = list(high_genes)
    adata.obs["leiden"] = "0"
    return adata


# Blood atlas marker sets used by the calibrated doublet tests
_T_MARKERS = ["CD3D", "CD3E", "CD2", "TRAC", "CD7"]
_CD4_MARKERS = ["CD4", "IL7R", "MAL", "TRAC"]
_B_MARKERS = ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5"]


class TestCheckEvidenceSufficiency:
    def test_no_markers(self, synthetic_pbmc):
        result = _check_evidence_sufficiency(synthetic_pbmc, "0", "leiden", [])
        assert result["flag"] == "NO_MARKERS"

    def test_good_coverage(self, synthetic_pbmc):
        # Use genes that exist in the synthetic data
        markers = ["CD3D", "CD3E", "CD2"]
        result = _check_evidence_sufficiency(
            synthetic_pbmc, synthetic_pbmc.obs["leiden"].iloc[0], "leiden", markers
        )
        assert result["flag"] in ("", "PARTIAL_EVIDENCE", "LOW_EVIDENCE")
        assert "Coverage:" in result["evidence"]

    def test_nonexistent_genes(self, synthetic_pbmc):
        markers = ["FAKE1", "FAKE2", "FAKE3"]
        cluster = synthetic_pbmc.obs["leiden"].iloc[0]
        result = _check_evidence_sufficiency(synthetic_pbmc, cluster, "leiden", markers)
        # All genes missing → 0 coverage → LOW_EVIDENCE
        assert result["flag"] in ("LOW_EVIDENCE", "NO_MARKERS")


class TestCheckNegativeMarkers:
    def test_no_negative_markers(self, synthetic_pbmc):
        cluster = synthetic_pbmc.obs["leiden"].iloc[0]
        result = _check_negative_markers(synthetic_pbmc, cluster, "leiden", [])
        assert result["flag"] == ""

    def test_with_negative_markers(self, synthetic_pbmc):
        cluster = synthetic_pbmc.obs["leiden"].iloc[0]
        result = _check_negative_markers(synthetic_pbmc, cluster, "leiden", ["GENE_0", "GENE_1"])
        assert isinstance(result["flag"], str)
        assert "evidence" in result


class TestCheckDoubletSignal:
    def test_no_doublet(self, synthetic_pbmc):
        from celltypepilot.data_adapter import get_all_markers_for_tissue, load_marker_atlas

        atlas = load_marker_atlas("human")
        markers = get_all_markers_for_tissue(atlas, "blood")
        cluster = synthetic_pbmc.obs["leiden"].iloc[0]
        result = _check_doublet_signal(synthetic_pbmc, cluster, "leiden", markers)
        assert isinstance(result["flag"], str)
        assert "evidence" in result

    def test_empty_markers(self, synthetic_pbmc):
        cluster = synthetic_pbmc.obs["leiden"].iloc[0]
        result = _check_doublet_signal(synthetic_pbmc, cluster, "leiden", {})
        assert result["flag"] == ""

    def test_same_lineage_coexpression_not_flagged(self, blood_markers):
        """T cell + CD4+ T cell co-expression is subtype refinement, not a doublet."""
        from celltypepilot.data_adapter import build_lineage_groups, load_marker_atlas

        atlas = load_marker_atlas("human")
        lineages = build_lineage_groups(atlas, "blood")
        adata = _make_expression_adata(_T_MARKERS + _CD4_MARKERS)
        result = _check_doublet_signal(
            adata, "0", "leiden", blood_markers, lineages
        )
        assert result["flag"] == ""

    def test_cross_lineage_coexpression_flagged(self, blood_markers):
        """T + B lineage markers co-expressed → genuine doublet signal."""
        from celltypepilot.data_adapter import build_lineage_groups, load_marker_atlas

        atlas = load_marker_atlas("human")
        lineages = build_lineage_groups(atlas, "blood")
        adata = _make_expression_adata(_T_MARKERS + _B_MARKERS)
        result = _check_doublet_signal(
            adata, "0", "leiden", blood_markers, lineages
        )
        assert result["flag"] == "POSSIBLE_DOUBLET"
        assert "Cross-lineage" in result["evidence"]

    def test_cross_lineage_without_lineage_map_still_flagged(self, blood_markers):
        """Backward compat: no lineage map falls back to pairwise comparison."""
        adata = _make_expression_adata(_T_MARKERS + _B_MARKERS)
        result = _check_doublet_signal(adata, "0", "leiden", blood_markers)
        assert result["flag"] == "POSSIBLE_DOUBLET"

    def test_shared_gene_only_expression_not_flagged(self):
        """Shared genes (e.g. GNLY in NK + CD8 panels) are weak evidence.

        A cluster expressing only the shared gene of a second panel must not
        count that panel as an active lineage signature.
        """
        markers = {
            "NK cell": {
                "positive_markers": ["NCAM1", "NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"]
            },
            "CD8+ T cell": {
                "positive_markers": ["CD8A", "CD8B", "GZMB", "PRF1", "GNLY"]
            },
        }
        # NK program fully on; from the CD8 panel only the shared genes are on
        adata = _make_expression_adata(["NCAM1", "NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"])
        result = _check_doublet_signal(adata, "0", "leiden", markers)
        assert result["flag"] == ""

    def test_specific_genes_override_shared_weighting(self):
        """When panel-specific genes are also expressed the doublet still fires."""
        markers = {
            "NK cell": {
                "positive_markers": ["NCAM1", "NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"]
            },
            "CD8+ T cell": {
                "positive_markers": ["CD8A", "CD8B", "GZMB", "PRF1", "GNLY"]
            },
        }
        adata = _make_expression_adata(
            ["NCAM1", "NKG7", "GNLY", "KLRD1", "PRF1", "GZMB", "CD8A", "CD8B"]
        )
        result = _check_doublet_signal(adata, "0", "leiden", markers)
        assert result["flag"] == "POSSIBLE_DOUBLET"


class TestBuildLineageGroups:
    def test_subtype_maps_to_parent(self):
        from celltypepilot.data_adapter import build_lineage_groups, load_marker_atlas

        atlas = load_marker_atlas("human")
        groups = build_lineage_groups(atlas, "blood")
        assert groups["T cell"] == "T cell"
        assert groups["CD4+ T cell"] == "T cell"
        assert groups["CD8+ T cell"] == "T cell"
        assert groups["B cell"] == "B cell"

    def test_unknown_tissue_falls_back_to_general(self):
        from celltypepilot.data_adapter import build_lineage_groups, load_marker_atlas

        atlas = load_marker_atlas("human")
        groups = build_lineage_groups(atlas, "nonexistent_tissue_xyz")
        # Unknown tissues fall back to the general atlas, matching
        # get_all_markers_for_tissue behaviour.
        assert groups == build_lineage_groups(atlas, "general")
        assert "Macrophage" in groups

    def test_empty_atlas(self):
        from celltypepilot.data_adapter import build_lineage_groups

        assert build_lineage_groups({}, "blood") == {}


class TestCheckOntologyConsistency:
    def test_valid_cl_id(self):
        result = _check_ontology_consistency("T cell", "CL:0000084")
        assert result["flag"] == ""

    def test_no_cl_id(self):
        result = _check_ontology_consistency("T cell", "")
        assert result["flag"] == "NO_CL_ID"

    def test_invalid_format(self):
        result = _check_ontology_consistency("T cell", "INVALID:123")
        assert result["flag"] == "INVALID_CL_FORMAT"

    def test_none_cl_id(self):
        result = _check_ontology_consistency("T cell", None)
        assert result["flag"] == "NO_CL_ID"


class TestRecalibrateConfidence:
    def test_no_flags_preserves(self):
        result = _recalibrate_confidence(CONFIDENCE_HIGH, [], {"flag": ""}, {"flag": ""})
        assert result == CONFIDENCE_HIGH

    def test_low_evidence_downgrades(self):
        result = _recalibrate_confidence(
            CONFIDENCE_HIGH, ["LOW_EVIDENCE"], {"flag": "LOW_EVIDENCE"}, {"flag": ""}
        )
        assert result == CONFIDENCE_REVIEW

    def test_neg_marker_conflict_downgrades(self):
        result = _recalibrate_confidence(
            CONFIDENCE_HIGH, ["NEG_MARKER_CONFLICT"], {"flag": ""}, {"flag": "NEG_MARKER_CONFLICT"}
        )
        assert result == CONFIDENCE_REVIEW

    def test_partial_evidence_mild_downgrade(self):
        result = _recalibrate_confidence(
            CONFIDENCE_HIGH, ["PARTIAL_EVIDENCE"], {"flag": "PARTIAL_EVIDENCE"}, {"flag": ""}
        )
        assert result == CONFIDENCE_LOW

    def test_no_markers_severe(self):
        result = _recalibrate_confidence(
            CONFIDENCE_HIGH, ["NO_MARKERS"], {"flag": "NO_MARKERS"}, {"flag": ""}
        )
        assert result == CONFIDENCE_REVIEW

    def test_unknown_confidence_input(self):
        result = _recalibrate_confidence("unknown_conf", [], {"flag": ""}, {"flag": ""})
        assert result in (CONFIDENCE_HIGH, CONFIDENCE_MEDIUM, CONFIDENCE_LOW, CONFIDENCE_REVIEW)


class TestCheckEnsembleAgreement:
    def test_no_ensemble_data(self):
        result = _check_ensemble_agreement("0", "T cell", {})
        assert result["flag"] == ""

    def test_agreement(self):
        ensemble = {
            "0": {"marker_score": 0.8, "ref_score": 0.7, "agreement": True, "source": "both"}
        }
        result = _check_ensemble_agreement("0", "T cell", ensemble)
        assert result["flag"] == ""

    def test_strong_disagreement(self):
        ensemble = {
            "0": {"marker_score": 0.9, "ref_score": 0.1, "agreement": False, "source": "both"}
        }
        result = _check_ensemble_agreement("0", "T cell", ensemble)
        assert result["flag"] == "ENSEMBLE_DISAGREEMENT"

    def test_mild_disagreement(self):
        ensemble = {
            "0": {"marker_score": 0.5, "ref_score": 0.35, "agreement": False, "source": "both"}
        }
        result = _check_ensemble_agreement("0", "T cell", ensemble)
        assert result["flag"] == "ENSEMBLE_MILD_DISAGREEMENT"

    def test_weak_reference_only(self):
        ensemble = {
            "0": {"marker_score": 0.0, "ref_score": 0.3, "agreement": True, "source": "reference"}
        }
        result = _check_ensemble_agreement("0", "T cell", ensemble)
        assert result["flag"] == "WEAK_REFERENCE_ONLY"


class TestGenerateCriticSummary:
    def test_all_pass(self):
        results = pd.DataFrame(
            {
                "cluster": ["0", "1", "2"],
                "critic_flags": ["PASS", "PASS", "PASS"],
                "critic_confidence": ["high", "medium", "low"],
            }
        )
        summary = generate_critic_summary(results)
        assert summary["total_clusters"] == 3
        assert summary["pass"] == 3
        assert summary["flagged"] == 0
        assert summary["confidence_distribution"]["high"] == 1

    def test_with_flags(self):
        results = pd.DataFrame(
            {
                "cluster": ["0", "1"],
                "critic_flags": ["PASS", "LOW_EVIDENCE; NEG_MARKER_CONFLICT"],
                "critic_confidence": ["high", "needs_review"],
            }
        )
        summary = generate_critic_summary(results)
        assert summary["flagged"] == 1
        assert "LOW_EVIDENCE" in summary["flag_types"]
        assert "NEG_MARKER_CONFLICT" in summary["flag_types"]
        assert "1" in summary["clusters_needing_review"]

    def test_narrative_present(self):
        results = pd.DataFrame(
            {
                "cluster": ["0", "1"],
                "critic_flags": ["PASS", "POSSIBLE_DOUBLET"],
                "critic_confidence": ["high", "needs_review"],
            }
        )
        summary = generate_critic_summary(results)
        assert "2 cluster(s) reviewed" in summary["narrative"]
        assert "POSSIBLE_DOUBLET" in summary["narrative"]

    def test_narrative_empty_run(self):
        assert format_run_narrative({"total_clusters": 0}) == "No clusters were reviewed."


class TestEvidenceSummary:
    def test_pass_summary(self):
        row = pd.Series(
            {
                "cluster": "3",
                "cell_type": "T cell",
                "critic_flags": "PASS",
                "critic_confidence": "high",
                "pct_overlap": 0.8,
                "combined_score": 0.75,
            }
        )
        summary = format_evidence_summary(row)
        assert "Cluster 3" in summary
        assert "T cell" in summary
        assert "[HIGH]" in summary
        assert "PASS" in summary
        assert "marker overlap 80%" in summary
        assert "Accept annotation" in summary

    def test_flagged_summary_carries_action(self):
        row = pd.Series(
            {
                "cluster": "5",
                "cell_type": "B cell",
                "critic_flags": "POSSIBLE_DOUBLET",
                "critic_confidence": "needs_review",
                "pct_overlap": 0.4,
                "combined_score": 0.5,
            }
        )
        summary = format_evidence_summary(row)
        assert "FLAGGED (POSSIBLE_DOUBLET)" in summary
        assert "Sub-cluster or mark as doublet" in summary

    def test_run_critic_adds_evidence_summary_column(
        self, synthetic_pbmc, annotation_summary
    ):
        from celltypepilot.data_adapter import load_marker_atlas

        atlas = load_marker_atlas("human")
        results = run_critic(synthetic_pbmc, "leiden", annotation_summary, atlas, "blood")
        assert "evidence_summary" in results.columns
        for _, row in results.iterrows():
            assert row["evidence_summary"].startswith("Cluster ")
            assert "|" in row["evidence_summary"]
