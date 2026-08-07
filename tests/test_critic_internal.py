"""Tests for critic internal functions — evidence, negative markers, doublet, ontology."""

import numpy as np
import pandas as pd
import anndata as ad
import pytest

from celltypepilot.critic import (
    _check_evidence_sufficiency,
    _check_negative_markers,
    _check_doublet_signal,
    _check_ontology_consistency,
    _recalibrate_confidence,
    _check_ensemble_agreement,
    generate_critic_summary,
)
from celltypepilot.constants import (
    CONFIDENCE_HIGH, CONFIDENCE_MEDIUM, CONFIDENCE_LOW, CONFIDENCE_REVIEW,
)


class TestCheckEvidenceSufficiency:
    def test_no_markers(self, synthetic_pbmc):
        result = _check_evidence_sufficiency(
            synthetic_pbmc, "0", "leiden", []
        )
        assert result["flag"] == "NO_MARKERS"

    def test_good_coverage(self, synthetic_pbmc):
        # Use genes that exist in the synthetic data
        markers = ["CD3D", "CD3E", "CD2"]
        result = _check_evidence_sufficiency(
            synthetic_pbmc,
            synthetic_pbmc.obs["leiden"].iloc[0],
            "leiden", markers
        )
        assert result["flag"] in ("", "PARTIAL_EVIDENCE", "LOW_EVIDENCE")
        assert "Coverage:" in result["evidence"]

    def test_nonexistent_genes(self, synthetic_pbmc):
        markers = ["FAKE1", "FAKE2", "FAKE3"]
        cluster = synthetic_pbmc.obs["leiden"].iloc[0]
        result = _check_evidence_sufficiency(
            synthetic_pbmc, cluster, "leiden", markers
        )
        # All genes missing → 0 coverage → LOW_EVIDENCE
        assert result["flag"] in ("LOW_EVIDENCE", "NO_MARKERS")


class TestCheckNegativeMarkers:
    def test_no_negative_markers(self, synthetic_pbmc):
        cluster = synthetic_pbmc.obs["leiden"].iloc[0]
        result = _check_negative_markers(
            synthetic_pbmc, cluster, "leiden", []
        )
        assert result["flag"] == ""

    def test_with_negative_markers(self, synthetic_pbmc):
        cluster = synthetic_pbmc.obs["leiden"].iloc[0]
        result = _check_negative_markers(
            synthetic_pbmc, cluster, "leiden", ["GENE_0", "GENE_1"]
        )
        assert isinstance(result["flag"], str)
        assert "evidence" in result


class TestCheckDoubletSignal:
    def test_no_doublet(self, synthetic_pbmc):
        from celltypepilot.data_adapter import load_marker_atlas, get_all_markers_for_tissue
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
        result = _recalibrate_confidence(
            CONFIDENCE_HIGH, [], {"flag": ""}, {"flag": ""}
        )
        assert result == CONFIDENCE_HIGH

    def test_low_evidence_downgrades(self):
        result = _recalibrate_confidence(
            CONFIDENCE_HIGH, ["LOW_EVIDENCE"],
            {"flag": "LOW_EVIDENCE"}, {"flag": ""}
        )
        assert result == CONFIDENCE_REVIEW

    def test_neg_marker_conflict_downgrades(self):
        result = _recalibrate_confidence(
            CONFIDENCE_HIGH, ["NEG_MARKER_CONFLICT"],
            {"flag": ""}, {"flag": "NEG_MARKER_CONFLICT"}
        )
        assert result == CONFIDENCE_REVIEW

    def test_partial_evidence_mild_downgrade(self):
        result = _recalibrate_confidence(
            CONFIDENCE_HIGH, ["PARTIAL_EVIDENCE"],
            {"flag": "PARTIAL_EVIDENCE"}, {"flag": ""}
        )
        assert result == CONFIDENCE_LOW

    def test_no_markers_severe(self):
        result = _recalibrate_confidence(
            CONFIDENCE_HIGH, ["NO_MARKERS"],
            {"flag": "NO_MARKERS"}, {"flag": ""}
        )
        assert result == CONFIDENCE_REVIEW

    def test_unknown_confidence_input(self):
        result = _recalibrate_confidence(
            "unknown_conf", [], {"flag": ""}, {"flag": ""}
        )
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
        results = pd.DataFrame({
            "cluster": ["0", "1", "2"],
            "critic_flags": ["PASS", "PASS", "PASS"],
            "critic_confidence": ["high", "medium", "low"],
        })
        summary = generate_critic_summary(results)
        assert summary["total_clusters"] == 3
        assert summary["pass"] == 3
        assert summary["flagged"] == 0
        assert summary["confidence_distribution"]["high"] == 1

    def test_with_flags(self):
        results = pd.DataFrame({
            "cluster": ["0", "1"],
            "critic_flags": ["PASS", "LOW_EVIDENCE; NEG_MARKER_CONFLICT"],
            "critic_confidence": ["high", "needs_review"],
        })
        summary = generate_critic_summary(results)
        assert summary["flagged"] == 1
        assert "LOW_EVIDENCE" in summary["flag_types"]
        assert "NEG_MARKER_CONFLICT" in summary["flag_types"]
        assert "1" in summary["clusters_needing_review"]
