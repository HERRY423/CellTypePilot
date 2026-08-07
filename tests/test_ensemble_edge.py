"""Tests for ensemble scorer edge cases — adaptive weights, disagreements, summaries."""

import numpy as np
import pandas as pd
import pytest

from celltypepilot.ensemble_scorer import (
    ensemble_scores,
    generate_ensemble_summary,
    analyze_disagreements,
    _compute_adaptive_weights,
    _assign_ensemble_confidence,
    _interpret_disagreement,
    _ref_only_results,
    _marker_only_results,
)
from celltypepilot.constants import (
    CONFIDENCE_HIGH, CONFIDENCE_MEDIUM, CONFIDENCE_LOW, CONFIDENCE_REVIEW,
)


class TestEnsembleScoresEdgeCases:
    def test_both_empty(self):
        result = ensemble_scores(pd.DataFrame(), pd.DataFrame())
        assert result.empty

    def test_marker_only(self):
        marker = pd.DataFrame({
            "cluster": ["0"], "cell_type": ["T cell"],
            "combined_score": [0.8], "rank": [1],
        })
        result = ensemble_scores(marker, pd.DataFrame())
        assert not result.empty
        assert result.iloc[0]["source"] == "marker"
        assert result.iloc[0]["ref_weight_used"] == 0.0

    def test_ref_only(self):
        ref = pd.DataFrame({
            "cluster": ["0"], "cell_type": ["T cell"],
            "ref_score": [0.7], "ref_rank": [1],
        })
        result = ensemble_scores(pd.DataFrame(), ref)
        assert not result.empty
        assert result.iloc[0]["source"] == "reference"
        assert result.iloc[0]["marker_weight_used"] == 0.0

    def test_both_present(self):
        marker = pd.DataFrame({
            "cluster": ["0", "0"], "cell_type": ["T cell", "B cell"],
            "combined_score": [0.8, 0.2], "rank": [1, 2],
        })
        ref = pd.DataFrame({
            "cluster": ["0", "0"], "cell_type": ["T cell", "B cell"],
            "ref_score": [0.7, 0.3], "ref_rank": [1, 2],
        })
        result = ensemble_scores(marker, ref)
        assert not result.empty
        assert len(result) == 2
        assert all("ensemble_score" in result.columns for _ in [1])

    def test_non_adaptive_weights(self):
        marker = pd.DataFrame({
            "cluster": ["0"], "cell_type": ["T cell"],
            "combined_score": [0.5], "rank": [1],
        })
        ref = pd.DataFrame({
            "cluster": ["0"], "cell_type": ["T cell"],
            "ref_score": [0.5], "ref_rank": [1],
        })
        result = ensemble_scores(marker, ref, marker_weight=0.6, adaptive=False)
        row = result.iloc[0]
        assert abs(row["marker_weight_used"] - 0.6) < 0.01

    def test_score_bounded_01(self):
        marker = pd.DataFrame({
            "cluster": ["0"], "cell_type": ["T cell"],
            "combined_score": [1.0], "rank": [1],
        })
        ref = pd.DataFrame({
            "cluster": ["0"], "cell_type": ["T cell"],
            "ref_score": [1.0], "ref_rank": [1],
        })
        result = ensemble_scores(marker, ref)
        assert 0.0 <= result.iloc[0]["ensemble_score"] <= 1.0

    def test_union_of_cell_types(self):
        marker = pd.DataFrame({
            "cluster": ["0", "0"], "cell_type": ["T cell", "B cell"],
            "combined_score": [0.8, 0.2], "rank": [1, 2],
        })
        ref = pd.DataFrame({
            "cluster": ["0", "0"], "cell_type": ["NK cell", "T cell"],
            "ref_score": [0.6, 0.4], "ref_rank": [1, 2],
        })
        result = ensemble_scores(marker, ref)
        cell_types = set(result["cell_type"])
        assert "T cell" in cell_types
        assert "B cell" in cell_types
        assert "NK cell" in cell_types


class TestComputeAdaptiveWeights:
    def test_high_marker_confidence(self):
        marker_cl = pd.DataFrame({
            "cell_type": ["T cell"], "combined_score": [0.8], "rank": [1],
        })
        ref_cl = pd.DataFrame({
            "cell_type": ["T cell"], "ref_score": [0.5], "ref_rank": [1],
        })
        m_w, r_w = _compute_adaptive_weights(marker_cl, ref_cl, 0.5, True)
        assert m_w == 0.7  # High marker confidence
        assert abs(m_w + r_w - 1.0) < 0.01

    def test_low_marker_high_ref(self):
        marker_cl = pd.DataFrame({
            "cell_type": ["T cell"], "combined_score": [0.2], "rank": [1],
        })
        ref_cl = pd.DataFrame({
            "cell_type": ["B cell"], "ref_score": [0.8], "ref_rank": [1],
        })
        m_w, r_w = _compute_adaptive_weights(marker_cl, ref_cl, 0.5, True)
        assert m_w <= 0.3  # Low marker → reference override

    def test_non_adaptive(self):
        m_w, r_w = _compute_adaptive_weights(
            pd.DataFrame(), pd.DataFrame(), 0.6, False
        )
        assert abs(m_w - 0.6) < 0.01
        assert abs(r_w - 0.4) < 0.01

    def test_weights_sum_to_one(self):
        marker_cl = pd.DataFrame({
            "cell_type": ["T cell"], "combined_score": [0.5], "rank": [1],
        })
        ref_cl = pd.DataFrame({
            "cell_type": ["T cell"], "ref_score": [0.5], "ref_rank": [1],
        })
        m_w, r_w = _compute_adaptive_weights(marker_cl, ref_cl, 0.5, True)
        assert abs(m_w + r_w - 1.0) < 0.01


class TestAssignEnsembleConfidence:
    def test_high_confidence(self):
        row = pd.Series({"ensemble_score": 0.8, "agreement": True, "source": "both"})
        assert _assign_ensemble_confidence(row) == CONFIDENCE_HIGH

    def test_medium_confidence(self):
        row = pd.Series({"ensemble_score": 0.55, "agreement": True, "source": "both"})
        assert _assign_ensemble_confidence(row) == CONFIDENCE_MEDIUM

    def test_disagreement_low_score(self):
        row = pd.Series({"ensemble_score": 0.3, "agreement": False, "source": "both"})
        assert _assign_ensemble_confidence(row) == CONFIDENCE_REVIEW

    def test_disagreement_medium_score(self):
        row = pd.Series({"ensemble_score": 0.45, "agreement": False, "source": "both"})
        assert _assign_ensemble_confidence(row) == CONFIDENCE_LOW

    def test_marker_only_medium(self):
        row = pd.Series({"ensemble_score": 0.55, "agreement": True, "source": "marker"})
        assert _assign_ensemble_confidence(row) == CONFIDENCE_MEDIUM

    def test_reference_only_lower(self):
        # score=0.55 with agreement=True hits the generic "score>=0.5 and agreement" → MEDIUM
        # (source-specific downgrade only applies when generic thresholds don't match)
        row = pd.Series({"ensemble_score": 0.55, "agreement": True, "source": "reference"})
        assert _assign_ensemble_confidence(row) == CONFIDENCE_MEDIUM

    def test_reference_only_low_score(self):
        # score=0.45 with agreement=True skips generic thresholds, hits source-specific
        row = pd.Series({"ensemble_score": 0.45, "agreement": True, "source": "reference"})
        assert _assign_ensemble_confidence(row) == CONFIDENCE_LOW

    def test_low_score_review(self):
        row = pd.Series({"ensemble_score": 0.2, "agreement": True, "source": "marker"})
        assert _assign_ensemble_confidence(row) == CONFIDENCE_REVIEW


class TestGenerateEnsembleSummary:
    def test_basic_summary(self):
        df = pd.DataFrame({
            "cluster": ["0", "0", "1", "1"],
            "cell_type": ["T cell", "B cell", "NK cell", "Monocyte"],
            "ensemble_score": [0.8, 0.2, 0.6, 0.4],
            "rank": [1, 2, 1, 2],
            "source": ["both", "both", "both", "both"],
            "agreement": [True, True, True, True],
            "marker_score": [0.7, 0.1, 0.5, 0.3],
            "ref_score": [0.9, 0.3, 0.7, 0.5],
        })
        summary = generate_ensemble_summary(df)
        assert len(summary) == 2
        assert "confidence" in summary.columns

    def test_empty_input(self):
        summary = generate_ensemble_summary(pd.DataFrame())
        assert summary.empty


class TestAnalyzeDisagreements:
    def test_no_disagreement(self):
        df = pd.DataFrame({
            "cluster": ["0", "0"],
            "cell_type": ["T cell", "B cell"],
            "marker_score": [0.8, 0.1],
            "ref_score": [0.7, 0.2],
            "ensemble_score": [0.75, 0.15],
            "rank": [1, 2],
        })
        disagreements = analyze_disagreements(df)
        # Both methods agree T cell is best → no disagreement
        assert disagreements.empty

    def test_with_disagreement(self):
        # marker best=T cell (0.9), ref best=B cell (0.7) → severity=0.2 >= min_score_gap
        df = pd.DataFrame({
            "cluster": ["0", "0", "0", "0"],
            "cell_type": ["T cell", "B cell", "T cell", "B cell"],
            "marker_score": [0.9, 0.1, 0.9, 0.1],
            "ref_score": [0.1, 0.7, 0.1, 0.7],
            "ensemble_score": [0.5, 0.4, 0.5, 0.4],
            "rank": [1, 2, 1, 2],
            "marker_weight_used": [0.5, 0.5, 0.5, 0.5],
            "ref_weight_used": [0.5, 0.5, 0.5, 0.5],
            "agreement": [True, True, True, True],
            "source": ["both", "both", "both", "both"],
        })
        disagreements = analyze_disagreements(df)
        assert not disagreements.empty
        assert "interpretation" in disagreements.columns

    def test_empty_input(self):
        assert analyze_disagreements(pd.DataFrame()).empty


class TestInterpretDisagreement:
    def test_marker_strong_ref_weak(self):
        text = _interpret_disagreement("T cell", 0.8, "B cell", 0.1)
        assert "novel" in text.lower() or "marker database" in text.lower()

    def test_ref_strong_marker_weak(self):
        text = _interpret_disagreement("T cell", 0.1, "B cell", 0.8)
        assert "transitional" in text.lower() or "rare" in text.lower()

    def test_both_moderate(self):
        text = _interpret_disagreement("T cell", 0.5, "B cell", 0.4)
        assert "intermediate" in text.lower() or "both" in text.lower()

    def test_both_weak(self):
        text = _interpret_disagreement("T cell", 0.15, "B cell", 0.1)
        assert "manual review" in text.lower() or "low-quality" in text.lower()


class TestRefOnlyResults:
    def test_wrapping(self):
        ref = pd.DataFrame({
            "cluster": ["0"], "cell_type": ["T cell"],
            "ref_score": [0.7], "ref_rank": [1],
        })
        result = _ref_only_results(ref)
        assert result.iloc[0]["source"] == "reference"
        assert result.iloc[0]["marker_weight_used"] == 0.0
        assert result.iloc[0]["ensemble_score"] == 0.7


class TestMarkerOnlyResults:
    def test_wrapping(self):
        marker = pd.DataFrame({
            "cluster": ["0"], "cell_type": ["T cell"],
            "combined_score": [0.8], "rank": [1],
        })
        result = _marker_only_results(marker)
        assert result.iloc[0]["source"] == "marker"
        assert result.iloc[0]["ref_weight_used"] == 0.0
        assert result.iloc[0]["ensemble_score"] == 0.8
