"""Tests for constants module — verify thresholds, paths, and configuration integrity."""

from celltypepilot.constants import (
    ATLAS_PATH,
    CB_PALETTE,
    CONFIDENCE_HIGH,
    CONFIDENCE_LOW,
    CONFIDENCE_MEDIUM,
    CONFIDENCE_REVIEW,
    CRITIC_DOUBLET_COEXPR_THRESHOLD,
    CRITIC_LOW_COVERAGE_THRESHOLD,
    CRITIC_NEG_MARKER_PCT_THRESHOLD,
    ENSEMBLE_AGREEMENT_THRESHOLD,
    ENSEMBLE_MARKER_HIGH,
    ENSEMBLE_MARKER_LOW,
    ENSEMBLE_REF_OVERRIDE,
    MARKER_FC_THRESHOLD,
    MARKER_PCT_THRESHOLD,
    MARKER_SPECIFICITY_THRESHOLD,
    MIN_CLUSTER_SIZE,
    OUTPUT_ANNOTATED,
    OUTPUT_EVIDENCE,
    OUTPUT_FIGURES_DIR,
    OUTPUT_MANIFEST,
    OUTPUT_REPORT,
    PKG_ROOT,
    REF_CORR_MIN_GENES,
    REF_KNN_DEFAULT_K,
    REF_KNN_MAX_K,
    REF_MIN_SHARED_GENES,
    SPECIES_HUMAN,
    SPECIES_MOUSE,
)


class TestConfidenceConstants:
    """Confidence level string constants."""

    def test_confidence_values(self):
        assert CONFIDENCE_HIGH == "high"
        assert CONFIDENCE_MEDIUM == "medium"
        assert CONFIDENCE_LOW == "low"
        assert CONFIDENCE_REVIEW == "needs_review"

    def test_confidence_levels_are_distinct(self):
        levels = {CONFIDENCE_HIGH, CONFIDENCE_MEDIUM, CONFIDENCE_LOW, CONFIDENCE_REVIEW}
        assert len(levels) == 4


class TestSpeciesConstants:
    def test_species_values(self):
        assert SPECIES_HUMAN == "human"
        assert SPECIES_MOUSE == "mouse"


class TestThresholds:
    """Verify threshold values are in valid ranges."""

    def test_marker_thresholds_in_range(self):
        assert 0.0 < MARKER_PCT_THRESHOLD <= 1.0
        assert MARKER_FC_THRESHOLD > 0
        assert 0.0 < MARKER_SPECIFICITY_THRESHOLD <= 1.0

    def test_critic_thresholds_in_range(self):
        assert 0.0 < CRITIC_NEG_MARKER_PCT_THRESHOLD <= 1.0
        assert 0.0 < CRITIC_DOUBLET_COEXPR_THRESHOLD <= 1.0
        assert 0.0 < CRITIC_LOW_COVERAGE_THRESHOLD <= 1.0

    def test_ensemble_thresholds_in_range(self):
        assert 0.0 < ENSEMBLE_AGREEMENT_THRESHOLD <= 1.0
        assert 0.0 < ENSEMBLE_MARKER_HIGH <= 1.0
        assert 0.0 < ENSEMBLE_MARKER_LOW <= 1.0
        assert 0.0 < ENSEMBLE_REF_OVERRIDE <= 1.0

    def test_ensemble_marker_high_gt_low(self):
        assert ENSEMBLE_MARKER_HIGH > ENSEMBLE_MARKER_LOW

    def test_reference_thresholds_positive(self):
        assert REF_MIN_SHARED_GENES > 0
        assert REF_CORR_MIN_GENES > 0
        assert REF_KNN_DEFAULT_K > 0
        assert REF_KNN_MAX_K >= REF_KNN_DEFAULT_K

    def test_min_cluster_size_positive(self):
        assert MIN_CLUSTER_SIZE > 0


class TestOutputConstants:
    def test_output_filenames(self):
        assert OUTPUT_ANNOTATED == "data.annotated.h5ad"
        assert OUTPUT_EVIDENCE == "evidence_table.csv"
        assert OUTPUT_MANIFEST == "manifest.json"
        assert OUTPUT_FIGURES_DIR == "figures"
        assert OUTPUT_REPORT == "report_draft.html"

    def test_output_filenames_have_extensions(self):
        for name in [OUTPUT_ANNOTATED, OUTPUT_EVIDENCE, OUTPUT_MANIFEST, OUTPUT_REPORT]:
            assert "." in name


class TestPaths:
    def test_pkg_root_exists(self):
        assert PKG_ROOT.exists()
        assert PKG_ROOT.is_dir()

    def test_atlas_path_exists(self):
        assert ATLAS_PATH.exists()
        assert ATLAS_PATH.suffix == ".json"


class TestColorPalette:
    def test_palette_not_empty(self):
        assert len(CB_PALETTE) >= 8

    def test_palette_hex_format(self):
        for color in CB_PALETTE:
            assert color.startswith("#")
            assert len(color) == 7

    def test_palette_unique_colors(self):
        assert len(set(CB_PALETTE)) == len(CB_PALETTE)
