"""Tests for the Scanpy-native Python API.

Validates that annotate(), inspect(), and critic_review() work correctly
with in-memory AnnData objects, without requiring file I/O.
"""

import tempfile
from pathlib import Path

import anndata as ad
import pandas as pd
import pytest

# ──────────────────────────────────────────────
# Fixtures
# ──────────────────────────────────────────────


@pytest.fixture
def pbmc_adata(synthetic_pbmc):
    """Alias for the shared synthetic PBMC fixture."""
    return synthetic_pbmc


# ──────────────────────────────────────────────
# inspect() API tests
# ──────────────────────────────────────────────


class TestInspectAPI:
    """Tests for the Scanpy-native inspect() function."""

    def test_inspect_returns_dict(self, pbmc_adata):
        from celltypepilot import inspect

        report = inspect(pbmc_adata)
        assert isinstance(report, dict)

    def test_inspect_detects_species(self, pbmc_adata):
        from celltypepilot import inspect

        report = inspect(pbmc_adata)
        assert report["species"] == "human"

    def test_inspect_reports_dimensions(self, pbmc_adata):
        from celltypepilot import inspect

        report = inspect(pbmc_adata)
        assert report["n_obs"] == 500
        assert report["n_vars"] == pbmc_adata.n_vars

    def test_inspect_finds_cluster_keys(self, pbmc_adata):
        from celltypepilot import inspect

        report = inspect(pbmc_adata)
        assert "leiden" in report["cluster_keys"]

    def test_inspect_finds_embedding_keys(self, pbmc_adata):
        from celltypepilot import inspect

        report = inspect(pbmc_adata)
        assert any("umap" in k for k in report["embedding_keys"])

    def test_inspect_with_cluster_key(self, pbmc_adata):
        from celltypepilot import inspect

        report = inspect(pbmc_adata, cluster_key="leiden")
        assert report["cluster_sizes"]  # non-empty
        assert not report["fatal"]

    def test_inspect_without_cluster_key_warns(self, pbmc_adata):
        from celltypepilot import inspect

        report = inspect(pbmc_adata)
        # Should warn about no cluster_key specified
        assert any("cluster_key" in w for w in report["warnings"])

    def test_inspect_no_file_fields(self, pbmc_adata):
        """inspect() should NOT have file-specific fields."""
        from celltypepilot import inspect

        report = inspect(pbmc_adata)
        assert "path" not in report
        assert "sha256" not in report


# ──────────────────────────────────────────────
# annotate() API tests
# ──────────────────────────────────────────────


class TestAnnotateAPI:
    """Tests for the Scanpy-native annotate() function."""

    def test_annotate_returns_anndata(self, pbmc_adata):
        from celltypepilot import annotate

        result = annotate(pbmc_adata, "leiden", species="human", tissue="blood")
        assert isinstance(result, ad.AnnData)

    def test_annotate_modifies_obs_in_place(self, pbmc_adata):
        from celltypepilot import annotate

        annotate(pbmc_adata, "leiden", species="human", tissue="blood")
        assert "ctp_cell_type" in pbmc_adata.obs.columns
        assert "ctp_confidence" in pbmc_adata.obs.columns
        assert "ctp_cl_id" in pbmc_adata.obs.columns
        assert "ctp_decision" in pbmc_adata.obs.columns

    def test_annotate_adds_all_ctp_columns(self, pbmc_adata):
        from celltypepilot import annotate

        annotate(pbmc_adata, "leiden", species="human", tissue="blood")
        expected_columns = [
            "ctp_cell_type",
            "ctp_cl_id",
            "ctp_confidence",
            "ctp_candidate_cell_type",
            "ctp_decision",
            "ctp_abstain_reason",
            "ctp_cell_state_candidate",
            "ctp_state_decision",
            "ctp_cell_state",
            "ctp_state_score",
            "ctp_state_confidence",
            "ctp_state_evidence",
            "ctp_display_label",
        ]
        for col in expected_columns:
            assert col in pbmc_adata.obs.columns, f"Missing column: {col}"

    def test_annotate_auto_detects_species(self, pbmc_adata):
        from celltypepilot import annotate

        result = annotate(pbmc_adata, "leiden", tissue="blood")
        # Should auto-detect human and still produce annotations
        assert "ctp_cell_type" in result.obs.columns

    def test_annotate_invalid_cluster_key_raises(self, pbmc_adata):
        from celltypepilot import annotate
        from celltypepilot.orchestrator import PipelineError

        with pytest.raises(PipelineError, match="cluster key"):
            annotate(pbmc_adata, "nonexistent_key", species="human", tissue="blood")

    def test_annotate_no_file_output(self, pbmc_adata):
        """annotate() without output_dir should NOT write files."""
        from celltypepilot import annotate

        annotate(pbmc_adata, "leiden", species="human", tissue="blood")
        # No output_dir → no files written. Just verify obs is populated.
        assert pbmc_adata.obs["ctp_cell_type"].notna().all()

    def test_annotate_with_output_dir(self, pbmc_adata):
        """annotate() with output_dir should write artifacts."""
        from celltypepilot import annotate

        with tempfile.TemporaryDirectory() as tmpdir:
            annotate(
                pbmc_adata,
                "leiden",
                species="human",
                tissue="blood",
                output_dir=tmpdir,
                no_figures=True,
            )
            assert (Path(tmpdir) / "evidence_table.csv").exists()
            assert (Path(tmpdir) / "manifest.json").exists()
            assert (Path(tmpdir) / "data.annotated.h5ad").exists()

    def test_annotate_with_progress_callback(self, pbmc_adata):
        from celltypepilot import annotate

        progress_calls = []

        def progress(step, total, message):
            progress_calls.append((step, total, message))

        annotate(
            pbmc_adata,
            "leiden",
            species="human",
            tissue="blood",
            progress=progress,
        )
        assert len(progress_calls) > 0
        assert all(step <= total for step, total, _ in progress_calls)

    def test_annotate_returns_same_object(self, pbmc_adata):
        """annotate() should return the same AnnData object (mutated in-place)."""
        from celltypepilot import annotate

        result = annotate(pbmc_adata, "leiden", species="human", tissue="blood")
        assert result is pbmc_adata

    def test_annotate_disable_states(self, pbmc_adata):
        from celltypepilot import annotate

        annotate(
            pbmc_adata,
            "leiden",
            species="human",
            tissue="blood",
            enable_states=False,
        )
        assert "ctp_cell_type" in pbmc_adata.obs.columns
        # State columns should still exist but with abstain defaults
        assert "ctp_state_decision" in pbmc_adata.obs.columns


# ──────────────────────────────────────────────
# critic_review() API tests
# ──────────────────────────────────────────────


class TestCriticReviewAPI:
    """Tests for the Scanpy-native critic_review() function."""

    def test_critic_review_returns_dict(self, pbmc_adata):
        from celltypepilot import critic_review

        # Get a valid cluster ID
        cluster = str(pbmc_adata.obs["leiden"].iloc[0])
        result = critic_review(pbmc_adata, "leiden", cluster, species="human", tissue="blood")
        assert isinstance(result, dict)

    def test_critic_review_has_expected_keys(self, pbmc_adata):
        from celltypepilot import critic_review

        cluster = str(pbmc_adata.obs["leiden"].iloc[0])
        result = critic_review(pbmc_adata, "leiden", cluster, species="human", tissue="blood")
        assert "cluster" in result
        assert "candidates" in result
        assert "critic_results" in result
        assert "critic_summary" in result

    def test_critic_review_candidates_format(self, pbmc_adata):
        from celltypepilot import critic_review

        cluster = str(pbmc_adata.obs["leiden"].iloc[0])
        result = critic_review(pbmc_adata, "leiden", cluster, species="human", tissue="blood")
        assert len(result["candidates"]) > 0
        for candidate in result["candidates"]:
            assert "cell_type" in candidate
            assert "score" in candidate
            assert "overlap" in candidate

    def test_critic_review_invalid_cluster_raises(self, pbmc_adata):
        from celltypepilot import critic_review
        from celltypepilot.orchestrator import PipelineError

        with pytest.raises(PipelineError, match="not found"):
            critic_review(
                pbmc_adata, "leiden", "nonexistent_cluster", species="human", tissue="blood"
            )

    def test_critic_review_invalid_cluster_key_raises(self, pbmc_adata):
        from celltypepilot import critic_review
        from celltypepilot.orchestrator import PipelineError

        with pytest.raises(PipelineError, match="cluster key"):
            critic_review(pbmc_adata, "nonexistent_key", "0", species="human", tissue="blood")

    def test_critic_review_auto_detects_species(self, pbmc_adata):
        from celltypepilot import critic_review

        cluster = str(pbmc_adata.obs["leiden"].iloc[0])
        result = critic_review(pbmc_adata, "leiden", cluster, tissue="blood")
        assert result["cluster"] == cluster


# ──────────────────────────────────────────────
# Lazy import tests
# ──────────────────────────────────────────────


class TestLazyImports:
    """Verify that the new API functions are lazily importable."""

    def test_annotate_is_lazy_importable(self):
        import celltypepilot

        assert hasattr(celltypepilot, "annotate")
        assert callable(celltypepilot.annotate)

    def test_inspect_is_lazy_importable(self):
        import celltypepilot

        assert hasattr(celltypepilot, "inspect")
        assert callable(celltypepilot.inspect)

    def test_critic_review_is_lazy_importable(self):
        import celltypepilot

        assert hasattr(celltypepilot, "critic_review")
        assert callable(celltypepilot.critic_review)

    def test_all_exports_present(self):
        import celltypepilot

        for name in celltypepilot.__all__:
            assert hasattr(celltypepilot, name), f"Missing export: {name}"


# ──────────────────────────────────────────────
# Integration: Scanpy workflow
# ──────────────────────────────────────────────


class TestScanpyWorkflowIntegration:
    """Test that the API fits naturally into a Scanpy workflow."""

    def test_typical_scanpy_workflow(self, pbmc_adata):
        """Simulate a typical Scanpy → CellTypePilot workflow."""
        import celltypepilot as ctp

        # Step 1: Inspect
        report = ctp.inspect(pbmc_adata, cluster_key="leiden")
        assert report["species"] == "human"
        assert "leiden" in report["cluster_keys"]

        # Step 2: Annotate
        adata = ctp.annotate(
            pbmc_adata,
            "leiden",
            species="human",
            tissue="blood",
        )
        assert "ctp_cell_type" in adata.obs.columns

        # Step 3: Use results in Scanpy plotting
        cell_type_dtype = adata.obs["ctp_cell_type"].dtype
        assert pd.api.types.is_string_dtype(cell_type_dtype) or isinstance(
            cell_type_dtype, pd.CategoricalDtype
        )

        # Step 4: Deep-review a flagged cluster
        flagged = adata.obs[adata.obs["ctp_confidence"] == "needs_review"]
        if not flagged.empty:
            cluster = str(flagged["leiden"].iloc[0])
            review = ctp.critic_review(adata, "leiden", cluster, species="human", tissue="blood")
            assert review["cluster"] == cluster

    def test_annotate_then_visualize(self, pbmc_adata):
        """Verify annotated data can be used with Scanpy plotting."""
        import celltypepilot as ctp

        ctp.annotate(pbmc_adata, "leiden", species="human", tissue="blood")

        # ctp_cell_type should be plottable via sc.pl.umap
        assert "ctp_cell_type" in pbmc_adata.obs.columns
        assert "X_umap" in pbmc_adata.obsm
