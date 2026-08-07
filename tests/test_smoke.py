"""Smoke tests for CellTypePilot using synthetic data."""

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import pytest


def make_synthetic_pbmc(n_cells=500, seed=42) -> ad.AnnData:
    """Create a synthetic PBMC-like AnnData for testing.

    Generates data with 5 cell types and known marker genes.
    """
    rng = np.random.RandomState(seed)

    # Define cell types and their markers
    cell_types = {
        "T cell": ["CD3D", "CD3E", "CD2", "TRAC", "IL7R"],
        "B cell": ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5"],
        "NK cell": ["NCAM1", "NKG7", "GNLY", "KLRD1", "PRF1"],
        "Monocyte": ["CD14", "LYZ", "S100A8", "S100A9", "FCN1"],
        "Dendritic cell": ["FCER1A", "CD1C", "CLEC10A", "ITGAX", "HLA-DRA"],
    }

    n_per_type = n_cells // len(cell_types)
    all_genes = set()
    for markers in cell_types.values():
        all_genes.update(markers)
    # Add background genes
    bg_genes = [f"GENE_{i}" for i in range(50)]
    all_genes = sorted(all_genes) + bg_genes
    all_genes = list(all_genes)

    X_list = []
    labels = []
    for ct_name, markers in cell_types.items():
        # Base expression: low for all genes
        expr = rng.exponential(0.5, size=(n_per_type, len(all_genes)))
        # Boost marker genes for this cell type
        for marker in markers:
            if marker in all_genes:
                idx = all_genes.index(marker)
                expr[:, idx] += rng.exponential(3.0, size=n_per_type)
        X_list.append(expr)
        labels.extend([ct_name] * n_per_type)

    X = np.vstack(X_list)
    # Log-transform
    X_log = np.log1p(X)

    var = pd.DataFrame(index=all_genes)
    obs = pd.DataFrame({"true_cell_type": labels})

    adata = ad.AnnData(X=X_log.astype(np.float32))
    adata.var_names = all_genes
    adata.obs = obs
    adata.layers["counts"] = X.astype(np.float32)

    # Add PCA and UMAP
    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Assign cluster labels based on true cell types (simulating clustering)
    # Each true cell type gets 1-2 clusters
    cluster_map = {}
    for i, ct in enumerate(labels):
        # Add some sub-clusters
        sub = 0 if i % 3 != 0 else 1
        cluster_map[i] = f"{ct}_{sub}"
    adata.obs["leiden"] = pd.Categorical(list(cluster_map.values()))

    return adata


class TestDataAdapter:
    """Tests for the data adapter module."""

    def test_detect_species_human(self):
        from celltypepilot.data_adapter import detect_species
        adata = make_synthetic_pbmc()
        assert detect_species(adata) == "human"

    def test_detect_species_mouse(self):
        from celltypepilot.data_adapter import detect_species
        adata = make_synthetic_pbmc()
        # Convert to mouse naming
        mouse_genes = [g[0].upper() + g[1:].lower() for g in adata.var_names]
        adata.var_names = mouse_genes
        assert detect_species(adata) == "mouse"

    def test_find_cluster_keys(self):
        from celltypepilot.data_adapter import find_cluster_keys
        adata = make_synthetic_pbmc()
        keys = find_cluster_keys(adata)
        assert "leiden" in keys

    def test_find_embedding_keys(self):
        from celltypepilot.data_adapter import find_embedding_keys
        adata = make_synthetic_pbmc()
        keys = find_embedding_keys(adata)
        assert any("umap" in k for k in keys)

    def test_inspect_adata(self):
        from celltypepilot.data_adapter import inspect_adata
        import tempfile
        tmpdir = tempfile.mkdtemp()
        h5ad_path = Path(tmpdir) / "test.h5ad"
        adata = make_synthetic_pbmc()
        adata.write(str(h5ad_path))
        report = inspect_adata(str(h5ad_path), cluster_key="leiden")
        assert report["n_obs"] == 500
        assert report["species"] == "human"
        assert "leiden" in report["cluster_keys"]

    def test_load_marker_atlas(self):
        from celltypepilot.data_adapter import load_marker_atlas
        atlas = load_marker_atlas("human")
        assert "tissues" in atlas
        assert "blood" in atlas["tissues"]

    def test_load_marker_atlas_mouse(self):
        from celltypepilot.data_adapter import load_marker_atlas
        atlas = load_marker_atlas("mouse")
        assert "tissues" in atlas
        # Check mouse gene conversion
        blood_types = atlas["tissues"]["blood"]["cell_types"]
        t_cell_markers = blood_types["T cell"]["positive_markers"]
        assert all(m[0].isupper() and (len(m) < 2 or m[1].islower()) for m in t_cell_markers if m.isalpha())


class TestMarkerScorer:
    """Tests for the marker scorer."""

    def test_compute_marker_scores(self):
        from celltypepilot.marker_scorer import compute_marker_scores
        from celltypepilot.data_adapter import load_marker_atlas, get_all_markers_for_tissue

        adata = make_synthetic_pbmc()
        atlas = load_marker_atlas("human")
        markers = get_all_markers_for_tissue(atlas, "blood")

        scores = compute_marker_scores(adata, "leiden", markers)
        assert not scores.empty
        assert "combined_score" in scores.columns
        assert "cluster" in scores.columns
        assert "cell_type" in scores.columns

    def test_generate_annotation_summary(self):
        from celltypepilot.marker_scorer import compute_marker_scores, generate_annotation_summary
        from celltypepilot.data_adapter import load_marker_atlas, get_all_markers_for_tissue

        adata = make_synthetic_pbmc()
        atlas = load_marker_atlas("human")
        markers = get_all_markers_for_tissue(atlas, "blood")

        scores = compute_marker_scores(adata, "leiden", markers)
        summary = generate_annotation_summary(scores, "leiden")
        assert not summary.empty
        assert "cell_type" in summary.columns
        assert "confidence" in summary.columns


class TestCritic:
    """Tests for the annotation critic."""

    def test_run_critic(self):
        from celltypepilot.marker_scorer import compute_marker_scores, generate_annotation_summary
        from celltypepilot.critic import run_critic, generate_critic_summary
        from celltypepilot.data_adapter import load_marker_atlas, get_all_markers_for_tissue

        adata = make_synthetic_pbmc()
        atlas = load_marker_atlas("human")
        markers = get_all_markers_for_tissue(atlas, "blood")

        scores = compute_marker_scores(adata, "leiden", markers)
        summary = generate_annotation_summary(scores, "leiden")
        critic_results = run_critic(adata, "leiden", summary, atlas, "blood")

        assert "critic_flags" in critic_results.columns
        assert "critic_confidence" in critic_results.columns

        critic_summary = generate_critic_summary(critic_results)
        assert critic_summary["total_clusters"] == len(summary)


class TestDoctor:
    """Tests for the doctor module."""

    def test_run_doctor(self):
        from celltypepilot.doctor import run_doctor
        report = run_doctor()
        assert report.python_ok is True
        assert len(report.dependencies) > 0


class TestProvenance:
    """Tests for provenance tracking."""

    def test_create_and_save_manifest(self):
        from celltypepilot.provenance import create_manifest, save_manifest, load_manifest

        with tempfile.TemporaryDirectory() as tmpdir:
            manifest = create_manifest(
                input_path="test.h5ad",
                data_hash="abc123",
                cluster_key="leiden",
                species="human",
                tissue="blood",
                parameters={"embedding_key": "X_umap"},
                output_dir=tmpdir,
            )
            assert manifest["celltypepilot_version"] == "0.1.0"
            assert manifest["mkg_version"] == "mkg-2026.08"

            path = save_manifest(manifest, tmpdir)
            loaded = load_manifest(path)
            assert loaded["input"]["sha256"] == "abc123"


class TestCLI:
    """Tests for the CLI interface."""

    def test_doctor_command(self):
        from typer.testing import CliRunner
        from celltypepilot.cli import app

        runner = CliRunner()
        result = runner.invoke(app, ["doctor"])
        assert result.exit_code == 0

    def test_version(self):
        from typer.testing import CliRunner
        from celltypepilot.cli import app

        runner = CliRunner()
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert "0.1.0" in result.output

    def test_markers_command(self):
        from typer.testing import CliRunner
        from celltypepilot.cli import app

        runner = CliRunner()
        result = runner.invoke(app, ["markers"])
        assert result.exit_code == 0
        assert "blood" in result.output

    def test_markers_tissue(self):
        from typer.testing import CliRunner
        from celltypepilot.cli import app

        runner = CliRunner()
        result = runner.invoke(app, ["markers", "--tissue", "blood"])
        assert result.exit_code == 0
        assert "T cell" in result.output

    def test_inspect_command(self):
        from typer.testing import CliRunner
        from celltypepilot.cli import app

        tmpdir = tempfile.mkdtemp()
        h5ad_path = Path(tmpdir) / "test_inspect.h5ad"
        adata = make_synthetic_pbmc()
        adata.write(str(h5ad_path))

        runner = CliRunner()
        result = runner.invoke(app, ["inspect", "--input", str(h5ad_path), "--cluster-key", "leiden"])
        assert result.exit_code == 0

    def test_annotate_command(self):
        from typer.testing import CliRunner
        from celltypepilot.cli import app

        with tempfile.TemporaryDirectory() as tmpdir:
            h5ad_path = Path(tmpdir) / "test.h5ad"
            output_dir = Path(tmpdir) / "output"
            adata = make_synthetic_pbmc()
            adata.write(str(h5ad_path))

            runner = CliRunner()
            result = runner.invoke(app, [
                "annotate",
                "--input", str(h5ad_path),
                "--cluster-key", "leiden",
                "--output", str(output_dir),
                "--species", "human",
                "--tissue", "blood",
                "--embedding-key", "X_umap",
            ])
            assert result.exit_code == 0
            assert (output_dir / "evidence_table.csv").exists()
            assert (output_dir / "manifest.json").exists()
            assert (output_dir / "data.annotated.h5ad").exists()
