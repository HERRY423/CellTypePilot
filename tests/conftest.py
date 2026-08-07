"""Shared test fixtures for CellTypePilot test suite."""

import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scanpy as sc


@pytest.fixture
def rng():
    """Reproducible random number generator."""
    return np.random.RandomState(42)


@pytest.fixture
def synthetic_pbmc(rng):
    """Create a synthetic PBMC-like AnnData for testing.

    Generates data with 5 cell types and known marker genes.
    """
    cell_types = {
        "T cell": ["CD3D", "CD3E", "CD2", "TRAC", "IL7R"],
        "B cell": ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5"],
        "NK cell": ["NCAM1", "NKG7", "GNLY", "KLRD1", "PRF1"],
        "Monocyte": ["CD14", "LYZ", "S100A8", "S100A9", "FCN1"],
        "Dendritic cell": ["FCER1A", "CD1C", "CLEC10A", "ITGAX", "HLA-DRA"],
    }

    n_per_type = 100
    all_genes = set()
    for markers in cell_types.values():
        all_genes.update(markers)
    bg_genes = [f"GENE_{i}" for i in range(50)]
    all_genes = sorted(all_genes) + bg_genes

    X_list = []
    labels = []
    for ct_name, markers in cell_types.items():
        expr = rng.exponential(0.5, size=(n_per_type, len(all_genes)))
        for marker in markers:
            if marker in all_genes:
                idx = all_genes.index(marker)
                expr[:, idx] += rng.exponential(3.0, size=n_per_type)
        X_list.append(expr)
        labels.extend([ct_name] * n_per_type)

    X = np.vstack(X_list)
    X_log = np.log1p(X)

    obs = pd.DataFrame({"true_cell_type": labels})

    adata = ad.AnnData(X=X_log.astype(np.float32))
    adata.var_names = all_genes
    adata.obs = obs
    adata.layers["counts"] = X.astype(np.float32)

    sc.pp.pca(adata, n_comps=20)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    cluster_map = {}
    for i, ct in enumerate(labels):
        sub = 0 if i % 3 != 0 else 1
        cluster_map[i] = f"{ct}_{sub}"
    adata.obs["leiden"] = pd.Categorical(list(cluster_map.values()))

    return adata


@pytest.fixture
def tmp_output_dir():
    """Temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def h5ad_path(synthetic_pbmc, tmp_output_dir):
    """Write synthetic data to a temporary .h5ad file."""
    path = tmp_output_dir / "test_data.h5ad"
    synthetic_pbmc.write(str(path))
    return path


@pytest.fixture
def blood_atlas():
    """Load the human blood marker atlas."""
    from celltypepilot.data_adapter import load_marker_atlas

    return load_marker_atlas("human")


@pytest.fixture
def blood_markers(blood_atlas):
    """Get blood tissue markers from atlas."""
    from celltypepilot.data_adapter import get_all_markers_for_tissue

    return get_all_markers_for_tissue(blood_atlas, "blood")


@pytest.fixture
def marker_scores(synthetic_pbmc, blood_markers):
    """Pre-computed marker scores for synthetic data."""
    from celltypepilot.marker_scorer import compute_marker_scores

    return compute_marker_scores(synthetic_pbmc, "leiden", blood_markers)


@pytest.fixture
def annotation_summary(marker_scores):
    """Pre-computed annotation summary."""
    from celltypepilot.marker_scorer import generate_annotation_summary

    return generate_annotation_summary(marker_scores, "leiden")


@pytest.fixture
def sample_critic_results(annotation_summary):
    """Mock critic results DataFrame for report/visualizer tests."""
    results = annotation_summary.copy()
    results["critic_flags"] = "PASS"
    results["critic_evidence"] = "Coverage: 5/5 (100%) positive markers expressed"
    results["critic_confidence"] = "high"
    results["critic_notes"] = ""
    return results
