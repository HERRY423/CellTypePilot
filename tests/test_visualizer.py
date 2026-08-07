"""Tests for the visualizer module — figure generation and output."""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

from celltypepilot.visualizer import (
    generate_all_figures,
    plot_umap_clusters,
    plot_umap_celltype,
    plot_umap_confidence,
    plot_confidence_bar,
    _generate_colors,
)


class TestGenerateColors:
    def test_small_palette(self):
        colors = _generate_colors(5)
        assert len(colors) == 5

    def test_large_palette(self):
        colors = _generate_colors(50)
        assert len(colors) == 50

    def test_single_color(self):
        colors = _generate_colors(1)
        assert len(colors) == 1

    def test_colors_are_hex(self):
        colors = _generate_colors(10)
        for c in colors:
            assert c.startswith("#")
            assert len(c) == 7


class TestPlotUmapClusters:
    def test_basic_plot(self, synthetic_pbmc, tmp_output_dir):
        path = plot_umap_clusters(
            synthetic_pbmc, "leiden", "X_umap", tmp_output_dir
        )
        assert path is not None
        assert Path(path).exists()
        assert Path(path).stat().st_size > 0

    def test_missing_embedding_returns_none(self, synthetic_pbmc, tmp_output_dir):
        path = plot_umap_clusters(
            synthetic_pbmc, "leiden", "X_nonexistent", tmp_output_dir
        )
        assert path is None


class TestPlotUmapCelltype:
    def test_basic_plot(self, synthetic_pbmc, sample_critic_results, tmp_output_dir):
        path = plot_umap_celltype(
            synthetic_pbmc, sample_critic_results, "leiden", "X_umap", tmp_output_dir
        )
        assert path is not None
        assert Path(path).exists()

    def test_missing_embedding(self, synthetic_pbmc, sample_critic_results, tmp_output_dir):
        path = plot_umap_celltype(
            synthetic_pbmc, sample_critic_results, "leiden", "X_missing", tmp_output_dir
        )
        assert path is None


class TestPlotUmapConfidence:
    def test_basic_plot(self, synthetic_pbmc, sample_critic_results, tmp_output_dir):
        path = plot_umap_confidence(
            synthetic_pbmc, sample_critic_results, "leiden", "X_umap", tmp_output_dir
        )
        assert path is not None
        assert Path(path).exists()

    def test_with_mixed_confidence(self, synthetic_pbmc, tmp_output_dir):
        annotations = pd.DataFrame({
            "cluster": ["0", "1", "2", "3"],
            "cell_type": ["T cell", "B cell", "NK cell", "Monocyte"],
            "critic_confidence": ["high", "medium", "low", "needs_review"],
        })
        path = plot_umap_confidence(
            synthetic_pbmc, annotations, "leiden", "X_umap", tmp_output_dir
        )
        assert path is not None


class TestPlotConfidenceBar:
    def test_basic_plot(self, sample_critic_results, tmp_output_dir):
        path = plot_confidence_bar(sample_critic_results, tmp_output_dir)
        assert path is not None
        assert Path(path).exists()

    def test_no_confidence_column(self, tmp_output_dir):
        df = pd.DataFrame({"cluster": ["0"], "cell_type": ["T cell"]})
        path = plot_confidence_bar(df, tmp_output_dir)
        assert path is None


class TestGenerateAllFigures:
    def test_generates_figures(self, synthetic_pbmc, sample_critic_results, tmp_output_dir):
        paths = generate_all_figures(
            synthetic_pbmc, "leiden", "X_umap",
            sample_critic_results, tmp_output_dir, "blood"
        )
        assert len(paths) >= 2  # At least umap_clusters + confidence bar
        for p in paths:
            assert Path(p).exists()

    def test_no_embedding_key(self, synthetic_pbmc, sample_critic_results, tmp_output_dir):
        paths = generate_all_figures(
            synthetic_pbmc, "leiden", "X_missing",
            sample_critic_results, tmp_output_dir, "blood"
        )
        # Should still generate confidence bar
        assert isinstance(paths, list)
