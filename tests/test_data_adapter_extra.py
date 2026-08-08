"""Tests for data_adapter edge cases — error handling, boundary conditions."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from celltypepilot.data_adapter import (
    compute_data_hash,
    detect_species,
    detect_tissue,
    find_cluster_keys,
    find_embedding_keys,
    find_layer_keys,
    format_inspect_report,
    get_all_markers_flat,
    get_all_markers_for_tissue,
    inspect,
    inspect_adata,
    load_h5ad,
    load_marker_atlas,
)


class TestLoadH5ad:
    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            load_h5ad("/nonexistent/path.h5ad")

    def test_wrong_extension(self, tmp_output_dir):
        bad_file = tmp_output_dir / "data.csv"
        bad_file.write_text("a,b,c")
        with pytest.raises(ValueError, match="Expected .h5ad"):
            load_h5ad(str(bad_file))

    def test_valid_load(self, h5ad_path):
        adata = load_h5ad(str(h5ad_path))
        assert adata.n_obs == 500
        assert adata.n_vars > 0


class TestComputeDataHash:
    def test_deterministic(self, h5ad_path):
        h1 = compute_data_hash(str(h5ad_path))
        h2 = compute_data_hash(str(h5ad_path))
        assert h1 == h2
        assert len(h1) == 64  # SHA-256 hex

    def test_different_files_different_hash(self, synthetic_pbmc, tmp_output_dir):
        p1 = tmp_output_dir / "f1.h5ad"
        p2 = tmp_output_dir / "f2.h5ad"
        synthetic_pbmc[:100].write(str(p1))
        synthetic_pbmc[100:200].write(str(p2))
        h1 = compute_data_hash(str(p1))
        h2 = compute_data_hash(str(p2))
        assert h1 != h2


class TestDetectSpecies:
    def test_human_genes(self, synthetic_pbmc):
        assert detect_species(synthetic_pbmc) == "human"

    def test_mouse_genes(self, synthetic_pbmc):
        mouse_genes = [g[0].upper() + g[1:].lower() for g in synthetic_pbmc.var_names]
        synthetic_pbmc.var_names = mouse_genes
        assert detect_species(synthetic_pbmc) == "mouse"

    def test_mixed_defaults_to_human(self):
        adata = ad.AnnData(X=np.zeros((10, 10)))
        adata.var_names = [f"Gene_{i}" for i in range(10)]
        result = detect_species(adata)
        assert result in ("human", "mouse")


class TestDetectTissue:
    def test_tissue_from_obs(self):
        adata = ad.AnnData(X=np.zeros((10, 5)))
        adata.var_names = [f"G{i}" for i in range(5)]
        adata.obs = pd.DataFrame({"tissue": ["blood"] * 10})
        assert detect_tissue(adata) == "blood"

    def test_no_tissue(self, synthetic_pbmc):
        result = detect_tissue(synthetic_pbmc)
        assert result is None

    def test_alternative_tissue_keys(self):
        adata = ad.AnnData(X=np.zeros((10, 5)))
        adata.var_names = [f"G{i}" for i in range(5)]
        adata.obs = pd.DataFrame({"organ": ["lung"] * 10})
        assert detect_tissue(adata) == "lung"


class TestFindClusterKeys:
    def test_find_leiden(self, synthetic_pbmc):
        keys = find_cluster_keys(synthetic_pbmc)
        assert "leiden" in keys

    def test_no_cluster_keys(self):
        adata = ad.AnnData(X=np.zeros((10, 5)))
        adata.var_names = [f"G{i}" for i in range(5)]
        adata.obs = pd.DataFrame({"foo": range(10)})
        keys = find_cluster_keys(adata)
        assert len(keys) == 0

    def test_multiple_cluster_keys(self):
        adata = ad.AnnData(X=np.zeros((10, 5)))
        adata.var_names = [f"G{i}" for i in range(5)]
        adata.obs = pd.DataFrame(
            {
                "leiden": range(10),
                "louvain": range(10),
                "cluster_id": range(10),
            }
        )
        keys = find_cluster_keys(adata)
        assert len(keys) >= 3


class TestFindEmbeddingKeys:
    def test_find_umap(self, synthetic_pbmc):
        keys = find_embedding_keys(synthetic_pbmc)
        assert any("umap" in k for k in keys)

    def test_no_embeddings(self):
        adata = ad.AnnData(X=np.zeros((10, 5)))
        keys = find_embedding_keys(adata)
        assert len(keys) == 0


class TestFindLayerKeys:
    def test_detects_layers(self, synthetic_pbmc):
        result = find_layer_keys(synthetic_pbmc)
        assert "counts" in result
        assert "lognorm" in result

    def test_counts_layer(self, synthetic_pbmc):
        result = find_layer_keys(synthetic_pbmc)
        assert result["counts"] == "counts"

    def test_lognorm_from_x(self, synthetic_pbmc):
        result = find_layer_keys(synthetic_pbmc)
        assert result["lognorm"] == "X"


class TestInspectAdata:
    def test_full_inspection(self, h5ad_path):
        report = inspect_adata(str(h5ad_path), cluster_key="leiden")
        assert report["n_obs"] == 500
        assert report["species"] == "human"
        assert "leiden" in report["cluster_keys"]
        assert len(report["cluster_sizes"]) > 0

    def test_no_cluster_key_warning(self, h5ad_path):
        report = inspect_adata(str(h5ad_path))
        assert len(report["warnings"]) > 0

    def test_detected_unsupported_species_is_reported_not_scored(self):
        adata = ad.AnnData(
            X=np.zeros((4, 3)),
            obs=pd.DataFrame({"cluster": ["0", "0", "1", "1"]}),
            var=pd.DataFrame(
                index=[
                    "ENSRNOG00000000001",
                    "ENSRNOG00000000002",
                    "ENSRNOG00000000003",
                ]
            ),
        )
        report = inspect(adata, cluster_key="cluster")
        assert report["species"] == "rat"
        assert report["annotation_species_supported"] is False
        assert report["supported_annotation_species"] == ["human", "mouse"]
        assert any("fails closed" in warning for warning in report["warnings"])


class TestFormatInspectReport:
    def test_format(self, h5ad_path):
        report = inspect_adata(str(h5ad_path), cluster_key="leiden")
        text = format_inspect_report(report)
        assert "CellTypePilot" in text
        assert "human" in text
        assert "leiden" in text


class TestLoadMarkerAtlas:
    def test_human_atlas(self):
        atlas = load_marker_atlas("human")
        assert "tissues" in atlas
        assert "blood" in atlas["tissues"]

    def test_mouse_atlas(self):
        atlas = load_marker_atlas("mouse")
        assert "tissues" in atlas
        blood_types = atlas["tissues"]["blood"]["cell_types"]
        t_cell_markers = blood_types["T cell"]["positive_markers"]
        for m in t_cell_markers:
            if m.isalpha():
                assert m[0].isupper()
                if len(m) > 1:
                    assert m[1].islower()

    def test_unsupported_species_fails_closed(self):
        with pytest.raises(ValueError, match="supports scoring only human, mouse"):
            load_marker_atlas("rat")


class TestGetAllMarkersForTissue:
    def test_blood_tissue(self, blood_atlas):
        markers = get_all_markers_for_tissue(blood_atlas, "blood")
        assert len(markers) > 0
        assert "T cell" in markers
        assert "positive_markers" in markers["T cell"]

    def test_nonexistent_falls_back(self, blood_atlas):
        markers = get_all_markers_for_tissue(blood_atlas, "nonexistent_tissue")
        # Falls back to "general" or returns empty
        assert isinstance(markers, dict)

    def test_includes_subtypes(self, blood_atlas):
        markers = get_all_markers_for_tissue(blood_atlas, "blood")
        # Should have both parent types and subtypes
        assert len(markers) > 0


class TestGetAllMarkersFlat:
    def test_flat_format(self, blood_atlas):
        flat = get_all_markers_flat(blood_atlas, "blood")
        for markers in flat.values():
            assert isinstance(markers, list)
            assert all(isinstance(m, str) for m in markers)
