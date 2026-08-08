"""Tests for the enhanced detection logic and the orchestrator module.

Covers:
- Species detection: Ensembl ID prefixes (multi-species), mixed naming
- Tissue detection: case-insensitive synonyms + keyword substring scan
- Orchestrator: pipeline execution, override application, helper functions
"""

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from celltypepilot.data_adapter import detect_species, detect_tissue, match_ensembl_species
from celltypepilot.orchestrator import (
    PipelineError,
    apply_overrides_to_h5ad,
    find_cluster_column,
    run_annotation_pipeline,
    write_annotations_to_adata,
)


def _make_adata(genes, obs_dict=None, n_cells=10):
    """Build a minimal AnnData with given var_names and obs columns."""
    adata = ad.AnnData(X=np.zeros((n_cells, len(genes)), dtype=np.float32))
    adata.var_names = list(genes)
    if obs_dict:
        adata.obs = pd.DataFrame(obs_dict, index=adata.obs_names)
    return adata


# ──────────────────────────────────────────────
# Species detection — Ensembl IDs
# ──────────────────────────────────────────────
class TestSpeciesEnsembl:
    @pytest.mark.parametrize(
        "prefix,expected",
        [
            ("ENSG", "human"),
            ("ENSMUSG", "mouse"),
            ("ENSRNOG", "rat"),
            ("ENSDARG", "zebrafish"),
            ("ENSGALG", "chicken"),
            ("ENSSSCG", "pig"),
            ("ENSBTAG", "cow"),
            ("ENSMMUG", "macaque"),
            ("ENSCAFG", "dog"),
        ],
    )
    def test_ensembl_prefix_species(self, prefix, expected):
        genes = [f"{prefix}{i:011d}" for i in range(100)]
        adata = _make_adata(genes)
        assert detect_species(adata) == expected

    def test_chicken_not_confused_with_human_prefix(self):
        # ENSGALG starts with ENSG — longest-prefix matching must win
        assert match_ensembl_species("ENSGALG00000000001") == "chicken"
        assert match_ensembl_species("ENSG00000000001") == "human"

    def test_match_ensembl_species_non_ensembl(self):
        assert match_ensembl_species("CD3D") is None
        assert match_ensembl_species("Cd3d") is None

    def test_mixed_ensembl_falls_back_to_symbols(self):
        # No single Ensembl prefix dominates → fall back to symbol conventions
        genes = (
            [f"ENSG{i:011d}" for i in range(30)]
            + [f"ENSMUSG{i:011d}" for i in range(30)]
            + [f"GENE{i}" for i in range(40)]
        )
        adata = _make_adata(genes)
        # Ambiguous symbols too → deterministic fallback to human
        assert detect_species(adata) == "human"

    def test_mixed_ensembl_with_dominant_species_symbols(self):
        # Ensembl prefixes mixed, but human symbols dominate the remainder
        genes = (
            [f"ENSG{i:011d}" for i in range(40)]
            + [f"ENSMUSG{i:011d}" for i in range(10)]
            + ["CD3D", "CD3E", "CD4", "CD8A", "MS4A1"] * 10
        )
        adata = _make_adata(genes)
        assert detect_species(adata) == "human"


class TestSpeciesSymbols:
    def test_ambiguous_mixed_symbols_defaults_human(self):
        # Balanced mix: human ALL-CAPS (incl. digits) vs mouse Title-case.
        # human_pattern = mouse_pattern = 5 → ambiguous → human fallback
        genes = ["CD3D", "IL7R", "FCER1A", "S100A8", "HLA-DRA"] + [
            "Cd3d",
            "Il7r",
            "Fcer1a",
            "Clec10a",
            "Itgax",
        ]
        adata = _make_adata(genes)
        assert detect_species(adata) == "human"

    def test_empty_var_names_defaults_human(self):
        adata = ad.AnnData(X=np.zeros((5, 0), dtype=np.float32))
        assert detect_species(adata) == "human"


# ──────────────────────────────────────────────
# Tissue detection — synonyms and case handling
# ──────────────────────────────────────────────
class TestTissueSynonyms:
    @pytest.mark.parametrize(
        "col,value",
        [
            ("Tissue", "blood"),  # capitalized
            ("TISSUE_TYPE", "lung"),  # upper snake case
            ("organ_system", "brain"),  # extended synonym
            ("Source", "pbmc"),  # capitalized variant synonym
            ("anatomy", "liver"),
            ("body_site", "skin"),
            ("anatomical_site", "gut"),
            ("sample_source", "kidney"),
        ],
    )
    def test_synonym_variants(self, col, value):
        adata = _make_adata(["CD3D", "CD3E"], obs_dict={col: [value] * 10})
        assert detect_tissue(adata) == value

    def test_keyword_substring_scan(self):
        # Column not in the synonym list but contains a keyword
        adata = _make_adata(["CD3D"], obs_dict={"tissue_of_origin": ["pancreas"] * 10})
        assert detect_tissue(adata) == "pancreas"

    def test_synonym_priority_over_keyword(self):
        # Exact synonym 'organ' should win over a keyword-matched column
        adata = _make_adata(
            ["CD3D"],
            obs_dict={"custom_tissue_info": ["lung"] * 10, "organ": ["heart"] * 10},
        )
        assert detect_tissue(adata) == "heart"

    def test_empty_and_nan_values_skipped(self):
        adata = _make_adata(
            ["CD3D"],
            obs_dict={"tissue": ["", np.nan, "blood", "blood"] + ["blood"] * 6},
        )
        assert detect_tissue(adata) == "blood"

    def test_all_empty_returns_none(self):
        adata = _make_adata(["CD3D"], obs_dict={"tissue": [""] * 10})
        assert detect_tissue(adata) is None

    def test_no_tissue_columns_returns_none(self):
        adata = _make_adata(["CD3D"], obs_dict={"condition": ["ctrl"] * 10})
        assert detect_tissue(adata) is None


# ──────────────────────────────────────────────
# Orchestrator — helpers
# ──────────────────────────────────────────────
class TestFindClusterColumn:
    def test_prefers_ctp_cl_id(self):
        obs = pd.DataFrame({"leiden": [0], "ctp_cl_id": ["0"]})
        assert find_cluster_column(obs) == "ctp_cl_id"

    def test_heuristic_cluster_name(self):
        obs = pd.DataFrame({"my_clusters": [0]})
        assert find_cluster_column(obs) == "my_clusters"

    def test_heuristic_cl_id(self):
        obs = pd.DataFrame({"sample_cl_id": [0]})
        assert find_cluster_column(obs) == "sample_cl_id"

    def test_none_when_absent(self):
        obs = pd.DataFrame({"condition": ["ctrl"]})
        assert find_cluster_column(obs) is None


# ──────────────────────────────────────────────
# Orchestrator — pipeline
# ──────────────────────────────────────────────
class TestRunAnnotationPipeline:
    def test_full_run_no_figures(self, synthetic_pbmc, tmp_output_dir):
        input_path = tmp_output_dir / "input.h5ad"
        synthetic_pbmc.write(str(input_path))
        out_dir = tmp_output_dir / "out"

        steps = []
        result = run_annotation_pipeline(
            input_path,
            "leiden",
            out_dir,
            no_figures=True,
            progress=lambda s, t, m: steps.append(s),
        )

        assert result["species"] == "human"
        assert result["tissue"] == "general"
        assert not result["critic_results"].empty
        assert steps == [1, 2, 3, 4, 5, 6, 8]  # step 7 is figures and was disabled
        for name in [
            "data.annotated.h5ad",
            "evidence_table.csv",
            "report_draft.html",
            "manifest.json",
            "methodology_draft.txt",
        ]:
            assert (out_dir / name).exists(), f"missing {name}"

        # Annotated obs columns written by write_annotations_to_adata
        adata = result["adata"]
        for col in ["ctp_cell_type", "ctp_cl_id", "ctp_confidence"]:
            assert col in adata.obs.columns
        assert {
            "ctp_cell_state_candidate",
            "ctp_state_decision",
            "ctp_cell_state",
            "ctp_display_label",
        } <= set(adata.obs.columns)

    def test_context_state_writeback_report_and_manifest_are_one_pipeline(
        self, synthetic_pbmc, tmp_output_dir
    ):
        input_path = tmp_output_dir / "input_context.h5ad"
        synthetic_pbmc.write(str(input_path))
        context_path = tmp_output_dir / "context.json"
        context_path.write_text(
            json.dumps(
                {
                    "schema_version": "celltypepilot.context.v1",
                    "species": "human",
                    "tissue": "general",
                    "condition": "synthetic acceptance condition",
                    "free_text": "Recorded for provenance, never expression evidence.",
                }
            ),
            encoding="utf-8",
        )

        result = run_annotation_pipeline(
            input_path,
            "leiden",
            tmp_output_dir / "context_out",
            species="human",
            tissue="general",
            context_file_path=context_path,
            no_figures=True,
        )
        params = result["manifest"]["parameters"]
        assert params["context_enabled"] is True
        assert len(params["context_sha256"]) == 64
        assert params["context_source_hashes"]["context_file_sha256"]
        assert params["state_contract"] == "identity_invariant_independent_axis_v1"
        assert "context_pack.normalized.json" in result["manifest"]["outputs"]
        assert "state_results.csv" in result["manifest"]["outputs"]
        assert result["paths"]["report"].exists()
        assert result["paths"]["annotated"].exists()

    def test_detected_tissue_from_obs(self, synthetic_pbmc, tmp_output_dir):
        synthetic_pbmc.obs["Tissue"] = "blood"
        input_path = tmp_output_dir / "input.h5ad"
        synthetic_pbmc.write(str(input_path))
        out_dir = tmp_output_dir / "out"

        result = run_annotation_pipeline(input_path, "leiden", out_dir, no_figures=True)
        assert result["tissue"] == "blood"

    def test_invalid_cluster_key_raises(self, synthetic_pbmc, tmp_output_dir):
        input_path = tmp_output_dir / "input.h5ad"
        synthetic_pbmc.write(str(input_path))
        with pytest.raises(PipelineError, match="cluster key"):
            run_annotation_pipeline(input_path, "no_such_key", tmp_output_dir / "out")

    def test_missing_input_raises(self, tmp_output_dir):
        with pytest.raises(FileNotFoundError):
            run_annotation_pipeline(
                tmp_output_dir / "missing.h5ad", "leiden", tmp_output_dir / "out"
            )


class TestWriteAnnotationsToAdata:
    def test_maps_clusters_and_fills_unknown(self, synthetic_pbmc, tmp_output_dir):
        critic_results = pd.DataFrame(
            {
                "cluster": ["T cell_0"],
                "cell_type": ["T cell"],
                "cl_id": ["CL:0000084"],
                "critic_confidence": ["high"],
            }
        )
        path = write_annotations_to_adata(synthetic_pbmc, critic_results, "leiden", tmp_output_dir)
        assert path.exists()
        mapped = synthetic_pbmc.obs[synthetic_pbmc.obs["leiden"] == "T cell_0"]
        assert (mapped["ctp_cell_type"] == "T cell").all()
        unmapped = synthetic_pbmc.obs[synthetic_pbmc.obs["leiden"] != "T cell_0"]
        assert (unmapped["ctp_cell_type"] == "Unknown").all()


# ──────────────────────────────────────────────
# Orchestrator — overrides
# ──────────────────────────────────────────────
def _make_annotated_h5ad(path):
    """Write a minimal annotated h5ad as produced by the pipeline."""
    n = 30
    adata = ad.AnnData(X=np.zeros((n, 5), dtype=np.float32))
    adata.var_names = [f"GENE_{i}" for i in range(5)]
    adata.obs = pd.DataFrame(
        {
            "ctp_cl_id": ["0"] * 20 + ["1"] * 10,
            "ctp_cell_type": ["T cell"] * 20 + ["B cell"] * 10,
            "ctp_confidence": ["high"] * 20 + ["medium"] * 10,
        },
        index=[f"cell_{i}" for i in range(n)],
    )
    adata.write(path)
    return adata


class TestApplyOverridesToH5ad:
    def test_apply_and_backup(self, tmp_output_dir):
        h5ad_path = tmp_output_dir / "data.annotated.h5ad"
        _make_annotated_h5ad(h5ad_path)

        overrides = {
            "0": {"new_type": "CD4 naive T cell", "reason": "manual review"},
            "1": {"new_type": "", "reason": "empty should be skipped"},
            "9": {"new_type": "Something", "reason": "unknown cluster"},
        }
        result = apply_overrides_to_h5ad(h5ad_path, overrides)

        assert result["applied"] == 1
        assert result["skipped"] == 2
        assert result["total"] == 3
        import os

        assert os.path.exists(result["backup"])

        updated = ad.read_h5ad(h5ad_path)
        mask = updated.obs["ctp_cl_id"] == "0"
        assert (updated.obs.loc[mask, "ctp_cell_type"] == "CD4 naive T cell").all()
        assert (updated.obs.loc[mask, "ctp_overridden"]).all()
        assert (updated.obs.loc[mask, "ctp_override_reason"] == "manual review").all()
        # Untouched cluster keeps its original label
        mask1 = updated.obs["ctp_cl_id"] == "1"
        assert (updated.obs.loc[mask1, "ctp_cell_type"] == "B cell").all()

    def test_detail_statuses(self, tmp_output_dir):
        h5ad_path = tmp_output_dir / "data.annotated.h5ad"
        _make_annotated_h5ad(h5ad_path)

        overrides = {
            "0": {"new_type": "X", "reason": ""},
            "1": {"new_type": "", "reason": ""},
            "42": {"new_type": "Y", "reason": ""},
        }
        result = apply_overrides_to_h5ad(h5ad_path, overrides)
        by_cluster = {d["cluster"]: d for d in result["details"]}
        assert by_cluster["0"]["status"] == "applied"
        assert by_cluster["1"]["status"] == "skipped"
        assert by_cluster["42"]["status"] == "skipped"
        assert by_cluster["42"]["reason"] == "No cells found"

    def test_missing_file_raises(self, tmp_output_dir):
        with pytest.raises(FileNotFoundError):
            apply_overrides_to_h5ad(tmp_output_dir / "nope.h5ad", {"0": {"new_type": "X"}})
