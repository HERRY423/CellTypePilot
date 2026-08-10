"""Unit and integration tests for the Atlas Governance System (Phase A)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from celltypepilot.atlas_conflict import (
    detect_marker_conflicts,
)
from celltypepilot.atlas_curation import (
    build_curation_queue,
)
from celltypepilot.atlas_diff import diff_atlases, format_diff_json, format_diff_report
from celltypepilot.atlas_governance import build_atlas_governance_report
from celltypepilot.atlas_lifecycle import (
    apply_deprecation_filter,
    compute_marker_weights,
    deprecate_cell_type,
    deprecate_marker_edge,
    sunset_check,
)
from celltypepilot.data_adapter import get_all_markers_for_tissue
from celltypepilot.marker_scorer import compute_marker_scores


@pytest.fixture
def sample_atlas():
    return {
        "version": "mkg-2026.08.1",
        "description": "Test marker atlas",
        "sources": ["PanglaoDB"],
        "tissues": {
            "blood": {
                "name": "Blood",
                "cell_types": {
                    "T cell": {
                        "cl_id": "CL:0000084",
                        "positive_markers": ["CD3D", "CD3E", "CD3G"],
                        "negative_markers": ["CD19"],
                        "marker_evidence": [
                            {
                                "gene": "CD3D",
                                "polarity": "positive",
                                "verification_status": "primary_source_verified",
                            },
                            {
                                "gene": "CD3E",
                                "polarity": "positive",
                                "verification_status": "literature_cooccurrence_supported",
                            },
                            {
                                "gene": "CD3G",
                                "polarity": "positive",
                                "verification_status": "aggregate_source_only_not_edge_verified",
                            },
                            {
                                "gene": "CD19",
                                "polarity": "negative",
                                "verification_status": "primary_source_verified",
                            },
                        ],
                        "subtypes": {
                            "CD4+ T cell": {
                                "cl_id": "CL:0000624",
                                "positive_markers": ["CD3D", "CD3E", "CD4"],
                                "negative_markers": ["CD8A"],
                                "marker_evidence": [
                                    {
                                        "gene": "CD4",
                                        "polarity": "positive",
                                        "verification_status": "database_record_verified",
                                    }
                                ],
                            }
                        },
                    },
                    "B cell": {
                        "cl_id": "CL:0000236",
                        "positive_markers": ["CD19", "MS4A1"],
                        "negative_markers": ["CD3D"],
                        "marker_evidence": [
                            {
                                "gene": "CD19",
                                "polarity": "positive",
                                "verification_status": "primary_source_verified",
                            },
                            {
                                "gene": "MS4A1",
                                "polarity": "positive",
                                "verification_status": "primary_source_verified",
                            },
                        ],
                    },
                },
            }
        },
    }


def test_compute_marker_weights(sample_atlas):
    records = sample_atlas["tissues"]["blood"]["cell_types"]["T cell"]["marker_evidence"]
    weights = compute_marker_weights(records)
    assert weights["CD3D"] == 1.0  # primary = rank 3 -> 1.0
    assert weights["CD3E"] == 0.7  # literature = rank 1 -> 0.7
    assert weights["CD3G"] == 0.5  # aggregate = rank 0 -> 0.5
    assert weights["CD19"] == 1.0  # primary = rank 3 -> 1.0


def test_deprecate_cell_type(sample_atlas):
    deprecate_cell_type(
        sample_atlas,
        cell_type="B cell",
        tissue="blood",
        reason="Obsolete nomenclature",
        superseded_by="B lymphocyte",
    )
    b_cell = sample_atlas["tissues"]["blood"]["cell_types"]["B cell"]
    assert b_cell["deprecated"] is True
    assert b_cell["deprecation_reason"] == "Obsolete nomenclature"
    assert b_cell["superseded_by"] == "B lymphocyte"


def test_deprecate_marker_edge(sample_atlas):
    deprecate_marker_edge(
        sample_atlas,
        gene="CD3G",
        cell_type="T cell",
        tissue="blood",
        reason="Non-specific expression in myeloid lineage",
    )
    records = sample_atlas["tissues"]["blood"]["cell_types"]["T cell"]["marker_evidence"]
    cd3g = next(r for r in records if r["gene"] == "CD3G")
    assert cd3g["deprecated"] is True
    assert cd3g["deprecation_reason"] == "Non-specific expression in myeloid lineage"


def test_sunset_check(sample_atlas):
    deprecate_cell_type(
        sample_atlas,
        cell_type="B cell",
        tissue="blood",
        reason="Obsolete",
    )
    sample_atlas["tissues"]["blood"]["cell_types"]["B cell"]["sunset_version"] = 2.0

    issues = sunset_check(sample_atlas, current_version=1.0)
    assert len(issues) == 0  # not sunsetted yet at 1.0

    issues_sunset = sunset_check(sample_atlas, current_version=2.0)
    assert len(issues_sunset) == 1
    assert "sunset version" in issues_sunset[0]


def test_apply_deprecation_filter(sample_atlas):
    deprecate_cell_type(
        sample_atlas,
        cell_type="B cell",
        tissue="blood",
        reason="Obsolete",
    )
    markers = get_all_markers_for_tissue(sample_atlas, "blood", include_deprecated=True)

    filtered = apply_deprecation_filter(markers, include_deprecated=False)
    assert "B cell" not in filtered
    assert "T cell" in filtered

    with_deprecated = apply_deprecation_filter(markers, include_deprecated=True)
    assert "B cell" in with_deprecated


def test_atlas_diff(sample_atlas):
    import copy

    new_atlas = copy.deepcopy(sample_atlas)
    deprecate_cell_type(new_atlas, "B cell", "blood", "Obsolete")

    diff = diff_atlases(sample_atlas, new_atlas)
    assert "blood/B cell" in diff.deprecated_cell_types

    report_text = format_diff_report(diff)
    assert "Deprecated Cell Types" in report_text

    json_data = format_diff_json(diff)
    assert "added_cell_types" in json_data


def test_atlas_conflict_detection(sample_atlas):
    import copy

    bad_atlas = copy.deepcopy(sample_atlas)
    # Inject polarity conflict: CD3D positive in T cell and negative in CD4+ T cell
    bad_atlas["tissues"]["blood"]["cell_types"]["T cell"]["subtypes"]["CD4+ T cell"][
        "negative_markers"
    ].append("CD3D")

    conflicts = detect_marker_conflicts(bad_atlas)
    assert len(conflicts) > 0
    polarity_conflicts = [c for c in conflicts if c.conflict_type == "polarity"]
    assert len(polarity_conflicts) > 0
    assert polarity_conflicts[0].gene == "CD3D"


def test_build_curation_queue(sample_atlas):
    queue = build_curation_queue(sample_atlas)
    assert isinstance(queue, pd.DataFrame)
    assert not queue.empty
    assert "gene" in queue.columns
    assert "priority" in queue.columns
    # Rank 0 (aggregate) should have highest priority (100)
    top_row = queue.iloc[0]
    assert top_row["priority"] == 100


def test_governance_report_integration():
    report = build_atlas_governance_report(include_packs=False)
    assert report["schema_version"] == "celltypepilot.atlas-governance.v1"
    assert "governance_health_score" in report
    assert 0.0 <= report["governance_health_score"] <= 100.0
    assert "needs_edge_curation" in report["aggregate"]


def test_evidence_weighted_scoring():
    # Build synthetic AnnData for test
    import anndata as ad

    X = np.array(
        [
            [10.0, 8.0, 0.0, 0.0],
            [9.0, 7.0, 0.0, 0.0],
            [0.0, 0.0, 9.0, 8.0],
            [0.0, 0.0, 8.0, 7.0],
        ]
    )
    obs = pd.DataFrame({"cluster": ["0", "0", "1", "1"]})
    var = pd.DataFrame(index=["CD3D", "CD3E", "CD19", "MS4A1"])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    markers = {
        "T cell": {
            "positive_markers": ["CD3D", "CD3E"],
            "negative_markers": ["CD19"],
            "marker_evidence": [
                {"gene": "CD3D", "verification_status": "primary_source_verified"},
                {"gene": "CD3E", "verification_status": "aggregate_source_only_not_edge_verified"},
            ],
        },
        "B cell": {
            "positive_markers": ["CD19", "MS4A1"],
            "negative_markers": ["CD3D"],
            "marker_evidence": [
                {"gene": "CD19", "verification_status": "primary_source_verified"},
                {"gene": "MS4A1", "verification_status": "primary_source_verified"},
            ],
        },
    }

    scores_unweighted = compute_marker_scores(adata, "cluster", markers, evidence_weighted=False)
    scores_weighted = compute_marker_scores(adata, "cluster", markers, evidence_weighted=True)

    assert not scores_unweighted.empty
    assert not scores_weighted.empty
    assert "combined_score" in scores_weighted.columns
