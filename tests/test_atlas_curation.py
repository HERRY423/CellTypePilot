"""Tests for the literature co-occurrence curation sweep (offline mocks only)."""

from __future__ import annotations

import copy
import json

import pytest

from celltypepilot.atlas_curation import (
    AGGREGATE_STATUS,
    CURATOR_ID,
    CurationError,
    apply_sweep_results,
    build_edge_query,
    sweep_edges,
    write_sweep_report,
)
from celltypepilot.data_adapter import (
    LITERATURE_COOCCURRENCE_STATUS,
    validate_atlas_provenance,
)


def _source() -> dict:
    return {
        "source_id": "cellmarker_2_0",
        "name": "CellMarker 2.0",
        "pmid": "36300619",
        "doi": "10.1093/nar/gkac947",
        "url": "http://bioinfo.life.hust.edu.cn/CellMarker/",
        "source_type": "aggregate_database",
    }


def _record(gene: str, polarity: str, tissue: str = "blood", status: str = AGGREGATE_STATUS) -> dict:
    return {
        "gene": gene,
        "polarity": polarity,
        "species": ["human"],
        "tissue": tissue,
        "state": "baseline_or_unspecified",
        "atlas_version": "mkg-test",
        "sources": [_source()],
        "evidence_scope": "identity",
        "verification_status": status,
    }


def _atlas() -> dict:
    return {
        "version": "mkg-test",
        "tissues": {
            "blood": {
                "cell_types": {
                    "T cell": {
                        "cl_id": "CL:0000084",
                        "positive_markers": ["CD3E"],
                        "negative_markers": ["CD19"],
                        "marker_evidence": [
                            _record("CD3E", "positive"),
                            _record("CD19", "negative"),
                        ],
                        "subtypes": {},
                    }
                }
            },
            "lung": {
                "cell_types": {
                    "Ciliated cell": {
                        "cl_id": "CL:0000064",
                        "positive_markers": ["TPPP3"],
                        "negative_markers": [],
                        "marker_evidence": [_record("TPPP3", "positive", tissue="lung")],
                        "subtypes": {},
                    }
                }
            },
        },
    }


# ──────────────────────────────────────────────
# Query building
# ──────────────────────────────────────────────


def test_build_edge_query_positive():
    query = build_edge_query("CD3E", "positive", ("T cell", "CD4+ T cell"))
    assert '"CD3E"[tiab]' in query
    assert '"CD4+ T cell"[tiab]' in query
    assert "marker[tiab]" in query


def test_build_edge_query_negative():
    query = build_edge_query("CD19", "negative", ("T cell",))
    assert "negative OR absence" in query
    assert "marker[tiab]" not in query


# ──────────────────────────────────────────────
# sweep_edges with an injected searcher
# ──────────────────────────────────────────────


def test_sweep_edges_supports_and_skips_verified(monkeypatch):
    monkeypatch.setattr("celltypepilot.atlas_curation.time.sleep", lambda _s: None)
    calls: list[str] = []

    def searcher(query: str) -> list[str]:
        calls.append(query)
        if "CD3E" in query:
            return ["1", "2", "3"]
        return []

    sweep = sweep_edges(_atlas(), searcher=searcher, delay_seconds=0)
    # all three edges start at aggregate status
    assert sweep["swept"] == 3
    assert sweep["supported"] == 1
    assert sweep["errors"] == 0
    supported_rows = [row for row in sweep["results"] if row["supported"]]
    assert len(supported_rows) == 1
    row = supported_rows[0]
    assert row["gene"] == "CD3E"
    assert row["tissue"] == "blood"
    assert row["cell_path"] == "T cell"
    assert row["pmids"] == ["1", "2", "3"]


def test_sweep_edges_tissue_filter_and_limit(monkeypatch):
    monkeypatch.setattr("celltypepilot.atlas_curation.time.sleep", lambda _s: None)
    sweep = sweep_edges(_atlas(), tissue="blood", searcher=lambda q: ["1"], delay_seconds=0)
    assert sweep["swept"] == 2
    sweep = sweep_edges(_atlas(), searcher=lambda q: ["1"], delay_seconds=0, limit=1)
    assert sweep["swept"] == 1


def test_sweep_edges_tolerates_searcher_errors(monkeypatch):
    monkeypatch.setattr("celltypepilot.atlas_curation.time.sleep", lambda _s: None)

    def searcher(query: str) -> list[str]:
        if "CD19" in query:
            raise RuntimeError("boom")
        return ["1", "2"]

    sweep = sweep_edges(_atlas(), searcher=searcher, delay_seconds=0)
    assert sweep["swept"] == 3
    assert sweep["errors"] == 1
    failed = [row for row in sweep["results"] if row["error"]]
    assert len(failed) == 1
    assert failed[0]["gene"] == "CD19"
    assert not failed[0]["supported"]


def test_sweep_edges_min_hits_threshold(monkeypatch):
    monkeypatch.setattr("celltypepilot.atlas_curation.time.sleep", lambda _s: None)
    sweep = sweep_edges(_atlas(), searcher=lambda q: ["1"], min_hits=2, delay_seconds=0)
    assert sweep["supported"] == 0


# ──────────────────────────────────────────────
# apply_sweep_results
# ──────────────────────────────────────────────


def _sweep_rows(atlas: dict, searcher) -> list[dict]:
    return sweep_edges(atlas, searcher=searcher, delay_seconds=0)["results"]


def test_apply_upgrades_supported_edges():
    atlas = _atlas()
    rows = _sweep_rows(atlas, lambda q: ["11", "22"] if "CD3E" in q else [])
    updated, applied = apply_sweep_results(atlas, rows, "mkg-test.1", verified_at="2026-01-01T00:00:00Z")

    assert applied == 1
    assert validate_atlas_provenance(updated) == []
    records = updated["tissues"]["blood"]["cell_types"]["T cell"]["marker_evidence"]
    upgraded = next(r for r in records if r["gene"] == "CD3E")
    assert upgraded["verification_status"] == LITERATURE_COOCCURRENCE_STATUS
    assert upgraded["curator"] == CURATOR_ID
    assert upgraded["verified_at"] == "2026-01-01T00:00:00Z"
    assert "11,22" in upgraded["evidence_locator"]
    assert upgraded["atlas_version"] == "mkg-test.1"
    # untouched edges keep their aggregate status
    untouched = next(r for r in records if r["gene"] == "CD19")
    assert untouched["verification_status"] == AGGREGATE_STATUS
    assert untouched["atlas_version"] == "mkg-test.1"
    assert updated["version"] == "mkg-test.1"


def test_apply_does_not_mutate_input_atlas():
    atlas = _atlas()
    snapshot = copy.deepcopy(atlas)
    rows = _sweep_rows(atlas, lambda q: ["11", "22"])
    apply_sweep_results(atlas, rows, "mkg-test.1")
    assert atlas == snapshot


def test_apply_requires_new_version():
    atlas = _atlas()
    rows = _sweep_rows(atlas, lambda q: ["11", "22"])
    with pytest.raises(CurationError):
        apply_sweep_results(atlas, rows, "  ")


def test_apply_fails_closed_when_atlas_invalid():
    atlas = _atlas()
    # Introduce a pre-existing invalid record: evidence list mismatch
    atlas["tissues"]["lung"]["cell_types"]["Ciliated cell"]["positive_markers"].append("GHOST")
    rows = _sweep_rows(atlas, lambda q: ["11", "22"])
    with pytest.raises(CurationError, match="provenance validation"):
        apply_sweep_results(atlas, rows, "mkg-test.1")


def test_apply_with_no_supported_rows_still_bumps_version():
    atlas = _atlas()
    rows = _sweep_rows(atlas, lambda q: [])
    updated, applied = apply_sweep_results(atlas, rows, "mkg-test.1")
    assert applied == 0
    assert updated["version"] == "mkg-test.1"
    assert validate_atlas_provenance(updated) == []


# ──────────────────────────────────────────────
# Report writing
# ──────────────────────────────────────────────


def test_write_sweep_report(tmp_path):
    atlas = _atlas()
    sweep = sweep_edges(atlas, searcher=lambda q: ["1", "2"], delay_seconds=0)
    output = write_sweep_report(sweep, tmp_path / "nested" / "sweep.json")
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["curator"] == CURATOR_ID
    assert payload["swept"] == 3
    assert payload["supported"] == 3
    assert len(payload["results"]) == 3
