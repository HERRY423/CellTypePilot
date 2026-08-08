"""Tests for the live Cell Ontology service (offline fixtures only)."""

from __future__ import annotations

import pytest

from celltypepilot.ontology import (
    ONTOLOGY_ENV_VAR,
    OntologyError,
    OntologyService,
    check_atlas_ontology,
    load_ontology,
    ontology_cache_path,
    ontology_cache_status,
    parse_obo,
    summarize_findings,
)

OBO_FIXTURE = """format-version: 1.2
ontology: cl

[Term]
id: CL:0000084
name: T cell
synonym: "T lymphocyte" EXACT []
is_a: CL:0000542 ! lymphocyte

[Term]
id: CL:0000063
name: obsolete cell by histology
is_obsolete: true
consider: CL:0000000

[Term]
id: CL:1000488
name: cholangiocyte
synonym: "bile duct epithelial cell" EXACT []

[Typedef]
id: part_of
name: part of
"""


@pytest.fixture()
def obo_file(tmp_path):
    path = tmp_path / "cl.obo"
    path.write_text(OBO_FIXTURE, encoding="utf-8")
    return path


@pytest.fixture()
def service(obo_file) -> OntologyService:
    return OntologyService(terms=parse_obo(obo_file), source=str(obo_file))


def _atlas(cell_types: dict) -> dict:
    return {"version": "mkg-test", "tissues": {"blood": {"cell_types": cell_types}}}


# ──────────────────────────────────────────────
# OBO parsing
# ──────────────────────────────────────────────


def test_parse_obo_terms(obo_file):
    terms = parse_obo(obo_file)
    assert set(terms) == {"CL:0000084", "CL:0000063", "CL:1000488"}
    t_cell = terms["CL:0000084"]
    assert t_cell.name == "T cell"
    assert t_cell.synonyms == ["T lymphocyte"]
    assert t_cell.parents == ["CL:0000542"]
    obsolete = terms["CL:0000063"]
    assert obsolete.is_obsolete
    assert obsolete.consider == ["CL:0000000"]


def test_parse_obo_empty_raises(tmp_path):
    path = tmp_path / "empty.obo"
    path.write_text("format-version: 1.2\n\n[Typedef]\nid: part_of\n", encoding="utf-8")
    with pytest.raises(OntologyError):
        parse_obo(path)


# ──────────────────────────────────────────────
# Atlas checks
# ──────────────────────────────────────────────


def test_check_clean_atlas(service):
    atlas = _atlas({"T cell": {"cl_id": "CL:0000084"}})
    assert check_atlas_ontology(service, atlas) == []


def test_check_matches_synonym(service):
    atlas = _atlas({"T lymphocyte": {"cl_id": "CL:0000084"}})
    assert check_atlas_ontology(service, atlas) == []


def test_check_malformed_cl_id(service):
    atlas = _atlas({"T cell": {"cl_id": "CL:123"}})
    findings = check_atlas_ontology(service, atlas)
    assert len(findings) == 1
    assert findings[0]["severity"] == "error"
    assert "malformed" in findings[0]["issue"]


def test_check_unknown_cl_id(service):
    atlas = _atlas({"Mystery": {"cl_id": "CL:9999999"}})
    findings = check_atlas_ontology(service, atlas)
    assert len(findings) == 1
    assert findings[0]["severity"] == "error"
    assert "not found" in findings[0]["issue"]


def test_check_obsolete_cl_id(service):
    atlas = _atlas({"Cholangiocyte": {"cl_id": "CL:0000063"}})
    findings = check_atlas_ontology(service, atlas)
    assert len(findings) == 1
    assert findings[0]["severity"] == "error"
    assert "obsolete" in findings[0]["issue"]
    assert "CL:0000000" in findings[0]["issue"]


def test_check_label_mismatch_is_warning(service):
    atlas = _atlas({"Hepatic cholangiocyte": {"cl_id": "CL:1000488"}})
    findings = check_atlas_ontology(service, atlas)
    assert len(findings) == 1
    assert findings[0]["severity"] == "warning"


def test_check_walks_subtypes(service):
    atlas = _atlas(
        {
            "T cell": {
                "cl_id": "CL:0000084",
                "subtypes": {"Bad child": {"cl_id": "CL:9999999"}},
            }
        }
    )
    findings = check_atlas_ontology(service, atlas)
    assert len(findings) == 1
    assert findings[0]["path"] == "blood/T cell/Bad child"


def test_summarize_findings():
    findings = [
        {"severity": "error", "path": "a", "cl_id": "CL:1", "issue": "x"},
        {"severity": "warning", "path": "b", "cl_id": "CL:2", "issue": "y"},
    ]
    summary = summarize_findings(findings, checked_nodes=10)
    assert summary == {"checked_nodes": 10, "errors": 1, "warnings": 1, "ok": False}


# ──────────────────────────────────────────────
# Cache behaviour (env-var isolated, no network)
# ──────────────────────────────────────────────


def test_cache_status_when_missing(tmp_path, monkeypatch):
    monkeypatch.setenv(ONTOLOGY_ENV_VAR, str(tmp_path))
    status = ontology_cache_status()
    assert status["cached"] is False
    assert "ontology update" in status["detail"]


def test_load_ontology_fails_closed(tmp_path, monkeypatch):
    monkeypatch.setenv(ONTOLOGY_ENV_VAR, str(tmp_path))
    with pytest.raises(OntologyError, match="ontology update"):
        load_ontology()


def test_load_ontology_from_cache(tmp_path, monkeypatch, obo_file):
    monkeypatch.setenv(ONTOLOGY_ENV_VAR, str(tmp_path))
    target = ontology_cache_path()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(OBO_FIXTURE, encoding="utf-8")
    service = load_ontology()
    assert service.resolve("CL:0000084").name == "T cell"
    assert service.label_of("CL:1000488") == "cholangiocyte"
    assert service.label_of("CL:0000000") == ""
