"""Tests for the literature module — data classes, query generation, formatting."""

from unittest.mock import patch, MagicMock
import json

import pytest

from celltypepilot.literature import (
    LiteratureHit,
    MarkerLiteratureEvidence,
    search_pubmed,
    search_biorxiv,
    validate_marker_in_literature,
    validate_annotation_with_literature,
    generate_mcp_search_queries,
    format_literature_for_report,
    check_mcp_availability,
    PUBMED_BASE,
    BIORXIV_BASE,
)


class TestLiteratureHit:
    def test_basic_creation(self):
        hit = LiteratureHit(
            title="Test Paper",
            authors="Smith J",
            journal="Nature",
            year=2025,
        )
        assert hit.title == "Test Paper"
        assert hit.year == 2025
        assert hit.source == "pubmed"
        assert hit.pmid is None

    def test_to_dict(self):
        hit = LiteratureHit(
            title="Paper", authors="Doe", journal="Science",
            year=2024, pmid="12345", relevance_score=0.95,
            abstract_snippet="A" * 300,
        )
        d = hit.to_dict()
        assert d["title"] == "Paper"
        assert d["pmid"] == "12345"
        assert d["relevance"] == 0.95
        assert len(d["snippet"]) == 200  # Truncated


class TestMarkerLiteratureEvidence:
    def test_basic_creation(self):
        ev = MarkerLiteratureEvidence(gene="CD3E", cell_type="T cell")
        assert ev.gene == "CD3E"
        assert ev.consensus == "unknown"
        assert ev.total_refs == 0

    def test_to_dict(self):
        hit = LiteratureHit(title="P", authors="A", journal="J", year=2024)
        ev = MarkerLiteratureEvidence(
            gene="CD3E", cell_type="T cell",
            hits=[hit], total_refs=1, consensus="supported",
        )
        d = ev.to_dict()
        assert d["gene"] == "CD3E"
        assert d["consensus"] == "supported"
        assert len(d["top_hits"]) == 1


class TestSearchPubmed:
    def test_network_error_returns_empty(self):
        import urllib.error
        with patch("celltypepilot.literature.urllib.request.urlopen") as mock:
            mock.side_effect = urllib.error.URLError("No network")
            hits = search_pubmed("CD3E T cell")
            assert hits == []

    def test_empty_results(self):
        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps({
            "esearchresult": {"idlist": []}
        }).encode("utf-8")
        mock_response.__enter__ = lambda s: s
        mock_response.__exit__ = MagicMock(return_value=False)

        with patch("celltypepilot.literature.urllib.request.urlopen",
                    return_value=mock_response):
            hits = search_pubmed("nonexistent_gene_xyz")
            assert hits == []


class TestSearchBiorxiv:
    def test_returns_list(self):
        hits = search_biorxiv("test query")
        assert isinstance(hits, list)

    def test_search_biorxiv_parsing(self):
        mock_response_data = {
            "resultList": {
                "result": [
                    {
                        "title": "Single-cell atlas of T cells.",
                        "authorString": "Alice A, Bob B, Charlie C, David D",
                        "bookOrReportDetails": {"publisher": "bioRxiv"},
                        "pubYear": "2025",
                        "doi": "10.1101/2025.01.01.123456",
                        "abstractText": "This is a preprint about T cells.",
                    }
                ]
            }
        }
        mock_resp = MagicMock()
        mock_resp.read.return_value = json.dumps(mock_response_data).encode("utf-8")
        mock_resp.__enter__ = lambda s: s
        mock_resp.__exit__ = MagicMock(return_value=False)

        with patch("celltypepilot.literature.urllib.request.urlopen", return_value=mock_resp):
            hits = search_biorxiv("CD3E T cell")
            assert len(hits) == 1
            hit = hits[0]
            assert hit.title == "Single-cell atlas of T cells"
            assert "Alice A" in hit.authors
            assert "et al." in hit.authors
            assert hit.journal == "bioRxiv"
            assert hit.year == 2025
            assert hit.doi == "10.1101/2025.01.01.123456"
            assert hit.source == "biorxiv"

    def test_network_error_returns_empty(self):
        import urllib.error
        with patch("celltypepilot.literature.urllib.request.urlopen") as mock:
            mock.side_effect = urllib.error.URLError("No network")
            hits = search_biorxiv("CD3E T cell")
            assert hits == []


class TestValidateMarkerInLiterature:
    def test_no_network_returns_unknown(self):
        with patch("celltypepilot.literature.search_pubmed", return_value=[]), \
             patch("celltypepilot.literature.search_biorxiv", return_value=[]):
            ev = validate_marker_in_literature("CD3E", "T cell")
            assert ev.consensus == "unknown"
            assert ev.total_refs == 0

    def test_single_hit_partially_supported(self):
        hit = LiteratureHit(title="P", authors="A", journal="J", year=2024)
        with patch("celltypepilot.literature.search_pubmed", return_value=[hit]), \
             patch("celltypepilot.literature.search_biorxiv", return_value=[]):
            ev = validate_marker_in_literature("CD3E", "T cell")
            assert ev.consensus == "partially_supported"
            assert ev.total_refs == 1

    def test_multiple_hits_supported(self):
        hits = [
            LiteratureHit(title="P1", authors="A", journal="J", year=2024),
            LiteratureHit(title="P2", authors="B", journal="N", year=2023),
        ]
        with patch("celltypepilot.literature.search_pubmed", return_value=hits), \
             patch("celltypepilot.literature.search_biorxiv", return_value=[]):
            ev = validate_marker_in_literature("CD3E", "T cell")
            assert ev.consensus == "supported"

    def test_biorxiv_fallback_supported(self):
        hit_pm = LiteratureHit(title="P1", authors="A", journal="J", year=2024)
        hit_bio = LiteratureHit(title="P2", authors="B", journal="bioRxiv", year=2025, source="biorxiv")
        with patch("celltypepilot.literature.search_pubmed", return_value=[hit_pm]), \
             patch("celltypepilot.literature.search_biorxiv", return_value=[hit_bio]):
            ev = validate_marker_in_literature("CD3E", "T cell", max_refs=2, include_biorxiv=True)
            assert ev.consensus == "supported"
            assert ev.total_refs == 2


class TestValidateAnnotationWithLiterature:
    def test_basic_validation(self):
        with patch("celltypepilot.literature.search_pubmed", return_value=[]):
            result = validate_annotation_with_literature(
                "T cell", ["CD3E", "CD3D"], ["MS4A1"]
            )
            assert "cell_type" in result
            assert result["cell_type"] == "T cell"
            assert result["positive_markers_checked"] == 2
            assert "overall_assessment" in result


class TestGenerateMcpQueries:
    def test_generates_queries(self):
        queries = generate_mcp_search_queries("T cell", ["CD3E", "CD3D"])
        assert len(queries) >= 2
        assert any("T cell" in q for q in queries)
        assert any("CD3E" in q for q in queries)

    def test_empty_markers(self):
        queries = generate_mcp_search_queries("T cell", [])
        assert len(queries) >= 1


class TestFormatLiteratureForReport:
    def test_empty_results(self):
        html = format_literature_for_report({"total_literature_refs": 0})
        assert "No literature evidence" in html

    def test_with_results(self):
        results = {
            "total_literature_refs": 5,
            "positive_markers_supported": 3,
            "positive_markers_checked": 4,
            "overall_assessment": "well_supported",
            "positive_evidence": [
                {
                    "gene": "CD3E",
                    "total_refs": 3,
                    "consensus": "supported",
                    "top_hits": [
                        {"authors": "Smith", "year": 2024, "title": "Study",
                         "journal": "Nature", "pmid": "123"}
                    ],
                }
            ],
        }
        html = format_literature_for_report(results)
        assert "CD3E" in html
        assert "well_supported" in html


class TestCheckMcpAvailability:
    def test_returns_dict(self):
        with patch("celltypepilot.literature.urllib.request.urlopen") as mock:
            mock.side_effect = Exception("No network")
            status = check_mcp_availability()
            assert isinstance(status, dict)
            assert "pubmed_direct" in status
