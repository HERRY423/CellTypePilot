"""CellTypePilot — Literature integration via MCP or direct API.

Optional module for literature-backed marker validation.
Supports PubMed, bioRxiv, and Cell Ontology via MCP tools or direct HTTP.

This module is OPTIONAL — CellTypePilot works without it.
When MCP tools are available, they provide richer evidence.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

# ──────────────────────────────────────────────
# Data classes
# ──────────────────────────────────────────────


@dataclass
class LiteratureHit:
    """A literature search result."""

    title: str
    authors: str
    journal: str
    year: int
    pmid: str | None = None
    doi: str | None = None
    abstract_snippet: str = ""
    relevance_score: float = 0.0
    source: str = "pubmed"  # pubmed / biorxiv / cell_ontology

    def to_dict(self) -> dict:
        return {
            "title": self.title,
            "authors": self.authors,
            "journal": self.journal,
            "year": self.year,
            "pmid": self.pmid,
            "doi": self.doi,
            "snippet": self.abstract_snippet[:200] if self.abstract_snippet else "",
            "relevance": round(self.relevance_score, 3),
            "source": self.source,
        }


@dataclass
class MarkerLiteratureEvidence:
    """Literature evidence for a specific marker gene."""

    gene: str
    cell_type: str
    hits: list[LiteratureHit] = field(default_factory=list)
    total_refs: int = 0
    consensus: str = "unknown"  # supported / disputed / unknown

    def to_dict(self) -> dict:
        return {
            "gene": self.gene,
            "cell_type": self.cell_type,
            "total_refs": self.total_refs,
            "consensus": self.consensus,
            "top_hits": [h.to_dict() for h in self.hits[:3]],
        }


# ──────────────────────────────────────────────
# PubMed E--utilities (direct HTTP, no MCP needed)
# ──────────────────────────────────────────────

PUBMED_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BIORXIV_BASE = "https://api.biorxiv.org/details/biorxiv"


def search_pubmed(
    query: str,
    max_results: int = 5,
    email: str = "celltypepilot@example.com",
) -> list[LiteratureHit]:
    """Search PubMed via E-utilities (no API key needed for basic use).

    Args:
        query: Search query (e.g., "CD3E T cell marker")
        max_results: Max results to return
        email: Contact email (required by NCBI, can be any valid email)

    Returns:
        List of LiteratureHit objects
    """
    try:
        # Search
        search_params = {
            "db": "pubmed",
            "term": query,
            "retmax": str(max_results),
            "retmode": "json",
            "email": email,
        }
        search_url = f"{PUBMED_BASE}/esearch.fcgi?{urllib.parse.urlencode(search_params)}"

        req = urllib.request.Request(search_url, headers={"User-Agent": "CellTypePilot/0.3.0"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            search_data = json.loads(resp.read().decode("utf-8"))

        pmids = search_data.get("esearchresult", {}).get("idlist", [])
        if not pmids:
            return []

        # Fetch details
        fetch_params = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "json",
            "email": email,
        }
        fetch_url = f"{PUBMED_BASE}/efetch.fcgi?{urllib.parse.urlencode(fetch_params)}"

        req = urllib.request.Request(fetch_url, headers={"User-Agent": "CellTypePilot/0.3.0"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            # efetch returns XML, but we can parse basic info from esummary
            pass

        # Use esummary for simpler JSON response
        summary_params = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "json",
            "email": email,
        }
        summary_url = f"{PUBMED_BASE}/esummary.fcgi?{urllib.parse.urlencode(summary_params)}"

        req = urllib.request.Request(summary_url, headers={"User-Agent": "CellTypePilot/0.3.0"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            summary_data = json.loads(resp.read().decode("utf-8"))

        hits = []
        result = summary_data.get("result", {})
        for pmid in pmids:
            article = result.get(pmid, {})
            if not article:
                continue

            # Extract authors
            authors_list = article.get("authors", [])
            authors = ", ".join(a.get("name", "") for a in authors_list[:3])
            if len(authors_list) > 3:
                authors += " et al."

            # Extract title
            title = article.get("title", "Unknown title")

            # Extract journal
            source = article.get("source", "")
            pub_date = article.get("pubdate", "")
            year = int(pub_date[:4]) if pub_date and pub_date[:4].isdigit() else 0

            hits.append(
                LiteratureHit(
                    title=title,
                    authors=authors,
                    journal=source,
                    year=year,
                    pmid=pmid,
                    source="pubmed",
                    relevance_score=1.0,
                )
            )

        return hits

    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError, KeyError):
        return []


def search_biorxiv(
    query: str,
    max_results: int = 5,
) -> list[LiteratureHit]:
    """Search bioRxiv and medRxiv preprints via Europe PMC REST API.

    Args:
        query: Search query (e.g., "CD3E T cell marker")
        max_results: Max results to return

    Returns:
        List of LiteratureHit objects with source="biorxiv"
    """
    try:
        encoded_query = urllib.parse.quote(
            f'(PUBLISHER:"bioRxiv" OR PUBLISHER:"medRxiv") AND ({query})'
        )
        url = (
            f"https://www.ebi.ac.uk/europepmc/webservices/rest/search"
            f"?query=SRC:PPR%20AND%20{encoded_query}&format=json&pageSize={max_results}"
        )

        req = urllib.request.Request(
            url,
            headers={
                "User-Agent": "CellTypePilot/0.3.0 (https://github.com/HERRY423/CellTypePilot)"
            },
        )
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode("utf-8"))

        hits = []
        results = data.get("resultList", {}).get("result", [])
        for r in results:
            publisher = r.get("bookOrReportDetails", {}).get("publisher", "bioRxiv")
            authors_raw = r.get("authorString", "Unknown authors")
            authors_parts = [a.strip() for a in authors_raw.split(",") if a.strip()]
            if len(authors_parts) > 3:
                authors = ", ".join(authors_parts[:3]) + " et al."
            else:
                authors = authors_raw

            title = r.get("title", "Unknown title").rstrip(".")
            doi = r.get("doi")
            pub_year = r.get("pubYear")
            year = int(pub_year) if pub_year and str(pub_year).isdigit() else 0
            snippet = r.get("abstractText", "")

            hits.append(
                LiteratureHit(
                    title=title,
                    authors=authors,
                    journal=publisher,
                    year=year,
                    doi=doi,
                    abstract_snippet=snippet,
                    source="biorxiv",
                    relevance_score=1.0,
                )
            )

        return hits

    except Exception:
        return []


# ──────────────────────────────────────────────
# Marker-literature validation
# ──────────────────────────────────────────────


def validate_marker_in_literature(
    gene: str,
    cell_type: str,
    max_refs: int = 3,
    include_biorxiv: bool = True,
) -> MarkerLiteratureEvidence:
    """Search literature for evidence supporting a marker-cell_type association.

    Args:
        gene: Gene symbol (e.g., "CD3E")
        cell_type: Cell type name (e.g., "T cells")
        max_refs: Max references to retrieve
        include_biorxiv: Whether to include bioRxiv preprints if PubMed results are sparse

    Returns:
        MarkerLiteratureEvidence with search results
    """
    # Build query
    query = f"{gene} {cell_type} marker expression"

    hits = search_pubmed(query, max_results=max_refs)

    # If pubmed hits are sparse and biorxiv is enabled, query biorxiv
    if include_biorxiv and len(hits) < max_refs:
        biorxiv_hits = search_biorxiv(query, max_results=max_refs - len(hits))
        hits.extend(biorxiv_hits)

    evidence = MarkerLiteratureEvidence(
        gene=gene,
        cell_type=cell_type,
        hits=hits,
        total_refs=len(hits),
    )

    # Determine consensus
    if len(hits) >= 2:
        evidence.consensus = "supported"
    elif len(hits) == 1:
        evidence.consensus = "partially_supported"
    else:
        evidence.consensus = "unknown"

    return evidence


def validate_annotation_with_literature(
    cell_type: str,
    positive_markers: list[str],
    negative_markers: list[str] | None = None,
    max_refs_per_marker: int = 2,
) -> dict:
    """Validate an entire annotation against literature.

    Args:
        cell_type: Assigned cell type
        positive_markers: Markers that should be expressed
        negative_markers: Markers that should NOT be expressed (optional)
        max_refs_per_marker: Max refs per marker

    Returns:
        Dict with validation summary
    """
    if negative_markers is None:
        negative_markers = []

    positive_evidence = []
    for gene in positive_markers[:5]:  # Limit to top 5
        ev = validate_marker_in_literature(gene, cell_type, max_refs=max_refs_per_marker)
        positive_evidence.append(ev.to_dict())

    negative_evidence = []
    for gene in negative_markers[:3]:  # Limit to top 3
        ev = validate_marker_in_literature(gene, cell_type, max_refs=max_refs_per_marker)
        negative_evidence.append(ev.to_dict())

    # Summary
    total_pos_refs = sum(e["total_refs"] for e in positive_evidence)
    supported_markers = sum(1 for e in positive_evidence if e["consensus"] == "supported")

    return {
        "cell_type": cell_type,
        "positive_markers_checked": len(positive_evidence),
        "positive_markers_supported": supported_markers,
        "total_literature_refs": total_pos_refs,
        "positive_evidence": positive_evidence,
        "negative_evidence": negative_evidence,
        "overall_assessment": (
            "well_supported"
            if supported_markers >= len(positive_evidence) * 0.6
            else "partially_supported"
            if supported_markers > 0
            else "not_validated"
        ),
    }


# ──────────────────────────────────────────────
# MCP tool integration (for Claude Code / Codex)
# ──────────────────────────────────────────────


def generate_mcp_search_queries(
    cell_type: str,
    markers: list[str],
) -> list[str]:
    """Generate search queries for MCP tools (PubMed, bioRxiv).

    These queries can be passed to MCP search tools if available.

    Args:
        cell_type: Cell type name
        markers: List of marker genes

    Returns:
        List of search query strings
    """
    queries = []

    # General cell type marker query
    queries.append(f'"{cell_type}" marker genes single-cell RNA-seq')

    # Individual marker queries
    for gene in markers[:3]:
        queries.append(f'"{gene}" "{cell_type}" expression marker')

    # Validation query
    queries.append(f'"{cell_type}" identification scRNA-seq validation')

    return queries


def format_literature_for_report(
    literature_results: dict,
) -> str:
    """Format literature evidence for inclusion in HTML report.

    Args:
        literature_results: Output from validate_annotation_with_literature()

    Returns:
        HTML string for report section
    """
    if not literature_results or literature_results.get("total_literature_refs", 0) == 0:
        return "<p>No literature evidence found. Manual validation recommended.</p>"

    html_parts = []
    html_parts.append('<div class="literature-section">')
    html_parts.append("<h4>Literature Validation</h4>")
    html_parts.append(
        f"<p>Found <strong>{literature_results['total_literature_refs']}</strong> "
        f"references supporting {literature_results['positive_markers_supported']}/"
        f"{literature_results['positive_markers_checked']} positive markers.</p>"
    )

    assessment = literature_results.get("overall_assessment", "unknown")
    color = {
        "well_supported": "#2ecc71",
        "partially_supported": "#f39c12",
        "not_validated": "#e74c3c",
    }
    html_parts.append(
        f'<p>Assessment: <span style="color: {color.get(assessment, "#95a5a6")}; '
        f'font-weight: bold;">{assessment}</span></p>'
    )

    # Top references
    for ev in literature_results.get("positive_evidence", []):
        if ev["total_refs"] > 0:
            html_parts.append(
                f"<p><strong>{ev['gene']}</strong>: {ev['total_refs']} refs "
                f"(consensus: {ev['consensus']})</p>"
            )
            for hit in ev.get("top_hits", []):
                html_parts.append(
                    f'<p class="ref">- {hit["authors"]} ({hit["year"]}). '
                    f"<em>{hit['title']}</em>. {hit['journal']}. "
                    f"PMID: {hit.get('pmid', 'N/A')}</p>"
                )

    html_parts.append("</div>")
    return "\n".join(html_parts)


# ──────────────────────────────────────────────
# Check MCP availability
# ──────────────────────────────────────────────


def check_mcp_availability() -> dict:
    """Check which MCP tools are available.

    Returns:
        Dict with availability status
    """
    mcp_status = {
        "pubmed_direct": False,
        "pubmed_mcp": False,
        "biorxiv_mcp": False,
        "cell_ontology_mcp": False,
    }

    # Check direct PubMed access (always works if network available)
    try:
        test_query = "test"
        params = {"db": "pubmed", "term": test_query, "retmax": "1", "retmode": "json"}
        url = f"{PUBMED_BASE}/esearch.fcgi?{urllib.parse.urlencode(params)}"
        req = urllib.request.Request(url, headers={"User-Agent": "CellTypePilot/0.3.0"})
        with urllib.request.urlopen(req, timeout=5) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            if "esearchresult" in data:
                mcp_status["pubmed_direct"] = True
    except Exception:
        pass

    return mcp_status
