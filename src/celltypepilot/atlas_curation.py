"""Automated literature co-occurrence curation for bundled marker atlases.

This module turns the published curation backlog into an executable,
auditable sweep: for every marker relationship still at aggregate-source
provenance, it queries PubMed for the gene + cell-type pair and, when
independent co-occurrence evidence exists, upgrades the relationship to
``literature_cooccurrence_supported`` with a full evidence locator.

Honesty contract (what this sweep does and does NOT claim):

- Co-occurrence in PubMed titles/abstracts is *supporting evidence that the
  association is discussed in the literature*. It is NOT a claim of
  marker specificity, functional validation, or primary-experiment review.
  The dedicated status name makes that boundary explicit.
- Edges without sufficient co-occurrence are left untouched (fail closed).
- Every upgrade records the exact query, hit count, top PMIDs, timestamp,
  and curator identity, so the sweep is replayable and auditable.
"""

from __future__ import annotations

import json
import time
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

from .data_adapter import (
    LITERATURE_COOCCURRENCE_STATUS,
    validate_atlas_provenance,
    EVIDENCE_STATUS_RANK
)
from .constants import CURATION_STALE_MONTHS
import pandas as pd

AGGREGATE_STATUS = "aggregate_source_only_not_edge_verified"
CURATOR_ID = "celltypepilot-pubmed-sweep-v1"
MIN_HITS_FOR_UPGRADE = 2
MAX_LOCATOR_PMIDS = 3
DEFAULT_RATE_LIMIT_SECONDS = 0.34  # NCBI: <= 3 requests/s without an API key


class CurationError(ValueError):
    """Raised when a curation sweep or apply step is malformed."""


def _iter_edges(atlas: dict):
    """Yield (tissue, cell_path, cell_info, record) for every evidence record."""

    def walk(cell_types: dict, tissue: str, parents: tuple[str, ...]):
        for name, info in cell_types.items():
            path = (*parents, name)
            yield tissue, path, info
            yield from walk(info.get("subtypes", {}), tissue, path)

    for tissue, tissue_info in atlas.get("tissues", {}).items():
        yield from walk(tissue_info.get("cell_types", {}), tissue, ())


def _display_cell_type(cell_path: tuple[str, ...]) -> str:
    return cell_path[-1].replace("_", " ")


def build_edge_query(gene: str, polarity: str, cell_path: tuple[str, ...]) -> str:
    """Build a conservative PubMed co-occurrence query for one edge."""
    cell_type = _display_cell_type(cell_path)
    if polarity == "negative":
        return f'"{gene}"[tiab] AND "{cell_type}"[tiab] AND (negative OR absence)'
    return f'"{gene}"[tiab] AND "{cell_type}"[tiab] AND marker[tiab]'


def default_searcher(query: str) -> list[str]:
    """PubMed search returning PMIDs, using the shared literature client."""
    from .literature import search_pubmed

    hits = search_pubmed(query, max_results=20)
    return [hit.pmid for hit in hits if hit.pmid]


def sweep_edges(
    atlas: dict,
    tissue: str | None = None,
    searcher=None,
    min_hits: int = MIN_HITS_FOR_UPGRADE,
    delay_seconds: float = DEFAULT_RATE_LIMIT_SECONDS,
    limit: int | None = None,
    progress=None,
) -> dict:
    """Query literature co-occurrence for every aggregate-status edge.

    Returns ``{"results": [...], "swept", "supported", "errors"}`` where each
    result row carries the atlas path coordinates plus ``supported`` and the
    locator fields. Never mutates the atlas.
    """
    searcher = searcher or default_searcher
    results: list[dict] = []
    errors = 0
    swept = 0
    edges = [
        (t, path, info, record)
        for t, path, info, record in _all_records(atlas)
        if record.get("verification_status") == AGGREGATE_STATUS
        and (tissue is None or t == tissue)
    ]
    if limit is not None:
        edges = edges[:limit]
    total = len(edges)
    for index, (edge_tissue, cell_path, info, record) in enumerate(edges, start=1):
        if progress is not None:
            progress(index, total, f"{record.get('gene')} × {cell_path[-1]} ({edge_tissue})")
        query = build_edge_query(str(record.get("gene", "")), str(record.get("polarity", "")), cell_path)
        row = {
            "tissue": edge_tissue,
            "cell_path": " > ".join(cell_path),
            "gene": record.get("gene", ""),
            "polarity": record.get("polarity", ""),
            "query": query,
            "supported": False,
            "total_hits": 0,
            "pmids": [],
            "error": "",
        }
        try:
            pmids = searcher(query)
            row["total_hits"] = len(pmids)
            row["pmids"] = pmids[:MAX_LOCATOR_PMIDS]
            row["supported"] = len(pmids) >= min_hits
        except Exception as exc:  # keep sweeping through transient failures
            errors += 1
            row["error"] = str(exc)
        results.append(row)
        swept += 1
        if delay_seconds > 0 and index < total:
            time.sleep(delay_seconds)
    supported = sum(1 for row in results if row["supported"])
    return {"results": results, "swept": swept, "supported": supported, "errors": errors}


def _all_records(atlas: dict):
    for edge_tissue, cell_path, info in _iter_edges(atlas):
        for record in info.get("marker_evidence", []):
            yield edge_tissue, cell_path, info, record


def apply_sweep_results(
    atlas: dict,
    results: list[dict],
    new_version: str,
    verified_at: str | None = None,
) -> tuple[dict, int]:
    """Return (updated_atlas, upgrades_applied) for supported sweep rows.

    Fails closed: the merged atlas must pass validate_atlas_provenance,
    otherwise a CurationError is raised and nothing is returned.
    """
    if not str(new_version).strip():
        raise CurationError("new_version is required to apply sweep results")
    supported_index = {
        (row["tissue"], row["cell_path"], row["gene"], row["polarity"]): row
        for row in results
        if row.get("supported")
    }
    updated = deepcopy(atlas)
    stamp = verified_at or datetime.now(timezone.utc).isoformat()
    applied = 0
    for edge_tissue, cell_path, info in _iter_edges(updated):
        path_key = " > ".join(cell_path)
        for record in info.get("marker_evidence", []):
            key = (edge_tissue, path_key, record.get("gene"), record.get("polarity"))
            row = supported_index.get(key)
            if row is None or record.get("verification_status") != AGGREGATE_STATUS:
                continue
            record["verification_status"] = LITERATURE_COOCCURRENCE_STATUS
            record["evidence_locator"] = (
                f"PubMed query {row['query']!r}: {row['total_hits']} co-occurrence hits; "
                f"top PMIDs {','.join(row['pmids'])}"
            )
            record["curator"] = CURATOR_ID
            record["verified_at"] = stamp
            applied += 1
    updated["version"] = new_version
    for _t, _path, info in _iter_edges(updated):
        for record in info.get("marker_evidence", []):
            record["atlas_version"] = new_version
    issues = validate_atlas_provenance(updated)
    if issues:
        raise CurationError(
            "Atlas failed provenance validation after applying sweep: " + "; ".join(issues[:5])
        )
    return updated, applied


def write_sweep_report(sweep: dict, output_path: str | Path) -> Path:
    """Persist sweep results as a JSON report for audit and re-application."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "curator": CURATOR_ID,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "swept": sweep["swept"],
        "supported": sweep["supported"],
        "errors": sweep["errors"],
        "min_hits_for_upgrade": MIN_HITS_FOR_UPGRADE,
        "results": sweep["results"],
    }
    output_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return output_path


def build_curation_queue(atlas: dict) -> pd.DataFrame:
    """Build a prioritized queue of edges needing curation."""
    from datetime import datetime, timezone
    
    rows = []
    now = datetime.now(timezone.utc)
    
    for edge_tissue, cell_path, info, record in _all_records(atlas):
        status = record.get("verification_status", AGGREGATE_STATUS)
        rank = EVIDENCE_STATUS_RANK.get(status, 0)
        
        # Priority logic: unverified edges are highest priority, stale edges next
        priority = 0
        if rank == 0:
            priority = 100
        else:
            verified_at = record.get("verified_at")
            if verified_at:
                try:
                    v_time = datetime.fromisoformat(verified_at.replace("Z", "+00:00"))
                    months_old = (now - v_time).days / 30.0
                    if months_old > CURATION_STALE_MONTHS:
                        priority = 50 + min(months_old, 49)
                except Exception:
                    priority = 75  # Unparseable date

        rows.append({
            "tissue": edge_tissue,
            "cell_type": " > ".join(cell_path),
            "gene": record.get("gene"),
            "polarity": record.get("polarity"),
            "status": status,
            "rank": rank,
            "priority": priority,
            "verified_at": record.get("verified_at")
        })
        
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("priority", ascending=False)
    return df


def promote_edge(atlas: dict, gene: str, cell_type: str, tissue: str, new_status: str, evidence: dict) -> dict:
    """Upgrade an edge's status and add audit trails."""
    from datetime import datetime, timezone
    
    for edge_tissue, cell_path, info, record in _all_records(atlas):
        if edge_tissue == tissue and cell_path[-1] == cell_type and record.get("gene") == gene:
            old_status = record.get("verification_status", AGGREGATE_STATUS)
            if EVIDENCE_STATUS_RANK.get(new_status, 0) <= EVIDENCE_STATUS_RANK.get(old_status, 0):
                raise CurationError("New status must be higher rank than old status")
                
            record["verification_status"] = new_status
            record["curator"] = evidence.get("curator", "manual")
            record["verified_at"] = datetime.now(timezone.utc).isoformat()
            
            if new_status == "primary_source_verified":
                record["evidence_locator"] = evidence.get("evidence_locator", "")
                if "sources" in evidence:
                    record["sources"] = evidence["sources"]
            elif new_status == "database_record_verified":
                record["source_record_id"] = evidence.get("source_record_id", "")
                record["source_record_url"] = evidence.get("source_record_url", "")
                
            return record
            
    raise CurationError(f"Edge not found: {gene} in {cell_type} ({tissue})")


def demote_edge(atlas: dict, gene: str, cell_type: str, tissue: str, new_status: str, reason: str) -> dict:
    """Downgrade an edge's status and add audit trails."""
    from datetime import datetime, timezone
    
    for edge_tissue, cell_path, info, record in _all_records(atlas):
        if edge_tissue == tissue and cell_path[-1] == cell_type and record.get("gene") == gene:
            old_status = record.get("verification_status", AGGREGATE_STATUS)
            if EVIDENCE_STATUS_RANK.get(new_status, 0) >= EVIDENCE_STATUS_RANK.get(old_status, 0):
                raise CurationError("New status must be lower rank than old status")
                
            record["verification_status"] = new_status
            record["curator"] = "manual_demotion"
            record["verified_at"] = datetime.now(timezone.utc).isoformat()
            record["demotion_reason"] = reason
            
            return record
            
    raise CurationError(f"Edge not found: {gene} in {cell_type} ({tissue})")
