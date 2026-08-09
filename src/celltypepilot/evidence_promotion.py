"""Human-gated marker-edge evidence promotion.

Automated searches and data-derived stability analyses may create proposals,
but they cannot promote atlas evidence. Promotion to a database-record or
primary-source status requires two distinct approving reviewers and a locked
evidence locator. Application returns a new atlas object and never mutates the
source atlas in place.
"""

from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .data_adapter import EVIDENCE_STATUS_RANK, validate_atlas_provenance

PROMOTION_SCHEMA = "celltypepilot.evidence-promotion.v1"
REVIEW_DECISIONS = {"approve", "reject", "needs_revision"}
HUMAN_GATED_STATUSES = {"database_record_verified", "primary_source_verified"}


class EvidencePromotionError(ValueError):
    """Raised when a promotion proposal violates the evidence contract."""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _canonical_hash(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _iter_records(atlas: dict):
    def walk(tissue: str, nodes: dict, parents: tuple[str, ...] = ()):
        for name, info in nodes.items():
            path = (*parents, name)
            for record in info.get("marker_evidence", []):
                yield tissue, path, record
            yield from walk(tissue, info.get("subtypes", {}), path)

    for tissue, info in atlas.get("tissues", {}).items():
        yield from walk(tissue, info.get("cell_types", {}))


def _find_record(atlas: dict, tissue: str, cell_path: str, gene: str, polarity: str):
    for edge_tissue, path, record in _iter_records(atlas):
        if (
            edge_tissue == tissue
            and " > ".join(path) == cell_path
            and str(record.get("gene")) == gene
            and str(record.get("polarity")) == polarity
        ):
            return record
    raise EvidencePromotionError(
        f"Marker edge not found: {tissue}/{cell_path}/{gene}/{polarity}"
    )


def build_promotion_proposal(
    atlas: dict,
    *,
    tissue: str,
    cell_path: str,
    gene: str,
    polarity: str,
    target_status: str,
    evidence: dict,
    requested_by: str,
    proposal_origin: str = "human_curator",
) -> dict:
    """Create a locked proposal without changing the atlas."""
    record = _find_record(atlas, tissue, cell_path, gene, polarity)
    current = str(record.get("verification_status", ""))
    if target_status not in EVIDENCE_STATUS_RANK:
        raise EvidencePromotionError(f"Unknown target_status {target_status!r}")
    if EVIDENCE_STATUS_RANK[target_status] <= EVIDENCE_STATUS_RANK.get(current, -1):
        raise EvidencePromotionError("target_status must be a strict evidence upgrade")
    if target_status in HUMAN_GATED_STATUSES and proposal_origin == "automated_search":
        # Automated systems may draft the proposal, but their output remains a
        # candidate until human review below. This field is retained for audit.
        evidence = dict(evidence)
        evidence["automated_candidate_only"] = True
    proposal = {
        "schema_version": PROMOTION_SCHEMA,
        "proposal_id": "ep-" + _canonical_hash(
            [atlas.get("version"), tissue, cell_path, gene, polarity, target_status, evidence]
        )[:16],
        "created_at": _utc_now(),
        "requested_by": str(requested_by),
        "proposal_origin": str(proposal_origin),
        "atlas_version": str(atlas.get("version", "")),
        "edge": {
            "tissue": tissue,
            "cell_path": cell_path,
            "gene": gene,
            "polarity": polarity,
        },
        "current_status": current,
        "target_status": target_status,
        "current_record_sha256": _canonical_hash(record),
        "evidence": deepcopy(evidence),
        "reviews": [],
        "status": "pending_human_review",
        "claim_boundary": (
            "A proposal is not atlas evidence. Only an applied, versioned proposal with "
            "the required independent human approvals changes runtime evidence status."
        ),
    }
    validate_promotion_proposal(proposal)
    return proposal


def validate_promotion_proposal(proposal: dict) -> None:
    if proposal.get("schema_version") != PROMOTION_SCHEMA:
        raise EvidencePromotionError(f"schema_version must be {PROMOTION_SCHEMA}")
    required = {
        "proposal_id",
        "requested_by",
        "atlas_version",
        "edge",
        "current_status",
        "target_status",
        "current_record_sha256",
        "evidence",
        "reviews",
    }
    missing = required - set(proposal)
    if missing:
        raise EvidencePromotionError(f"Promotion proposal missing fields: {sorted(missing)}")
    edge_missing = {"tissue", "cell_path", "gene", "polarity"} - set(proposal["edge"])
    if edge_missing:
        raise EvidencePromotionError(f"Promotion edge missing fields: {sorted(edge_missing)}")
    target = proposal["target_status"]
    evidence = proposal.get("evidence") or {}
    if target == "database_record_verified":
        for field in ("source_record_id", "source_record_url", "curator_notes"):
            if not str(evidence.get(field, "")).strip():
                raise EvidencePromotionError(f"database promotion requires evidence.{field}")
    if target == "primary_source_verified":
        for field in ("evidence_locator", "curator_notes"):
            if not str(evidence.get(field, "")).strip():
                raise EvidencePromotionError(f"primary promotion requires evidence.{field}")
        sources = evidence.get("sources") or []
        if not any(source.get("source_type") == "primary" for source in sources):
            raise EvidencePromotionError(
                "primary promotion requires a source with source_type='primary'"
            )


def add_promotion_review(
    proposal: dict,
    *,
    reviewer: str,
    decision: str,
    notes: str,
) -> dict:
    """Return a proposal with one append-only reviewer decision."""
    validate_promotion_proposal(proposal)
    reviewer = str(reviewer).strip()
    if not reviewer:
        raise EvidencePromotionError("reviewer is required")
    if reviewer == str(proposal.get("requested_by", "")):
        raise EvidencePromotionError("requester cannot approve their own promotion")
    if decision not in REVIEW_DECISIONS:
        raise EvidencePromotionError(f"decision must be one of {sorted(REVIEW_DECISIONS)}")
    if any(str(item.get("reviewer")) == reviewer for item in proposal.get("reviews", [])):
        raise EvidencePromotionError(f"reviewer {reviewer!r} already reviewed this proposal")
    updated = deepcopy(proposal)
    updated["reviews"].append(
        {
            "reviewer": reviewer,
            "decision": decision,
            "notes": str(notes),
            "reviewed_at": _utc_now(),
        }
    )
    approvals = {item["reviewer"] for item in updated["reviews"] if item["decision"] == "approve"}
    if any(item["decision"] == "reject" for item in updated["reviews"]):
        updated["status"] = "rejected"
    elif len(approvals) >= 2:
        updated["status"] = "approved_pending_apply"
    else:
        updated["status"] = "pending_human_review"
    return updated


def apply_approved_promotion(atlas: dict, proposal: dict, *, new_version: str) -> dict:
    """Apply a two-reviewer-approved proposal to a new atlas version."""
    validate_promotion_proposal(proposal)
    approvals = {
        item["reviewer"]
        for item in proposal.get("reviews", [])
        if item.get("decision") == "approve"
    }
    if len(approvals) < 2 or proposal.get("status") != "approved_pending_apply":
        raise EvidencePromotionError("Promotion requires two distinct approving reviewers")
    if any(item.get("decision") == "reject" for item in proposal.get("reviews", [])):
        raise EvidencePromotionError("Rejected proposal cannot be applied")
    if not str(new_version).strip() or new_version == atlas.get("version"):
        raise EvidencePromotionError("new_version must differ from the source atlas version")
    if str(atlas.get("version", "")) != str(proposal["atlas_version"]):
        raise EvidencePromotionError("Proposal atlas_version does not match source atlas")

    edge = proposal["edge"]
    source_record = _find_record(atlas, **edge)
    if _canonical_hash(source_record) != proposal["current_record_sha256"]:
        raise EvidencePromotionError("Marker edge changed after proposal creation")

    updated = deepcopy(atlas)
    record = _find_record(updated, **edge)
    target = proposal["target_status"]
    evidence = deepcopy(proposal["evidence"])
    record["verification_status"] = target
    record["curator"] = ";".join(sorted(approvals))
    record["verified_at"] = _utc_now()
    record["promotion_proposal_id"] = proposal["proposal_id"]
    record["curator_notes"] = evidence.get("curator_notes", "")
    if target == "database_record_verified":
        record["source_record_id"] = evidence["source_record_id"]
        record["source_record_url"] = evidence["source_record_url"]
    elif target == "primary_source_verified":
        record["evidence_locator"] = evidence["evidence_locator"]
        record["sources"] = evidence["sources"]

    updated["version"] = new_version
    for _tissue, _path, item in _iter_records(updated):
        item["atlas_version"] = new_version
    issues = validate_atlas_provenance(updated)
    if issues:
        raise EvidencePromotionError(
            "Promoted atlas failed provenance validation: " + "; ".join(issues[:5])
        )
    return updated


def write_promotion_proposal(proposal: dict, path: str | Path) -> Path:
    validate_promotion_proposal(proposal)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(proposal, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return target
