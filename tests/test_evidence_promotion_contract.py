import json
from pathlib import Path

import pytest

from celltypepilot.evidence_promotion import (
    EvidencePromotionError,
    add_promotion_review,
    apply_approved_promotion,
    build_promotion_proposal,
)

ATLAS_PATH = (
    Path(__file__).parents[1]
    / "src"
    / "celltypepilot"
    / "data"
    / "packs"
    / "lung_evidence_v0_1"
    / "marker_atlas.json"
)


def _proposal(atlas):
    return build_promotion_proposal(
        atlas,
        tissue="lung",
        cell_path="Capillary endothelial cell",
        gene="PECAM1",
        polarity="positive",
        target_status="database_record_verified",
        evidence={
            "source_record_id": "curated-record-123",
            "source_record_url": "https://example.org/records/123",
            "curator_notes": "Marker-to-cell edge checked against the named database record.",
        },
        requested_by="requester@example.org",
        proposal_origin="automated_search",
    )


def test_promotion_requires_two_independent_humans_and_versions_output():
    atlas = json.loads(ATLAS_PATH.read_text(encoding="utf-8"))
    proposal = _proposal(atlas)
    assert proposal["evidence"]["automated_candidate_only"] is True

    with pytest.raises(EvidencePromotionError, match="requester cannot approve"):
        add_promotion_review(
            proposal,
            reviewer="requester@example.org",
            decision="approve",
            notes="self approval",
        )

    proposal = add_promotion_review(
        proposal, reviewer="reviewer-a", decision="approve", notes="record checked"
    )
    assert proposal["status"] == "pending_human_review"
    proposal = add_promotion_review(
        proposal, reviewer="reviewer-b", decision="approve", notes="independent check"
    )
    assert proposal["status"] == "approved_pending_apply"

    promoted = apply_approved_promotion(atlas, proposal, new_version="lung-evidence-0.1.1-reviewed")
    assert atlas["version"] == "lung-evidence-0.1.0-beta"
    record = promoted["tissues"]["lung"]["cell_types"]["Capillary endothelial cell"][
        "marker_evidence"
    ][0]
    assert promoted["version"] == "lung-evidence-0.1.1-reviewed"
    assert record["verification_status"] == "database_record_verified"
    assert set(record["curator"].split(";")) == {"reviewer-a", "reviewer-b"}


def test_promotion_fails_if_locked_edge_changed_after_review():
    atlas = json.loads(ATLAS_PATH.read_text(encoding="utf-8"))
    proposal = _proposal(atlas)
    proposal = add_promotion_review(
        proposal, reviewer="reviewer-a", decision="approve", notes="checked"
    )
    proposal = add_promotion_review(
        proposal, reviewer="reviewer-b", decision="approve", notes="checked"
    )
    atlas["tissues"]["lung"]["cell_types"]["Capillary endothelial cell"]["marker_evidence"][0][
        "evidence_scope"
    ] = "tampered"
    with pytest.raises(EvidencePromotionError, match="changed after proposal"):
        apply_approved_promotion(atlas, proposal, new_version="tamper-must-fail")
