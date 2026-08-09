from __future__ import annotations

import json

import pandas as pd

from celltypepilot.agent_evidence import (
    GAP_SCHEMA,
    attach_contrastive_evidence,
    build_actionable_evidence_gaps,
    build_contrastive_evidence,
    write_agent_evidence_artifacts,
)
from celltypepilot.agent_protocol import (
    AGENT_DECISION_SCHEMA,
    agent_decision,
    validate_agent_decision,
)
from celltypepilot.mcp_server import tool_evidence_gap_queue, tool_read_evidence


def _marker_scores() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "cluster": "0",
                "cell_type": "Candidate A",
                "cl_id": "CL:1",
                "combined_score": 0.42,
                "rank": 1,
                "pos_supporting_markers": "G1;G2",
                "pos_missing_markers": "G3",
                "pos_silent_markers": "G4",
                "neg_expressed_markers": "N1",
                "marker_provenance_status": "aggregate_source_only_not_edge_verified",
            },
            {
                "cluster": "0",
                "cell_type": "Candidate B",
                "cl_id": "CL:2",
                "combined_score": 0.39,
                "rank": 2,
                "pos_supporting_markers": "G1;G5",
                "pos_missing_markers": "G6",
                "pos_silent_markers": "",
                "neg_expressed_markers": "",
                "marker_provenance_status": "literature_cooccurrence_supported",
            },
        ]
    )


def test_contrastive_evidence_explains_existing_ranking_without_reranking():
    scores = _marker_scores()
    contrast = build_contrastive_evidence(scores)
    row = contrast.iloc[0]

    assert row["selected_candidate"] == "Candidate A"
    assert row["alternative_candidate"] == "Candidate B"
    assert row["selected_score"] == 0.42
    assert row["alternative_score"] == 0.39
    assert row["score_margin"] == 0.03
    assert row["shared_supporting_markers"] == "G1"
    assert row["selected_only_supporting_markers"] == "G2"
    assert row["alternative_only_supporting_markers"] == "G5"
    assert row["ranking_semantics"] == "relative_evidence_signal_not_probability"


def test_unknown_becomes_bounded_actionable_gap_not_a_replacement_label():
    contrast = build_contrastive_evidence(_marker_scores())
    critic = pd.DataFrame(
        [
            {
                "cluster": "0",
                "cell_type": "Unknown",
                "candidate_cell_type": "Candidate A",
                "candidate_cl_id": "CL:1",
                "decision": "abstain",
                "critic_confidence": "needs_review",
                "critic_flags": (
                    "LOW_DE_SUPPORT; NEG_MARKER_CONFLICT; AGGREGATE_PROVENANCE_ONLY"
                ),
                "abstain_reason": "LOW_DE_SUPPORT; NEG_MARKER_CONFLICT",
                "n_expected_markers": 4,
                "n_missing_markers": 1,
                "n_silent_markers": 1,
                "pct_overlap": 0.25,
                "neg_expressed_markers": "N1",
            }
        ]
    )
    enriched = attach_contrastive_evidence(critic, contrast)
    gaps = build_actionable_evidence_gaps(enriched)
    gap_types = {item["gap_type"] for item in gaps["clusters"][0]["gaps"]}

    assert gaps["schema_version"] == GAP_SCHEMA
    assert gaps["n_unknown_clusters"] == 1
    assert {
        "marker_not_measured",
        "marker_present_but_silent",
        "directional_support_gap",
        "negative_marker_conflict",
        "aggregate_provenance_gap",
    } <= gap_types
    assert "keep_unknown" in gaps["clusters"][0]["allowed_next_actions"]
    assert gaps["clusters"][0]["candidate_cell_type"] == "Candidate A"
    assert "replacement_label" not in gaps["clusters"][0]
    assert gaps["clusters"][0]["human_action_required"] is True


def test_agent_decision_protocol_has_stable_fields_and_protects_reserved_keys():
    result = agent_decision(
        operation="prepare_annotation",
        status="ready",
        decision_scope="plan",
        allowed_next_actions=["annotate_from_plan"],
        artifact_paths={"plan": "plan.json"},
        human_action_required=False,
        claim_boundary="No annotation claim.",
        schema_version="legacy-schema-must-not-overwrite",
    )

    validate_agent_decision(result)
    assert result["schema_version"] == AGENT_DECISION_SCHEMA
    assert result["status"] == "ready"
    assert "do_not_use_cell_state_to_rescue_identity" in result["forbidden_claims"]


def test_contrast_and_gap_artifacts_round_trip(tmp_path):
    contrast = build_contrastive_evidence(_marker_scores())
    critic = pd.DataFrame(
        [
            {
                "cluster": "0",
                "cell_type": "Unknown",
                "candidate_cell_type": "Candidate A",
                "decision": "abstain",
                "critic_flags": "LOW_DE_SUPPORT",
                "n_expected_markers": 4,
                "n_missing_markers": 1,
                "n_silent_markers": 1,
                "pct_overlap": 0.25,
            }
        ]
    )
    enriched = attach_contrastive_evidence(critic, contrast)
    gaps = build_actionable_evidence_gaps(enriched)
    enriched.to_csv(tmp_path / "evidence_table.csv", index=False)
    paths = write_agent_evidence_artifacts(tmp_path, contrast, gaps)

    loaded_contrast = pd.read_csv(paths["contrastive_evidence"])
    loaded_gaps = json.loads(paths["evidence_gaps"].read_text(encoding="utf-8"))
    assert loaded_contrast.iloc[0]["alternative_candidate"] == "Candidate B"
    assert loaded_gaps["schema_version"] == GAP_SCHEMA
    assert loaded_gaps["n_unknown_clusters"] == 1

    evidence_response = tool_read_evidence(str(tmp_path), cluster="0")
    gap_response = tool_evidence_gap_queue(str(tmp_path))
    assert evidence_response["status"] == "available"
    assert evidence_response["contrastive_evidence"][0]["alternative_candidate"] == "Candidate B"
    assert evidence_response["actionable_evidence_gaps"][0]["human_action_required"] is True
    assert gap_response["schema_version"] == GAP_SCHEMA
    assert gap_response["mutation_policy"] == "read_only_bounded_actions_no_label_selection"
