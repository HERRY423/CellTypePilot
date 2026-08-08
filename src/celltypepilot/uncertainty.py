"""Product language for scores, uncertainty, and statistical claim boundaries.

CellTypePilot intentionally separates three ideas that are easy to conflate:

1. evidence scores produced during one annotation run,
2. rule-based critic confidence categories, and
3. calibrated statistical risk/probability claims that require separate data.

This module makes those semantics machine-readable so the CLI, reports, MCP
facade, and downstream review tools can use the same language.
"""

from __future__ import annotations

from typing import Any

import pandas as pd

UNCERTAINTY_LANGUAGE_SCHEMA = "celltypepilot.uncertainty-language.v1"

EVIDENCE_SCORE_SEMANTICS = "heuristic_evidence_score_not_probability"
CRITIC_CONFIDENCE_SEMANTICS = "rule_based_review_category_not_probability"
CALIBRATED_PROBABILITY_SEMANTICS = "not_available_from_annotation_run"
UNKNOWN_LABEL_SEMANTICS = "safety_abstention_not_biological_class"
OOD_NOVELTY_SEMANTICS = "not_assessed_in_annotation_run"
RISK_POLICY_NOT_ASSESSED = "not_assessed_in_annotation_run"
RISK_POLICY_APPLIED = "downgrade_only_threshold_policy_applied"

PRODUCT_CLAIM_BOUNDARY = (
    "Annotation-run scores are evidence-ranking signals, not calibrated probabilities. "
    "Unknown is a fail-closed abstention decision, not a biological cell class. "
    "Per-dataset batch robustness, complex-sample robustness, OOD/novelty detection, and "
    "selective-risk guarantees require separate benchmark/calibration artifacts."
)


def _policy_applied(calibration_policy: dict[str, Any] | None) -> bool:
    return bool(
        calibration_policy
        and calibration_policy.get("schema_version") == "celltypepilot.abstention-policy.v1"
    )


def attach_uncertainty_language(
    annotations: pd.DataFrame,
    calibration_policy: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Attach stable uncertainty-language columns to annotation evidence rows.

    The function is deliberately additive and backward-compatible. It keeps
    ``combined_score`` and ``critic_confidence`` unchanged, then adds explicit
    semantics so external agents and web reviewers do not treat those fields as
    probability estimates.
    """
    output = annotations.copy()
    policy_applied = _policy_applied(calibration_policy)
    threshold = (
        float(calibration_policy["threshold"])
        if policy_applied and calibration_policy is not None and "threshold" in calibration_policy
        else pd.NA
    )

    if "combined_score" in output:
        output["evidence_score"] = pd.to_numeric(output["combined_score"], errors="coerce")
    elif "evidence_score" not in output:
        output["evidence_score"] = pd.NA

    output["evidence_score_source"] = "combined_score"
    output["evidence_score_semantics"] = EVIDENCE_SCORE_SEMANTICS
    output["critic_confidence_semantics"] = CRITIC_CONFIDENCE_SEMANTICS
    output["calibrated_probability"] = pd.NA
    output["calibrated_probability_semantics"] = CALIBRATED_PROBABILITY_SEMANTICS
    output["selective_risk_policy_applied"] = policy_applied
    output["selective_risk_policy_threshold"] = threshold
    output["selective_risk_policy_semantics"] = (
        RISK_POLICY_APPLIED if policy_applied else RISK_POLICY_NOT_ASSESSED
    )
    output["ood_novelty_signal"] = OOD_NOVELTY_SEMANTICS
    output["unknown_label_semantics"] = UNKNOWN_LABEL_SEMANTICS
    output["uncertainty_language_schema"] = UNCERTAINTY_LANGUAGE_SCHEMA
    if "decision" in output:
        output["abstention_decision"] = output["decision"].fillna("accepted").astype(str)
    elif "abstention_decision" not in output:
        output["abstention_decision"] = "accepted"
    return output


def build_uncertainty_language_manifest(
    calibration_policy: dict[str, Any] | None = None,
    uses_reference: bool = False,
) -> dict[str, Any]:
    """Return the canonical uncertainty-language block for manifest.json/MCP."""
    policy_applied = _policy_applied(calibration_policy)
    risk_policy: dict[str, Any] = {
        "applied": policy_applied,
        "semantics": RISK_POLICY_APPLIED if policy_applied else RISK_POLICY_NOT_ASSESSED,
    }
    if policy_applied and calibration_policy is not None:
        risk_policy.update(
            {
                "threshold": calibration_policy.get("threshold"),
                "dataset_role": calibration_policy.get("dataset_role"),
                "empirical_selective_error": calibration_policy.get("empirical_selective_error"),
                "empirical_coverage": calibration_policy.get("empirical_coverage"),
                "claim_boundary": (
                    "This policy can only downgrade annotation calls to Unknown. Its empirical "
                    "risk/coverage values describe the separate calibration artifact and are not "
                    "a proof of robustness on the current dataset."
                ),
            }
        )

    return {
        "schema_version": UNCERTAINTY_LANGUAGE_SCHEMA,
        "score_columns": {
            "combined_score": EVIDENCE_SCORE_SEMANTICS,
            "evidence_score": EVIDENCE_SCORE_SEMANTICS,
        },
        "confidence_columns": {
            "critic_confidence": CRITIC_CONFIDENCE_SEMANTICS,
        },
        "probability_columns": {
            "calibrated_probability": CALIBRATED_PROBABILITY_SEMANTICS,
        },
        "abstention": {
            "unknown_label": UNKNOWN_LABEL_SEMANTICS,
            "decision_column": "abstention_decision",
            "fail_closed": True,
        },
        "selective_risk_policy": risk_policy,
        "ood_and_novelty": {
            "status": OOD_NOVELTY_SEMANTICS,
            "claim_boundary": (
                "Novel cell types, species outside the supported atlas scope, and out-of-domain "
                "samples must be handled as review/abstention cases unless a separate OOD "
                "detector is validated."
            ),
        },
        "reference_score_note": (
            "Reference backends may expose model-native confidence or probability-like outputs, "
            "but the exported CellTypePilot ensemble score remains an evidence-ranking score "
            "unless independently calibrated."
            if uses_reference
            else "No reference backend was used in this run."
        ),
        "product_claim_boundary": PRODUCT_CLAIM_BOUNDARY,
    }
