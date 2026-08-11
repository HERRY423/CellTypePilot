"""Ontology-aware selective identity decisions over backend candidates."""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from typing import Any

import pandas as pd

from .constants import (
    CONFIDENCE_HIGH,
    CONFIDENCE_LOW,
    CONFIDENCE_MEDIUM,
    CONFIDENCE_REVIEW,
)

SELECTOR_SCHEMA = "celltypepilot.hierarchical-selective-decision.v1"

DEFAULT_POLICY: dict[str, Any] = {
    "policy_id": "hierarchical-selective-default-v1",
    "min_independent_backends_leaf": 2,
    "min_independent_backends_ancestor": 2,
    "max_ancestor_hops": 2,
    "allow_llm_decision_vote": False,
    "calibrated": False,
}


class SelectiveDecisionError(ValueError):
    """Raised when a selective-decision policy is unsafe or malformed."""


def validate_policy(policy: dict | None) -> dict:
    if policy is not None and not isinstance(policy, dict):
        raise SelectiveDecisionError("Selective decision policy must be a JSON object")
    merged = {**DEFAULT_POLICY, **(policy or {})}
    if not isinstance(merged.get("policy_id"), str) or not merged["policy_id"].strip():
        raise SelectiveDecisionError("policy_id must be a non-empty string")
    for key in ("min_independent_backends_leaf", "min_independent_backends_ancestor"):
        value = merged.get(key)
        if not isinstance(value, int) or value < 1:
            raise SelectiveDecisionError(f"{key} must be a positive integer")
    max_hops = merged.get("max_ancestor_hops")
    if not isinstance(max_hops, int) or max_hops < 0:
        raise SelectiveDecisionError("max_ancestor_hops must be a non-negative integer")
    if not isinstance(merged.get("allow_llm_decision_vote"), bool):
        raise SelectiveDecisionError("allow_llm_decision_vote must be boolean")
    if merged.get("calibrated") is not False:
        raise SelectiveDecisionError(
            "This selector policy cannot claim calibration; use a separate held-out calibration artifact"
        )
    return merged


def _ancestor_distances(label: str, parent_by_name: dict[str, str]) -> dict[str, int]:
    distances = {label: 0}
    current = label
    while current in parent_by_name:
        parent = parent_by_name[current]
        if parent in distances:
            break
        distances[parent] = distances[current] + 1
        current = parent
    return distances


def _safe_common_ancestor(
    labels: list[str], parent_by_name: dict[str, str], max_hops: int
) -> tuple[str | None, int | None]:
    if not labels:
        return None, None
    paths = [_ancestor_distances(label, parent_by_name) for label in labels]
    common = set(paths[0]).intersection(*(set(path) for path in paths[1:]))
    eligible = [label for label in common if max(path[label] for path in paths) <= max_hops]
    if not eligible:
        return None, None
    chosen = min(
        eligible,
        key=lambda label: (
            max(path[label] for path in paths),
            sum(path[label] for path in paths),
            label,
        ),
    )
    return chosen, max(path[chosen] for path in paths)


def _candidate_set(rows: pd.DataFrame) -> str:
    labels = [str(value) for value in rows["canonical_cell_type"] if str(value)]
    return ";".join(dict.fromkeys(labels))


def select_hierarchical_identities(
    candidates: pd.DataFrame,
    resolver: dict,
    clusters: Iterable[str],
    policy: dict | None = None,
) -> pd.DataFrame:
    """Select leaf, governed ancestor, or abstention without score blending."""
    active_policy = validate_policy(policy)
    parent_by_name = resolver.get("parent_by_name", {})
    cl_by_name = resolver.get("cl_by_name", {})
    records: list[dict] = []

    for cluster in [str(value) for value in clusters]:
        cluster_rows = candidates[candidates["cluster"].astype(str) == cluster].copy()
        top = cluster_rows[
            (cluster_rows["rank"] == 1)
            & cluster_rows["identity_resolved"].astype(bool)
            & (cluster_rows["decision_role"] == "decision_candidate")
        ].copy()
        if active_policy["allow_llm_decision_vote"]:
            top = pd.concat(
                [
                    top,
                    cluster_rows[
                        (cluster_rows["rank"] == 1)
                        & cluster_rows["identity_resolved"].astype(bool)
                        & (cluster_rows["decision_role"] == "hypothesis_only")
                    ],
                ],
                ignore_index=True,
            )

        # Correlated adapters count once.  This is deterministic and avoids
        # manufacturing consensus by running multiple wrappers over one family.
        top = top.sort_values(["independence_group", "backend"])
        independent = top.drop_duplicates("independence_group", keep="first")
        labels = independent["canonical_cell_type"].astype(str).tolist()
        backend_names = independent["backend"].astype(str).tolist()
        groups = independent["independence_group"].astype(str).tolist()
        fallback_rows = cluster_rows[
            (cluster_rows["rank"] == 1) & cluster_rows["identity_resolved"].astype(bool)
        ]
        candidate_set = _candidate_set(
            pd.concat([independent, fallback_rows], ignore_index=True).drop_duplicates(
                "canonical_cell_type"
            )
        )
        counts = Counter(labels)
        n_independent = len(groups)
        max_agreement = max(counts.values(), default=0)
        agreement_fraction = max_agreement / n_independent if n_independent else 0.0

        selected = "Unknown"
        selected_cl = ""
        decision = "abstain"
        reason = "NO_DECISION_ELIGIBLE_BACKEND"
        ancestor_hops: int | None = None

        if n_independent < active_policy["min_independent_backends_leaf"]:
            reason = "INSUFFICIENT_INDEPENDENT_BACKENDS"
        elif len(counts) == 1:
            selected = labels[0]
            selected_cl = cl_by_name.get(selected, "")
            decision = "accepted_leaf"
            reason = ""
            candidate_set = selected
        else:
            ancestor, ancestor_hops = _safe_common_ancestor(
                labels, parent_by_name, active_policy["max_ancestor_hops"]
            )
            if ancestor and n_independent >= active_policy["min_independent_backends_ancestor"]:
                selected = ancestor
                selected_cl = cl_by_name.get(selected, "")
                decision = "accepted_ancestor"
                reason = "BACKEND_DISAGREEMENT_COLLAPSED_TO_GOVERNED_ANCESTOR"
            else:
                reason = "BACKEND_DISAGREEMENT_NO_SAFE_ANCESTOR"

        display_candidate = selected
        if decision == "abstain":
            display_candidate = candidate_set.split(";", 1)[0] if candidate_set else "Unknown"
            selected_cl = cl_by_name.get(display_candidate, "")

        records.append(
            {
                "schema_version": SELECTOR_SCHEMA,
                "cluster": cluster,
                "selective_decision": decision,
                "selected_cell_type": selected,
                "selected_cl_id": cl_by_name.get(selected, "") if selected != "Unknown" else "",
                "selective_candidate_cell_type": display_candidate,
                "selective_candidate_cl_id": selected_cl,
                "candidate_set": candidate_set,
                "n_independent_backends": n_independent,
                "supporting_backends": ";".join(backend_names),
                "independence_groups": ";".join(groups),
                "backend_agreement_fraction": round(agreement_fraction, 4),
                "ancestor_hops": ancestor_hops,
                "selective_abstain_reason": reason if decision == "abstain" else "",
                "selection_note": reason if decision == "accepted_ancestor" else "",
                "decision_policy_id": active_policy["policy_id"],
                "calibrated_probability": None,
                "claim_boundary": (
                    "agreement_is_descriptive_not_calibrated_probability; identity_requires_critic_review"
                ),
            }
        )
    return pd.DataFrame(records)


def attach_marker_evidence(
    decisions: pd.DataFrame, marker_scores: pd.DataFrame, resolver: dict
) -> pd.DataFrame:
    """Build critic input using marker evidence for the selected candidate."""
    rows: list[dict] = []
    cl_by_name = resolver.get("cl_by_name", {})
    for _, decision in decisions.iterrows():
        label = str(decision["selective_candidate_cell_type"])
        match = marker_scores[
            (marker_scores["cluster"].astype(str) == str(decision["cluster"]))
            & (marker_scores["cell_type"].astype(str) == label)
        ]
        evidence = match.iloc[0].to_dict() if not match.empty else {}
        evidence.update(decision.to_dict())
        evidence["cell_type"] = label
        evidence["cl_id"] = cl_by_name.get(label, evidence.get("cl_id", ""))
        evidence.setdefault("combined_score", 0.0)
        evidence.setdefault("pct_overlap", 0.0)
        evidence.setdefault("neg_conflict", 0.0)
        evidence.setdefault("specificity", 0.0)
        score = float(evidence.get("combined_score", 0.0) or 0.0)
        overlap = float(evidence.get("pct_overlap", 0.0) or 0.0)
        conflict = float(evidence.get("neg_conflict", 0.0) or 0.0)
        specificity = float(evidence.get("specificity", 0.0) or 0.0)
        if score >= 0.7 and overlap >= 0.5 and conflict < 0.1 and specificity >= 0.5:
            evidence["confidence"] = CONFIDENCE_HIGH
        elif score >= 0.5 and overlap >= 0.3 and conflict < 0.2:
            evidence["confidence"] = CONFIDENCE_MEDIUM
        elif score >= 0.3:
            evidence["confidence"] = CONFIDENCE_LOW
        else:
            evidence["confidence"] = CONFIDENCE_REVIEW
        rows.append(evidence)
    return pd.DataFrame(rows)


def enforce_selective_decisions(
    critic_results: pd.DataFrame, decisions: pd.DataFrame
) -> pd.DataFrame:
    """Ensure marker/critic logic can downgrade but never rescue the selector."""
    output = critic_results.copy()
    decision_map = decisions.set_index("cluster").to_dict(orient="index")
    for index, row in output.iterrows():
        selective = decision_map[str(row["cluster"])]
        if selective["selective_decision"] != "abstain":
            continue
        reason = selective["selective_abstain_reason"]
        output.at[index, "candidate_cell_type"] = selective["selective_candidate_cell_type"]
        output.at[index, "candidate_cl_id"] = selective["selective_candidate_cl_id"]
        output.at[index, "cell_type"] = "Unknown"
        if "cl_id" in output.columns:
            output.at[index, "cl_id"] = ""
        output.at[index, "decision"] = "abstain"
        prior_reason = str(row.get("abstain_reason", "") or "")
        output.at[index, "abstain_reason"] = "; ".join(
            value for value in (reason, prior_reason) if value
        )
        prior_flags = str(row.get("critic_flags", "") or "")
        if prior_flags == "PASS":
            prior_flags = ""
        flag = f"SELECTIVE_{reason}"
        output.at[index, "critic_flags"] = "; ".join(
            value for value in (prior_flags, flag) if value
        )
    return output


def build_selector_manifest(decisions: pd.DataFrame, policy: dict | None = None) -> dict:
    active_policy = validate_policy(policy)
    counts = decisions["selective_decision"].value_counts().to_dict() if not decisions.empty else {}
    return {
        "schema_version": SELECTOR_SCHEMA,
        "policy": active_policy,
        "decision_counts": counts,
        "marker_role": "evidence_only_not_candidate_vote",
        "llm_role": (
            "decision_candidate" if active_policy["allow_llm_decision_vote"] else "hypothesis_only"
        ),
        "score_fusion": "prohibited_across_backend_semantics",
        "claim_boundary": (
            "agreement is descriptive and not a calibrated probability or selective-risk guarantee"
        ),
    }
