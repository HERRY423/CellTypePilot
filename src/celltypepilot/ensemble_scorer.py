"""Ensemble Scorer — adaptive fusion of marker + reference scores.

Combines the deterministic marker-based scoring with reference embedding
scores to produce more robust annotations, especially for:
- Continuous differentiation trajectories
- Rare transitional states
- Cases where one method fails but the other succeeds

The ensemble uses adaptive weighting: when marker scoring is confident,
it dominates; when markers are ambiguous, reference embedding takes over.
Disagreements between the two methods are flagged as potential novel
or transitional states.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

from .constants import (
    CONFIDENCE_HIGH, CONFIDENCE_MEDIUM, CONFIDENCE_LOW, CONFIDENCE_REVIEW,
)

logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────
# Adaptive weight schedule
# ──────────────────────────────────────────────

# When marker score is above this threshold, trust markers more
MARKER_HIGH_CONFIDENCE = 0.6
# When marker score is below this, trust reference more
MARKER_LOW_CONFIDENCE = 0.3
# Minimum reference score to trust reference override
REF_OVERRIDE_THRESHOLD = 0.5


def ensemble_scores(
    marker_scores: pd.DataFrame,
    ref_scores: pd.DataFrame,
    marker_weight: float = 0.5,
    adaptive: bool = True,
) -> pd.DataFrame:
    """Combine marker and reference scores with adaptive weighting.

    For each cluster, produces a fused score for each cell type that
    appears in either scorer's results.

    Args:
        marker_scores: Output from marker_scorer.compute_marker_scores()
        ref_scores: Output from reference_scorer.score_by_reference()
        marker_weight: Base weight for marker score (0-1).
            Only used when adaptive=False.
        adaptive: If True, dynamically adjust weights based on
            marker confidence. High marker confidence → marker-heavy;
            low marker confidence → reference-heavy.

    Returns:
        DataFrame with columns:
            cluster, cell_type, ensemble_score, marker_score, ref_score,
            marker_weight_used, ref_weight_used, agreement,
            source (marker/ref/both), rank
    """
    if marker_scores.empty and ref_scores.empty:
        return pd.DataFrame()

    if marker_scores.empty:
        return _ref_only_results(ref_scores)

    if ref_scores.empty:
        return _marker_only_results(marker_scores)

    # Build per-cluster ensemble
    clusters = sorted(
        set(marker_scores["cluster"].unique()) |
        set(ref_scores["cluster"].unique())
    )

    all_results = []

    for cl in clusters:
        m_cl = marker_scores[marker_scores["cluster"] == cl].copy()
        r_cl = ref_scores[ref_scores["cluster"] == cl].copy()

        # Determine adaptive weights for this cluster
        m_weight, r_weight = _compute_adaptive_weights(
            m_cl, r_cl, base_marker_weight=marker_weight, adaptive=adaptive
        )

        # Build unified cell type list (union of both scorers)
        m_types = dict(zip(m_cl["cell_type"], m_cl["combined_score"])) if not m_cl.empty else {}
        r_types = dict(zip(r_cl["cell_type"], r_cl["ref_score"])) if not r_cl.empty else {}

        all_types = sorted(set(list(m_types.keys()) + list(r_types.keys())))

        for ct in all_types:
            m_score = m_types.get(ct, 0.0)
            r_score = r_types.get(ct, 0.0)

            # Fused score
            ensemble = m_weight * m_score + r_weight * r_score
            ensemble = max(0.0, min(1.0, ensemble))

            # Source tracking
            has_marker = ct in m_types and m_types[ct] > 0.01
            has_ref = ct in r_types and r_types[ct] > 0.01

            if has_marker and has_ref:
                source = "both"
                agreement = abs(m_score - r_score) < 0.2
            elif has_marker:
                source = "marker"
                agreement = True  # No disagreement possible
            else:
                source = "reference"
                agreement = True

            all_results.append({
                "cluster": str(cl),
                "cell_type": ct,
                "ensemble_score": round(ensemble, 4),
                "marker_score": round(m_score, 4),
                "ref_score": round(r_score, 4),
                "marker_weight_used": round(m_weight, 3),
                "ref_weight_used": round(r_weight, 3),
                "agreement": agreement,
                "source": source,
            })

    df = pd.DataFrame(all_results)
    if df.empty:
        return df

    # Rank within each cluster
    df["rank"] = df.groupby("cluster")["ensemble_score"].rank(
        ascending=False, method="first"
    ).astype(int)
    df = df.sort_values(["cluster", "rank"])

    return df


def _compute_adaptive_weights(
    marker_cl: pd.DataFrame,
    ref_cl: pd.DataFrame,
    base_marker_weight: float,
    adaptive: bool,
) -> tuple[float, float]:
    """Compute adaptive weights for a single cluster.

    Strategy:
    - If marker top-1 has high confidence → weight markers more
    - If marker top-1 has low confidence → weight reference more
    - If reference top-1 is much higher than marker → trust reference
    """
    if not adaptive:
        return base_marker_weight, 1.0 - base_marker_weight

    # Get marker top-1 score
    if not marker_cl.empty:
        m_best = marker_cl.nsmallest(1, "rank")
        m_top_score = float(m_best.iloc[0]["combined_score"])
    else:
        m_top_score = 0.0

    # Get reference top-1 score
    if not ref_cl.empty:
        r_best = ref_cl.nsmallest(1, "ref_rank")
        r_top_score = float(r_best.iloc[0]["ref_score"])
    else:
        r_top_score = 0.0

    # Adaptive schedule
    if m_top_score >= MARKER_HIGH_CONFIDENCE:
        # Marker is confident → trust it
        m_weight = 0.7
    elif m_top_score <= MARKER_LOW_CONFIDENCE:
        # Marker is uncertain → lean on reference
        if r_top_score >= REF_OVERRIDE_THRESHOLD:
            m_weight = 0.2  # Reference override
        else:
            m_weight = 0.3  # Both uncertain → slight marker preference
    else:
        # Medium confidence → balanced
        m_weight = 0.5

    # Special case: strong disagreement with confident reference
    if not ref_cl.empty and not marker_cl.empty:
        m_top_type = marker_cl.nsmallest(1, "rank").iloc[0]["cell_type"]
        r_top_type = ref_cl.nsmallest(1, "ref_rank").iloc[0]["cell_type"]

        if m_top_type != r_top_type and r_top_score > 0.7 and m_top_score < 0.3:
            # Reference strongly disagrees and is confident → trust reference
            m_weight = 0.15

    r_weight = 1.0 - m_weight
    return m_weight, r_weight


def _ref_only_results(ref_scores: pd.DataFrame) -> pd.DataFrame:
    """Wrap reference-only results in ensemble format."""
    df = ref_scores.copy()
    df["ensemble_score"] = df["ref_score"]
    df["marker_score"] = 0.0
    df["marker_weight_used"] = 0.0
    df["ref_weight_used"] = 1.0
    df["agreement"] = True
    df["source"] = "reference"
    df["rank"] = df["ref_rank"]
    return df


def _marker_only_results(marker_scores: pd.DataFrame) -> pd.DataFrame:
    """Wrap marker-only results in ensemble format."""
    df = marker_scores.copy()
    df["ensemble_score"] = df["combined_score"]
    df["ref_score"] = 0.0
    df["marker_weight_used"] = 1.0
    df["ref_weight_used"] = 0.0
    df["agreement"] = True
    df["source"] = "marker"
    df["rank"] = df.get("rank", range(1, len(df) + 1))
    return df


# ──────────────────────────────────────────────
# Ensemble summary
# ──────────────────────────────────────────────

def generate_ensemble_summary(
    ensemble_df: pd.DataFrame,
) -> pd.DataFrame:
    """Generate per-cluster annotation summary from ensemble scores.

    Returns DataFrame with: cluster, cell_type, ensemble_score,
    confidence, source, agreement, marker_score, ref_score
    """
    if ensemble_df.empty:
        return pd.DataFrame()

    # Take top-1 for each cluster
    top1 = ensemble_df[ensemble_df["rank"] == 1].copy()
    top1["confidence"] = top1.apply(_assign_ensemble_confidence, axis=1)

    return top1[[
        "cluster", "cell_type", "ensemble_score", "confidence",
        "source", "agreement", "marker_score", "ref_score",
    ]].reset_index(drop=True)


def _assign_ensemble_confidence(row: pd.Series) -> str:
    """Assign confidence based on ensemble score and agreement."""
    score = row.get("ensemble_score", 0)
    agreement = row.get("agreement", True)
    source = row.get("source", "marker")
    m_score = row.get("marker_score", 0)
    r_score = row.get("ref_score", 0)

    # Both methods agree and score is high → high confidence
    if score >= 0.6 and agreement:
        return CONFIDENCE_HIGH

    # Good score but some disagreement → medium
    if score >= 0.5 and agreement:
        return CONFIDENCE_MEDIUM

    # Disagreement between methods → needs review
    if not agreement:
        if score >= 0.4:
            return CONFIDENCE_LOW
        return CONFIDENCE_REVIEW

    # Single-source results
    if source == "marker" and score >= 0.5:
        return CONFIDENCE_MEDIUM
    if source == "reference" and score >= 0.5:
        return CONFIDENCE_LOW  # Reference-only → slightly lower confidence

    if score >= 0.3:
        return CONFIDENCE_LOW

    return CONFIDENCE_REVIEW


# ──────────────────────────────────────────────
# Disagreement analysis
# ──────────────────────────────────────────────

def analyze_disagreements(
    ensemble_df: pd.DataFrame,
    min_score_gap: float = 0.2,
) -> pd.DataFrame:
    """Find clusters where marker and reference strongly disagree.

    These are candidates for:
    - Transitional / intermediate states
    - Novel cell types not in the reference
    - Marker database gaps
    - Technical artifacts

    Returns:
        DataFrame with disagreement details, sorted by severity.
    """
    if ensemble_df.empty:
        return pd.DataFrame()

    # Find top-1 from each source per cluster
    disagreements = []
    clusters = ensemble_df["cluster"].unique()

    for cl in clusters:
        cl_data = ensemble_df[ensemble_df["cluster"] == cl]

        # Best marker cell type
        m_best = cl_data.nsmallest(1, "rank") if "rank" in cl_data.columns else pd.DataFrame()
        # Best reference cell type
        r_best = cl_data.nsmallest(1, "ref_rank") if "ref_rank" in cl_data.columns else pd.DataFrame()

        if m_best.empty or r_best.empty:
            continue

        m_type = m_best.iloc[0]["cell_type"]
        r_type = r_best.iloc[0]["cell_type"]
        m_score = float(m_best.iloc[0]["marker_score"])
        r_score = float(r_best.iloc[0]["ref_score"])

        if m_type != r_type:
            severity = abs(m_score - r_score)
            if severity >= min_score_gap:
                disagreements.append({
                    "cluster": cl,
                    "marker_type": m_type,
                    "marker_score": round(m_score, 4),
                    "ref_type": r_type,
                    "ref_score": round(r_score, 4),
                    "severity": round(severity, 4),
                    "interpretation": _interpret_disagreement(
                        m_type, m_score, r_type, r_score
                    ),
                })

    df = pd.DataFrame(disagreements)
    if not df.empty:
        df = df.sort_values("severity", ascending=False)

    return df


def _interpret_disagreement(
    m_type: str, m_score: float,
    r_type: str, r_score: float,
) -> str:
    """Interpret the biological meaning of a disagreement."""
    if m_score > 0.5 and r_score < 0.2:
        return (
            f"Marker scoring strongly supports '{m_type}' but reference "
            f"embedding does not recognize '{m_type}'. Possible novel "
            f"subtype or marker database gap."
        )
    elif r_score > 0.5 and m_score < 0.2:
        return (
            f"Reference embedding strongly supports '{r_type}' but marker "
            f"scoring finds weak evidence. Possible transitional state or "
            f"rare cell type with atypical marker expression."
        )
    elif m_score > 0.3 and r_score > 0.3:
        return (
            f"Both methods find evidence but disagree: markers→'{m_type}' "
            f"vs reference→'{r_type}'. Likely a differentiation intermediate "
            f"sharing features of both types. Consider trajectory analysis."
        )
    else:
        return (
            f"Both methods weak: markers→'{m_type}' ({m_score:.2f}) "
            f"vs reference→'{r_type}' ({r_score:.2f}). "
            f"Cluster may be low-quality, doublet-enriched, or novel. "
            f"Manual review recommended."
        )
