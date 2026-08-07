"""CellTypePilot — Single-cell annotation intelligence layer."""

from .critic import generate_critic_summary, run_critic
from .ensemble_scorer import analyze_disagreements, ensemble_scores, generate_ensemble_summary
from .marker_scorer import compute_marker_scores, generate_annotation_summary
from .reference_scorer import (
    check_reference_backends,
    detect_transitional_states,
    score_by_reference,
)

__version__ = "0.1.0"
MKG_VERSION = "mkg-2026.08"  # Marker Knowledge Graph version

__all__ = [
    "__version__",
    "MKG_VERSION",
    "compute_marker_scores",
    "generate_annotation_summary",
    "score_by_reference",
    "detect_transitional_states",
    "check_reference_backends",
    "ensemble_scores",
    "generate_ensemble_summary",
    "analyze_disagreements",
    "run_critic",
    "generate_critic_summary",
]
