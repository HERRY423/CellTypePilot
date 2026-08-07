"""CellTypePilot — Single-cell annotation intelligence layer."""

__version__ = "0.1.0"
MKG_VERSION = "mkg-2026.08"  # Marker Knowledge Graph version

# Public API
from .marker_scorer import compute_marker_scores, generate_annotation_summary
from .reference_scorer import score_by_reference, detect_transitional_states, check_reference_backends
from .ensemble_scorer import ensemble_scores, generate_ensemble_summary, analyze_disagreements
from .critic import run_critic, generate_critic_summary
