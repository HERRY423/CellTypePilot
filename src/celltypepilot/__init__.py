"""CellTypePilot — local-first single-cell annotation review plugin.

Public functions are imported lazily so lightweight commands such as ``doctor``,
``benchmark``, and ``calibrate`` do not initialize Scanpy or optional backends.
"""

from __future__ import annotations

from importlib import import_module

__version__ = "0.3.0"
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

_LAZY_EXPORTS = {
    "compute_marker_scores": (".marker_scorer", "compute_marker_scores"),
    "generate_annotation_summary": (".marker_scorer", "generate_annotation_summary"),
    "score_by_reference": (".reference_scorer", "score_by_reference"),
    "detect_transitional_states": (".reference_scorer", "detect_transitional_states"),
    "check_reference_backends": (".reference_scorer", "check_reference_backends"),
    "ensemble_scores": (".ensemble_scorer", "ensemble_scores"),
    "generate_ensemble_summary": (".ensemble_scorer", "generate_ensemble_summary"),
    "analyze_disagreements": (".ensemble_scorer", "analyze_disagreements"),
    "run_critic": (".critic", "run_critic"),
    "generate_critic_summary": (".critic", "generate_critic_summary"),
}


def __getattr__(name: str):
    if name not in _LAZY_EXPORTS:
        raise AttributeError(name)
    module_name, attribute = _LAZY_EXPORTS[name]
    value = getattr(import_module(module_name, __name__), attribute)
    globals()[name] = value
    return value
