"""Lifecycle management for the marker atlas.

Handles deprecation, sunsetting, and evidence-weighted scoring support.
"""

from __future__ import annotations

import logging
from typing import Any

from .constants import DEPRECATION_SUNSET_DEFAULT_VERSIONS, EVIDENCE_WEIGHT_MAP
from .data_adapter import EVIDENCE_STATUS_RANK

logger = logging.getLogger(__name__)


def deprecate_cell_type(
    atlas: dict[str, Any],
    cell_type: str,
    tissue: str,
    reason: str,
    superseded_by: str | None = None
) -> None:
    """Set deprecated=True on a cell type node."""
    tissue_data = atlas.get("tissues", {}).get(tissue)
    if not tissue_data:
        return

    def _deprecate(cell_types: dict[str, Any]) -> bool:
        if cell_type in cell_types:
            info = cell_types[cell_type]
            info["deprecated"] = True
            info["deprecation_reason"] = reason
            if superseded_by:
                info["superseded_by"] = superseded_by
            return True
        for ct_info in cell_types.values():
            if _deprecate(ct_info.get("subtypes", {})):
                return True
        return False

    if not _deprecate(tissue_data.get("cell_types", {})):
        logger.warning(f"Cell type {cell_type} not found in tissue {tissue} for deprecation.")


def deprecate_marker_edge(
    atlas: dict[str, Any],
    gene: str,
    cell_type: str,
    tissue: str,
    reason: str
) -> None:
    """Set deprecated=True on a marker edge."""
    tissue_data = atlas.get("tissues", {}).get(tissue)
    if not tissue_data:
        return

    def _deprecate_edge(cell_types: dict[str, Any]) -> bool:
        if cell_type in cell_types:
            info = cell_types[cell_type]
            records = info.get("marker_evidence", [])
            found = False
            for record in records:
                if record.get("gene") == gene:
                    record["deprecated"] = True
                    record["deprecation_reason"] = reason
                    found = True
            return found
        for ct_info in cell_types.values():
            if _deprecate_edge(ct_info.get("subtypes", {})):
                return True
        return False

    if not _deprecate_edge(tissue_data.get("cell_types", {})):
        logger.warning(f"Marker {gene} for cell type {cell_type} not found in tissue {tissue} for deprecation.")


def _to_float_version(v: str | float) -> float:
    if isinstance(v, (int, float)):
        return float(v)
    import re
    m = re.search(r"(\d+(?:\.\d+)?)", str(v))
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    return 1.0


def sunset_check(atlas: dict[str, Any], current_version: str | float = 1.0) -> list[str]:
    """Find nodes/edges past their sunset_version."""
    issues = []
    curr_v = _to_float_version(current_version)

    def _check(cell_types: dict[str, Any], path: str):
        for ct_name, ct_info in cell_types.items():
            current_path = f"{path}/{ct_name}"
            if ct_info.get("deprecated"):
                s_ver = _to_float_version(
                    ct_info.get("sunset_version", curr_v + DEPRECATION_SUNSET_DEFAULT_VERSIONS)
                )
                if curr_v >= s_ver:
                    issues.append(f"Cell type {current_path} is past its sunset version {s_ver}.")

            for record in ct_info.get("marker_evidence", []):
                if record.get("deprecated"):
                    s_ver = _to_float_version(
                        record.get("sunset_version", curr_v + DEPRECATION_SUNSET_DEFAULT_VERSIONS)
                    )
                    if curr_v >= s_ver:
                        issues.append(
                            f"Marker edge {record.get('gene')} in {current_path} is past its sunset version {s_ver}."
                        )

            _check(ct_info.get("subtypes", {}), current_path)

    for tissue, tissue_data in atlas.get("tissues", {}).items():
        _check(tissue_data.get("cell_types", {}), tissue)

    return issues


def apply_deprecation_filter(markers_dict: dict[str, dict], include_deprecated: bool = False) -> dict[str, dict]:
    """Filter deprecated entries from markers_dict if not included."""
    if include_deprecated:
        return markers_dict
        
    filtered_dict = {}
    for ct, info in markers_dict.items():
        if info.get("deprecated"):
            continue
            
        filtered_info = dict(info)
        filtered_evidence = []
        deprecated_genes = set()
        
        for record in info.get("marker_evidence", []):
            if record.get("deprecated"):
                deprecated_genes.add(record.get("gene"))
            else:
                filtered_evidence.append(record)
                
        filtered_info["marker_evidence"] = filtered_evidence
        filtered_info["positive_markers"] = [g for g in info.get("positive_markers", []) if g not in deprecated_genes]
        filtered_info["negative_markers"] = [g for g in info.get("negative_markers", []) if g not in deprecated_genes]
        
        filtered_dict[ct] = filtered_info
        
    return filtered_dict


def compute_marker_weights(evidence_records: list[dict[str, Any]]) -> dict[str, float]:
    """Derive weights from evidence rank."""
    weights = {}
    for record in evidence_records:
        gene = record.get("gene")
        status = record.get("verification_status", "aggregate_source_only_not_edge_verified")
        rank = EVIDENCE_STATUS_RANK.get(status, 0)
        weight = EVIDENCE_WEIGHT_MAP.get(rank, 0.5)
        # Handle duplicate genes by taking the max weight
        if gene:
            weights[gene] = max(weights.get(gene, 0.0), weight)
    return weights


def validate_lifecycle_fields(atlas: dict[str, Any]) -> list[str]:
    """Validate v2 fields, return warnings."""
    warnings = []
    
    def _validate(cell_types: dict[str, Any], path: str):
        for ct_name, ct_info in cell_types.items():
            current_path = f"{path}/{ct_name}"
            if ct_info.get("deprecated") and not ct_info.get("deprecation_reason"):
                warnings.append(f"Cell type {current_path} is deprecated but missing 'deprecation_reason'.")
                
            for record in ct_info.get("marker_evidence", []):
                if record.get("deprecated") and not record.get("deprecation_reason"):
                    warnings.append(f"Marker edge {record.get('gene')} in {current_path} is deprecated but missing 'deprecation_reason'.")
                    
            _validate(ct_info.get("subtypes", {}), current_path)

    for tissue, tissue_data in atlas.get("tissues", {}).items():
        _validate(tissue_data.get("cell_types", {}), tissue)
        
    return warnings
