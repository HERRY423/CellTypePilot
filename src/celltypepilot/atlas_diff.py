"""Atlas version comparison engine.

Computes structural differences between two marker atlases.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List


@dataclass
class AtlasDiff:
    """Tracks changes between two atlas versions."""
    added_cell_types: List[str] = field(default_factory=list)
    removed_cell_types: List[str] = field(default_factory=list)
    added_markers: List[Dict[str, str]] = field(default_factory=list)
    removed_markers: List[Dict[str, str]] = field(default_factory=list)
    deprecated_cell_types: List[str] = field(default_factory=list)
    deprecated_markers: List[Dict[str, str]] = field(default_factory=list)
    status_changes: List[Dict[str, Any]] = field(default_factory=list)


def _flatten_atlas(atlas: dict[str, Any]) -> dict[str, Any]:
    """Flatten atlas into easily comparable dictionary of node info and edge info."""
    nodes = {}
    edges = {}

    def _walk(cell_types: dict[str, Any], tissue: str, path: str):
        for ct_name, ct_info in cell_types.items():
            current_path = f"{tissue}/{path}/{ct_name}" if path else f"{tissue}/{ct_name}"
            nodes[current_path] = {
                "deprecated": ct_info.get("deprecated", False)
            }
            
            for record in ct_info.get("marker_evidence", []):
                gene = record.get("gene")
                polarity = record.get("polarity")
                if gene and polarity:
                    edge_key = f"{current_path}|{gene}|{polarity}"
                    edges[edge_key] = {
                        "tissue": tissue,
                        "cell_type": current_path,
                        "gene": gene,
                        "polarity": polarity,
                        "verification_status": record.get("verification_status"),
                        "deprecated": record.get("deprecated", False)
                    }
                    
            _walk(ct_info.get("subtypes", {}), tissue, current_path)

    for tissue, tissue_data in atlas.get("tissues", {}).items():
        _walk(tissue_data.get("cell_types", {}), tissue, "")
        
    return {"nodes": nodes, "edges": edges}


def diff_atlases(old_atlas: dict[str, Any], new_atlas: dict[str, Any]) -> AtlasDiff:
    """Compute structural differences between two atlases."""
    old_flat = _flatten_atlas(old_atlas)
    new_flat = _flatten_atlas(new_atlas)
    
    diff = AtlasDiff()
    
    old_nodes = old_flat["nodes"]
    new_nodes = new_flat["nodes"]
    
    diff.added_cell_types = list(set(new_nodes.keys()) - set(old_nodes.keys()))
    diff.removed_cell_types = list(set(old_nodes.keys()) - set(new_nodes.keys()))
    
    for node, info in new_nodes.items():
        if node in old_nodes:
            if info["deprecated"] and not old_nodes[node]["deprecated"]:
                diff.deprecated_cell_types.append(node)
                
    old_edges = old_flat["edges"]
    new_edges = new_flat["edges"]
    
    for edge in set(new_edges.keys()) - set(old_edges.keys()):
        diff.added_markers.append(new_edges[edge])
        
    for edge in set(old_edges.keys()) - set(new_edges.keys()):
        diff.removed_markers.append(old_edges[edge])
        
    for edge, info in new_edges.items():
        if edge in old_edges:
            old_info = old_edges[edge]
            if info["deprecated"] and not old_info["deprecated"]:
                diff.deprecated_markers.append(info)
            if info["verification_status"] != old_info["verification_status"]:
                diff.status_changes.append({
                    "cell_type": info["cell_type"],
                    "gene": info["gene"],
                    "old_status": old_info["verification_status"],
                    "new_status": info["verification_status"]
                })
                
    return diff


def format_diff_report(diff: AtlasDiff) -> str:
    """Format differences into a human-readable report."""
    lines = ["Atlas Diff Report", "=" * 20]
    
    if diff.added_cell_types:
        lines.append(f"\nAdded Cell Types ({len(diff.added_cell_types)}):")
        for ct in diff.added_cell_types:
            lines.append(f"  + {ct}")
            
    if diff.removed_cell_types:
        lines.append(f"\nRemoved Cell Types ({len(diff.removed_cell_types)}):")
        for ct in diff.removed_cell_types:
            lines.append(f"  - {ct}")
            
    if diff.deprecated_cell_types:
        lines.append(f"\nDeprecated Cell Types ({len(diff.deprecated_cell_types)}):")
        for ct in diff.deprecated_cell_types:
            lines.append(f"  ~ {ct}")
            
    if diff.added_markers:
        lines.append(f"\nAdded Markers ({len(diff.added_markers)}):")
        for m in diff.added_markers:
            lines.append(f"  + {m['gene']} ({m['polarity']}) in {m['cell_type']}")
            
    if diff.removed_markers:
        lines.append(f"\nRemoved Markers ({len(diff.removed_markers)}):")
        for m in diff.removed_markers:
            lines.append(f"  - {m['gene']} ({m['polarity']}) in {m['cell_type']}")
            
    if diff.deprecated_markers:
        lines.append(f"\nDeprecated Markers ({len(diff.deprecated_markers)}):")
        for m in diff.deprecated_markers:
            lines.append(f"  ~ {m['gene']} ({m['polarity']}) in {m['cell_type']}")
            
    if diff.status_changes:
        lines.append(f"\nStatus Changes ({len(diff.status_changes)}):")
        for sc in diff.status_changes:
            lines.append(f"  ~ {sc['gene']} in {sc['cell_type']}: {sc['old_status']} -> {sc['new_status']}")
            
    if len(lines) == 2:
        lines.append("\nNo changes detected.")
        
    return "\n".join(lines)


def format_diff_json(diff: AtlasDiff) -> dict[str, Any]:
    """Format differences into a machine-readable JSON dictionary."""
    return {
        "added_cell_types": diff.added_cell_types,
        "removed_cell_types": diff.removed_cell_types,
        "added_markers": diff.added_markers,
        "removed_markers": diff.removed_markers,
        "deprecated_cell_types": diff.deprecated_cell_types,
        "deprecated_markers": diff.deprecated_markers,
        "status_changes": diff.status_changes
    }
