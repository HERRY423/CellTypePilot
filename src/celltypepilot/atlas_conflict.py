"""Intra-atlas consistency checker.

Detects marker conflicts within the knowledge graph.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .constants import CONFLICT_JACCARD_THRESHOLD


@dataclass
class ConflictRecord:
    """Record of a detected conflict in the marker atlas."""

    conflict_type: str
    gene: str
    cell_type_a: str
    cell_type_b: str
    severity: str
    resolution_hint: str


def _jaccard(set1: set, set2: set) -> float:
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union > 0 else 0.0


def detect_marker_conflicts(atlas: dict[str, Any]) -> list[ConflictRecord]:
    """Detect conflicts in the marker atlas."""
    conflicts = []

    # We will compute flat representations of cell types
    flat_nodes = {}

    def _walk(cell_types: dict[str, Any], tissue: str, parent: str | None):
        for ct_name, ct_info in cell_types.items():
            pos = set(ct_info.get("positive_markers", []))
            neg = set(ct_info.get("negative_markers", []))
            flat_nodes[f"{tissue}/{ct_name}"] = {
                "name": ct_name,
                "tissue": tissue,
                "parent": parent,
                "positive": pos,
                "negative": neg,
                "subtypes": list(ct_info.get("subtypes", {}).keys()),
            }
            _walk(ct_info.get("subtypes", {}), tissue, f"{tissue}/{ct_name}")

    for tissue, tissue_data in atlas.get("tissues", {}).items():
        _walk(tissue_data.get("cell_types", {}), tissue, None)

    node_names = list(flat_nodes.keys())

    # 1. Polarity and Uniqueness
    for i in range(len(node_names)):
        for j in range(i + 1, len(node_names)):
            node_a = node_names[i]
            node_b = node_names[j]
            info_a = flat_nodes[node_a]
            info_b = flat_nodes[node_b]

            # Polarity conflicts (same gene +/- in related types - naive check here just looks at any overlap)
            # A more robust check might only check lineage, but here we flag direct contradictions
            pos_neg_overlap_ab = info_a["positive"] & info_b["negative"]
            for gene in pos_neg_overlap_ab:
                conflicts.append(
                    ConflictRecord(
                        conflict_type="polarity",
                        gene=gene,
                        cell_type_a=node_a,
                        cell_type_b=node_b,
                        severity="high"
                        if info_a["parent"] == node_b or info_b["parent"] == node_a
                        else "medium",
                        resolution_hint="Gene cannot be a positive marker for one and negative for another related type.",
                    )
                )
            pos_neg_overlap_ba = info_b["positive"] & info_a["negative"]
            for gene in pos_neg_overlap_ba:
                conflicts.append(
                    ConflictRecord(
                        conflict_type="polarity",
                        gene=gene,
                        cell_type_a=node_b,
                        cell_type_b=node_a,
                        severity="high"
                        if info_a["parent"] == node_b or info_b["parent"] == node_a
                        else "medium",
                        resolution_hint="Gene cannot be a positive marker for one and negative for another related type.",
                    )
                )

            # Uniqueness violations
            if (
                info_a["parent"] != node_b
                and info_b["parent"] != node_a
                and info_a["parent"] != info_b["parent"]
            ):
                sim = _jaccard(info_a["positive"], info_b["positive"])
                if sim > CONFLICT_JACCARD_THRESHOLD:
                    conflicts.append(
                        ConflictRecord(
                            conflict_type="uniqueness",
                            gene="multiple",
                            cell_type_a=node_a,
                            cell_type_b=node_b,
                            severity="high",
                            resolution_hint=f"Jaccard similarity {sim:.2f} exceeds threshold {CONFLICT_JACCARD_THRESHOLD}.",
                        )
                    )

    # 2. Hierarchy contradictions (child missing parent's core markers)
    for node, info in flat_nodes.items():
        if info["parent"]:
            parent_info = flat_nodes[info["parent"]]
            missing_core = parent_info["positive"] - info["positive"]
            for gene in missing_core:
                conflicts.append(
                    ConflictRecord(
                        conflict_type="hierarchy",
                        gene=gene,
                        cell_type_a=info["parent"],
                        cell_type_b=node,
                        severity="medium",
                        resolution_hint="Child should inherit core positive markers from parent.",
                    )
                )

    return conflicts


def validate_no_blocking_conflicts(atlas: dict[str, Any]) -> bool:
    """Return True if there are no high-severity conflicts blocking usage."""
    conflicts = detect_marker_conflicts(atlas)
    return not any(c.severity == "high" for c in conflicts)
