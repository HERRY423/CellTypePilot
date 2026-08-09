"""Freeze broad-class label maps before public benchmark predictions are scored."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import anndata as ad
import pandas as pd

from celltypepilot.data_adapter import get_all_markers_for_tissue, load_marker_atlas


def broad_class(label: str) -> str:
    value = label.strip().lower()
    if value in {"", "unknown", "unassigned", "abstain", "nan", "na"}:
        return "Unknown"
    if "plasma" in value or "plasmablast" in value:
        return "plasma_cell"
    if "natural killer" in value or value.startswith("nk ") or " nk cell" in value:
        return "nk_cell"
    if "mast cell" in value:
        return "mast_cell"
    if (
        re.search(r"\bb cell\b", value)
        or "b lymphocyte" in value
        or value.startswith("b-cell")
    ):
        return "b_cell"
    if (
        re.search(r"\bt cell\b", value)
        or "t-cell" in value
        or "t lymphocyte" in value
        or value.startswith("cd4")
        or value.startswith("cd8")
        or "thymocyte" in value
    ):
        return "t_cell"
    if any(
        token in value
        for token in ("myeloid", "monocyte", "macrophage", "dendritic", "neutrophil")
    ):
        return "myeloid"
    if "endothelial" in value:
        return "endothelial"
    if "pericyte" in value:
        return "pericyte"
    if "smooth muscle" in value:
        return "smooth_muscle"
    if "fibroblast" in value or "stromal" in value:
        return "fibroblast_stromal"
    if any(
        token in value
        for token in (
            "epithelial",
            "enterocyte",
            "goblet",
            "paneth",
            "tuft",
            "stem cell",
            "transit amplifying",
            "enteroendocrine",
            "alveolar type",
            "ciliated",
            "club cell",
            "basal cell",
            "ionocyte",
            "m cell of gut",
            "mucus secreting",
            "neuroendocrine",
            "serous cell",
        )
    ):
        return "epithelial"
    return "other"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    registry_path = Path(args.registry).resolve()
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    atlas = load_marker_atlas("human")
    summaries = []

    for cohort in registry["cohorts"]:
        data_path = (registry_path.parent / cohort["local_path"]).resolve()
        output_path = (registry_path.parent / cohort["label_map_path"]).resolve()
        if output_path.exists() and not args.force:
            summaries.append({"cohort_id": cohort["cohort_id"], "status": "already_present"})
            continue
        dataset = ad.read_h5ad(data_path, backed="r")
        truth_key = cohort["metadata"]["truth_key"]
        truth_labels = sorted(dataset.obs[truth_key].astype(str).unique())
        if dataset.file is not None:
            dataset.file.close()
        atlas_labels = sorted(get_all_markers_for_tissue(atlas, cohort["tissue"]).keys())
        atlas_labels.append("Unknown")
        rows = [
            {
                "method": "__truth__",
                "raw_label": label,
                "canonical_label": broad_class(label),
            }
            for label in truth_labels
        ]
        for method in registry["required_methods"]:
            source_labels = (
                sorted(set(atlas_labels) | set(truth_labels))
                if method == "celltypepilot"
                else sorted(set(truth_labels) | {"Unknown"})
            )
            rows.extend(
                {
                    "method": method,
                    "raw_label": label,
                    "canonical_label": broad_class(label),
                }
                for label in source_labels
            )
        frame = pd.DataFrame(rows).drop_duplicates(["method", "raw_label"])
        output_path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(output_path, index=False, lineterminator="\n")
        audit_path = output_path.with_suffix(".audit.json")
        audit = {
            "schema_version": "celltypepilot.label-map-audit.v1",
            "cohort_id": cohort["cohort_id"],
            "evaluation_level": "predeclared_broad_cell_class",
            "truth_labels": len(truth_labels),
            "celltypepilot_atlas_labels": len(atlas_labels),
            "canonical_classes": sorted(frame["canonical_label"].unique()),
            "other_truth_labels": sorted(
                label for label in truth_labels if broad_class(label) == "other"
            ),
            "mapping_policy": "ordered_token_aware_rules_in_scripts/build_public_label_maps.py",
            "label_map_sha256": sha256(output_path),
            "frozen_before_prediction_scoring": True,
        }
        audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
        summaries.append(
            {
                "cohort_id": cohort["cohort_id"],
                "status": "written",
                "label_map": str(output_path),
                "audit": str(audit_path),
                "other_truth_labels": len(audit["other_truth_labels"]),
            }
        )
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
