"""Annotation orchestrator — pipeline business logic separated from the CLI.

The CLI layer (``cli.py``) only parses arguments and renders output; all
file I/O orchestration, scoring coordination, and override application
lives here so it can be reused by the Web Inspector and tested directly.
"""

from __future__ import annotations

import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional

import pandas as pd

from .constants import OUTPUT_ANNOTATED, SPECIES_HUMAN, SPECIES_MOUSE

# progress(step, total, message)
ProgressFn = Optional[Callable[[int, int, str], None]]


class PipelineError(ValueError):
    """Raised when the annotation pipeline cannot proceed."""


def find_cluster_column(obs: pd.DataFrame) -> Optional[str]:
    """Locate the cluster label column in obs (ctp_cl_id first, then heuristics)."""
    if "ctp_cl_id" in obs.columns:
        return "ctp_cl_id"
    for col in obs.columns:
        lower = str(col).lower()
        if "cluster" in lower or "cl_id" in lower:
            return col
    return None


def write_annotations_to_adata(
    adata,
    critic_results: pd.DataFrame,
    cluster_key: str,
    output_dir: Path,
) -> Path:
    """Write annotation results back into adata obs and save."""
    cluster_to_ct = dict(zip(critic_results["cluster"], critic_results["cell_type"]))
    cluster_to_cl = dict(zip(
        critic_results["cluster"],
        critic_results.get("cl_id", [""] * len(critic_results)),
    ))
    cluster_to_conf = dict(zip(
        critic_results["cluster"],
        critic_results.get("critic_confidence", [""] * len(critic_results)),
    ))

    adata.obs["ctp_cell_type"] = adata.obs[cluster_key].map(cluster_to_ct).fillna("Unknown")
    adata.obs["ctp_cl_id"] = adata.obs[cluster_key].map(cluster_to_cl).fillna("")
    adata.obs["ctp_confidence"] = adata.obs[cluster_key].map(cluster_to_conf).fillna("unknown")

    output_path = Path(output_dir) / OUTPUT_ANNOTATED
    adata.write(output_path)
    return output_path


def run_annotation_pipeline(
    input_path: str | Path,
    cluster_key: str,
    output_dir: str | Path,
    species: Optional[str] = None,
    tissue: Optional[str] = None,
    embedding_key: Optional[str] = None,
    layer: Optional[str] = None,
    no_figures: bool = False,
    progress: ProgressFn = None,
) -> dict:
    """Run the full marker-scoring annotation pipeline.

    Returns a result dict with: adata, critic_results, critic_summary,
    manifest, figure_paths, output_path, species, tissue, embedding_key,
    data_hash.

    Raises:
        PipelineError: invalid cluster key or no annotations produced.
        FileNotFoundError: input file missing.
    """
    from .data_adapter import (
        load_h5ad, compute_data_hash, detect_species, detect_tissue,
        find_embedding_keys, load_marker_atlas, get_all_markers_for_tissue,
    )
    from .marker_scorer import compute_marker_scores, generate_annotation_summary
    from .critic import run_critic, generate_critic_summary
    from .visualizer import generate_all_figures
    from .reporter import save_evidence_table, generate_html_report, generate_methodology_text
    from .provenance import create_manifest, update_manifest_outputs, save_manifest

    def _emit(step: int, msg: str):
        if progress is not None:
            progress(step, 6, msg)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Step 1: Load data
    _emit(1, "Loading data...")
    adata = load_h5ad(input_path)
    data_hash = compute_data_hash(input_path)

    # Step 2: Detect/auto-set parameters
    _emit(2, "Detecting parameters...")
    if species is None:
        species = detect_species(adata)
    if species not in (SPECIES_HUMAN, SPECIES_MOUSE):
        # Marker atlas only covers human/mouse; continue with human symbols
        # but surface the situation to the caller via the result dict.
        pass
    if tissue is None:
        tissue = detect_tissue(adata) or "general"
    if embedding_key is None:
        candidates = find_embedding_keys(adata)
        if candidates:
            embedding_key = candidates[0]

    if cluster_key not in adata.obs.columns:
        raise PipelineError(
            f"cluster key '{cluster_key}' not found in obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # Step 3: Marker scoring
    _emit(3, "Computing marker scores...")
    atlas = load_marker_atlas(species)
    markers = get_all_markers_for_tissue(atlas, tissue)

    scores = compute_marker_scores(adata, cluster_key, markers, layer=layer)
    summary = generate_annotation_summary(scores, cluster_key)

    if summary.empty:
        raise PipelineError(
            "No annotations generated. Check marker gene overlap with your data."
        )

    # Step 4: Critic
    _emit(4, "Running Annotation Critic...")
    critic_results = run_critic(adata, cluster_key, summary, atlas, tissue)
    critic_summary = generate_critic_summary(critic_results)

    # Step 5: Figures
    figure_paths = []
    if not no_figures and embedding_key:
        _emit(5, "Generating figures...")
        figure_paths = generate_all_figures(
            adata, cluster_key, embedding_key, critic_results, output_path, tissue
        )

    # Step 6: Save outputs
    _emit(6, "Saving outputs...")

    evidence_path = save_evidence_table(critic_results, output_path)

    manifest = create_manifest(
        input_path=str(input_path),
        data_hash=data_hash,
        cluster_key=cluster_key,
        species=species,
        tissue=tissue,
        parameters={
            "embedding_key": embedding_key,
            "layer": layer,
        },
        output_dir=output_path,
    )

    report_path = generate_html_report(
        critic_results, critic_results, critic_summary, manifest, figure_paths, output_path
    )

    method_text = generate_methodology_text(manifest, critic_summary, critic_results)
    method_path = output_path / "methodology_draft.txt"
    with open(method_path, "w", encoding="utf-8") as f:
        f.write(method_text)

    manifest = update_manifest_outputs(manifest, output_path)
    manifest_path = save_manifest(manifest, output_path)

    annotated_path = write_annotations_to_adata(adata, critic_results, cluster_key, output_path)

    return {
        "adata": adata,
        "critic_results": critic_results,
        "critic_summary": critic_summary,
        "manifest": manifest,
        "figure_paths": figure_paths,
        "output_path": output_path,
        "species": species,
        "tissue": tissue,
        "embedding_key": embedding_key,
        "data_hash": data_hash,
        "paths": {
            "evidence": evidence_path,
            "report": report_path,
            "methodology": method_path,
            "manifest": manifest_path,
            "annotated": annotated_path,
        },
    }


def apply_overrides_to_h5ad(
    h5ad_path: str | Path,
    overrides: dict,
    backup_dir: Optional[str | Path] = None,
) -> dict:
    """Apply annotation overrides to an annotated .h5ad file.

    Creates a timestamped backup of the original file before modifying.

    Returns:
        Summary dict: applied, skipped, total, backup, details.
    """
    import scanpy as sc

    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"No annotated data at {h5ad_path}")

    backup_root = Path(backup_dir) if backup_dir else h5ad_path.parent
    backup_name = (
        f"{h5ad_path.stem}.backup_"
        f"{datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')}{h5ad_path.suffix}"
    )
    backup_path = backup_root / backup_name
    shutil.copy2(h5ad_path, backup_path)

    adata = sc.read_h5ad(h5ad_path)
    obs = adata.obs

    # After an h5ad round-trip string columns come back as Categorical;
    # convert so new labels can be assigned without adding categories first.
    for col in ("ctp_cell_type", "ctp_override_reason"):
        if col in obs.columns and isinstance(obs[col].dtype, pd.CategoricalDtype):
            adata.obs[col] = obs[col].astype(object)
    if "ctp_overridden" in obs.columns and isinstance(obs["ctp_overridden"].dtype, pd.CategoricalDtype):
        adata.obs["ctp_overridden"] = obs["ctp_overridden"].astype(bool)

    cluster_col = find_cluster_column(obs)

    applied = 0
    skipped = 0
    details = []

    for cluster_id, override in overrides.items():
        new_type = override.get("new_type", "")
        reason = override.get("reason", "")
        if not new_type:
            skipped += 1
            details.append({"cluster": cluster_id, "status": "skipped", "reason": "Empty new_type"})
            continue

        if cluster_col is None:
            skipped += 1
            details.append({"cluster": cluster_id, "status": "error", "reason": "No cluster column found"})
            continue

        mask = obs[cluster_col].astype(str) == str(cluster_id)
        n_cells = mask.sum()

        if n_cells == 0:
            skipped += 1
            details.append({"cluster": cluster_id, "status": "skipped", "reason": "No cells found"})
            continue

        old_type = obs.loc[mask, "ctp_cell_type"].iloc[0] if "ctp_cell_type" in obs.columns else "Unknown"

        if "ctp_cell_type" in obs.columns:
            adata.obs.loc[mask, "ctp_cell_type"] = new_type
        if "ctp_override_reason" not in obs.columns:
            adata.obs["ctp_override_reason"] = ""
        adata.obs.loc[mask, "ctp_override_reason"] = reason
        if "ctp_overridden" not in obs.columns:
            adata.obs["ctp_overridden"] = False
        adata.obs.loc[mask, "ctp_overridden"] = True

        applied += 1
        details.append({
            "cluster": cluster_id,
            "old_type": old_type,
            "new_type": new_type,
            "n_cells": int(n_cells),
            "reason": reason,
            "status": "applied",
        })

    adata.write(h5ad_path)

    return {
        "applied": applied,
        "skipped": skipped,
        "total": len(overrides),
        "backup": str(backup_path),
        "details": details,
    }


def regenerate_figures_after_override(
    output_dir: str | Path,
    adata,
    cluster_col: str,
) -> list:
    """Regenerate figures after overrides were applied. Returns figure paths."""
    from .data_adapter import find_embedding_keys
    from .visualizer import generate_all_figures

    candidates = find_embedding_keys(adata)
    if not candidates:
        return []
    return generate_all_figures(
        adata, cluster_col, candidates[0], None, Path(output_dir), "general"
    )
