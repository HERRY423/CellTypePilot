"""Annotation orchestrator — pipeline business logic separated from the CLI.

The CLI layer (``cli.py``) only parses arguments and renders output; all
file I/O orchestration, scoring coordination, and override application
lives here so it can be reused by the Web Inspector and tested directly.
"""

from __future__ import annotations

import shutil
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .constants import (
    MARKER_FC_THRESHOLD,
    MARKER_FDR_THRESHOLD,
    MARKER_PCT_THRESHOLD,
    OUTPUT_ANNOTATED,
    SPECIES_HUMAN,
    SPECIES_MOUSE,
)

# progress(step, total, message)
ProgressFn = Callable[[int, int, str], None] | None


class PipelineError(ValueError):
    """Raised when the annotation pipeline cannot proceed."""


def find_cluster_column(obs: pd.DataFrame) -> str | None:
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
    cluster_to_ct = dict(zip(critic_results["cluster"], critic_results["cell_type"], strict=True))
    cluster_to_cl = dict(
        zip(
            critic_results["cluster"],
            critic_results.get("cl_id", [""] * len(critic_results)),
            strict=True,
        )
    )
    cluster_to_conf = dict(
        zip(
            critic_results["cluster"],
            critic_results.get("critic_confidence", [""] * len(critic_results)),
            strict=True,
        )
    )
    cluster_to_candidate = dict(
        zip(
            critic_results["cluster"],
            critic_results.get("candidate_cell_type", critic_results["cell_type"]),
            strict=True,
        )
    )
    cluster_to_decision = dict(
        zip(
            critic_results["cluster"],
            critic_results.get("decision", ["accepted"] * len(critic_results)),
            strict=True,
        )
    )
    cluster_to_reason = dict(
        zip(
            critic_results["cluster"],
            critic_results.get("abstain_reason", [""] * len(critic_results)),
            strict=True,
        )
    )
    cluster_to_state_candidate = dict(
        zip(
            critic_results["cluster"],
            critic_results.get("cell_state_candidate", ["Unknown"] * len(critic_results)),
            strict=True,
        )
    )
    cluster_to_state_decision = dict(
        zip(
            critic_results["cluster"],
            critic_results.get("state_decision", ["abstain"] * len(critic_results)),
            strict=True,
        )
    )
    cluster_to_state_score = dict(
        zip(
            critic_results["cluster"],
            critic_results.get("state_score", [0.0] * len(critic_results)),
            strict=True,
        )
    )
    cluster_to_state_confidence = dict(
        zip(
            critic_results["cluster"],
            critic_results.get("state_confidence", ["needs_review"] * len(critic_results)),
            strict=True,
        )
    )
    cluster_to_state_evidence = dict(
        zip(
            critic_results["cluster"],
            critic_results.get("state_evidence", ["state_not_scored"] * len(critic_results)),
            strict=True,
        )
    )
    cluster_to_display = dict(
        zip(
            critic_results["cluster"],
            critic_results.get("display_label", critic_results["cell_type"]),
            strict=True,
        )
    )

    cluster_series = adata.obs[cluster_key].astype(str)
    adata.obs["ctp_cell_type"] = cluster_series.map(cluster_to_ct).fillna("Unknown")
    adata.obs["ctp_cl_id"] = cluster_series.map(cluster_to_cl).fillna("")
    adata.obs["ctp_confidence"] = cluster_series.map(cluster_to_conf).fillna("unknown")
    adata.obs["ctp_candidate_cell_type"] = cluster_series.map(cluster_to_candidate).fillna(
        "Unknown"
    )
    adata.obs["ctp_decision"] = cluster_series.map(cluster_to_decision).fillna("abstain")
    adata.obs["ctp_abstain_reason"] = cluster_series.map(cluster_to_reason).fillna(
        "cluster_not_scored"
    )
    adata.obs["ctp_cell_state_candidate"] = cluster_series.map(cluster_to_state_candidate).fillna(
        "Unknown"
    )
    adata.obs["ctp_state_decision"] = cluster_series.map(cluster_to_state_decision).fillna(
        "abstain"
    )
    adata.obs["ctp_cell_state"] = adata.obs["ctp_cell_state_candidate"].where(
        adata.obs["ctp_state_decision"] == "supported", "Unknown"
    )
    adata.obs["ctp_state_score"] = cluster_series.map(cluster_to_state_score).fillna(0.0)
    adata.obs["ctp_state_confidence"] = cluster_series.map(cluster_to_state_confidence).fillna(
        "needs_review"
    )
    adata.obs["ctp_state_evidence"] = cluster_series.map(cluster_to_state_evidence).fillna(
        "state_not_scored"
    )
    adata.obs["ctp_display_label"] = cluster_series.map(cluster_to_display).fillna("Unknown")

    output_path = Path(output_dir) / OUTPUT_ANNOTATED
    adata.write(output_path)
    return output_path


def run_annotation_pipeline(
    input_path: str | Path,
    cluster_key: str,
    output_dir: str | Path,
    species: str | None = None,
    tissue: str | None = None,
    embedding_key: str | None = None,
    layer: str | None = None,
    no_figures: bool = False,
    reference_path: str | Path | None = None,
    ref_label_key: str = "cell_type",
    model_path: str | None = None,
    reference_backend: str = "auto",
    marker_weight: float = 0.5,
    no_ensemble: bool = False,
    allow_unverified_reference: bool = False,
    marker_evidence_policy: str = "database",
    calibration_policy_path: str | Path | None = None,
    context_text: str | None = None,
    context_file_path: str | Path | None = None,
    custom_markers_path: str | Path | None = None,
    enable_states: bool = True,
    progress: ProgressFn = None,
) -> dict:
    """Run marker, optional reference, ensemble, critic, and artifact generation.

    Returns a result dict with: adata, critic_results, critic_summary,
    manifest, figure_paths, output_path, species, tissue, embedding_key,
    data_hash.

    Raises:
        PipelineError: invalid cluster key or no annotations produced.
        FileNotFoundError: input file missing.
    """
    from .critic import generate_critic_summary, run_critic
    from .data_adapter import (
        compute_data_hash,
        detect_species,
        detect_tissue,
        find_embedding_keys,
        get_all_markers_for_tissue,
        load_h5ad,
        load_marker_atlas,
        summarize_atlas_evidence,
        validate_atlas_provenance,
    )
    from .ensemble_scorer import (
        analyze_disagreements,
        ensemble_scores,
        generate_ensemble_summary,
    )
    from .marker_scorer import compute_marker_scores, generate_annotation_summary
    from .provenance import create_manifest, save_manifest, update_manifest_outputs
    from .reporter import generate_html_report, generate_methodology_text, save_evidence_table
    from .visualizer import generate_all_figures

    def _emit(step: int, msg: str):
        if progress is not None:
            progress(step, 8, msg)

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
    if layer is not None and layer not in adata.layers:
        raise PipelineError(f"layer '{layer}' not found. Available layers: {list(adata.layers)}")

    # Resolve governed context before selecting candidate marker panels. Free
    # text is retained for provenance only; only structured hypotheses merge.
    from .context_pack import (
        ContextPackError,
        context_manifest_parameters,
        load_context_pack,
        merge_identity_hypotheses,
        resolve_atlas_tissue,
    )

    try:
        context_pack = load_context_pack(
            context_text=context_text,
            context_file=context_file_path,
            custom_markers_file=custom_markers_path,
            species=species,
            tissue=tissue,
        )
    except ContextPackError as exc:
        raise PipelineError(f"Context safety gate failed: {exc}") from exc
    context_enabled = bool(context_text or context_file_path or custom_markers_path)

    # Step 3: Marker scoring
    _emit(3, "Computing marker scores...")
    atlas = load_marker_atlas(species)
    try:
        atlas_tissue = resolve_atlas_tissue(tissue, atlas, context_pack)
    except ContextPackError as exc:
        raise PipelineError(f"Context safety gate failed: {exc}") from exc
    provenance_issues = validate_atlas_provenance(atlas)
    if provenance_issues:
        raise PipelineError(
            "Marker atlas provenance validation failed: " + "; ".join(provenance_issues[:5])
        )
    markers = get_all_markers_for_tissue(
        atlas,
        atlas_tissue,
        evidence_policy=marker_evidence_policy,
    )
    markers = merge_identity_hypotheses(markers, context_pack)
    atlas_evidence_summary = summarize_atlas_evidence(atlas, atlas_tissue)
    n_marker_relationships = sum(
        len(info.get("positive_markers", [])) + len(info.get("negative_markers", []))
        for info in markers.values()
    )

    scores = compute_marker_scores(adata, cluster_key, markers, layer=layer)
    summary = generate_annotation_summary(scores, cluster_key)

    if summary.empty and reference_path is None and model_path is None:
        raise PipelineError("No annotations generated. Check marker gene overlap with your data.")

    # Step 4: Optional reference scoring and ensemble fusion. This is part of
    # the same artifact-producing pipeline, not a side command.
    _emit(4, "Computing optional reference and ensemble scores...")
    ref_scores = pd.DataFrame()
    ensemble_df = pd.DataFrame()
    transitions = pd.DataFrame()
    disagreements = pd.DataFrame()
    uses_reference = reference_path is not None or model_path is not None
    reference_hash = None
    model_hash = None
    reference_contract = None
    if uses_reference:
        from .reference_registry import ReferenceContractError
        from .reference_scorer import detect_transitional_states, score_by_reference

        reference = load_h5ad(reference_path) if reference_path is not None else None
        if reference_path is not None:
            reference_hash = compute_data_hash(reference_path)
        if model_path is not None and Path(model_path).is_file():
            model_hash = compute_data_hash(model_path)
        try:
            ref_scores = score_by_reference(
                adata,
                cluster_key,
                reference=reference,
                ref_label_key=ref_label_key,
                model_path=model_path,
                backend=reference_backend,
                species=species,
                tissue=tissue,
                allow_unverified_reference=allow_unverified_reference,
            )
        except ReferenceContractError as exc:
            raise PipelineError(f"Reference safety gate failed: {exc}") from exc
        reference_contract = dict(ref_scores.attrs.get("reference_contract", {}))
        ensemble_df = ensemble_scores(
            pd.DataFrame() if no_ensemble else scores,
            ref_scores,
            marker_weight=marker_weight,
            adaptive=True,
        )
        ensemble_summary = generate_ensemble_summary(ensemble_df)
        summary = _merge_ensemble_annotation_evidence(ensemble_summary, scores, markers)
        transitions = detect_transitional_states(ref_scores, scores)
        disagreements = analyze_disagreements(ensemble_df)

    # Step 5: Critic
    _emit(5, "Running Annotation Critic...")
    critic_results = run_critic(
        adata,
        cluster_key,
        summary,
        atlas,
        atlas_tissue,
        ensemble_info=ensemble_df if not ensemble_df.empty else None,
        layer=layer,
        evidence_policy=marker_evidence_policy,
        marker_definitions=markers,
    )
    calibration_policy = None
    calibration_policy_hash = None
    if calibration_policy_path is not None:
        import json

        from .calibration import CalibrationError, apply_policy_to_annotations

        policy_path = Path(calibration_policy_path)
        if not policy_path.exists():
            raise PipelineError(f"Calibration policy not found: {policy_path}")
        try:
            calibration_policy = json.loads(policy_path.read_text(encoding="utf-8"))
            critic_results = apply_policy_to_annotations(critic_results, calibration_policy)
        except (json.JSONDecodeError, CalibrationError, KeyError, ValueError) as exc:
            raise PipelineError(f"Invalid calibration policy: {exc}") from exc
        calibration_policy_hash = compute_data_hash(policy_path)
    critic_summary = generate_critic_summary(critic_results)

    # Step 6: State Lens. This is an independent output axis; the attach
    # function asserts that canonical identity columns are unchanged.
    _emit(6, "Scoring independent cell states...")
    state_results = pd.DataFrame()
    if enable_states:
        from .state_scorer import (
            StateScoringError,
            attach_state_results,
            load_state_definitions,
            score_cell_states,
        )

        try:
            state_definitions = load_state_definitions(species, tissue, context_pack)
            state_results = score_cell_states(
                adata,
                cluster_key,
                critic_results,
                state_definitions,
                layer=layer,
            )
            critic_results = attach_state_results(critic_results, state_results)
        except StateScoringError as exc:
            raise PipelineError(f"State safety gate failed: {exc}") from exc

    # Step 7: Figures
    figure_paths = []
    if not no_figures and embedding_key:
        _emit(7, "Generating figures...")
        figure_paths = generate_all_figures(
            adata, cluster_key, embedding_key, critic_results, output_path, tissue
        )

    # Step 8: Save all outputs and hash them into one manifest.
    _emit(8, "Saving outputs...")

    evidence_path = save_evidence_table(critic_results, output_path)

    auxiliary_paths = {}
    for name, frame in (
        ("reference_scores", ref_scores),
        ("ensemble_scores", ensemble_df),
        ("transitional_states", transitions),
        ("disagreements", disagreements),
        ("state_results", state_results),
    ):
        if not frame.empty:
            frame_path = output_path / f"{name}.csv"
            frame.to_csv(frame_path, index=False)
            auxiliary_paths[name] = frame_path

    manifest = create_manifest(
        input_path=str(input_path),
        data_hash=data_hash,
        cluster_key=cluster_key,
        species=species,
        tissue=tissue,
        parameters={
            "embedding_key": embedding_key,
            "layer": layer,
            "reference_path": str(reference_path) if reference_path is not None else None,
            "ref_label_key": ref_label_key if uses_reference else None,
            "model_path": model_path,
            "reference_backend": reference_backend if uses_reference else None,
            "marker_weight": marker_weight if uses_reference and not no_ensemble else None,
            "no_ensemble": no_ensemble if uses_reference else None,
            "reference_sha256": reference_hash,
            "model_sha256": model_hash,
            "reference_contract": reference_contract,
            "allow_unverified_reference": allow_unverified_reference if uses_reference else None,
            "de_method": "wilcoxon",
            "de_direction": "positive",
            "de_log2fc_min": MARKER_FC_THRESHOLD,
            "de_fdr_max": MARKER_FDR_THRESHOLD,
            "marker_expression_fraction_min": MARKER_PCT_THRESHOLD,
            "pipeline_stages": [
                "context",
                "marker",
                "reference" if uses_reference else "reference_skipped",
                "ensemble" if uses_reference and not no_ensemble else "ensemble_skipped",
                "critic",
                "state" if enable_states else "state_skipped",
                "writeback",
                "report",
                "manifest",
            ],
            "marker_atlas_schema": atlas.get("schema_version"),
            "marker_relationships_in_scope": n_marker_relationships,
            "marker_relationships_total_tissue": atlas_evidence_summary["total_relationships"],
            "marker_provenance_validation": "passed",
            "marker_evidence_policy": marker_evidence_policy,
            "marker_evidence_summary": atlas_evidence_summary,
            "atlas_tissue": atlas_tissue,
            **context_manifest_parameters(context_pack, context_enabled),
            "state_lens_enabled": enable_states,
            "state_contract": "identity_invariant_independent_axis_v1",
            "calibration_policy_path": str(calibration_policy_path)
            if calibration_policy_path is not None
            else None,
            "calibration_policy_sha256": calibration_policy_hash,
            "calibration_policy": calibration_policy,
        },
        output_dir=output_path,
    )

    if context_enabled:
        import json

        context_output = output_path / "context_pack.normalized.json"
        context_output.write_text(
            json.dumps(context_pack, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

    report_path = generate_html_report(
        critic_results, critic_results, critic_summary, manifest, figure_paths, output_path
    )

    method_text = generate_methodology_text(manifest, critic_summary, critic_results)
    method_path = output_path / "methodology_draft.txt"
    with open(method_path, "w", encoding="utf-8") as f:
        f.write(method_text)

    annotated_path = write_annotations_to_adata(adata, critic_results, cluster_key, output_path)
    manifest = update_manifest_outputs(manifest, output_path)
    manifest_path = save_manifest(manifest, output_path)

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
            **auxiliary_paths,
        },
    }


def _merge_ensemble_annotation_evidence(
    ensemble_summary: pd.DataFrame,
    marker_scores: pd.DataFrame,
    markers: dict[str, dict],
) -> pd.DataFrame:
    """Attach marker evidence and ontology identifiers to ensemble top calls."""
    if ensemble_summary.empty:
        return ensemble_summary

    result = ensemble_summary.copy()
    result["combined_score"] = result["ensemble_score"]
    evidence_columns = [
        column
        for column in marker_scores.columns
        if column not in {"cl_id", "combined_score", "rank"}
    ]
    if not marker_scores.empty:
        result = result.merge(
            marker_scores[evidence_columns],
            on=["cluster", "cell_type"],
            how="left",
        )
    result["cl_id"] = result["cell_type"].map(
        {cell_type: info.get("cl_id", "") for cell_type, info in markers.items()}
    )
    for column in ("pct_overlap", "mean_log2fc", "specificity", "neg_conflict"):
        if column not in result:
            result[column] = 0.0
        result[column] = result[column].fillna(0.0)
    return result


def apply_overrides_to_h5ad(
    h5ad_path: str | Path,
    overrides: dict,
    backup_dir: str | Path | None = None,
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
    if "ctp_overridden" in obs.columns and isinstance(
        obs["ctp_overridden"].dtype, pd.CategoricalDtype
    ):
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
            details.append(
                {"cluster": cluster_id, "status": "error", "reason": "No cluster column found"}
            )
            continue

        mask = obs[cluster_col].astype(str) == str(cluster_id)
        n_cells = mask.sum()

        if n_cells == 0:
            skipped += 1
            details.append({"cluster": cluster_id, "status": "skipped", "reason": "No cells found"})
            continue

        old_type = (
            obs.loc[mask, "ctp_cell_type"].iloc[0] if "ctp_cell_type" in obs.columns else "Unknown"
        )

        if "ctp_cell_type" in obs.columns:
            adata.obs.loc[mask, "ctp_cell_type"] = new_type
        if "ctp_override_reason" not in obs.columns:
            adata.obs["ctp_override_reason"] = ""
        adata.obs.loc[mask, "ctp_override_reason"] = reason
        if "ctp_overridden" not in obs.columns:
            adata.obs["ctp_overridden"] = False
        adata.obs.loc[mask, "ctp_overridden"] = True

        applied += 1
        details.append(
            {
                "cluster": cluster_id,
                "old_type": old_type,
                "new_type": new_type,
                "n_cells": int(n_cells),
                "reason": reason,
                "status": "applied",
            }
        )

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
