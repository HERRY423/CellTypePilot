"""Truth-free lineage-coverage audits for backend-neutral candidate selection."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from .candidate_backends import (
    aggregate_cell_candidates,
    concatenate_candidates,
    normalize_candidate_table,
)
from .data_adapter import load_marker_atlas
from .hierarchical_selector import select_hierarchical_identities
from .identity_contract import (
    build_identity_resolver,
    collect_pack_identity_contract,
    compose_marker_definitions,
)
from .pack_manager import merge_marker_atlas, resolve_extension_packs
from .validation_domains import load_validation_domains

LINEAGE_AUDIT_SCHEMA = "celltypepilot.lineage-coverage-audit.v1"


class LineageCoverageError(ValueError):
    """Raised when a lineage audit cannot be constructed without outcome access."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fine_lineage(label: str) -> str:
    value = " ".join(str(label).casefold().replace("+", " positive ").split())
    if any(
        token in value
        for token in (
            "t cell",
            "b cell",
            "nk cell",
            "lymphocyte",
            "plasma cell",
            "plasmablast",
            "innate lymphoid",
        )
    ):
        return "lymphoid"
    if any(
        token in value
        for token in (
            "myeloid",
            "monocyte",
            "macrophage",
            "dendritic",
            "neutrophil",
            "mast cell",
            "leukocyte",
        )
    ):
        return "myeloid"
    if any(token in value for token in ("endothelial", "pericyte", "vascular")):
        return "vascular"
    if any(
        token in value
        for token in ("fibroblast", "stromal", "smooth muscle", "schwann", "mesenchymal")
    ):
        return "stromal"
    if any(
        token in value
        for token in (
            "epithelial",
            "alveolar",
            "ciliated",
            "club cell",
            "goblet",
            "basal cell",
            "enterocyte",
            "stem cell",
            "paneth",
            "tuft",
            "neuroendocrine",
            "malignant",
        )
    ):
        return "epithelial"
    return "other"


def _domain_lineage(label: str, domain_id: str) -> str:
    fine = _fine_lineage(label)
    if domain_id == "lung" and fine in {"myeloid", "lymphoid"}:
        return "immune"
    if domain_id == "tumor_microenvironment":
        if fine in {"myeloid", "lymphoid"}:
            return "immune"
        if fine == "epithelial":
            return "malignant_or_epithelial"
    return fine


def build_selector_lineage_audit(
    predictions_path: str | Path,
    cluster_map_path: str | Path,
    output_dir: str | Path,
    *,
    domain_id: str,
    methods: tuple[str, ...] = ("celltypist", "popv", "singler"),
) -> dict:
    """Re-select OOF backend predictions and audit lineages without reading truth.

    ``predictions_path`` must contain fold-isolated cell predictions. The audit
    deliberately ignores any truth table or accuracy metric; it establishes only
    that candidate generation no longer collapses to a single lineage.
    """
    domains = load_validation_domains()
    if domain_id not in domains["domains"]:
        raise LineageCoverageError(f"Unknown validation domain: {domain_id}")
    contract = domains["domains"][domain_id]["atlas_contract"]
    predictions_path = Path(predictions_path).resolve()
    cluster_map_path = Path(cluster_map_path).resolve()
    if not predictions_path.is_file() or not cluster_map_path.is_file():
        raise LineageCoverageError("Predictions and cluster-map artifacts are required")
    predictions = pd.read_csv(predictions_path, dtype=str)
    required = {"cell_id", "method", "predicted_label"}
    if not required.issubset(predictions.columns):
        raise LineageCoverageError(
            f"Predictions lack columns: {sorted(required - set(predictions))}"
        )
    selected = predictions[predictions["method"].isin(methods)].copy()
    missing_methods = sorted(set(methods) - set(selected["method"]))
    if missing_methods:
        raise LineageCoverageError(f"Required backend predictions are absent: {missing_methods}")
    if "raw_predicted_label" in selected:
        selected["predicted_label"] = selected["raw_predicted_label"].fillna(
            selected["predicted_label"]
        )
    cluster_map = pd.read_csv(cluster_map_path, dtype=str)
    if not {"cell_id", "cluster"}.issubset(cluster_map.columns):
        raise LineageCoverageError("Cluster map requires cell_id and cluster")
    assignments = cluster_map.set_index("cell_id")["cluster"]
    assignments.index = assignments.index.astype(str)

    records, warnings = resolve_extension_packs(contract["required_packs"], "human")
    atlas = load_marker_atlas("human")
    atlas, merge_warnings = merge_marker_atlas(atlas, records, "human")
    pack_contract = collect_pack_identity_contract(records)
    tissue = (
        "lung"
        if domain_id == "lung"
        else "gut"
        if domain_id == "gut_ibd"
        else "tumor_microenvironment"
    )
    _, scope = compose_marker_definitions(atlas, tissue, pack_contract=pack_contract)
    resolver = build_identity_resolver(atlas, scope["active_tissues"], pack_contract)

    tables = []
    for method in methods:
        frame = selected[selected["method"] == method].copy()
        if frame.duplicated("cell_id").any():
            raise LineageCoverageError(f"Backend {method} contains duplicate cell predictions")
        frame["backend"] = method
        aggregated = aggregate_cell_candidates(frame, assignments)
        tables.append(
            normalize_candidate_table(
                aggregated,
                resolver,
                default_backend=method,
                source_artifact=f"oof_predictions_sha256:{_sha256(predictions_path)}",
            )
        )
    candidates = concatenate_candidates(tables)
    decisions = select_hierarchical_identities(
        candidates, resolver, sorted(assignments.astype(str).unique())
    )
    decisions["fine_lineage"] = decisions["selective_candidate_cell_type"].map(_fine_lineage)
    decisions["domain_lineage"] = decisions["selective_candidate_cell_type"].map(
        lambda label: _domain_lineage(label, domain_id)
    )
    observed = sorted(set(decisions.loc[decisions["domain_lineage"] != "other", "domain_lineage"]))
    required_lineages = list(contract["required_lineages"])
    missing = sorted(set(required_lineages) - set(observed))
    n_accepted = int(decisions["selective_decision"].str.startswith("accepted").sum())
    status = "passed" if not missing and n_accepted > 0 else "failed"

    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    decisions_path = output / "selector_lineage_decisions.csv"
    candidates_path = output / "backend_candidates.normalized.csv"
    cell_predictions_path = output / "selector_cell_predictions.csv"
    decisions.to_csv(decisions_path, index=False, lineterminator="\n")
    candidates.to_csv(candidates_path, index=False, lineterminator="\n")
    decision_by_cluster = decisions.set_index("cluster")
    cell_predictions = pd.DataFrame(
        {
            "cell_id": assignments.index.astype(str),
            "cluster": assignments.astype(str).to_numpy(),
        }
    )
    selected_label = decision_by_cluster["selected_cell_type"].to_dict()
    agreement = decision_by_cluster["backend_agreement_fraction"].to_dict()
    cell_predictions["predicted_label"] = (
        cell_predictions["cluster"].map(selected_label).fillna("Unknown")
    )
    cell_predictions["confidence"] = pd.to_numeric(
        cell_predictions["cluster"].map(agreement), errors="coerce"
    ).fillna(0.0)
    cell_predictions["method"] = "celltypepilot"
    cell_predictions.to_csv(cell_predictions_path, index=False, lineterminator="\n")
    manifest = {
        "schema_version": LINEAGE_AUDIT_SCHEMA,
        "domain_id": domain_id,
        "status": status,
        "source_predictions_sha256": _sha256(predictions_path),
        "source_cluster_map_sha256": _sha256(cluster_map_path),
        "methods": list(methods),
        "required_packs": contract["required_packs"],
        "pack_warnings": [*warnings, *merge_warnings],
        "active_tissues": scope["active_tissues"],
        "required_lineages": required_lineages,
        "observed_lineages": observed,
        "missing_lineages": missing,
        "n_clusters": int(len(decisions)),
        "n_accepted_clusters": n_accepted,
        "decision_counts": decisions["selective_decision"].value_counts().to_dict(),
        "truth_access": "prohibited_not_read",
        "claim_boundary": (
            "Passing proves multi-lineage candidate addressability in these OOF artifacts only; "
            "it is not accuracy, calibration, biological validation, or a reason to relax abstention."
        ),
        "artifacts": {
            "decisions": decisions_path.name,
            "candidates": candidates_path.name,
            "cell_predictions": cell_predictions_path.name,
        },
    }
    manifest_path = output / "lineage_coverage_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return {
        "manifest": manifest,
        "manifest_path": manifest_path,
        "decisions": decisions,
        "candidates": candidates,
        "cell_predictions": cell_predictions,
    }
