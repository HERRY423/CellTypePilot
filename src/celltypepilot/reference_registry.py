"""Fail-closed reference contracts for plugin-facing annotation workflows."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import anndata as ad


class ReferenceContractError(ValueError):
    """Raised when a reference is incompatible or lacks required provenance."""


# Small by design: only models whose intended scope is explicit may be selected
# automatically. Installing CellTypist must never make an immune model the silent
# default for an unrelated tissue.
CELLTYPIST_MODEL_REGISTRY: dict[str, dict[str, Any]] = {
    "celltypist/immune_all_low": {
        "model_name": "Immune_All_Low.pkl",
        "species": ["human"],
        "tissues": ["blood"],
        "scope": "immune_cell_types_and_subtypes",
        "selection": "registry_approved",
    }
}


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def select_registered_celltypist_model(species: str, tissue: str) -> dict[str, Any]:
    """Return the one compatible default model, or fail instead of guessing."""
    matches = [
        {"registry_id": key, **value}
        for key, value in CELLTYPIST_MODEL_REGISTRY.items()
        if species in value["species"] and tissue in value["tissues"]
    ]
    if len(matches) != 1:
        raise ReferenceContractError(
            "No unique registry-approved CellTypist model matches "
            f"species={species!r}, tissue={tissue!r}. Supply an explicit model plus "
            "a verified .pkl.json sidecar, or provide a contracted reference .h5ad."
        )
    return matches[0]


def _model_sidecar_path(model_path: str | Path) -> Path:
    path = Path(model_path)
    return path.with_suffix(path.suffix + ".json")


def validate_model_sidecar(
    model_path: str | Path,
    species: str,
    tissue: str,
    allow_unverified: bool = False,
) -> dict[str, Any]:
    """Validate an explicit model against its adjacent provenance sidecar."""
    path = Path(model_path)
    if not path.is_file():
        raise ReferenceContractError(f"Reference model does not exist: {path}")
    sidecar = _model_sidecar_path(path)
    if not sidecar.exists():
        if allow_unverified:
            return {
                "status": "unverified_override",
                "model_path": str(path),
                "reason": f"missing sidecar {sidecar.name}",
            }
        raise ReferenceContractError(
            f"Reference model sidecar is required: {sidecar}. "
            "It must declare species, tissues, source, version, labels, and sha256."
        )
    metadata = json.loads(sidecar.read_text(encoding="utf-8"))
    required = {"species", "tissues", "source", "version", "labels", "sha256"}
    missing = required - set(metadata)
    if missing:
        raise ReferenceContractError(f"Model sidecar missing fields: {sorted(missing)}")
    if species not in metadata["species"] or tissue not in metadata["tissues"]:
        raise ReferenceContractError(
            f"Model scope does not include species={species!r}, tissue={tissue!r}"
        )
    observed = file_sha256(path)
    if observed.lower() != str(metadata["sha256"]).lower():
        raise ReferenceContractError("Model sha256 does not match its provenance sidecar")
    return {"status": "verified", "model_path": str(path), **metadata}


def validate_reference_adata(
    reference: ad.AnnData,
    species: str,
    tissue: str,
    label_key: str,
    allow_unverified: bool = False,
) -> dict[str, Any]:
    """Validate a custom AnnData reference contract and query compatibility."""
    if label_key not in reference.obs.columns:
        raise ReferenceContractError(f"Reference label column {label_key!r} is missing")
    if reference.obs[label_key].isna().any():
        raise ReferenceContractError("Reference labels must not contain missing values")

    metadata = reference.uns.get("celltypepilot_reference")
    if not isinstance(metadata, dict):
        if allow_unverified:
            return {
                "status": "unverified_override",
                "reason": "missing reference.uns['celltypepilot_reference'] contract",
                "n_reference_cells": int(reference.n_obs),
            }
        raise ReferenceContractError(
            "Custom reference is missing reference.uns['celltypepilot_reference']. "
            "Declare species, tissues, source, version, label_ontology, and training_studies."
        )

    required = {"species", "tissues", "source", "version", "label_ontology", "training_studies"}
    missing = required - set(metadata)
    if missing:
        raise ReferenceContractError(f"Reference contract missing fields: {sorted(missing)}")
    if species != metadata["species"]:
        raise ReferenceContractError(
            f"Reference species {metadata['species']!r} does not match query {species!r}"
        )
    if tissue not in metadata["tissues"]:
        raise ReferenceContractError(
            f"Reference tissues {metadata['tissues']!r} do not include query tissue {tissue!r}"
        )
    return {
        "status": "verified",
        "n_reference_cells": int(reference.n_obs),
        **metadata,
    }
