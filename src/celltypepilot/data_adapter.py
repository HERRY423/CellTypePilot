"""Data Adapter — load, inspect, and validate .h5ad input."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from rich.console import Console

from .constants import (
    SPECIES_HUMAN, SPECIES_MOUSE, MIN_CLUSTER_SIZE, ATLAS_PATH,
)

console = Console()


def load_h5ad(path: str | Path) -> ad.AnnData:
    """Load an AnnData object from .h5ad."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    if path.suffix not in (".h5ad", ".h5"):
        raise ValueError(f"Expected .h5ad file, got: {path.suffix}")
    console.print(f"[blue]Loading[/blue] {path} ...")
    adata = ad.read_h5ad(path)
    console.print(f"  → {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def compute_data_hash(path: str | Path) -> str:
    """Compute SHA-256 hash of input file for provenance."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def detect_species(adata: ad.AnnData) -> str:
    """Auto-detect species from gene naming conventions."""
    var_names = list(adata.var_names[:500])
    # Mouse genes: first letter uppercase, rest lowercase (e.g., Cd3d)
    mouse_pattern = sum(1 for g in var_names if g and g[0].isupper() and (len(g) < 2 or g[1].islower()))
    # Human genes: all uppercase (e.g., CD3D)
    human_pattern = sum(1 for g in var_names if g and g == g.upper() and g.isalpha())

    if human_pattern > mouse_pattern * 2:
        return SPECIES_HUMAN
    elif mouse_pattern > human_pattern * 2:
        return SPECIES_MOUSE
    else:
        # Default to human; let user confirm
        return SPECIES_HUMAN


def detect_tissue(adata: ad.AnnData) -> Optional[str]:
    """Try to detect tissue from obs metadata."""
    tissue_keys = ["tissue", "tissue_type", "organ", "sample_tissue"]
    for key in tissue_keys:
        if key in adata.obs.columns:
            vals = adata.obs[key].dropna().unique()
            if len(vals) > 0:
                return str(vals[0])
    return None


def find_cluster_keys(adata: ad.AnnData) -> list[str]:
    """Find candidate cluster keys in obs."""
    cluster_keywords = ["leiden", "louvain", "cluster", "group", "label", "celltype", "cell_type"]
    found = []
    for key in adata.obs.columns:
        key_lower = key.lower()
        for kw in cluster_keywords:
            if kw in key_lower:
                found.append(key)
                break
    return found


def find_embedding_keys(adata: ad.AnnData) -> list[str]:
    """Find candidate embedding keys in obsm."""
    embedding_keywords = ["umap", "tsne", "pca", "embedding", "x_", "spatial"]
    found = []
    for key in adata.obsm.keys():
        key_lower = key.lower()
        for kw in embedding_keywords:
            if kw in key_lower:
                found.append(key)
                break
    return found


def find_layer_keys(adata: ad.AnnData) -> dict:
    """Identify count and normalized layers."""
    result = {"counts": None, "lognorm": None}
    if adata.X is not None:
        # Check if X looks like raw counts (non-negative integers or sparse)
        x_sample = adata.X[:min(100, adata.n_obs), :min(100, adata.n_vars)]
        if hasattr(x_sample, "toarray"):
            x_sample = x_sample.toarray()
        x_arr = np.asarray(x_sample)
        if np.all(x_arr >= 0) and np.allclose(x_arr, x_arr.astype(int)):
            result["counts"] = "X"
        else:
            result["lognorm"] = "X"

    for layer_name in adata.layers:
        if layer_name is None:
            continue
        name_lower = str(layer_name).lower()
        if any(kw in name_lower for kw in ["count", "raw", "unnormalized"]):
            result["counts"] = layer_name
        elif any(kw in name_lower for kw in ["log", "norm", "normalized"]):
            result["lognorm"] = layer_name
        elif any(kw in name_lower for kw in ["scale", "scaled"]):
            pass  # scaled data is not suitable for marker scoring

    return result


def inspect_adata(
    path: str | Path,
    cluster_key: Optional[str] = None,
    embedding_key: Optional[str] = None,
) -> dict:
    """Full inspection report of an h5ad file.

    Returns a structured dict with all relevant metadata.
    """
    adata = load_h5ad(path)
    data_hash = compute_data_hash(path)

    report = {
        "path": str(path),
        "sha256": data_hash,
        "n_obs": adata.n_obs,
        "n_vars": adata.n_vars,
        "species": detect_species(adata),
        "tissue": detect_tissue(adata),
        "obs_columns": list(adata.obs.columns),
        "obsm_keys": list(adata.obsm.keys()),
        "layers": list(adata.layers.keys()),
        "layer_info": find_layer_keys(adata),
        "has_raw": adata.raw is not None,
        "cluster_keys": find_cluster_keys(adata),
        "embedding_keys": find_embedding_keys(adata),
        "cluster_sizes": {},
        "warnings": [],
        "fatal": [],
    }

    # Cluster info
    if cluster_key and cluster_key in adata.obs.columns:
        sizes = adata.obs[cluster_key].value_counts().to_dict()
        report["cluster_sizes"] = {str(k): int(v) for k, v in sizes.items()}
        small_clusters = {k: v for k, v in sizes.items() if v < MIN_CLUSTER_SIZE}
        if small_clusters:
            report["warnings"].append(
                f"{len(small_clusters)} cluster(s) have fewer than {MIN_CLUSTER_SIZE} cells: "
                f"{list(small_clusters.keys())}"
            )
    elif not cluster_key:
        candidates = find_cluster_keys(adata)
        if candidates:
            report["warnings"].append(
                f"No --cluster-key specified. Candidates found: {candidates}"
            )
        else:
            report["fatal"].append("No cluster key found in obs. Cannot proceed without clustering.")

    # Embedding info
    if embedding_key and embedding_key in adata.obsm.keys():
        report["embedding_shape"] = list(adata.obsm[embedding_key].shape)
    elif not embedding_key:
        candidates = find_embedding_keys(adata)
        if candidates:
            report["warnings"].append(
                f"No --embedding-key specified. Candidates found: {candidates}"
            )

    # Gene ID convention
    sample_genes = list(adata.var_names[:100])
    all_upper = all(g == g.upper() for g in sample_genes if g.isalpha())
    first_cap = all(g[0].isupper() and g[1:].islower() for g in sample_genes if g.isalpha() and len(g) > 1)
    if all_upper:
        report["gene_id_convention"] = "human_symbols"
    elif first_cap:
        report["gene_id_convention"] = "mouse_symbols"
    elif all(g.startswith("ENSG") for g in sample_genes):
        report["gene_id_convention"] = "ensembl_human"
    elif all(g.startswith("ENSMUSG") for g in sample_genes):
        report["gene_id_convention"] = "ensembl_mouse"
    else:
        report["gene_id_convention"] = "unknown"

    return report


def format_inspect_report(report: dict) -> str:
    """Format inspection report as human-readable string."""
    lines = [
        "=" * 60,
        "CellTypePilot — Data Inspection Report",
        "=" * 60,
        f"File:        {report['path']}",
        f"SHA-256:     {report['sha256'][:16]}...",
        f"Cells:       {report['n_obs']:,}",
        f"Genes:       {report['n_vars']:,}",
        f"Species:     {report['species']}",
        f"Tissue:      {report['tissue'] or '(not detected)'}",
        f"Gene IDs:    {report['gene_id_convention']}",
        f"Has .raw:    {report['has_raw']}",
        "",
        "Layers:      " + (", ".join(str(l) for l in report["layers"]) if report["layers"] else "(none)"),
        f"  Counts:    {report['layer_info']['counts'] or '(not found)'}",
        f"  LogNorm:   {report['layer_info']['lognorm'] or '(not found)'}",
        "",
        "Obs columns: " + ", ".join(report["obs_columns"][:20]),
        "Obsm keys:   " + ", ".join(str(k) for k in report["obsm_keys"]),
        "",
        ("Cluster keys:    " + ", ".join(report["cluster_keys"]) if report["cluster_keys"] else "Cluster keys:    (none found)"),
        ("Embedding keys:   " + ", ".join(report["embedding_keys"]) if report["embedding_keys"] else "Embedding keys:   (none found)"),
    ]

    if report["cluster_sizes"]:
        lines.append("")
        lines.append("Cluster sizes:")
        for k, v in sorted(report["cluster_sizes"].items(), key=lambda x: -x[1]):
            lines.append(f"  {k}: {v:,} cells")

    if report["warnings"]:
        lines.append("")
        lines.append("  [WARN] Warnings:")
        for w in report["warnings"]:
            lines.append(f"  - {w}")

    if report["fatal"]:
        lines.append("")
        lines.append("  [X] Fatal:")
        for f in report["fatal"]:
            lines.append(f"  - {f}")

    lines.append("=" * 60)
    return "\n".join(lines)


def load_marker_atlas(species: str = "human") -> dict:
    """Load the built-in marker knowledge graph."""
    with open(ATLAS_PATH, "r", encoding="utf-8") as f:
        atlas = json.load(f)

    if species == SPECIES_MOUSE:
        # Convert human gene symbols to mouse conventions
        atlas = _convert_atlas_to_mouse(atlas)

    return atlas


def _convert_atlas_to_mouse(atlas: dict) -> dict:
    """Convert human gene symbols in atlas to mouse conventions."""
    import copy
    mouse_atlas = copy.deepcopy(atlas)
    gene_map = atlas.get("mouse_gene_map", {}).get("exceptions", {})

    def convert_gene(gene: str) -> str:
        if gene in gene_map:
            return gene_map[gene]
        # Default: capitalize first letter, lowercase rest
        if gene and gene.isalpha():
            return gene[0].upper() + gene[1:].lower()
        return gene

    for tissue_key, tissue_data in mouse_atlas.get("tissues", {}).items():
        for ct_key, ct_data in tissue_data.get("cell_types", {}).items():
            if "positive_markers" in ct_data:
                ct_data["positive_markers"] = [convert_gene(g) for g in ct_data["positive_markers"]]
            if "negative_markers" in ct_data:
                ct_data["negative_markers"] = [convert_gene(g) for g in ct_data["negative_markers"]]
            # Recurse into subtypes
            for sub_key, sub_data in ct_data.get("subtypes", {}).items():
                if "positive_markers" in sub_data:
                    sub_data["positive_markers"] = [convert_gene(g) for g in sub_data["positive_markers"]]
                if "negative_markers" in sub_data:
                    sub_data["negative_markers"] = [convert_gene(g) for g in sub_data["negative_markers"]]

    return mouse_atlas


def get_all_markers_for_tissue(atlas: dict, tissue: str) -> dict[str, dict]:
    """Get all cell type markers for a given tissue.

    Returns: {cell_type_name: {positive_markers: [...], negative_markers: [...], cl_id: ...}}
    """
    tissue_data = atlas.get("tissues", {}).get(tissue)
    if not tissue_data:
        # Fall back to general tissue
        tissue_data = atlas.get("tissues", {}).get("general")
    if not tissue_data:
        return {}

    result = {}
    for ct_name, ct_info in tissue_data.get("cell_types", {}).items():
        pos = list(ct_info.get("positive_markers", []))
        neg = list(ct_info.get("negative_markers", []))
        cl_id = ct_info.get("cl_id", "")
        result[ct_name] = {"positive_markers": pos, "negative_markers": neg, "cl_id": cl_id}

        # Also include subtypes
        for sub_name, sub_info in ct_info.get("subtypes", {}).items():
            sub_pos = list(sub_info.get("positive_markers", []))
            sub_neg = list(sub_info.get("negative_markers", []))
            sub_cl = sub_info.get("cl_id", "")
            result[sub_name] = {"positive_markers": sub_pos, "negative_markers": sub_neg, "cl_id": sub_cl}

    return result


def get_all_markers_flat(atlas: dict, tissue: str) -> dict[str, list[str]]:
    """Get a flat mapping of cell_type → positive markers for quick scoring."""
    markers = get_all_markers_for_tissue(atlas, tissue)
    return {ct: info["positive_markers"] for ct, info in markers.items()}
