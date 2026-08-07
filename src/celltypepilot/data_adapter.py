"""Data Adapter — load, inspect, and validate .h5ad input."""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from rich.console import Console

from .constants import (
    ATLAS_PATH,
    ENSEMBL_PREFIX_SPECIES,
    MIN_CLUSTER_SIZE,
    SPECIES_DOMINANCE_RATIO,
    SPECIES_HUMAN,
    SPECIES_MOUSE,
    SPECIES_SYMBOL_RATIO,
    TISSUE_COLUMN_KEYWORDS,
    TISSUE_COLUMN_SYNONYMS,
)

logger = logging.getLogger(__name__)
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


def match_ensembl_species(gene: str) -> str | None:
    """Map a single gene ID to a species via Ensembl prefix, or None."""
    for prefix, species in ENSEMBL_PREFIX_SPECIES:
        if gene.startswith(prefix):
            return species
    return None


def detect_species(adata: ad.AnnData) -> str:
    """Auto-detect species from gene identifiers.

    Detection order:
    1. Ensembl ID prefixes (ENSG/ENSMUSG/ENSRNOG/ENSDARG/...) — the
       authoritative signal, supports human, mouse, rat, zebrafish,
       chicken, pig, cow, macaque, dog.
    2. Gene-symbol conventions: ALL-CAPS → human, Title-case → mouse.
    3. Ambiguous / mixed naming falls back to human with a warning.
    """
    var_names = [str(g) for g in adata.var_names[:1000]]
    if not var_names:
        logger.warning("Empty var_names; defaulting species to human")
        return SPECIES_HUMAN

    n_sampled = len(var_names)

    # 1. Ensembl prefix voting (longest-prefix match per gene)
    prefix_counts: dict[str, int] = {}
    for gene in var_names:
        sp = match_ensembl_species(gene)
        if sp is not None:
            prefix_counts[sp] = prefix_counts.get(sp, 0) + 1

    if prefix_counts:
        best_species, best_count = max(prefix_counts.items(), key=lambda kv: kv[1])
        if best_count >= n_sampled * SPECIES_DOMINANCE_RATIO:
            return best_species
        logger.warning(
            "Mixed Ensembl ID prefixes detected (%s); falling back to symbol conventions",
            prefix_counts,
        )

    # 2. Symbol conventions
    # Human symbols: all uppercase (e.g., CD3D, S100A8, HLA-DRA).
    # Do NOT require isalpha() — most real symbols contain digits.
    human_pattern = sum(
        1 for g in var_names if g and g == g.upper() and any(c.isalpha() for c in g)
    )
    # Mouse/rat symbols: first letter uppercase, rest lowercase (e.g., Cd3d).
    # Exclude ALL-CAPS genes — they match this shape too but are the human
    # convention, so counting them here would bias mixed datasets to mouse.
    mouse_pattern = sum(
        1
        for g in var_names
        if g and g[0].isupper() and (len(g) < 2 or g[1].islower()) and g != g.upper()
    )

    if human_pattern > mouse_pattern * SPECIES_SYMBOL_RATIO:
        return SPECIES_HUMAN
    if mouse_pattern > human_pattern * SPECIES_SYMBOL_RATIO:
        return SPECIES_MOUSE

    # 3. Ambiguous — default to human but flag it
    logger.warning(
        "Species detection ambiguous (human-like=%d, mouse-like=%d of %d genes); "
        "defaulting to human. Pass --species explicitly to override.",
        human_pattern,
        mouse_pattern,
        n_sampled,
    )
    return SPECIES_HUMAN


def _first_nonempty_value(series: pd.Series) -> str | None:
    """Return the first non-null, non-empty value of an obs column."""
    vals = series.dropna().unique()
    for v in vals:
        s = str(v).strip()
        if s and s.lower() != "nan":
            return s
    return None


def detect_tissue(adata: ad.AnnData) -> str | None:
    """Try to detect tissue from obs metadata.

    Matching is case-insensitive: first an exact synonym lookup
    (tissue, organ, organ_system, source, anatomy, body_site, ...),
    then a substring keyword scan of all obs columns.
    """
    col_map = {str(col).strip().lower(): col for col in adata.obs.columns}

    # 1. Exact (case-insensitive) synonym match, in priority order
    for synonym in TISSUE_COLUMN_SYNONYMS:
        if synonym in col_map:
            value = _first_nonempty_value(adata.obs[col_map[synonym]])
            if value:
                return value

    # 2. Substring keyword scan (e.g., "Tissue", "organ_system", "anatomy_region")
    for keyword in TISSUE_COLUMN_KEYWORDS:
        for col_lower, col in col_map.items():
            if keyword in col_lower:
                value = _first_nonempty_value(adata.obs[col])
                if value:
                    return value

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
    for key in adata.obsm:
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
        x_sample = adata.X[: min(100, adata.n_obs), : min(100, adata.n_vars)]
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
    cluster_key: str | None = None,
    embedding_key: str | None = None,
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
            report["warnings"].append(f"No --cluster-key specified. Candidates found: {candidates}")
        else:
            report["fatal"].append(
                "No cluster key found in obs. Cannot proceed without clustering."
            )

    # Embedding info
    if embedding_key and embedding_key in adata.obsm:
        report["embedding_shape"] = list(adata.obsm[embedding_key].shape)
    elif not embedding_key:
        candidates = find_embedding_keys(adata)
        if candidates:
            report["warnings"].append(
                f"No --embedding-key specified. Candidates found: {candidates}"
            )

    # Gene ID convention (majority vote over sampled genes)
    sample_genes = [str(g) for g in adata.var_names[:100]]
    ensembl_votes: dict[str, int] = {}
    for g in sample_genes:
        sp = match_ensembl_species(g)
        if sp is not None:
            ensembl_votes[sp] = ensembl_votes.get(sp, 0) + 1

    alpha_genes = [g for g in sample_genes if g.isalpha()]
    all_upper = bool(alpha_genes) and all(g == g.upper() for g in alpha_genes)
    first_cap = bool(alpha_genes) and all(
        g[0].isupper() and g[1:].islower() for g in alpha_genes if len(g) > 1
    )
    majority = len(sample_genes) // 2 + 1

    if ensembl_votes:
        best_species, best_count = max(ensembl_votes.items(), key=lambda kv: kv[1])
        if best_count >= majority:
            report["gene_id_convention"] = f"ensembl_{best_species}"
        else:
            report["gene_id_convention"] = "mixed_ensembl"
    elif all_upper:
        report["gene_id_convention"] = "human_symbols"
    elif first_cap:
        report["gene_id_convention"] = "mouse_symbols"
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
        "Layers:      "
        + (", ".join(str(name) for name in report["layers"]) if report["layers"] else "(none)"),
        f"  Counts:    {report['layer_info']['counts'] or '(not found)'}",
        f"  LogNorm:   {report['layer_info']['lognorm'] or '(not found)'}",
        "",
        "Obs columns: " + ", ".join(report["obs_columns"][:20]),
        "Obsm keys:   " + ", ".join(str(k) for k in report["obsm_keys"]),
        "",
        (
            "Cluster keys:    " + ", ".join(report["cluster_keys"])
            if report["cluster_keys"]
            else "Cluster keys:    (none found)"
        ),
        (
            "Embedding keys:   " + ", ".join(report["embedding_keys"])
            if report["embedding_keys"]
            else "Embedding keys:   (none found)"
        ),
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
    with open(ATLAS_PATH, encoding="utf-8") as f:
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

    for _tissue_key, tissue_data in mouse_atlas.get("tissues", {}).items():
        for _ct_key, ct_data in tissue_data.get("cell_types", {}).items():
            if "positive_markers" in ct_data:
                ct_data["positive_markers"] = [convert_gene(g) for g in ct_data["positive_markers"]]
            if "negative_markers" in ct_data:
                ct_data["negative_markers"] = [convert_gene(g) for g in ct_data["negative_markers"]]
            # Recurse into subtypes
            for _sub_key, sub_data in ct_data.get("subtypes", {}).items():
                if "positive_markers" in sub_data:
                    sub_data["positive_markers"] = [
                        convert_gene(g) for g in sub_data["positive_markers"]
                    ]
                if "negative_markers" in sub_data:
                    sub_data["negative_markers"] = [
                        convert_gene(g) for g in sub_data["negative_markers"]
                    ]

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
            result[sub_name] = {
                "positive_markers": sub_pos,
                "negative_markers": sub_neg,
                "cl_id": sub_cl,
            }

    return result


def get_all_markers_flat(atlas: dict, tissue: str) -> dict[str, list[str]]:
    """Get a flat mapping of cell_type → positive markers for quick scoring."""
    markers = get_all_markers_for_tissue(atlas, tissue)
    return {ct: info["positive_markers"] for ct, info in markers.items()}
