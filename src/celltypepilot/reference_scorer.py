"""Reference Embedding Scorer — deep learning reference mapping.

Complements the marker-based scorer by projecting query cells into a
reference embedding space and transferring labels via KNN or model
prediction. This captures continuous differentiation trajectories and
rare transitional states that pure marker overlap scoring misses.

Three backends (auto-selected by availability):

1. **CellTypist** (default) — pre-trained logistic regression models
   on CellxGene Census. Fast, no GPU needed, 500+ tissue models.
   Best for: standard human/mouse tissues.

2. **scANVI / scVI** — custom reference atlas with semi-supervised
   VAE. Projects query into shared latent space, KNN label transfer.
   Best for: custom atlases, cross-species, disease-specific refs.

3. **Correlation** — lightweight fallback using mean expression
   correlation with reference profiles. Works with any reference
   AnnData. No extra dependencies.

All backends output the same format: per-cluster probability vectors
over cell types, compatible with the ensemble scorer.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.stats import pearsonr

logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────
# Backend detection
# ──────────────────────────────────────────────

def _check_celltypist() -> bool:
    try:
        import celltypist
        return True
    except ImportError:
        return False


def _check_scvi() -> bool:
    try:
        import scvi
        return True
    except ImportError:
        return False


def check_reference_backends() -> dict:
    """Check which reference scoring backends are available.

    Returns:
        Dict with backend availability status.
    """
    return {
        "celltypist": _check_celltypist(),
        "scanvi": _check_scvi(),
        "correlation": True,  # Always available
    }


# ──────────────────────────────────────────────
# Core interface
# ──────────────────────────────────────────────

def score_by_reference(
    query: ad.AnnData,
    cluster_key: str,
    reference: Optional[ad.AnnData] = None,
    ref_label_key: str = "cell_type",
    model_path: Optional[str] = None,
    backend: str = "auto",
    n_neighbors: int = 15,
    gene_map: Optional[dict[str, str]] = None,
) -> pd.DataFrame:
    """Score clusters by reference embedding mapping.

    Projects query cells into a reference space and transfers labels.
    Returns per-cluster probability distributions over cell types.

    Args:
        query: Query AnnData (pre-clustered)
        cluster_key: Column in query.obs with cluster labels
        reference: Reference AnnData with cell type labels
        ref_label_key: Column in reference.obs with cell type labels
        model_path: Path to pre-trained model (CellTypist .pkl or scVI dir)
        backend: "celltypist", "scanvi", "correlation", or "auto"
        n_neighbors: Number of neighbors for KNN transfer
        gene_map: Optional gene name mapping {query_name: ref_name}

    Returns:
        DataFrame with columns:
            cluster, cell_type, ref_score, ref_rank,
            probability_vector (JSON), top5_types, top5_scores
    """
    if backend == "auto":
        backend = _auto_select_backend(reference, model_path)

    logger.info(f"Using reference backend: {backend}")

    if backend == "celltypist":
        return _score_celltypist(query, cluster_key, model_path, gene_map)
    elif backend == "scanvi":
        return _score_scanvi(query, cluster_key, reference, ref_label_key, gene_map)
    elif backend == "correlation":
        return _score_correlation(query, cluster_key, reference, ref_label_key, gene_map)
    else:
        raise ValueError(f"Unknown backend: {backend}. Use: celltypist, scanvi, correlation, auto")


def _auto_select_backend(
    reference: Optional[ad.AnnData],
    model_path: Optional[str],
) -> str:
    """Auto-select the best available backend."""
    # If model path provided → CellTypist
    if model_path and _check_celltypist():
        return "celltypist"

    # If reference AnnData + scvi available → scANVI
    if reference is not None and _check_scvi():
        return "scanvi"

    # If reference AnnData → correlation (always works)
    if reference is not None:
        return "correlation"

    # If CellTypist available → use default model
    if _check_celltypist():
        return "celltypist"

    raise RuntimeError(
        "No reference scoring backend available.\n"
        "Install one of:\n"
        "  pip install celltypist          # Pre-trained models (recommended)\n"
        "  pip install scvi-tools           # Custom reference atlas\n"
        "  Or provide a reference .h5ad     # Correlation-based (no extra deps)"
    )


# ──────────────────────────────────────────────
# Backend 1: CellTypist
# ──────────────────────────────────────────────

def _score_celltypist(
    query: ad.AnnData,
    cluster_key: str,
    model_path: Optional[str],
    gene_map: Optional[dict[str, str]],
) -> pd.DataFrame:
    """Score using CellTypist pre-trained models.

    CellTypist uses logistic regression trained on millions of cells.
    Outputs per-cell probability distributions over 500+ cell types.
    """
    import celltypist
    from celltypist import models

    # Apply gene mapping if needed
    query_work = _apply_gene_map(query, gene_map)

    # Load model
    if model_path:
        model = models.Model.load(model_path)
    else:
        # Use default immune or general model
        try:
            model = models.Model.load(model="Immune_All_Low.pkl")
        except Exception:
            try:
                model = models.Model.load(model="Pan_Immune_Low.pkl")
            except Exception:
                raise RuntimeError(
                    "No CellTypist model found. Download models:\n"
                    "  python -c 'import celltypist; celltypist.models.download_models()'\n"
                    "Or specify --model-path."
                )

    # Normalize query to match model expectations
    query_norm = query_work.copy()
    # CellTypist expects log1p-normalized, highly variable genes
    sc_t = __import__("scanpy")
    if "total_counts" not in query_norm.obs.columns:
        sc_t.pp.normalize_total(query_norm)
        sc_t.pp.log1p(query_norm)

    # Subset to model genes
    model_genes = [g for g in model.genes if g in query_norm.var_names]
    if len(model_genes) < 50:
        raise RuntimeError(
            f"Only {len(model_genes)} query genes overlap with model. "
            f"Check gene naming (HGNC symbols expected)."
        )
    query_norm = query_norm[:, model_genes]

    # Predict (chunked for large datasets)
    predictions = celltypist.annotate(
        query_norm,
        model=model,
        majority_voting=False,  # We want per-cell probabilities
    )

    # Extract probability matrix
    proba = predictions.probability_matrix  # cells × cell_types

    # Aggregate per cluster
    return _aggregate_cluster_probabilities(
        proba, query_work.obs[cluster_key], source="celltypist"
    )


# ──────────────────────────────────────────────
# Backend 2: scANVI / scVI
# ──────────────────────────────────────────────

def _score_scanvi(
    query: ad.AnnData,
    cluster_key: str,
    reference: ad.AnnData,
    ref_label_key: str,
    gene_map: Optional[dict[str, str]],
) -> pd.DataFrame:
    """Score using scANVI reference mapping.

    Trains scANVI on reference, projects query, transfers labels via
    soft assignment (posterior probability per cell type).
    """
    import scvi
    import scanpy as sc

    query_work = _apply_gene_map(query, gene_map)
    ref_work = _apply_gene_map(reference, gene_map)

    # Harmonize gene sets
    shared_genes = list(set(query_work.var_names) & set(ref_work.var_names))
    if len(shared_genes) < 200:
        raise RuntimeError(
            f"Only {len(shared_genes)} shared genes between query and reference. "
            f"Need at least 200 for scANVI mapping."
        )

    # Combine for joint training
    ref_sub = ref_work[:, shared_genes].copy()
    query_sub = query_work[:, shared_genes].copy()

    # Label the combined object
    ref_sub.obs["_ref_labels"] = ref_sub.obs[ref_label_key].values
    ref_sub.obs["_is_reference"] = True
    query_sub.obs["_ref_labels"] = "unknown"
    query_sub.obs["_is_reference"] = False

    combined = ad.concat([ref_sub, query_sub], merge="same")

    # Setup and train scANVI
    scvi.model.SCANVI.setup_anndata(
        combined,
        labels_key="_ref_labels",
        batch_key="_is_reference",
        unlabeled_category="unknown",
    )
    scvi_model = scvi.model.SCANVI(combined)
    scvi_model.train(max_epochs=20, check_val_every_n_epoch=5)

    # Get soft labels for query cells
    query_mask = combined.obs["_is_reference"] == False
    query_latent = scvi_model.get_latent_representation(combined[query_mask])
    soft_labels = scvi_model.predict(combined[query_mask])  # DataFrame: cells × cell_types

    # Aggregate per cluster
    return _aggregate_cluster_probabilities(
        soft_labels, query_work.obs[cluster_key], source="scanvi"
    )


# ──────────────────────────────────────────────
# Backend 3: Correlation (lightweight fallback)
# ──────────────────────────────────────────────

def _score_correlation(
    query: ad.AnnData,
    cluster_key: str,
    reference: ad.AnnData,
    ref_label_key: str,
    gene_map: Optional[dict[str, str]],
    n_top_genes: int = 2000,
) -> pd.DataFrame:
    """Score by correlation with reference mean expression profiles.

    For each cell type in the reference, computes the mean expression
    profile. Then for each query cluster, computes the Pearson
    correlation with each reference profile.

    This is a lightweight fallback that works with any reference
    AnnData without extra dependencies.
    """
    query_work = _apply_gene_map(query, gene_map)
    ref_work = _apply_gene_map(reference, gene_map)

    # Shared genes
    shared_genes = list(set(query_work.var_names) & set(ref_work.var_names))
    if len(shared_genes) < 100:
        raise RuntimeError(
            f"Only {len(shared_genes)} shared genes. Need at least 100."
        )

    # Compute reference profiles (mean expression per cell type)
    ref_sub = ref_work[:, shared_genes]
    ref_labels = ref_sub.obs[ref_label_key]
    unique_types = ref_labels.unique()

    ref_profiles = {}
    for ct in unique_types:
        mask = ref_labels == ct
        subset = ref_sub[mask]
        if sparse.issparse(subset.X):
            mean_expr = np.asarray(subset.X.mean(axis=0)).flatten()
        else:
            mean_expr = subset.X.mean(axis=0)
        ref_profiles[ct] = mean_expr

    # Compute query cluster profiles
    query_sub = query_work[:, shared_genes]
    clusters = query_sub.obs[cluster_key].unique()

    cluster_profiles = {}
    for cl in clusters:
        mask = query_sub.obs[cluster_key] == cl
        subset = query_sub[mask]
        if sparse.issparse(subset.X):
            mean_expr = np.asarray(subset.X.mean(axis=0)).flatten()
        else:
            mean_expr = subset.X.mean(axis=0)
        cluster_profiles[str(cl)] = mean_expr

    # Compute correlations
    results = []
    for cl_str, cl_profile in cluster_profiles.items():
        for ct, ref_profile in ref_profiles.items():
            # Pearson correlation
            corr, _ = pearsonr(cl_profile, ref_profile)
            corr = max(0.0, corr)  # Clip negative correlations to 0
            results.append({
                "cluster": cl_str,
                "cell_type": ct,
                "ref_score": round(corr, 4),
            })

    df = pd.DataFrame(results)
    if df.empty:
        return df

    # Rank within each cluster
    df["ref_rank"] = df.groupby("cluster")["ref_score"].rank(
        ascending=False, method="first"
    ).astype(int)
    df = df.sort_values(["cluster", "ref_rank"])

    # Add top-5 summary
    for cl_str in df["cluster"].unique():
        mask = df["cluster"] == cl_str
        top5 = df[mask].nsmallest(5, "ref_rank")
        top5_types = ";".join(top5["cell_type"].values)
        top5_scores = ";".join(f"{s:.3f}" for s in top5["ref_score"].values)
        df.loc[mask, "top5_types"] = top5_types
        df.loc[mask, "top5_scores"] = top5_scores

    return df


# ──────────────────────────────────────────────
# Shared utilities
# ──────────────────────────────────────────────

def _apply_gene_map(
    adata: ad.AnnData,
    gene_map: Optional[dict[str, str]],
) -> ad.AnnData:
    """Apply gene name mapping (e.g., mouse → human orthologs)."""
    if not gene_map:
        return adata

    new_var_names = [gene_map.get(g, g) for g in adata.var_names]
    adata = adata.copy()
    adata.var_names = new_var_names
    # Remove duplicates (keep first)
    adata = adata[:, ~adata.var_names.duplicated()]
    return adata


def _aggregate_cluster_probabilities(
    proba: pd.DataFrame,
    cluster_labels: pd.Series,
    source: str = "reference",
) -> pd.DataFrame:
    """Aggregate per-cell probability distributions to per-cluster.

    For each cluster, computes the mean probability for each cell type.

    Args:
        proba: DataFrame (cells × cell_types) with probabilities
        cluster_labels: Series with cluster assignments per cell
        source: Backend name for logging

    Returns:
        DataFrame in the same format as marker_scorer output
    """
    proba = proba.copy()
    proba.index = cluster_labels.values[:len(proba)]

    # Mean probability per cluster
    cluster_proba = proba.groupby(proba.index).mean()

    results = []
    for cl_str in cluster_proba.index:
        scores = cluster_proba.loc[cl_str].sort_values(ascending=False)
        for rank, (ct, score) in enumerate(scores.items(), 1):
            results.append({
                "cluster": str(cl_str),
                "cell_type": ct,
                "ref_score": round(float(score), 4),
                "ref_rank": rank,
            })

    df = pd.DataFrame(results)
    if df.empty:
        return df

    # Add top-5 summary
    for cl_str in df["cluster"].unique():
        mask = df["cluster"] == cl_str
        top5 = df[mask].nsmallest(5, "ref_rank")
        df.loc[mask, "top5_types"] = ";".join(top5["cell_type"].values)
        df.loc[mask, "top5_scores"] = ";".join(
            f"{s:.3f}" for s in top5["ref_score"].values
        )

    return df


# ──────────────────────────────────────────────
# Trajectory detection
# ──────────────────────────────────────────────

def detect_transitional_states(
    ref_scores: pd.DataFrame,
    marker_scores: pd.DataFrame,
    agreement_threshold: float = 0.3,
    top_n: int = 3,
) -> pd.DataFrame:
    """Detect clusters where marker and reference scorers disagree.

    Disagreement often indicates:
    - Continuous differentiation intermediates
    - Novel or rare cell states not in the reference
    - Technical artifacts requiring manual review

    Args:
        ref_scores: Output from score_by_reference()
        marker_scores: Output from compute_marker_scores()
        agreement_threshold: Min score difference to flag as disagreement
        top_n: Number of top candidates to compare

    Returns:
        DataFrame with disagreement analysis per cluster
    """
    if ref_scores.empty or marker_scores.empty:
        return pd.DataFrame()

    results = []
    clusters = marker_scores["cluster"].unique()

    for cl in clusters:
        # Get top-N from marker scorer
        m_top = marker_scores[marker_scores["cluster"] == cl].nsmallest(
            top_n, "rank"
        )
        # Get top-N from reference scorer
        r_top = ref_scores[ref_scores["cluster"] == cl].nsmallest(
            top_n, "ref_rank"
        )

        if m_top.empty or r_top.empty:
            continue

        m_best = m_top.iloc[0]
        r_best = r_top.iloc[0]

        m_score = float(m_top.iloc[0]["combined_score"])
        r_score = float(r_top.iloc[0]["ref_score"])
        m_type = m_top.iloc[0]["cell_type"]
        r_type = r_top.iloc[0]["cell_type"]

        # Agreement check
        same_type = m_type == r_type
        score_gap = abs(m_score - r_score)

        # Detect transitional signal:
        # Reference scorer gives high confidence to a DIFFERENT type
        # than marker scorer, or marker scorer has low confidence
        # but reference is confident
        is_transitional = False
        transition_type = ""

        if not same_type:
            # Check if reference top-2 includes marker's top-1 (and vice versa)
            r_top_types = set(r_top["cell_type"].values[:top_n])
            m_top_types = set(m_top["cell_type"].values[:top_n])

            if m_type in r_top_types and r_type in m_top_types:
                is_transitional = True
                transition_type = "mixed_signal"
            elif r_score > 0.5 and m_score < 0.4:
                is_transitional = True
                transition_type = "reference_override"
            elif m_score > 0.5 and r_score < 0.3:
                is_transitional = True
                transition_type = "marker_only"

        # Check for multi-lineage signal (broad distribution)
        if not r_top.empty and len(r_top) >= 2:
            top1_score = float(r_top.iloc[0]["ref_score"])
            top2_score = float(r_top.iloc[1]["ref_score"])
            if top1_score > 0.2 and top2_score > 0.2:
                entropy = -sum(
                    float(r_top.iloc[i]["ref_score"]) *
                    np.log2(float(r_top.iloc[i]["ref_score"]) + 1e-10)
                    for i in range(min(5, len(r_top)))
                )
                if entropy > 1.5:  # High entropy = diffuse distribution
                    is_transitional = True
                    transition_type = transition_type or "diffuse_distribution"

        results.append({
            "cluster": cl,
            "marker_type": m_type,
            "marker_score": round(m_score, 4),
            "ref_type": r_type,
            "ref_score": round(r_score, 4),
            "agreement": same_type,
            "score_gap": round(score_gap, 4),
            "is_transitional": is_transitional,
            "transition_type": transition_type,
            "ref_top5": r_top.iloc[0].get("top5_types", ""),
        })

    return pd.DataFrame(results)
