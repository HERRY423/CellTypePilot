"""Reference Embedding Scorer — deep learning reference mapping.

Complements the marker-based scorer by projecting query cells into a
reference embedding space and transferring labels via KNN or model
prediction. This captures continuous differentiation trajectories and
rare transitional states that pure marker overlap scoring misses.

Four backends (auto-selected by availability):

1. **CellTypist** (default) — pre-trained logistic regression models
   on CellxGene Census. Fast, no GPU needed, 500+ tissue models.
   Best for: standard human/mouse tissues.

2. **scANVI / scVI** — custom reference atlas with semi-supervised
   VAE. Projects query into shared latent space, KNN label transfer.
   Best for: custom atlases, cross-species, disease-specific refs.

3. **KNN** — direct KNN label transfer on PCA/cosine similarity
   between query and reference. Lightweight, no model training.
   Best for: quick reference mapping without deep learning.

4. **Correlation** — lightweight fallback using mean expression
   correlation with reference profiles. Works with any reference
   AnnData. No extra dependencies.

All backends output the same format: per-cluster probability vectors
over cell types, compatible with the ensemble scorer.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse
from scipy.stats import pearsonr
from sklearn.neighbors import NearestNeighbors

from .constants import (
    REF_MIN_SHARED_GENES, REF_CORR_MIN_GENES,
    REF_SCANVI_MIN_GENES, REF_CELLTYPIST_MIN_GENES,
    REF_KNN_DEFAULT_K, REF_KNN_MAX_K,
)

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
        "knn": True,       # Always available (sklearn)
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
    n_neighbors: int = REF_KNN_DEFAULT_K,
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
        backend: "celltypist", "scanvi", "knn", "correlation", or "auto"
        n_neighbors: Number of neighbors for KNN transfer
        gene_map: Optional gene name mapping {query_name: ref_name}

    Returns:
        DataFrame with columns:
            cluster, cell_type, ref_score, ref_rank, top5_types, top5_scores
    """
    if backend == "auto":
        backend = _auto_select_backend(reference, model_path)

    logger.info(f"Using reference backend: {backend}")

    if backend == "celltypist":
        return _score_celltypist(query, cluster_key, model_path, gene_map)
    elif backend == "scanvi":
        return _score_scanvi(query, cluster_key, reference, ref_label_key, gene_map)
    elif backend == "knn":
        return _score_knn(query, cluster_key, reference, ref_label_key,
                          gene_map, n_neighbors)
    elif backend == "correlation":
        return _score_correlation(query, cluster_key, reference, ref_label_key, gene_map)
    else:
        raise ValueError(
            f"Unknown backend: {backend}. "
            f"Use: celltypist, scanvi, knn, correlation, auto"
        )


def _auto_select_backend(
    reference: Optional[ad.AnnData],
    model_path: Optional[str],
) -> str:
    """Auto-select the best available backend.

    Priority: CellTypist > scANVI > KNN > correlation
    """
    # If model path provided → CellTypist
    if model_path and _check_celltypist():
        return "celltypist"

    # If reference AnnData + scvi available → scANVI
    if reference is not None and _check_scvi():
        return "scanvi"

    # If reference AnnData → KNN (lightweight, better than correlation)
    if reference is not None:
        return "knn"

    # If CellTypist available → use default model
    if _check_celltypist():
        return "celltypist"

    raise RuntimeError(
        "No reference scoring backend available.\n"
        "Install one of:\n"
        "  pip install celltypist          # Pre-trained models (recommended)\n"
        "  pip install scvi-tools           # Custom reference atlas\n"
        "  Or provide a reference .h5ad     # KNN/correlation (no extra deps)"
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

    query_work = _apply_gene_map(query, gene_map)

    # Load model
    model = _load_celltypist_model(models, model_path)

    # Normalize query to match model expectations
    query_norm = query_work.copy()
    if "total_counts" not in query_norm.obs.columns:
        sc.pp.normalize_total(query_norm)
        sc.pp.log1p(query_norm)

    # Subset to model genes
    model_genes = [g for g in model.genes if g in query_norm.var_names]
    if len(model_genes) < REF_CELLTYPIST_MIN_GENES:
        raise RuntimeError(
            f"Only {len(model_genes)} query genes overlap with model "
            f"(need ≥{REF_CELLTYPIST_MIN_GENES}). "
            f"Check gene naming (HGNC symbols expected)."
        )
    query_norm = query_norm[:, model_genes]

    # Predict
    predictions = celltypist.annotate(
        query_norm, model=model, majority_voting=False,
    )

    # Extract probability matrix (cells × cell_types)
    proba = predictions.probability_matrix

    return _aggregate_cluster_probabilities(
        proba, query_work.obs[cluster_key], source="celltypist"
    )


def _load_celltypist_model(models, model_path: Optional[str]):
    """Load CellTypist model with fallback chain."""
    if model_path:
        return models.Model.load(model_path)

    # Try common default models
    for name in ["Immune_All_Low.pkl", "Pan_Immune_Low.pkl"]:
        try:
            return models.Model.load(model=name)
        except Exception:
            continue

    raise RuntimeError(
        "No CellTypist model found. Download models:\n"
        "  python -c 'import celltypist; celltypist.models.download_models()'\n"
        "Or specify --model-path."
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

    query_work = _apply_gene_map(query, gene_map)
    ref_work = _apply_gene_map(reference, gene_map)

    # Harmonize gene sets
    shared_genes = sorted(set(query_work.var_names) & set(ref_work.var_names))
    if len(shared_genes) < REF_SCANVI_MIN_GENES:
        raise RuntimeError(
            f"Only {len(shared_genes)} shared genes between query and reference "
            f"(need ≥{REF_SCANVI_MIN_GENES})."
        )

    # Combine for joint training
    ref_sub = ref_work[:, shared_genes].copy()
    query_sub = query_work[:, shared_genes].copy()

    ref_sub.obs["_ref_labels"] = ref_sub.obs[ref_label_key].astype(str).values
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

    # Get soft labels for query cells only
    query_mask = ~combined.obs["_is_reference"].values
    soft_labels = scvi_model.predict(combined[query_mask])

    return _aggregate_cluster_probabilities(
        soft_labels, query_work.obs[cluster_key], source="scanvi"
    )


# ──────────────────────────────────────────────
# Backend 3: KNN label transfer
# ──────────────────────────────────────────────

def _score_knn(
    query: ad.AnnData,
    cluster_key: str,
    reference: ad.AnnData,
    ref_label_key: str,
    gene_map: Optional[dict[str, str]],
    n_neighbors: int = REF_KNN_DEFAULT_K,
    n_hvg: int = 2000,
) -> pd.DataFrame:
    """Score by KNN label transfer on PCA space.

    Computes PCA on the combined reference + query (using reference HVGs),
    then for each query cell finds its K nearest reference neighbors and
    transfers labels weighted by inverse distance.

    Memory-safe: operates on PCA-reduced space, not full expression.
    """
    query_work = _apply_gene_map(query, gene_map)
    ref_work = _apply_gene_map(reference, gene_map)

    # Shared genes
    shared_genes = sorted(set(query_work.var_names) & set(ref_work.var_names))
    if len(shared_genes) < REF_MIN_SHARED_GENES:
        raise RuntimeError(
            f"Only {len(shared_genes)} shared genes (need ≥{REF_MIN_SHARED_GENES})."
        )

    ref_sub = ref_work[:, shared_genes].copy()
    query_sub = query_work[:, shared_genes].copy()

    # Normalize both
    sc.pp.normalize_total(ref_sub)
    sc.pp.log1p(ref_sub)
    sc.pp.normalize_total(query_sub)
    sc.pp.log1p(query_sub)

    # Select HVGs from reference (or use shared genes if few)
    if len(shared_genes) > n_hvg:
        sc.pp.highly_variable_genes(ref_sub, n_top_genes=n_hvg)
        hvg = ref_sub.var_names[ref_sub.var["highly_variable"]].tolist()
    else:
        hvg = shared_genes

    ref_sub = ref_sub[:, hvg]
    query_sub = query_sub[:, hvg]

    # PCA on reference
    sc.pp.pca(ref_sub, n_comps=min(50, len(hvg) - 1))

    # Project query into reference PCA space
    # Use the reference PCA loadings to transform query
    query_pca = _project_to_pca(query_sub, ref_sub.uns["pca"]["variance_ratio"],
                                 ref_sub.varm["PCs"])

    # KNN: for each query cell, find nearest reference cells
    k = min(n_neighbors, ref_sub.n_obs)
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean", algorithm="auto")
    nn.fit(ref_sub.obsm["X_pca"])
    distances, indices = nn.kneighbors(query_pca)

    # Weighted vote by inverse distance
    ref_labels = ref_sub.obs[ref_label_key].values
    unique_types = sorted(set(ref_labels))
    type_to_idx = {ct: i for i, ct in enumerate(unique_types)}

    # Build per-cell probability vectors
    n_query = query_sub.n_obs
    n_types = len(unique_types)
    proba_matrix = np.zeros((n_query, n_types))

    for cell_i in range(n_query):
        for neighbor_j in range(k):
            dist = distances[cell_i, neighbor_j]
            weight = 1.0 / (dist + 1e-8)  # Inverse distance
            ct = ref_labels[indices[cell_i, neighbor_j]]
            proba_matrix[cell_i, type_to_idx[ct]] += weight

    # Normalize rows to sum to 1
    row_sums = proba_matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    proba_matrix /= row_sums

    # Convert to DataFrame
    proba_df = pd.DataFrame(proba_matrix, columns=unique_types, index=query_sub.obs_names)

    return _aggregate_cluster_probabilities(
        proba_df, query_work.obs[cluster_key], source="knn"
    )


def _project_to_pca(
    query_adata: ad.AnnData,
    variance_ratio: np.ndarray,
    pcs: np.ndarray,
) -> np.ndarray:
    """Project query data into reference PCA space.

    Args:
        query_adata: Normalized + log1p query (already subset to HVGs)
        variance_ratio: From reference .uns["pca"]["variance_ratio"]
        pcs: From reference .varm["PCs"] — genes × components

    Returns:
        PCA coordinates for query cells (n_cells × n_components)
    """
    X = query_adata.X
    if sparse.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float64)

    # Center using reference mean (approximate: use query mean as proxy)
    # In production, store reference means; here we center per-gene
    X_centered = X - X.mean(axis=0)

    # Project
    return X_centered @ pcs


# ──────────────────────────────────────────────
# Backend 4: Correlation (lightweight fallback)
# ──────────────────────────────────────────────

def _score_correlation(
    query: ad.AnnData,
    cluster_key: str,
    reference: ad.AnnData,
    ref_label_key: str,
    gene_map: Optional[dict[str, str]],
    n_top_genes: int = 2000,
) -> pd.DataFrame:
    """Score by Pearson correlation with reference mean expression profiles.

    For each cell type in the reference, computes the mean expression
    profile. Then for each query cluster, computes the Pearson
    correlation with each reference profile.

    Uses HVG selection from reference to focus on informative genes
    and reduce memory for large datasets.
    """
    query_work = _apply_gene_map(query, gene_map)
    ref_work = _apply_gene_map(reference, gene_map)

    # Shared genes
    shared_genes = sorted(set(query_work.var_names) & set(ref_work.var_names))
    if len(shared_genes) < REF_CORR_MIN_GENES:
        raise RuntimeError(
            f"Only {len(shared_genes)} shared genes (need ≥{REF_CORR_MIN_GENES})."
        )

    # Subset to shared genes
    ref_sub = ref_work[:, shared_genes].copy()
    query_sub = query_work[:, shared_genes].copy()

    # Normalize both
    sc.pp.normalize_total(ref_sub)
    sc.pp.log1p(ref_sub)
    sc.pp.normalize_total(query_sub)
    sc.pp.log1p(query_sub)

    # Select top genes for correlation (reduce noise + memory)
    use_genes = shared_genes
    if len(shared_genes) > n_top_genes:
        sc.pp.highly_variable_genes(ref_sub, n_top_genes=n_top_genes)
        hvg = ref_sub.var_names[ref_sub.var["highly_variable"]].tolist()
        if len(hvg) >= REF_CORR_MIN_GENES:
            use_genes = hvg
            ref_sub = ref_sub[:, use_genes]
            query_sub = query_sub[:, use_genes]

    # Compute reference profiles (mean expression per cell type)
    ref_labels = ref_sub.obs[ref_label_key]
    unique_types = ref_labels.unique()

    ref_profiles = {}
    for ct in unique_types:
        mask = (ref_labels == ct).values
        X = ref_sub[mask].X
        if sparse.issparse(X):
            ref_profiles[ct] = np.asarray(X.mean(axis=0)).flatten()
        else:
            ref_profiles[ct] = np.asarray(X.mean(axis=0)).flatten()

    # Compute query cluster profiles
    clusters = query_sub.obs[cluster_key].unique()
    cluster_profiles = {}
    for cl in clusters:
        mask = (query_sub.obs[cluster_key] == cl).values
        X = query_sub[mask].X
        if sparse.issparse(X):
            cluster_profiles[str(cl)] = np.asarray(X.mean(axis=0)).flatten()
        else:
            cluster_profiles[str(cl)] = np.asarray(X.mean(axis=0)).flatten()

    # Compute correlation matrix (vectorized)
    n_clusters = len(cluster_profiles)
    n_types = len(ref_profiles)
    cl_names = list(cluster_profiles.keys())
    ct_names = list(ref_profiles.keys())

    # Stack into matrices for vectorized correlation
    ref_mat = np.column_stack([ref_profiles[ct] for ct in ct_names])  # genes × types
    cl_mat = np.column_stack([cluster_profiles[cl] for cl in cl_names])  # genes × clusters

    # Vectorized Pearson correlation via z-scored dot product
    ref_z = (ref_mat - ref_mat.mean(axis=0)) / (ref_mat.std(axis=0) + 1e-10)
    cl_z = (cl_mat - cl_mat.mean(axis=0)) / (cl_mat.std(axis=0) + 1e-10)
    corr_mat = (cl_z.T @ ref_z) / len(use_genes)  # clusters × types

    # Clip negatives to 0
    corr_mat = np.maximum(corr_mat, 0.0)

    # Build results DataFrame
    results = []
    for i, cl_str in enumerate(cl_names):
        for j, ct in enumerate(ct_names):
            results.append({
                "cluster": cl_str,
                "cell_type": ct,
                "ref_score": round(float(corr_mat[i, j]), 4),
            })

    df = pd.DataFrame(results)
    if df.empty:
        return df

    # Rank within each cluster
    df["ref_rank"] = df.groupby("cluster")["ref_score"].rank(
        ascending=False, method="first"
    ).astype(int)
    df = df.sort_values(["cluster", "ref_rank"])

    # Vectorized top-5 summary
    df = _add_top5_summary(df)

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
        DataFrame with: cluster, cell_type, ref_score, ref_rank, top5_*
    """
    proba = proba.copy()
    proba.index = cluster_labels.values[:len(proba)].astype(str)

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

    # Vectorized top-5 summary
    df = _add_top5_summary(df)

    return df


def _add_top5_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Add top5_types and top5_scores columns (vectorized)."""
    # Get top-5 per cluster
    top5 = df.nsmallest(5, "ref_rank") if len(df) >= 5 else df.sort_values("ref_rank")

    top5_agg = top5.groupby("cluster").agg(
        top5_types=("cell_type", lambda x: ";".join(x.values)),
        top5_scores=("ref_score", lambda x: ";".join(f"{v:.3f}" for v in x.values)),
    )

    # Merge back
    df = df.merge(top5_agg, left_on="cluster", right_index=True, how="left")
    return df


# ──────────────────────────────────────────────
# Trajectory / transitional state detection
# ──────────────────────────────────────────────

def detect_transitional_states(
    ref_scores: pd.DataFrame,
    marker_scores: pd.DataFrame,
    top_n: int = 3,
) -> pd.DataFrame:
    """Detect clusters where marker and reference scorers disagree.

    Disagreement often indicates:
    - Continuous differentiation intermediates
    - Novel or rare cell states not in the reference
    - Technical artifacts requiring manual review

    Uses three signals:
    1. Cross-ranking: each method's top-1 appears in the other's top-N
    2. Confidence asymmetry: one method confident, other uncertain
    3. Distribution entropy: broad probability = diffuse identity

    Args:
        ref_scores: Output from score_by_reference()
        marker_scores: Output from compute_marker_scores()
        top_n: Number of top candidates to compare

    Returns:
        DataFrame with disagreement analysis per cluster
    """
    if ref_scores.empty or marker_scores.empty:
        return pd.DataFrame()

    results = []
    clusters = marker_scores["cluster"].unique()

    for cl in clusters:
        m_top = marker_scores[marker_scores["cluster"] == cl].nsmallest(top_n, "rank")
        r_top = ref_scores[ref_scores["cluster"] == cl].nsmallest(top_n, "ref_rank")

        if m_top.empty or r_top.empty:
            continue

        m_score = float(m_top.iloc[0]["combined_score"])
        r_score = float(r_top.iloc[0]["ref_score"])
        m_type = m_top.iloc[0]["cell_type"]
        r_type = r_top.iloc[0]["cell_type"]

        same_type = m_type == r_type
        score_gap = abs(m_score - r_score)

        # ── Signal 1: Cross-ranking ──
        is_transitional = False
        transition_type = ""

        if not same_type:
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

        # ── Signal 2: Distribution entropy ──
        if not r_top.empty and len(r_top) >= 2:
            top_scores = r_top["ref_score"].values[:min(5, len(r_top))].astype(float)
            top_scores = top_scores / (top_scores.sum() + 1e-10)  # Normalize
            entropy = -np.sum(top_scores * np.log2(top_scores + 1e-10))
            if entropy > 1.5:
                is_transitional = True
                transition_type = transition_type or "diffuse_distribution"

        # ── Signal 3: Confidence asymmetry ──
        if not same_type and not is_transitional:
            if r_score - m_score > 0.3 and r_score > 0.5:
                is_transitional = True
                transition_type = "asymmetric_confidence"

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
