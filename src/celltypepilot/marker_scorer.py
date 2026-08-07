"""Marker Scorer — deterministic marker-based scoring for cell type annotation.

This is the zero-cost, fully reproducible first path in the candidate generation engine.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from .constants import (
    CONFIDENCE_HIGH,
    CONFIDENCE_LOW,
    CONFIDENCE_MEDIUM,
    CONFIDENCE_REVIEW,
    MARKER_PCT_THRESHOLD,
)


def compute_marker_scores(
    adata: ad.AnnData,
    cluster_key: str,
    markers: dict[str, dict],
    layer: str | None = None,
) -> pd.DataFrame:
    """Score each cluster against each cell type using marker overlap.

    For each cluster:
    1. Compute DE genes (Wilcoxon rank-sum)
    2. For each candidate cell type, compute:
       - pct_overlap: fraction of positive markers that are DE and expressed
       - mean_expression: average expression of positive markers in the cluster
       - specificity: fraction of marker-expressing cells that are in this cluster
       - neg_conflict: fraction of negative markers that ARE expressed (penalty)

    Returns a DataFrame with columns:
        cluster, cell_type, pct_overlap, mean_expression, specificity,
        neg_conflict, combined_score, rank
    """
    # Step 1: Compute DE genes per cluster
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon", n_genes=adata.n_vars)

    # Extract DE results
    de_results = _extract_de_results(adata, cluster_key)

    # Step 2: Score each cluster × cell_type combination
    results = []
    clusters = adata.obs[cluster_key].unique()

    for cluster in clusters:
        cluster_de = de_results.get(str(cluster), pd.DataFrame())
        if cluster_de.empty:
            continue

        de_genes = set(cluster_de["gene"].values)
        de_genes_fc = dict(zip(cluster_de["gene"], cluster_de["logfoldchange"], strict=True))

        for ct_name, ct_info in markers.items():
            pos_markers = ct_info.get("positive_markers", [])
            neg_markers = ct_info.get("negative_markers", [])

            # Positive marker analysis
            pos_detected = [g for g in pos_markers if g in adata.var_names]
            pos_de = [g for g in pos_detected if g in de_genes]
            pos_expressed = _get_expressed_markers(adata, cluster, cluster_key, pos_detected)

            pct_overlap = len(pos_de) / max(len(pos_detected), 1)
            mean_fc = np.mean([de_genes_fc.get(g, 0) for g in pos_de]) if pos_de else 0.0
            pct_expressed = len(pos_expressed) / max(len(pos_detected), 1)

            # Specificity: how specific are these markers to this cluster?
            specificity = _compute_specificity(adata, cluster, cluster_key, pos_detected)

            # Negative marker analysis
            neg_detected = [g for g in neg_markers if g in adata.var_names]
            neg_expressed = _get_expressed_markers(adata, cluster, cluster_key, neg_detected)
            neg_conflict = len(neg_expressed) / max(len(neg_detected), 1) if neg_detected else 0.0

            # Combined score
            combined = (
                0.35 * pct_overlap
                + 0.25 * min(mean_fc / 2.0, 1.0)  # normalize FC to [0,1]
                + 0.20 * specificity
                + 0.20 * pct_expressed
                - 0.30 * neg_conflict  # penalty
            )
            combined = max(0.0, min(1.0, combined))

            results.append(
                {
                    "cluster": str(cluster),
                    "cell_type": ct_name,
                    "cl_id": ct_info.get("cl_id", ""),
                    "n_pos_markers": len(pos_detected),
                    "n_pos_de": len(pos_de),
                    "pct_overlap": round(pct_overlap, 4),
                    "mean_log2fc": round(mean_fc, 4),
                    "pct_expressed": round(pct_expressed, 4),
                    "specificity": round(specificity, 4),
                    "n_neg_markers": len(neg_detected),
                    "neg_conflict": round(neg_conflict, 4),
                    "combined_score": round(combined, 4),
                }
            )

    df = pd.DataFrame(results)
    if df.empty:
        return df

    # Rank cell types within each cluster
    df["rank"] = (
        df.groupby("cluster")["combined_score"].rank(ascending=False, method="first").astype(int)
    )
    df = df.sort_values(["cluster", "rank"])

    return df


def _extract_de_results(adata: ad.AnnData, cluster_key: str) -> dict[str, pd.DataFrame]:
    """Extract DE results from scanpy's rank_genes_groups into per-cluster DataFrames."""
    results = {}
    clusters = adata.obs[cluster_key].unique()

    for cluster in clusters:
        cluster_str = str(cluster)
        try:
            result = sc.get.rank_genes_groups(adata, group=cluster_str, key="rank_genes_groups")
            df = pd.DataFrame(
                {
                    "gene": result["names"],
                    "logfoldchange": result["logfoldchanges"],
                    "pval": result["pvals"],
                    "pval_adj": result["pvals_adj"],
                }
            )
            results[cluster_str] = df
        except Exception:
            # Fallback: extract from stored results
            try:
                names = adata.uns["rank_genes_groups"]["names"]
                lfc = adata.uns["rank_genes_groups"]["logfoldchanges"]
                if (
                    isinstance(names, np.ndarray)
                    and names.dtype.names
                    and cluster_str in names.dtype.names
                ):
                    df = pd.DataFrame(
                        {
                            "gene": names[cluster_str],
                            "logfoldchange": lfc[cluster_str],
                        }
                    )
                    results[cluster_str] = df
            except (KeyError, TypeError):
                pass

    return results


def _get_expressed_markers(
    adata: ad.AnnData,
    cluster: str,
    cluster_key: str,
    markers: list[str],
    threshold: float = 0.0,
) -> list[str]:
    """Get markers expressed in >MARKER_PCT_THRESHOLD of cells in a cluster."""
    mask = adata.obs[cluster_key] == cluster
    subset = adata[mask]
    expressed = []
    for gene in markers:
        if gene not in adata.var_names:
            continue
        idx = list(adata.var_names).index(gene)
        expr = subset.X[:, idx]
        expr = expr.toarray().flatten() if sparse.issparse(expr) else np.asarray(expr).flatten()
        pct = np.mean(expr > threshold)
        if pct >= MARKER_PCT_THRESHOLD:
            expressed.append(gene)
    return expressed


def _compute_specificity(
    adata: ad.AnnData,
    cluster: str,
    cluster_key: str,
    markers: list[str],
) -> float:
    """Compute how specific markers are to this cluster vs others."""
    if not markers:
        return 0.0

    mask = adata.obs[cluster_key] == cluster

    specificities = []
    for gene in markers:
        if gene not in adata.var_names:
            continue
        idx = list(adata.var_names).index(gene)
        expr = adata.X[:, idx]
        expr = expr.toarray().flatten() if sparse.issparse(expr) else np.asarray(expr).flatten()

        n_expressing = np.sum(expr > 0)
        n_in_cluster = np.sum((expr > 0) & mask.values)

        spec = n_in_cluster / n_expressing if n_expressing > 0 else 0.0
        specificities.append(spec)

    return np.mean(specificities) if specificities else 0.0


def assign_confidence(scores_row: pd.Series) -> str:
    """Assign confidence level based on scoring metrics."""
    score = scores_row.get("combined_score", 0)
    pct = scores_row.get("pct_overlap", 0)
    neg = scores_row.get("neg_conflict", 0)
    spec = scores_row.get("specificity", 0)

    if score >= 0.7 and pct >= 0.5 and neg < 0.1 and spec >= 0.5:
        return CONFIDENCE_HIGH
    elif score >= 0.5 and pct >= 0.3 and neg < 0.2:
        return CONFIDENCE_MEDIUM
    elif score >= 0.3:
        return CONFIDENCE_LOW
    else:
        return CONFIDENCE_REVIEW


def generate_annotation_summary(
    scores: pd.DataFrame,
    cluster_key: str,
) -> pd.DataFrame:
    """Generate per-cluster annotation summary from marker scores.

    Returns DataFrame with: cluster, cell_type, cl_id, combined_score, confidence, rank
    """
    if scores.empty:
        return pd.DataFrame()

    # Take top-1 for each cluster
    top1 = scores[scores["rank"] == 1].copy()
    top1["confidence"] = top1.apply(assign_confidence, axis=1)

    return top1[
        [
            "cluster",
            "cell_type",
            "cl_id",
            "combined_score",
            "confidence",
            "pct_overlap",
            "mean_log2fc",
            "specificity",
            "neg_conflict",
        ]
    ].reset_index(drop=True)
