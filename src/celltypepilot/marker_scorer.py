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
    MARKER_FC_THRESHOLD,
    MARKER_FDR_THRESHOLD,
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
    rank_kwargs = {
        "groupby": cluster_key,
        "method": "wilcoxon",
        "n_genes": adata.n_vars,
        "pts": True,
        "use_raw": False,
    }
    if layer is not None:
        rank_kwargs["layer"] = layer
    sc.tl.rank_genes_groups(adata, **rank_kwargs)

    # Extract DE results
    de_results = _extract_de_results(adata, cluster_key)

    # Step 2: Score each cluster × cell_type combination
    results = []
    clusters = adata.obs[cluster_key].unique()

    for cluster in clusters:
        cluster_de = de_results.get(str(cluster), pd.DataFrame())
        if cluster_de.empty:
            continue

        significant_de = cluster_de[
            (pd.to_numeric(cluster_de["pval_adj"], errors="coerce") <= MARKER_FDR_THRESHOLD)
            & (pd.to_numeric(cluster_de["logfoldchange"], errors="coerce") >= MARKER_FC_THRESHOLD)
        ].copy()
        de_genes = set(significant_de["gene"].astype(str))
        de_genes_fc = dict(
            zip(
                significant_de["gene"].astype(str),
                pd.to_numeric(significant_de["logfoldchange"], errors="coerce"),
                strict=True,
            )
        )

        for ct_name, ct_info in markers.items():
            pos_markers = ct_info.get("positive_markers", [])
            neg_markers = ct_info.get("negative_markers", [])
            provenance_records = ct_info.get("marker_evidence", [])
            provenance_sources = sorted(
                {
                    f"{source.get('source_id', '')}:PMID:{source.get('pmid', '')}"
                    for record in provenance_records
                    for source in record.get("sources", [])
                }
            )

            # Positive marker analysis
            pos_present = [g for g in pos_markers if g in adata.var_names]
            pos_missing = [g for g in pos_markers if g not in adata.var_names]
            expression_pcts = _get_marker_expression_pcts(
                adata, cluster, cluster_key, pos_present, layer=layer
            )
            pos_expressed = [g for g, pct in expression_pcts.items() if pct >= MARKER_PCT_THRESHOLD]
            pos_silent = [g for g in pos_present if g not in pos_expressed]
            # Supporting markers pass direction, logFC, FDR, and expression gates.
            pos_de = [g for g in pos_expressed if g in de_genes]

            expected_denominator = max(len(pos_markers), 1)
            pct_overlap = len(pos_de) / expected_denominator
            mean_fc = np.mean([de_genes_fc.get(g, 0) for g in pos_de]) if pos_de else 0.0
            pct_expressed = len(pos_expressed) / expected_denominator

            # Specificity: how specific are these markers to this cluster?
            specificity = _compute_specificity(
                adata, cluster, cluster_key, pos_present, layer=layer
            )

            # Negative marker analysis
            neg_present = [g for g in neg_markers if g in adata.var_names]
            neg_expressed = _get_expressed_markers(
                adata, cluster, cluster_key, neg_present, layer=layer
            )
            neg_conflict = len(neg_expressed) / max(len(neg_markers), 1) if neg_markers else 0.0

            # Combined score
            combined = (
                0.35 * pct_overlap
                + 0.25 * max(0.0, min(mean_fc / 2.0, 1.0))  # normalize FC to [0,1]
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
                    "n_pos_markers": len(pos_markers),
                    "n_pos_present": len(pos_present),
                    "n_pos_missing": len(pos_missing),
                    "n_pos_silent": len(pos_silent),
                    "n_pos_de": len(pos_de),
                    "pos_present_markers": ";".join(pos_present),
                    "pos_missing_markers": ";".join(pos_missing),
                    "pos_silent_markers": ";".join(pos_silent),
                    "pos_supporting_markers": ";".join(pos_de),
                    "pct_overlap": round(pct_overlap, 4),
                    "mean_log2fc": round(mean_fc, 4),
                    "pct_expressed": round(pct_expressed, 4),
                    "specificity": round(specificity, 4),
                    "n_neg_markers": len(neg_markers),
                    "n_neg_present": len(neg_present),
                    "neg_expressed_markers": ";".join(neg_expressed),
                    "neg_conflict": round(neg_conflict, 4),
                    "combined_score": round(combined, 4),
                    "n_provenance_relationships": len(provenance_records),
                    "marker_provenance_status": ";".join(
                        sorted(
                            {
                                record.get("verification_status", "missing")
                                for record in provenance_records
                            }
                        )
                    ),
                    "marker_provenance_sources": ";".join(provenance_sources),
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
            result = sc.get.rank_genes_groups_df(adata, group=cluster_str, key="rank_genes_groups")
            df = pd.DataFrame(
                {
                    "gene": result["names"].to_numpy(),
                    "logfoldchange": result["logfoldchanges"].to_numpy(),
                    "pval": result["pvals"].to_numpy(),
                    "pval_adj": result["pvals_adj"].to_numpy(),
                }
            )
            results[cluster_str] = df
        except Exception:
            # Fallback: extract from stored results
            try:
                names = adata.uns["rank_genes_groups"]["names"]
                lfc = adata.uns["rank_genes_groups"]["logfoldchanges"]
                pvals = adata.uns["rank_genes_groups"]["pvals"]
                pvals_adj = adata.uns["rank_genes_groups"]["pvals_adj"]
                if (
                    isinstance(names, np.ndarray)
                    and names.dtype.names
                    and cluster_str in names.dtype.names
                ):
                    df = pd.DataFrame(
                        {
                            "gene": names[cluster_str],
                            "logfoldchange": lfc[cluster_str],
                            "pval": pvals[cluster_str],
                            "pval_adj": pvals_adj[cluster_str],
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
    layer: str | None = None,
) -> list[str]:
    """Get markers expressed in >MARKER_PCT_THRESHOLD of cells in a cluster."""
    pcts = _get_marker_expression_pcts(
        adata, cluster, cluster_key, markers, threshold=threshold, layer=layer
    )
    return [gene for gene, pct in pcts.items() if pct >= MARKER_PCT_THRESHOLD]


def _get_marker_expression_pcts(
    adata: ad.AnnData,
    cluster: str,
    cluster_key: str,
    markers: list[str],
    threshold: float = 0.0,
    layer: str | None = None,
) -> dict[str, float]:
    """Return within-cluster detectable-expression fractions for present markers."""
    mask = adata.obs[cluster_key].astype(str) == str(cluster)
    matrix = adata.layers[layer] if layer is not None else adata.X
    gene_idx = {str(g): i for i, g in enumerate(adata.var_names)}
    pcts: dict[str, float] = {}
    for gene in markers:
        idx = gene_idx.get(gene)
        if idx is None:
            continue
        expr = matrix[mask.values, idx]
        expr = expr.toarray().ravel() if sparse.issparse(expr) else np.asarray(expr).ravel()
        pcts[gene] = float(np.mean(expr > threshold)) if expr.size else 0.0
    return pcts


def _compute_specificity(
    adata: ad.AnnData,
    cluster: str,
    cluster_key: str,
    markers: list[str],
    layer: str | None = None,
) -> float:
    """Compute how specific markers are to this cluster vs others."""
    if not markers:
        return 0.0

    mask = adata.obs[cluster_key].astype(str) == str(cluster)
    matrix = adata.layers[layer] if layer is not None else adata.X
    gene_idx = {str(g): i for i, g in enumerate(adata.var_names)}

    specificities = []
    for gene in markers:
        if gene not in adata.var_names:
            continue
        idx = gene_idx[gene]
        expr = matrix[:, idx]
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

    return top1.drop(columns=["rank"], errors="ignore").reset_index(drop=True)
