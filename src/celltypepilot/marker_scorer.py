"""Marker Scorer — deterministic marker-based scoring for cell type annotation.

This is the zero-cost, fully reproducible first path in the candidate generation engine.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from .atlas_lifecycle import compute_marker_weights
from .constants import (
    CONFIDENCE_HIGH,
    CONFIDENCE_LOW,
    CONFIDENCE_MEDIUM,
    CONFIDENCE_REVIEW,
    DISEASE_FC_RELAXATION,
    MARKER_FC_THRESHOLD,
    MARKER_FDR_THRESHOLD,
    MARKER_PCT_THRESHOLD,
)
from .gene_aliases import build_var_alias_index, resolve_marker_list


def compute_marker_scores(
    adata: ad.AnnData,
    cluster_key: str,
    markers: dict[str, dict],
    layer: str | None = None,
    evidence_weighted: bool = False,
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
    # Detect disease/tumor microenvironment context from adata
    is_disease_context = False
    obs_text = (
        " ".join([str(col) for col in adata.obs.columns])
        + " "
        + " ".join([str(v) for v in adata.obs.iloc[0].values if isinstance(v, str)])
    )
    obs_text = obs_text.lower()
    if any(
        k in obs_text for k in ["tumor", "gbm", "cancer", "glioma", "lesion", "inflamed", "disease"]
    ):
        is_disease_context = True

    effective_fc_threshold = (
        MARKER_FC_THRESHOLD * DISEASE_FC_RELAXATION if is_disease_context else MARKER_FC_THRESHOLD
    )

    # Pre-build alias index
    alias_index = build_var_alias_index(adata.var_names)

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
            & (
                pd.to_numeric(cluster_de["logfoldchange"], errors="coerce")
                >= effective_fc_threshold
            )
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
            context_positive = set(ct_info.get("context_positive_markers", []))
            atlas_positive = set(ct_info.get("atlas_positive_markers", pos_markers))
            provenance_records = ct_info.get("marker_evidence", [])
            provenance_sources = sorted(
                {
                    f"{source.get('source_id', '')}:PMID:{source.get('pmid', '')}"
                    for record in provenance_records
                    for source in record.get("sources", [])
                }
            )

            # Positive marker analysis with alias resolution
            pos_present, pos_missing = resolve_marker_list(pos_markers, alias_index)
            expression_pcts = _get_marker_expression_pcts(
                adata, cluster, cluster_key, pos_present, layer=layer
            )
            pos_expressed = [g for g, pct in expression_pcts.items() if pct >= MARKER_PCT_THRESHOLD]
            pos_silent = [g for g in pos_present if g not in pos_expressed]

            # Supporting markers pass direction, logFC, FDR, and expression gates.
            pos_de = [g for g in pos_expressed if g in de_genes]
            context_supporting = [g for g in pos_de if g in context_positive]
            atlas_supporting = [g for g in pos_de if g in atlas_positive]

            gene_weights = compute_marker_weights(provenance_records)

            expected_denominator = max(len(pos_markers), 1)

            if evidence_weighted:
                weighted_numerator = sum(gene_weights.get(g, 0.5) for g in pos_de)
                weighted_denominator = sum(gene_weights.get(g, 0.5) for g in pos_markers)
                pct_overlap = weighted_numerator / max(weighted_denominator, 0.001)
                evidence_weight_mean = (
                    np.mean([gene_weights.get(g, 0.5) for g in pos_de]) if pos_de else 0.0
                )
            else:
                pct_overlap = len(pos_de) / expected_denominator
                evidence_weight_mean = 0.0

            mean_fc = np.mean([de_genes_fc.get(g, 0) for g in pos_de]) if pos_de else 0.0
            pct_expressed = len(pos_expressed) / expected_denominator

            # Specificity: how specific are these markers to this cluster?
            specificity = _compute_specificity(
                adata, cluster, cluster_key, pos_present, layer=layer
            )

            # Negative marker analysis with alias resolution
            neg_present, _ = resolve_marker_list(neg_markers, alias_index)
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
                    "context_origin": bool(ct_info.get("context_origin", False)),
                    "context_review_status": ct_info.get("context_review_status", "not_applicable"),
                    "context_supporting_markers": ";".join(context_supporting),
                    "n_context_supporting_markers": len(context_supporting),
                    "n_atlas_supporting_markers": len(atlas_supporting),
                    "context_only_support": bool(pos_de) and not atlas_supporting,
                    "evidence_weight_mean": round(evidence_weight_mean, 4),
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


def compute_profile_correlation_scores(
    adata: ad.AnnData,
    cluster_key: str,
    markers: dict[str, dict],
    layer: str | None = None,
) -> pd.DataFrame:
    """Compute an independent secondary score based on cluster mean profile correlation.

    This provides an internal secondary scoring method (zero external dependencies)
    to enable ensemble consensus voting and detect score instability.

    Returns DataFrame with: cluster, cell_type, ref_score
    """
    matrix = adata.layers[layer] if layer is not None else adata.X
    clusters = sorted([str(c) for c in adata.obs[cluster_key].unique()])
    gene_idx_map = {str(g): i for i, g in enumerate(adata.var_names)}
    alias_index = build_var_alias_index(adata.var_names)

    # Compute mean expression profile per cluster
    cluster_means: dict[str, np.ndarray] = {}
    for c in clusters:
        mask = (adata.obs[cluster_key].astype(str) == str(c)).values
        sub = matrix[mask]
        mean_profile = sub.mean(axis=0)
        mean_profile = (
            mean_profile.A1 if sparse.issparse(mean_profile) else np.asarray(mean_profile).ravel()
        )
        cluster_means[c] = mean_profile

    results = []
    for c in clusters:
        c_profile = cluster_means[c]
        for ct_name, ct_info in markers.items():
            pos_markers = ct_info.get("positive_markers", [])
            pos_present, _ = resolve_marker_list(pos_markers, alias_index)

            indices = [gene_idx_map[g] for g in pos_present if g in gene_idx_map]
            if not indices:
                results.append({"cluster": c, "cell_type": ct_name, "ref_score": 0.0})
                continue

            expr_vals = c_profile[indices]
            # Higher average expression of positive markers = stronger profile score
            profile_score = float(np.mean(expr_vals))
            # Squashing to [0, 1] range using soft sigmoid
            ref_score = float(1.0 / (1.0 + np.exp(-profile_score + 1.0)))
            results.append({"cluster": c, "cell_type": ct_name, "ref_score": round(ref_score, 4)})

    return pd.DataFrame(results)
