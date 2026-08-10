"""Bootstrap resampling for confidence intervals and stability analysis.

Provides non-parametric bootstrap infrastructure for:
- Metric confidence intervals (accuracy, F1, etc.)
- Cluster assignment stability under cell resampling
- Evidence score confidence intervals per cluster
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from dataclasses import dataclass

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class BootstrapResult:
    """Result of a bootstrap confidence interval computation."""

    point_estimate: float
    ci_lower: float
    ci_upper: float
    se: float
    n_boot: int
    ci_level: float

    def __repr__(self) -> str:
        return (
            f"BootstrapResult({self.point_estimate:.4f} "
            f"[{self.ci_lower:.4f}, {self.ci_upper:.4f}], "
            f"SE={self.se:.4f}, n_boot={self.n_boot})"
        )


def bootstrap_metric_ci(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    metric_fn: Callable[[np.ndarray, np.ndarray], float],
    n_boot: int = 1000,
    ci: float = 0.95,
    seed: int = 42,
) -> BootstrapResult:
    """Bootstrap confidence intervals for a scalar metric.

    Resamples (y_true, y_pred) pairs with replacement, computes
    metric_fn on each resample, and returns percentile CI.

    Parameters
    ----------
    y_true : array-like of true labels
    y_pred : array-like of predicted labels
    metric_fn : callable (y_true, y_pred) -> float
    n_boot : number of bootstrap replicates
    ci : confidence level (e.g. 0.95 for 95% CI)
    seed : random seed for reproducibility
    """
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    n_samples = len(y_true)

    if n_samples == 0:
        return BootstrapResult(np.nan, np.nan, np.nan, np.nan, n_boot, ci)

    point_estimate = metric_fn(y_true, y_pred)
    boot_stats = []

    for _ in range(n_boot):
        indices = rng.choice(n_samples, size=n_samples, replace=True)
        boot_y_true = y_true[indices]
        boot_y_pred = y_pred[indices]
        try:
            val = metric_fn(boot_y_true, boot_y_pred)
            boot_stats.append(val)
        except Exception:
            pass

    if not boot_stats:
        return BootstrapResult(point_estimate, np.nan, np.nan, np.nan, n_boot, ci)

    boot_stats = np.array(boot_stats)
    alpha = 1.0 - ci
    ci_lower = float(np.percentile(boot_stats, 100 * (alpha / 2)))
    ci_upper = float(np.percentile(boot_stats, 100 * (1 - alpha / 2)))
    se = float(np.std(boot_stats, ddof=1))

    return BootstrapResult(point_estimate, ci_lower, ci_upper, se, n_boot, ci)


def grouped_bootstrap_metric_ci(
    values: np.ndarray,
    groups: np.ndarray,
    *,
    strata: np.ndarray | None = None,
    n_boot: int = 2000,
    ci: float = 0.95,
    seed: int = 42,
) -> BootstrapResult:
    """Bootstrap a mean at the biological-unit level.

    ``values`` must already contain one observation per independent unit (normally
    one donor).  Groups are resampled, never cells.  When ``strata`` is supplied,
    strata (normally studies) are resampled first and units are then resampled
    within each selected stratum.  This prevents a large donor or study from
    acquiring artificial precision simply because it contributed more cells.
    """
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups, dtype=str)
    if values.shape != groups.shape:
        raise ValueError("values and groups must have the same shape")
    if strata is not None and np.asarray(strata).shape != values.shape:
        raise ValueError("strata must have the same shape as values")
    if n_boot < 1:
        raise ValueError("n_boot must be at least 1")
    if not 0 < ci < 1:
        raise ValueError("ci must be between 0 and 1")

    frame = pd.DataFrame({"value": values, "group": groups}).dropna(subset=["value"])
    if strata is not None:
        frame["stratum"] = np.asarray(strata, dtype=str)[frame.index]
    if frame.empty:
        return BootstrapResult(np.nan, np.nan, np.nan, np.nan, n_boot, ci)
    if frame["group"].duplicated().any():
        raise ValueError("grouped bootstrap requires one row per independent group")

    point_estimate = float(frame["value"].mean())
    rng = np.random.default_rng(seed)
    boot_stats: list[float] = []
    if strata is None:
        group_values = frame["value"].to_numpy()
        for _ in range(n_boot):
            boot_stats.append(
                float(rng.choice(group_values, len(group_values), replace=True).mean())
            )
    else:
        stratum_names = frame["stratum"].drop_duplicates().to_numpy()
        by_stratum = {
            name: group["value"].to_numpy() for name, group in frame.groupby("stratum", sort=False)
        }
        for _ in range(n_boot):
            sampled_strata = rng.choice(stratum_names, len(stratum_names), replace=True)
            sampled_units = []
            for name in sampled_strata:
                unit_values = by_stratum[name]
                sampled_units.extend(
                    rng.choice(unit_values, len(unit_values), replace=True).tolist()
                )
            boot_stats.append(float(np.mean(sampled_units)))

    samples = np.asarray(boot_stats, dtype=float)
    alpha = 1.0 - ci
    return BootstrapResult(
        point_estimate=point_estimate,
        ci_lower=float(np.quantile(samples, alpha / 2)),
        ci_upper=float(np.quantile(samples, 1 - alpha / 2)),
        se=float(np.std(samples, ddof=1)) if len(samples) > 1 else 0.0,
        n_boot=n_boot,
        ci_level=ci,
    )


def bootstrap_cluster_stability(
    adata: ad.AnnData,
    cluster_key: str,
    n_boot: int = 200,
    subsample_frac: float = 0.8,
    seed: int = 42,
) -> pd.DataFrame:
    """Assess cluster stability via subsampling and re-clustering.

    For each bootstrap replicate:
    1. Subsample cells (without replacement) at subsample_frac
    2. Re-cluster with leiden at the same resolution
    3. Compute ARI of each original cluster vs best-matching boot cluster

    Requires >= 50 cells per cluster for reliable stability estimates.
    Clusters with fewer cells get NaN stability scores.
    """
    from sklearn.metrics import adjusted_rand_score

    rng = np.random.default_rng(seed)
    clusters = adata.obs[cluster_key].unique()
    results = []

    # Determine leiden resolution from uns if available
    resolution = 1.0
    if "leiden" in adata.uns and isinstance(adata.uns["leiden"], dict):
        params = adata.uns["leiden"].get("params", {})
        resolution = params.get("resolution", 1.0)

    n_cells = len(adata)
    n_subsample = int(n_cells * subsample_frac)

    if n_subsample < 50:
        logger.warning(
            "Fewer than 50 cells after subsampling; cluster stability analysis may be unreliable."
        )

    for cl in clusters:
        mask = adata.obs[cluster_key] == cl
        n_cluster_cells = int(mask.sum())

        if n_cluster_cells < 50:
            results.append(
                {
                    "cluster": cl,
                    "stability_score": np.nan,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                    "n_cells": n_cluster_cells,
                    "n_boot_successful": 0,
                }
            )
            continue

        boot_scores = []
        for _ in range(n_boot):
            try:
                indices = rng.choice(n_cells, size=n_subsample, replace=False)
                sub_adata = adata[indices].copy()

                # Require neighbors graph for leiden
                if "neighbors" not in sub_adata.uns:
                    import scanpy as sc

                    sc.pp.neighbors(
                        sub_adata, use_rep="X_pca" if "X_pca" in sub_adata.obsm else None
                    )

                import scanpy as sc

                sc.tl.leiden(sub_adata, resolution=resolution, key_added="boot_cluster")

                # Binary ARI: this-cluster vs rest in both original and boot labels
                orig_labels = (sub_adata.obs[cluster_key] == cl).astype(int).values

                best_ari = -1.0
                for boot_cl in sub_adata.obs["boot_cluster"].unique():
                    boot_labels = (sub_adata.obs["boot_cluster"] == boot_cl).astype(int).values
                    ari = adjusted_rand_score(orig_labels, boot_labels)
                    best_ari = max(best_ari, ari)

                boot_scores.append(best_ari)
            except Exception as e:
                logger.debug("Bootstrap re-clustering failed: %s", e)

        if boot_scores:
            scores = np.array(boot_scores)
            results.append(
                {
                    "cluster": cl,
                    "stability_score": float(np.mean(scores)),
                    "ci_lower": float(np.percentile(scores, 2.5)),
                    "ci_upper": float(np.percentile(scores, 97.5)),
                    "n_cells": n_cluster_cells,
                    "n_boot_successful": len(boot_scores),
                }
            )
        else:
            results.append(
                {
                    "cluster": cl,
                    "stability_score": np.nan,
                    "ci_lower": np.nan,
                    "ci_upper": np.nan,
                    "n_cells": n_cluster_cells,
                    "n_boot_successful": 0,
                }
            )

    return pd.DataFrame(results)


def bootstrap_evidence_score_ci(
    adata: ad.AnnData,
    cluster_key: str,
    marker_sets: dict[str, dict],
    n_boot: int = 500,
    seed: int = 42,
) -> pd.DataFrame:
    """Bootstrap CI for evidence scores by resampling cells within each cluster.

    For each cluster and each candidate cell type in marker_sets:
    1. Resample cells within the cluster with replacement
    2. Recompute mean expression overlap with positive markers
    3. Derive a simplified combined score analog
    4. Report percentile CI across bootstrap replicates

    Parameters
    ----------
    adata : AnnData with expression data
    cluster_key : obs column identifying clusters
    marker_sets : {cell_type: {positive_markers: [...], ...}} from atlas
    n_boot : number of bootstrap replicates
    seed : random seed

    Returns
    -------
    DataFrame with columns: cluster, cell_type, score_point,
        score_ci_lower, score_ci_upper, score_se
    """
    rng = np.random.default_rng(seed)
    results = []

    clusters = adata.obs[cluster_key].unique()
    gene_names = list(adata.var_names)
    gene_set = set(gene_names)

    for cl in clusters:
        mask = adata.obs[cluster_key] == cl
        n_cells = int(mask.sum())

        if n_cells < 10:
            continue

        # Get expression matrix for this cluster
        import scipy.sparse as sp

        X_cluster = adata[mask].X
        if sp.issparse(X_cluster):
            X_cluster = X_cluster.toarray()

        # Find best matching cell type
        best_ct = None
        best_score = -1.0

        for ct_name, ct_info in marker_sets.items():
            pos_markers = ct_info.get("positive_markers", [])
            pos_in_data = [g for g in pos_markers if g in gene_set]

            if not pos_in_data:
                continue

            pos_indices = [gene_names.index(g) for g in pos_in_data]

            # Point estimate: fraction of positive markers expressed
            mean_expr = X_cluster[:, pos_indices].mean(axis=0)
            pct_expressed = float(np.mean(mean_expr > 0))
            overlap = pct_expressed

            if overlap > best_score:
                best_score = overlap
                best_ct = ct_name

        if best_ct is None:
            results.append(
                {
                    "cluster": cl,
                    "cell_type": "Unknown",
                    "score_point": 0.0,
                    "score_ci_lower": 0.0,
                    "score_ci_upper": 0.0,
                    "score_se": 0.0,
                }
            )
            continue

        # Bootstrap CI for the best cell type
        pos_markers = marker_sets[best_ct].get("positive_markers", [])
        pos_in_data = [g for g in pos_markers if g in gene_set]
        pos_indices = [gene_names.index(g) for g in pos_in_data]

        boot_scores = []
        for _ in range(n_boot):
            boot_idx = rng.choice(n_cells, size=n_cells, replace=True)
            X_boot = X_cluster[boot_idx]
            mean_expr = X_boot[:, pos_indices].mean(axis=0)
            pct_expressed = float(np.mean(mean_expr > 0))
            boot_scores.append(pct_expressed)

        scores = np.array(boot_scores)
        results.append(
            {
                "cluster": cl,
                "cell_type": best_ct,
                "score_point": float(np.mean(scores)),
                "score_ci_lower": float(np.percentile(scores, 2.5)),
                "score_ci_upper": float(np.percentile(scores, 97.5)),
                "score_se": float(np.std(scores, ddof=1)),
            }
        )

    return (
        pd.DataFrame(results)
        if results
        else pd.DataFrame(
            columns=[
                "cluster",
                "cell_type",
                "score_point",
                "score_ci_lower",
                "score_ci_upper",
                "score_se",
            ]
        )
    )
