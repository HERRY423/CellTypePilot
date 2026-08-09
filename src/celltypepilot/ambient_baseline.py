"""Ambient RNA Baseline Decontamination Module.

Computes dataset-wide background gene expression fractions and identifies
ambient RNA contaminants (e.g., Gfap in brain, Hba-a1 in blood) to prevent
false-positive negative marker conflict flags and doublet calls.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
from scipy import sparse

# High-frequency ambient contaminants by tissue/context
KNOWN_HIGH_AMBIENT_GENES: dict[str, set[str]] = {
    "brain": {"Gfap", "Slc1a3", "Plp1", "Mbp", "Mog", "Aqp4", "S100b", "GFAP", "PLP1", "MBP"},
    "blood": {"HBA1", "HBA2", "HBB", "HBD", "Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt"},
    "pbmc": {"HBA1", "HBA2", "HBB", "HBD", "Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt"},
    "lung": {"SFTPB", "SFTPC", "Sftpb", "Sftpc"},
    "liver": {"ALB", "Alb", "HP", "Hp"},
}


def compute_ambient_baseline(
    adata: ad.AnnData,
    layer: str | None = None,
    threshold: float = 0.0,
) -> dict[str, float]:
    """Compute global background expression fraction for all genes across the dataset.

    Args:
        adata: AnnData object
        layer: Matrix layer name or None for X
        threshold: Expression cutoff for calling a gene expressed in a cell

    Returns:
        dict mapping gene symbol -> dataset-wide fraction of expressing cells
    """
    matrix = adata.layers[layer] if layer is not None else adata.X
    n_cells = adata.n_obs

    if n_cells == 0:
        return {}

    if sparse.issparse(matrix):
        # Count non-zero entries per column
        expr_counts = np.asarray((matrix > threshold).sum(axis=0)).ravel()
    else:
        expr_counts = np.asarray((matrix > threshold).sum(axis=0)).ravel()

    global_fractions = expr_counts / float(n_cells)
    return {str(gene): float(frac) for gene, frac in zip(adata.var_names, global_fractions)}


def is_ambient_contamination(
    gene: str,
    cluster_pct: float,
    global_pct: float,
    tissue: str | None = None,
    ambient_fold_threshold: float = 2.5,
    min_cluster_pct: float = 0.15,
) -> bool:
    """Determine if a gene's expression in a cluster is ambient RNA contamination.

    A gene expression signal is considered ambient contamination if:
    1. The cluster expression fraction is below min_cluster_pct (e.g. < 15%), OR
    2. The cluster expression fraction does NOT significantly exceed (>= fold_threshold)
       the dataset-wide global expression fraction (indicating widespread background).

    Known high-ambient genes (like Gfap in brain) receive an additional safety margin.

    Args:
        gene: Gene symbol
        cluster_pct: Fraction of cells in the cluster expressing the gene
        global_pct: Dataset-wide fraction of cells expressing the gene
        tissue: Tissue type if known
        ambient_fold_threshold: Ratio of cluster_pct / global_pct required to be non-ambient
        min_cluster_pct: Expression cutoff below which signal is considered noise

    Returns:
        True if the signal is consistent with ambient RNA contamination.
    """
    if cluster_pct < min_cluster_pct:
        return True

    # If dataset-wide global expression is very low, cluster_pct is likely specific
    if global_pct < 0.05:
        return False

    fold_over_global = cluster_pct / global_pct

    # Adjust threshold for known contaminant genes in specific tissue
    tissue_key = (tissue or "").lower()
    known_contaminants = KNOWN_HIGH_AMBIENT_GENES.get(tissue_key, set())
    effective_threshold = ambient_fold_threshold * 1.5 if gene in known_contaminants else ambient_fold_threshold

    return fold_over_global < effective_threshold
