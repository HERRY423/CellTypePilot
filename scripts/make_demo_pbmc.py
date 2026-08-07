"""Generate a realistic synthetic PBMC demo dataset for CellTypePilot.

Biologically plausible design:
- Sparse background (~10% dropout-style noise), markers expressed only in
  their true cell types.
- 5 clean populations: CD4+ T, CD8+ T, naive B, NK, classical monocyte.
- 1 engineered doublet cluster (T + B signatures mixed) to demonstrate the
  critic's POSSIBLE_DOUBLET detection.
- 1 weak-evidence cluster (few markers expressed) to demonstrate
  LOW_EVIDENCE / PARTIAL_EVIDENCE flags.

Run: python scripts/make_demo_pbmc.py
Writes: demo_pbmc.h5ad (overwrites the bundled demo file)
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

RNG = np.random.default_rng(20260807)

# Cell populations and their canonical markers (subset of the blood atlas)
POPULATIONS: dict[str, dict] = {
    "CD4+ T cell": {
        "markers": ["CD3D", "CD3E", "CD2", "TRAC", "CD4", "IL7R", "MAL"],
        "n_cells": 120,
    },
    "CD8+ T cell": {
        "markers": ["CD3D", "CD3E", "CD2", "TRAC", "CD8A", "CD8B", "GZMK"],
        "n_cells": 100,
    },
    "Naive B cell": {
        "markers": ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "TCL1A"],
        "n_cells": 100,
    },
    "NK cell": {
        "markers": ["NCAM1", "NKG7", "GNLY", "KLRD1", "PRF1", "GZMB"],
        "n_cells": 90,
    },
    "Classical monocyte": {
        "markers": ["CD14", "LYZ", "S100A8", "S100A9", "FCN1", "VCAN"],
        "n_cells": 110,
    },
    "Dendritic cell": {
        "markers": ["FCER1A", "CD1C", "CLEC10A", "ITGAX", "HLA-DRA"],
        "n_cells": 60,
    },
}

N_BG_GENES = 300  # background genes not in any marker panel


def _make_population(markers: list[str], n_cells: int, all_genes: list[str]) -> np.ndarray:
    """Sparse background + strong marker expression in the right columns."""
    gene_idx = {g: i for i, g in enumerate(all_genes)}
    expr = np.zeros((n_cells, len(all_genes)), dtype=np.float32)

    # Sparse transcriptional background: ~8% of entries, low magnitude
    bg_mask = RNG.random((n_cells, len(all_genes))) < 0.08
    expr[bg_mask] = RNG.exponential(0.15, size=int(bg_mask.sum())).astype(np.float32)

    # True markers: expressed in ~90% of cells with realistic magnitudes
    for g in markers:
        col = gene_idx[g]
        on = RNG.random(n_cells) < 0.90
        expr[on, col] = RNG.exponential(1.8, size=int(on.sum())).astype(np.float32) + 0.5

    return expr


def build_demo() -> ad.AnnData:
    all_markers: list[str] = []
    for pop in POPULATIONS.values():
        for g in pop["markers"]:
            if g not in all_markers:
                all_markers.append(g)
    bg_genes = [f"GENE_{i:03d}" for i in range(N_BG_GENES)]
    all_genes = all_markers + bg_genes

    blocks: list[np.ndarray] = []
    labels: list[str] = []
    for ct, pop in POPULATIONS.items():
        blocks.append(_make_population(pop["markers"], pop["n_cells"], all_genes))
        labels.extend([ct] * pop["n_cells"])

    # Engineered doublet cluster: half T-cell program, half B-cell program
    n_doublet = 50
    doublet = np.zeros((n_doublet, len(all_genes)), dtype=np.float32)
    t_block = _make_population(POPULATIONS["CD4+ T cell"]["markers"], n_doublet, all_genes)
    b_block = _make_population(POPULATIONS["Naive B cell"]["markers"], n_doublet, all_genes)
    doublet = np.maximum(t_block, b_block)  # co-expression of both programs
    blocks.append(doublet)
    labels.extend(["Doublet"] * n_doublet)

    # Weak-evidence cluster: DC program with only 1/5 markers expressed
    n_weak = 40
    weak = np.zeros((n_weak, len(all_genes)), dtype=np.float32)
    bg_mask = RNG.random((n_weak, len(all_genes))) < 0.08
    weak[bg_mask] = RNG.exponential(0.15, size=int(bg_mask.sum())).astype(np.float32)
    gene_idx = {g: i for i, g in enumerate(all_genes)}
    on = RNG.random(n_weak) < 0.85
    weak[on, gene_idx["FCER1A"]] = RNG.exponential(1.8, size=int(on.sum())) + 0.5
    blocks.append(weak)
    labels.extend(["LowEvidence"] * n_weak)

    X_counts = np.vstack(blocks)
    X_log = np.log1p(X_counts)

    obs = pd.DataFrame(
        {
            "ground_truth": labels,
            "tissue": "blood",
        },
        index=[f"cell_{i}" for i in range(len(labels))],
    )

    adata = ad.AnnData(X=X_log.astype(np.float32), obs=obs)
    adata.var_names = all_genes
    adata.layers["counts"] = X_counts.astype(np.float32)

    # Standard embedding + clustering so inspect/annotate work out of the box
    sc.pp.pca(adata, n_comps=20, random_state=42)
    sc.pp.neighbors(adata, random_state=42)
    sc.tl.umap(adata, random_state=42)
    sc.tl.leiden(adata, resolution=0.9, random_state=42)

    return adata


def main() -> None:
    adata = build_demo()
    out = Path("demo_pbmc.h5ad")
    adata.write(out)
    print(f"Wrote {out}: {adata.n_obs} cells x {adata.n_vars} genes")
    print("Clusters:", dict(adata.obs["leiden"].value_counts()))


if __name__ == "__main__":
    main()
