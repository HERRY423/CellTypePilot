"""CellTypePilot — Seurat .rds adapter.

Bridges Seurat R objects to AnnData for the CellTypePilot pipeline.
Supports two modes:
1. rpy2 (if installed) — direct in-process conversion
2. External R script — generates a temporary .h5ad via SeuratDisk

This allows Seurat users to stay in their familiar ecosystem while
using CellTypePilot's annotation capabilities.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd


# ──────────────────────────────────────────────
# rpy2-based conversion (preferred if available)
# ──────────────────────────────────────────────

def _check_rpy2() -> bool:
    """Check if rpy2 is available."""
    try:
        import rpy2
        return True
    except ImportError:
        return False


def _convert_seurat_via_rpy2(rds_path: str | Path, chunk_size: int = 10_000) -> ad.AnnData:
    """Convert Seurat .rds to AnnData using rpy2.

    Memory-safe: extracts sparse matrix components directly from R
    without materializing a dense intermediate. For very large datasets
    (>500k cells), falls back to chunked column-wise extraction.

    Args:
        rds_path: Path to .rds file
        chunk_size: Cells per chunk for large-dataset mode (default 10k)
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.conversion import localconverter
    import scipy.sparse as sp

    rds_path = Path(rds_path).resolve()

    # Load Seurat object
    ro.r(f'''
        library(Seurat)
        seurat_obj <- readRDS("{rds_path.as_posix()}")
    ''')

    # Extract count matrix slot (sparse-preserving)
    ro.r('''
        default_assay <- DefaultAssay(seurat_obj)
        counts_mat <- GetAssayData(seurat_obj, assay = default_assay, slot = "counts")
        if (ncol(counts_mat) == 0 || nrow(counts_mat) == 0) {
            counts_mat <- GetAssayData(seurat_obj, assay = default_assay, slot = "data")
        }
    ''')

    # ── Sparse-preserving extraction ──────────────────────────
    # Extract dgCMatrix components (i, p, x) directly from R
    # without ever calling as.matrix() — avoids OOM on large datasets.
    n_genes = int(ro.r("nrow(counts_mat)")[0])
    n_cells = int(ro.r("ncol(counts_mat)")[0])
    nnz_approx = int(ro.r("""
        if (is(counts_mat, "dgCMatrix") || is(counts_mat, "CsparseMatrix")) {
            length(counts_mat@x)
        } else if (is(counts_mat, "dgTMatrix") || is(counts_mat, "TsparseMatrix")) {
            length(counts_mat@i)
        } else {
            # Estimate: if we can't determine, assume ~10% nonzero
            as.integer(nrow(counts_mat) * ncol(counts_mat) * 0.1)
        }
    """)[0])

    # Decide strategy based on dataset size
    total_elements = n_genes * n_cells
    LARGE_THRESHOLD = 500_000_000  # 500M elements (~4GB dense float64)

    if total_elements > LARGE_THRESHOLD and nnz_approx > 100_000_000:
        # Very large & dense: chunked column-wise extraction
        X = _extract_sparse_chunked(ro, sp, n_genes, n_cells, chunk_size)
    else:
        # Standard path: extract CSC components directly
        X = _extract_sparse_direct(ro, sp, n_genes, n_cells)

    # Extract gene names and cell barcodes
    gene_names = list(ro.r("rownames(counts_mat)"))
    cell_names = list(ro.r("colnames(counts_mat)"))

    # Extract metadata
    with localconverter(ro.default_converter + pandas2ri.converter):
        metadata = ro.conversion.rpy2py(ro.r("seurat_obj@meta.data"))

    # Extract embeddings
    obsm = {}
    ro.r('''
        reductions <- Reductions(seurat_obj)
    ''')
    reduction_names = list(ro.r("reductions"))

    for red_name in reduction_names:
        ro.r(f'emb <- Embeddings(seurat_obj, reduction = "{red_name}")')
        emb = np.array(ro.r("emb"))
        key = f"X_{red_name.lower()}"
        obsm[key] = emb

    # Build AnnData
    var = pd.DataFrame(index=gene_names)
    obs = pd.DataFrame(index=cell_names)

    # Merge Seurat metadata
    for col in metadata.columns:
        obs[col] = metadata[col].values

    adata = ad.AnnData(X=X, obs=obs, var=var, obsm=obsm)

    return adata


def _extract_sparse_direct(ro, sp, n_genes: int, n_cells: int) -> sp.csr_matrix:
    """Extract sparse matrix from R without dense intermediate.

    Converts R's dgCMatrix (CSC) to scipy CSR directly via component
    extraction. The R matrix is genes × cells; we need cells × genes,
    so we transpose during construction.

    Memory: O(nnz) — never materializes the full matrix.
    """
    # Ensure CSC format in R, then extract i/p/x arrays
    ro.r('''
        library(Matrix)
        .ctp_dgc <- as(counts_mat, "dgCMatrix")
    ''')

    # Extract components individually (small memory footprint)
    with localconverter(ro.default_converter + numpy2ri.converter):
        i_indices = np.array(ro.r(".ctp_dgc@i"))      # 0-based row indices
        p_indices = np.array(ro.r(".ctp_dgc@p"))      # column pointers
        x_values = np.array(ro.r(".ctp_dgc@x"))       # non-zero values

    # Build scipy CSC matrix (genes × cells in R orientation)
    csc_genes_cells = sp.csc_matrix(
        (x_values, i_indices, p_indices),
        shape=(n_genes, n_cells),
    )

    # Transpose to cells × genes → becomes CSR (efficient transpose of CSC)
    return csc_genes_cells.T.tocsr()


def _extract_sparse_chunked(ro, sp, n_genes: int, n_cells: int, chunk_size: int) -> sp.csr_matrix:
    """Chunked column-wise extraction for extremely large datasets.

    Processes cells in chunks of `chunk_size` columns, building the
    final CSR matrix from sparse blocks. Peak memory is bounded by
    chunk_size × n_genes rather than n_cells × n_genes.

    Args:
        ro: rpy2 robjects module
        sp: scipy.sparse module
        n_genes: Number of genes (rows in R matrix)
        n_cells: Number of cells (columns in R matrix)
        chunk_size: Cells per chunk

    Returns:
        scipy.sparse.csr_matrix of shape (n_cells, n_genes)
    """
    import gc

    n_chunks = (n_cells + chunk_size - 1) // chunk_size
    chunks = []

    for chunk_idx in range(n_chunks):
        col_start = chunk_idx * chunk_size + 1  # R is 1-based
        col_end = min((chunk_idx + 1) * chunk_size, n_cells)

        ro.r(f'''
            .ctp_chunk <- counts_mat[, {col_start}:{col_end}, drop = FALSE]
            .ctp_chunk_dgc <- as(.ctp_chunk, "dgCMatrix")
        ''')

        with localconverter(ro.default_converter + numpy2ri.converter):
            i_idx = np.array(ro.r(".ctp_chunk_dgc@i"))
            p_idx = np.array(ro.r(".ctp_chunk_dgc@p"))
            x_val = np.array(ro.r(".ctp_chunk_dgc@x"))

        chunk_n_cells = col_end - col_start + 1
        chunk_csc = sp.csc_matrix(
            (x_val, i_idx, p_idx),
            shape=(n_genes, chunk_n_cells),
        )
        # Transpose chunk: chunk_n_cells × n_genes (CSR)
        chunks.append(chunk_csc.T.tocsr())

        # Free R memory
        ro.r("rm(.ctp_chunk, .ctp_chunk_dgc); gc()")
        gc.collect()

    # Stack all chunks vertically (cells axis)
    return sp.vstack(chunks, format="csr")


# ──────────────────────────────────────────────
# External R script fallback
# ──────────────────────────────────────────────

R_CONVERSION_SCRIPT = '''
#!/usr/bin/env Rscript
# CellTypePilot Seurat -> h5ad conversion script
# Usage: Rscript convert_seurat.rds <input.rds> <output.h5ad>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript convert_seurat.rds <input.rds> <output.h5ad>")
}

input_rds <- args[1]
output_h5ad <- args[2]

# Check required packages
if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat R package is required. Install with: install.packages('Seurat')")
}

if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
    stop("SeuratDisk R package is required. Install with: install.packages('SeuratDisk')")
}

library(Seurat)
library(SeuratDisk)

# Load Seurat object
cat("Loading Seurat object from:", input_rds, "\\n")
seurat_obj <- readRDS(input_rds)

# Save as temporary .h5Seurat
temp_h5seurat <- tempfile(fileext = ".h5Seurat")
SaveH5Seurat(seurat_obj, filename = temp_h5seurat)

# Convert to .h5ad
cat("Converting to h5ad...\\n")
Convert(temp_h5seurat, dest = output_h5ad, overwrite = TRUE)

cat("Conversion complete:", output_h5ad, "\\n")

# Clean up
unlink(temp_h5seurat)
'''


def _convert_seurat_via_rscript(rds_path: str | Path) -> ad.AnnData:
    """Convert Seurat .rds to AnnData using external R script.

    Requires: R with Seurat and SeuratDisk packages installed.
    """
    rds_path = Path(rds_path).resolve()

    # Check if R is available
    r_path = shutil.which("Rscript")
    if r_path is None:
        raise RuntimeError(
            "Rscript not found in PATH. Either install rpy2 (pip install rpy2) "
            "or install R with Seurat and SeuratDisk packages."
        )

    # Create temporary directory for conversion
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Write R script
        script_path = tmpdir / "convert_seurat.rds"
        script_path.write_text(R_CONVERSION_SCRIPT, encoding="utf-8")

        # Output path
        output_h5ad = tmpdir / "converted.h5ad"

        # Run R script
        result = subprocess.run(
            [r_path, str(script_path), str(rds_path), str(output_h5ad)],
            capture_output=True,
            text=True,
            timeout=300,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"R conversion failed:\n{result.stderr}\n\n"
                "Make sure Seurat and SeuratDisk R packages are installed:\n"
                "  install.packages('Seurat')\n"
                "  install.packages('SeuratDisk')"
            )

        if not output_h5ad.exists():
            raise FileNotFoundError(f"Conversion did not produce output: {output_h5ad}")

        # Load the converted h5ad
        adata = ad.read_h5ad(output_h5ad)

    return adata


# ──────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────

def load_seurat_rds(rds_path: str | Path) -> ad.AnnData:
    """Load a Seurat .rds file and convert to AnnData.

    Tries rpy2 first (faster, in-process), falls back to external R script.

    Args:
        rds_path: Path to the .rds file

    Returns:
        AnnData object ready for CellTypePilot pipeline

    Raises:
        RuntimeError: If neither rpy2 nor R is available
        FileNotFoundError: If the .rds file doesn't exist
    """
    rds_path = Path(rds_path)

    if not rds_path.exists():
        raise FileNotFoundError(f"Seurat .rds file not found: {rds_path}")

    if not rds_path.suffix.lower() == ".rds":
        raise ValueError(f"Expected .rds file, got: {rds_path.suffix}")

    # Try rpy2 first
    if _check_rpy2():
        try:
            return _convert_seurat_via_rpy2(rds_path)
        except Exception as e:
            # Fall back to R script
            pass

    # Fall back to external R script
    return _convert_seurat_via_rscript(rds_path)


def check_seurat_support() -> dict:
    """Check Seurat .rds support availability.

    Returns:
        Dict with support status and details
    """
    support = {
        "rpy2_available": False,
        "r_available": False,
        "seurat_rds_supported": False,
        "details": {},
    }

    # Check rpy2
    try:
        import rpy2
        support["rpy2_available"] = True
        support["details"]["rpy2_version"] = rpy2.__version__
    except ImportError:
        pass

    # Check R
    r_path = shutil.which("Rscript")
    if r_path:
        support["r_available"] = True
        support["details"]["r_path"] = r_path

        # Check Seurat/SeuratDisk
        try:
            result = subprocess.run(
                [r_path, "-e", "cat(packageVersion('Seurat'))"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            if result.returncode == 0:
                support["details"]["seurat_version"] = result.stdout.strip()
        except Exception:
            pass

    # Overall support
    support["seurat_rds_supported"] = support["rpy2_available"] or support["r_available"]

    return support


def seurat_to_h5ad(rds_path: str | Path, output_path: str | Path) -> Path:
    """Convert Seurat .rds to .h5ad file.

    Useful for users who want to convert their data once and reuse.

    Args:
        rds_path: Path to input .rds file
        output_path: Path for output .h5ad file

    Returns:
        Path to the created .h5ad file
    """
    adata = load_seurat_rds(rds_path)
    output_path = Path(output_path)
    adata.write_h5ad(output_path)
    return output_path
