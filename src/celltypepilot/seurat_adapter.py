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


def _convert_seurat_via_rpy2(rds_path: str | Path) -> ad.AnnData:
    """Convert Seurat .rds to AnnData using rpy2.

    Requires: rpy2, Seurat (R package), SeuratDisk (R package)
    """
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    import scipy.sparse as sp

    rds_path = Path(rds_path).resolve()

    # Load Seurat object
    ro.r(f'''
        library(Seurat)
        seurat_obj <- readRDS("{rds_path.as_posix()}")
    ''')

    # Extract count matrix
    ro.r('''
        # Get the default assay
        default_assay <- DefaultAssay(seurat_obj)
        counts_mat <- GetAssayData(seurat_obj, assay = default_assay, slot = "counts")
        if (ncol(counts_mat) == 0 || nrow(counts_mat) == 0) {
            counts_mat <- GetAssayData(seurat_obj, assay = default_assay, slot = "data")
        }
    ''')

    # Convert count matrix to scipy sparse
    counts_r = ro.r("counts_mat")
    try:
        # Try as sparse matrix
        from rpy2.robjects import numpy2ri
        with localconverter(ro.default_converter + numpy2ri.converter):
            counts_dense = np.array(ro.r("as.matrix(counts_mat)"))
        X = sp.csr_matrix(counts_dense.T)  # Cells x Genes
    except Exception:
        # Fallback: extract as dense
        counts_dense = np.array(ro.r("as.matrix(counts_mat)"))
        X = sp.csr_matrix(counts_dense.T)

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
