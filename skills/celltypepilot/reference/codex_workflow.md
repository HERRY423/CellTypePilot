# CellTypePilot — Codex Workflow Reference

> Detailed workflow guide for Codex CLI users.

## Quick start (one-liner)

```bash
# Install
pip install celltypepilot

# Run everything in one step
celltypepilot annotate --input data.h5ad --cluster-key leiden --tissue blood --output ./results
```

## Step-by-step workflow

### 1. Environment check

```bash
celltypepilot doctor
```

Output tells you:
- Python version (need >= 3.10)
- Core dependencies (anndata, scanpy, numpy, pandas, scipy, matplotlib, seaborn, typer, rich)
- Optional dependencies (python-docx, cupy, scvi-tools, decoupler)
- Capability tier: "full" = all core deps met, "degraded" = missing deps

### 2. Data inspection

```bash
celltypepilot inspect --input data.h5ad --json
```

JSON output includes:
```json
{
  "path": "data.h5ad",
  "n_obs": 10000,
  "n_vars": 500,
  "species": "human",
  "tissue": "blood",
  "cluster_keys": ["leiden", "louvain"],
  "embedding_keys": ["X_umap", "X_pca"],
  "cluster_sizes": {"0": 2000, "1": 1500, ...},
  "gene_id_convention": "human_symbols",
  "warnings": []
}
```

Use this to confirm parameters before annotation.

### 3. Annotation

```bash
celltypepilot annotate \
  --input data.h5ad \
  --cluster-key leiden \
  --output ./results \
  --species human \
  --tissue blood \
  --embedding-key X_umap
```

What happens internally:
1. Load h5ad, compute SHA-256 hash
2. Load Marker Knowledge Graph for species/tissue
3. For each cluster: compute DE genes (Wilcoxon), score against all cell type markers
4. Assign top-scoring cell type per cluster
5. Run Annotation Critic on each assignment
6. Generate figures (UMAP, dotplot, confidence)
7. Write HTML report, evidence table, manifest, methodology draft

### 4. Review results

Read the evidence table:
```bash
cat results/evidence_table.csv
```

Key columns:
- `cluster` — cluster ID
- `cell_type` — assigned cell type
- `combined_score` — marker score [0, 1]
- `critic_confidence` — high/medium/low/needs_review
- `critic_flags` — PASS or specific flags

Open the HTML report:
```bash
# macOS
open results/report_draft.html
# Linux
xdg-open results/report_draft.html
# Windows
start results\report_draft.html
```

### 5. Deep review (optional)

For clusters flagged by the critic:

```bash
celltypepilot critic --input data.h5ad --cluster-key leiden --focus 7
```

Output shows:
- Top-5 candidate cell types with scores
- Which markers support each candidate
- Negative marker conflicts
- Doublet signal assessment

### 6. Use results in Python

```python
import scanpy as sc

adata = sc.read_h5ad("results/data.annotated.h5ad")

# New columns in obs:
# - ctp_cell_type: assigned cell type
# - ctp_cl_id: Cell Ontology ID
# - ctp_confidence: critic confidence level

sc.pl.umap(adata, color=["ctp_cell_type", "ctp_confidence"])
```

## Codex quick-start script

For a one-step workflow, use the bundled script:

```bash
# Linux/macOS
bash scripts/codex_annotate.sh data.h5ad leiden ./output blood

# Windows PowerShell
.\scripts\codex_annotate.ps1 -Input data.h5ad -ClusterKey leiden -OutputDir .\output -Tissue blood
```

The script automatically:
1. Runs `celltypepilot doctor` to check environment
2. Runs `celltypepilot inspect --json` to auto-detect species, tissue, embedding
3. Runs `celltypepilot annotate` with detected parameters
4. Prints a summary of output files

## Marker atlas reference

List available tissues:
```bash
celltypepilot markers
```

List cell types and markers for a tissue:
```bash
celltypepilot markers --tissue blood
```

JSON output:
```bash
celltypepilot markers --tissue blood --json
```

## Troubleshooting

| Problem | Solution |
|---|---|
| `celltypepilot: command not found` | Run `pip install celltypepilot`, or `pip install -e .` from a development checkout |
| `No cluster key found` | Check `celltypepilot inspect` output for available keys |
| `No annotations generated` | Marker gene overlap too low — check species/tissue match |
| `GBK/Unicode encoding error` | Fixed in v0.1.0; update with `pip install --upgrade celltypepilot` |
| `leidenalg not found` | Use louvain clusters or assign cluster labels manually |
