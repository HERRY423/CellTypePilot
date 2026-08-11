---
description: Inspect single-cell data — auto-detect species, tissue, clusters, embeddings.
---

# CellTypePilot Inspect

Quickly understand the structure and content of a `.h5ad` file before annotation.

> **Naming note**: registered as `/ctp-inspect` (not `/inspect`) to avoid a future collision
> with a possible built-in `/inspect` (requested upstream) and to stay distinct from the
> `claude inspect` CLI. The underlying CLI check is `celltypepilot inspect`.

## Instructions

1. **Run inspection**:
   ```
   celltypepilot inspect -i <path> --json
   ```

2. **Review the report**:
   - Species auto-detection (human/mouse gene naming patterns)
   - Tissue context (from obs metadata)
   - Cluster keys found (leiden, louvain, etc.) with cluster counts
   - Embedding keys (UMAP, tSNE, PCA)
   - Layers available (counts, log-normalized, etc.)
   - Gene ID convention (symbols vs Ensembl)
   - Any warnings (missing embeddings, unusual gene counts, etc.)

3. **Report to user**: Summarize the data characteristics and confirm parameters before annotation.

## Notes

- This is a read-only operation — no files are modified
- Use `--json` for structured output the agent can parse
- If species/tissue cannot be auto-detected, ask the user before proceeding
