---
description: Annotate single-cell clusters with CellTypePilot — evidence-based cell-type labeling.
---

# Annotate

Run the full CellTypePilot annotation pipeline on pre-clustered single-cell data.

## Instructions

1. **Environment check** — Run `celltypepilot doctor` first. If it fails, stop and help the user install dependencies.

2. **Inspect data** — Run `celltypepilot inspect -i <path> --json` to auto-detect:
   - Species (human/mouse)
   - Tissue context
   - Cluster key (leiden, louvain, etc.)
   - Embedding key (UMAP, tSNE)
   - Layer info

3. **Confirm parameters** — Ask the user:
   - Confirm/correct species and tissue
   - Choose cluster key if multiple found
   - Choose embedding key for visualization

4. **Run annotation**:
   ```
   celltypepilot annotate \
     --input <path> \
     --cluster-key <key> \
     --output <output_dir> \
     --species <human|mouse> \
     --tissue <tissue> \
     --embedding-key <key>
   ```

5. **Present results**:
   - Show confidence distribution (high/medium/low/needs_review)
   - Highlight critic-flagged clusters
   - Point user to HTML report and methodology draft

## Notes

- Input must be a pre-clustered `.h5ad` file
- All outputs are written to `<output_dir>/`
- Use `/critic <cluster_id>` to deep-review any flagged cluster
- Use `/ctp-inspect <path>` to re-inspect data without annotating
