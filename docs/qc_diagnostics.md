# QC diagnostic contracts (plugin ecosystem)

QC axes are **diagnostic review signals**. They never rescue, upgrade, or rewrite
identity labels (`Unknown` stays `Unknown`).

## Axes (composable, optional)

| Axis | Typical source | Missing metadata |
|------|----------------|------------------|
| `low_rna` | `n_genes_by_counts` / `nGene` | `not_assessed_missing_metadata` |
| `high_mito` | `pct_counts_mt` / `percent_mito` | `not_assessed_missing_metadata` |
| `doublet` | obs score **or** external Scrublet/scDblFinder CSV | `not_assessed_missing_tool_output_and_metadata` |
| `ambient_rna` | obs score **or** SoupX/DecontX-style CSV | `not_assessed_missing_tool_output_and_metadata` |
| `sample_enrichment` | cluster × sample dominance | `not_assessed_missing_*_metadata` |
| `batch_sensitivity` | batch/platform column presence; performance later in benchmark | `not_assessed_missing_batch_metadata` |

## Fail-closed language

- Missing evidence → **`not_assessed_*`**
- **Never** return `clean`, `qc_clean`, `doublet_free`, `artifact_free`, etc.
- `NO_CELLS_FLAGGED` / `assessed_no_flags_among_available_axes` means “no positive flags among assessed cells”, **not** biological absence of artifacts.

## External tool tables

CSV contract (doublet / ambient):

```text
cell_id,doublet_score
cell-1,0.91
```

or with explicit flags:

```text
cell_id,predicted_doublet
cell-1,true
```

```bash
celltypepilot annotate -i data.h5ad -k leiden -o out/ \
  --doublet-table scrublet.csv --ambient-table soup.csv

celltypepilot qc-diagnostics -i data.h5ad -k leiden \
  --doublet-table scrublet.csv -o qc_out/ --json
```

Outputs: `qc_diagnostics.json`, `qc_diagnostics.csv` with
`can_rescue_identity=false` on every row.

## Benchmark path

`robustness.qc_stratified_performance` and `sample_enrichment_diagnostics` use the
same **not_assessed** policy for missing predeclared columns/values. Identity
metrics are never adjusted upward because QC looks good.

## Agent guidance

1. Surface QC flags for human review priority.
2. Do **not** change identity because doublet tool is clean or missing.
3. Treat missing ambient/doublet metadata as **not_assessed**, not “no doublets”.
