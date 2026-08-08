# CellTypePilot — Output Reference

> Detailed schemas for all output files.

## Output Directory Layout

```
<output_dir>/
├── data.annotated.h5ad          # AnnData with annotation results
├── evidence_table.csv           # Per-cluster evidence and critic results
├── state_results.csv            # Independent state evidence (when enabled)
├── context_pack.normalized.json # Governed context (when supplied)
├── manifest.json                # Run provenance manifest
├── report_draft.html            # Comprehensive HTML report
├── methodology_draft.txt        # Draft methods paragraph
└── figures/
    ├── umap_cluster.png         # UMAP by cluster
    ├── umap_celltype.png        # UMAP by cell type
    ├── umap_confidence.png      # UMAP by confidence level
    ├── marker_dotplot.png       # Marker gene dotplot
    └── confidence_distribution.png  # Confidence bar chart
```

## 1. data.annotated.h5ad

Standard AnnData with identity and state columns in `.obs`:

| Column | Type | Description |
|---|---|---|
| `ctp_cell_type` | str | Accepted identity, or `Unknown` after abstention |
| `ctp_cl_id` | str | Cell Ontology ID (e.g., CL:0000084) |
| `ctp_confidence` | str | Critic-calibrated confidence: high/medium/low/needs_review |
| `ctp_candidate_cell_type` | str | Best identity candidate, retained across abstention |
| `ctp_decision` | str | Identity decision: accepted or abstain |
| `ctp_abstain_reason` | str | Machine-readable fail-closed reason(s) |
| `ctp_cell_state_candidate` | str | Best independent state candidate |
| `ctp_state_decision` | str | supported, hypothesis, or abstain |
| `ctp_cell_state` | str | Supported state, otherwise `Unknown` |
| `ctp_state_score` | float | Heuristic state evidence score |
| `ctp_state_confidence` | str | Heuristic state confidence |
| `ctp_state_evidence` | str | Compact state evidence summary |
| `ctp_display_label` | str | Identity plus supported state for display only |

The original file is never modified. The annotated copy is written to the output directory.

## 2. evidence_table.csv

One row per cluster. Columns:

| Column | Type | Description |
|---|---|---|
| `cluster` | str | Cluster ID (from the cluster key) |
| `cell_type` | str | Assigned cell type name |
| `cl_id` | str | Cell Ontology ID |
| `combined_score` | float | Combined marker score [0, 1] |
| `confidence` | str | Pre-critic confidence level |
| `pct_overlap` | float | Fraction of markers that are DE |
| `mean_log2fc` | float | Mean log2 fold change of DE markers |
| `specificity` | float | Mean marker specificity to this cluster |
| `neg_conflict` | float | Fraction of negative markers expressed |
| `critic_flags` | str | Semicolon-separated critic flags (or "PASS") |
| `critic_evidence` | str | Detailed evidence text from critic |
| `critic_confidence` | str | Post-critic confidence level |
| `critic_notes` | str | Human-readable notes from critic |

Identity coverage includes every expected marker. Missing genes and present-but-silent genes are
reported separately. Context-derived support and its review status are also retained so that
context-only evidence cannot be mistaken for bundled atlas evidence.

## 3. state_results.csv

One row per cluster when the State Lens is enabled. It records the best state candidate, state
decision, score, confidence, supporting markers, negative conflicts, and separate missing/silent
marker counts. State attachment is identity-invariant and cannot alter an identity abstention.

## 4. context_pack.normalized.json

Written only when context or custom markers are provided. It contains the normalized
`celltypepilot.context.v1` pack and canonical/source hashes. Free text is preserved for audit but
does not contribute marker evidence.

## 5. manifest.json

```json
{
  "celltypepilot_version": "0.3.0",
  "mkg_version": "mkg-2026.08",
  "timestamp": "2026-08-06T12:00:00+00:00",
  "input": {
    "path": "/path/to/data.h5ad",
    "sha256": "abc123..."
  },
  "parameters": {
    "cluster_key": "leiden",
    "species": "human",
    "tissue": "blood",
    "embedding_key": "X_umap",
    "layer": null,
    "context_schema_version": "celltypepilot.context.v1",
    "context_sha256": null,
    "state_lens_enabled": true,
    "state_contract": "identity_invariant_independent_axis_v1"
  },
  "outputs": {
    "evidence_table.csv": {
      "sha256": "def456...",
      "size_bytes": 1234
    },
    "figures/umap_cluster.png": {
      "sha256": "ghi789...",
      "size_bytes": 56789
    }
  }
}
```

## 6. report_draft.html

Self-contained HTML report with embedded CSS. Sections:
1. **Overview**: Stats grid showing cluster count, confidence distribution, critic pass/fail
2. **Annotation Results**: Full table with cluster, cell type, CL ID, score, confidence, critic status
3. **Figures**: Grid of all generated figures with captions
4. **Critic Details**: Expanded cards for each flagged cluster showing evidence and notes

## 7. methodology_draft.txt

A plain-text paragraph suitable for adaptation into a paper's Methods section. Example:

> Cell type annotation was performed using CellTypePilot (v0.3.0), an evidence-driven
> annotation pipeline with built-in critic review. Marker gene evidence was sourced from
> the CellTypePilot Marker Knowledge Graph (MKG mkg-2026.08), a curated atlas integrating
> PanglaoDB, CellMarker, and Cell Ontology resources. For each of the N clusters identified
> by leiden clustering, marker gene overlap, expression specificity, fold-change magnitude,
> and negative marker conflict were scored to generate candidate annotations with confidence
> levels. A rules-based Annotation Critic module reviewed each assignment for evidence
> sufficiency, negative marker conflicts, potential doublet signatures, and ontology
> consistency. Of N clusters, X were assigned high confidence, Y medium confidence, and
> Z were flagged for manual review.

## JSON Output Mode

All commands support `--json` for structured output. The annotate command's JSON output:

```json
{
  "annotations": [
    {
      "cluster": "0",
      "cell_type": "CD4+ T cell",
      "cl_id": "CL:0000624",
      "combined_score": 0.82,
      "confidence": "high",
      "critic_flags": "PASS",
      "critic_confidence": "high"
    }
  ],
  "critic_summary": {
    "total_clusters": 15,
    "pass": 12,
    "flagged": 3,
    "confidence_distribution": {
      "high": 8,
      "medium": 4,
      "low": 2,
      "needs_review": 1
    }
  },
  "manifest": { ... }
}
```
