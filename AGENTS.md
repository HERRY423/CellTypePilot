# CellTypePilot — Codex Agent Instructions

> Single-cell cell-type annotation with evidence, critic review, and publication-ready output.

## When to use

Activate CellTypePilot when the user wants to:
- Annotate, label, or identify cell types in pre-clustered scRNA-seq / spatial data
- Work with `.h5ad` files containing clustered single-cell data
- Understand "what are these clusters?"
- Get publication-ready cell-type annotations with evidence

Trigger phrases: "annotate my clusters", "what cell types are these?", "label my scRNA-seq",
"figure out the cell types", "run celltypepilot", "annotate my h5ad".

## Prerequisites

Before running anything, verify the environment:

```bash
celltypepilot doctor
```

If core dependencies are missing, install them:

```bash
pip install -e .
```

CellTypePilot runs on CPU, requires Python >= 3.10, scanpy, anndata, and matplotlib.
No MCP servers, no pixi, no conda environment needed for the basic path.

## Workflow — Three Stages

### Stage 1: Inspect the data

Run inspection to understand the input:

```bash
celltypepilot inspect --input <path-to-h5ad> --json
```

This reports:
- Species (auto-detected from gene naming: human ALL-CAPS vs mouse capitalized)
- Tissue (from obs metadata if available)
- Cluster keys found (leiden, louvain, etc.)
- Embedding keys found (UMAP, tSNE, etc.)
- Layer info (counts vs log-normalized)
- Cluster sizes and warnings

After inspection, confirm with the user:
1. Species and tissue context
2. Which cluster key to use
3. Which embedding to use for visualization

### Stage 2: Annotate

Run the full annotation pipeline:

```bash
celltypepilot annotate \
  --input <path> \
  --cluster-key <key> \
  --output <output_dir> \
  --species <human|mouse> \
  --tissue <tissue> \
  --embedding-key <key>
```

This produces:
- `data.annotated.h5ad` — AnnData with `ctp_cell_type`, `ctp_cl_id`, `ctp_confidence` in obs
- `evidence_table.csv` — per-cluster evidence: scores, markers, critic flags, confidence
- `figures/` — UMAP (clusters, cell types, confidence), dotplot, confidence distribution
- `report_draft.html` — comprehensive HTML report with all figures embedded
- `methodology_draft.txt` — draft methods paragraph for papers
- `manifest.json` — run provenance: versions, parameters, data hash, output hashes

After annotation, present the results:
1. How many clusters annotated, confidence distribution
2. Which clusters were flagged by the critic and why
3. Point the user to the HTML report and figures

### Stage 3: Deep review (optional)

If the user wants to investigate a specific flagged cluster:

```bash
celltypepilot critic --input <path> --cluster-key <key> --focus <cluster_id>
```

This shows the top-5 candidate cell types with scores and detailed critic evidence.

### Stage 4: Literature validation (optional)

For literature-backed marker validation (requires network access to PubMed):

```bash
# Check if literature search is available
celltypepilot doctor

# Validate specific markers against literature
celltypepilot literature --cell-type "T cells" --markers "CD3E,CD4,CD8A" --json
```

This searches PubMed for evidence supporting each marker-cell_type association.
Useful when the critic flags a cluster and you want additional validation.

## Available commands

| Command | Purpose |
|---|---|
| `celltypepilot doctor` | Check environment, dependencies, MCP status |
| `celltypepilot inspect -i <path>` | Inspect h5ad: species, tissue, clusters, embeddings |
| `celltypepilot annotate -i <path> -k <key>` | Full annotation pipeline |
| `celltypepilot critic -i <path> -k <key> -f <cluster>` | Deep-review a specific cluster |
| `celltypepilot markers -t <tissue>` | List available cell types and markers |
| `celltypepilot literature -c <type> -m <markers>` | Literature validation via PubMed |

All commands support `--json` for structured output.

## Supported tissues

Built-in Marker Knowledge Graph (MKG mkg-2026.08) covers:
- **Blood/PBMC**: T cells (CD4/CD8/naive/memory/Treg/Th1/Th17), B cells, NK cells, monocytes, DCs, platelets
- **Lung**: alveolar macrophages, AT1/AT2, ciliated, club, goblet, basal
- **Liver**: hepatocytes, Kupffer, cholangiocytes, endothelial, stellate
- **Brain**: excitatory/inhibitory neurons, astrocytes, oligodendrocytes, OPCs, microglia
- **Kidney**: proximal tubule, loop of Henle, DCT, podocytes, collecting duct
- **Gut**: enterocytes, goblet, enteroendocrine, Paneth, stem cells
- **Skin**: keratinocytes, melanocytes, Langerhans, fibroblasts
- **Heart**: cardiomyocytes, fibroblasts, endothelial, smooth muscle
- **Pancreas**: alpha/beta/delta cells, ductal, acinar
- **Skeletal muscle**: myofibers, satellite cells, FAPs
- **General**: endothelial, pericytes, fibroblasts, macrophages, mast cells, epithelial

## Critic flags reference

| Flag | Meaning |
|---|---|
| `PASS` | All checks passed |
| `LOW_EVIDENCE` | <20% markers detected — manual review needed |
| `PARTIAL_EVIDENCE` | 20-50% coverage — consider review |
| `NEG_MARKER_CONFLICT` | Negative markers expressed — likely misannotation |
| `POSSIBLE_DOUBLET` | Two lineage signatures co-expressed — sub-cluster or mark doublet |

## Confidence levels

| Level | Criteria |
|---|---|
| `high` | score >= 0.7, overlap >= 0.5, no neg conflicts, high specificity |
| `medium` | score >= 0.5, overlap >= 0.3, mild neg conflicts |
| `low` | score >= 0.3 |
| `needs_review` | score < 0.3 or major conflicts |

## Key design principles

1. **Evidence over black-box** — every annotation has marker evidence behind it
2. **Critic is the soul** — it doubts, checks, and flags; a flagged cluster is a success
3. **Cost-aware** — deterministic scoring is free; LLM reasoning is optional
4. **Artifact-centric** — all compute writes durable files (CSV, PNG, JSON)
5. **Publication-ready** — figures, evidence table, and methodology paragraph included

## Project layout

```
celltypepilot/
├── AGENTS.md                          # ← you are here (Codex instructions)
├── skills/celltypepilot/SKILL.md      # Claude Code skill (alternative entry)
├── src/celltypepilot/                 # Python package (shared backend)
│   ├── cli.py                         # CLI entry point (5 commands)
│   ├── data_adapter.py                # h5ad loading, species/tissue detection
│   ├── data/marker_atlas.json         # Built-in marker knowledge graph
│   ├── marker_scorer.py               # DE + marker overlap scoring
│   ├── critic.py                      # Annotation Critic review
│   ├── visualizer.py                  # UMAP, dotplot, confidence figures
│   ├── provenance.py                  # manifest.json generation
│   ├── reporter.py                    # HTML report + methodology text
│   └── doctor.py                      # Environment check
└── tests/                             # 18 smoke tests
```
