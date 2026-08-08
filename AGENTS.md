# CellTypePilot — Codex Agent Instructions

> Single-cell annotation plugin with evidence, conservative abstention, and reviewable draft output.

## When to use

Activate CellTypePilot when the user wants to:
- Annotate, label, or identify cell types in pre-clustered scRNA-seq / spatial data
- Work with `.h5ad` files containing clustered single-cell data
- Convert Seurat `.rds` files to `.h5ad` for analysis
- Understand "what are these clusters?"
- Get auditable draft cell-type annotations with evidence for human review
- Review annotations interactively in a web panel

Trigger phrases: "annotate my clusters", "what cell types are these?", "label my scRNA-seq",
"figure out the cell types", "run celltypepilot", "annotate my h5ad", "convert my rds",
"review annotations in browser".

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

Optional extras:
- `pip install -e ".[web]"` — Web Inspector (Flask-based interactive review panel)
- `pip install -e ".[seurat]"` — Seurat .rds conversion (requires rpy2 or R)
- `pip install -e ".[all]"` — All optional features

## Workflow — Four Stages

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

If the user provides a Seurat `.rds` file, convert first:

```bash
celltypepilot convert-rds --input <path.rds> --output <path.h5ad>
```

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
- `data.annotated.h5ad` — AnnData with final label, candidate, decision, abstain reason, CL ID, and confidence in obs
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
| `celltypepilot doctor` | Check environment, dependencies, MCP status, license |
| `celltypepilot inspect -i <path>` | Inspect h5ad: species, tissue, clusters, embeddings |
| `celltypepilot annotate -i <path> -k <key>` | Full annotation pipeline |
| `celltypepilot benchmark -i <path> --truth-key <key> --study-key <key> --donor-key <key>` | Lock/evaluate independent holdouts |
| `celltypepilot benchmark-run ...` | Execute fold-isolated CellTypePilot/CellTypist and configured SingleR/Azimuth/popV adapters |
| `celltypepilot calibrate ...` | Fit a downgrade-only abstention policy on a separate calibration dataset |
| `celltypepilot critic -i <path> -k <key> -f <cluster>` | Deep-review a specific cluster |
| `celltypepilot markers -t <tissue>` | List available cell types and markers |
| `celltypepilot literature -c <type> -m <markers>` | Literature validation via PubMed |
| `celltypepilot inspect-web -o <dir>` | Launch Web Inspector (interactive review panel) |
| `celltypepilot convert-rds -i <path.rds>` | Convert Seurat .rds → .h5ad |
| `celltypepilot license status` | Check license tier and features |

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

Premium atlas (requires academic/commercial license):
- **Tumor microenvironment**: TAMs, CAFs, Tregs, exhausted T cells, MDSCs, malignant cells
- **Developing brain**: radial glia, intermediate progenitors, migrating neurons
- **Inflamed tissue**: activated fibroblasts, M1/M2 macrophages
- **Immune activation**: activated CD4/CD8, plasmablasts, activated DCs

## Critic flags reference

| Flag | Meaning |
|---|---|
| `PASS` | All checks passed |
| `LOW_EVIDENCE` | <20% of all expected markers expressed — abstain |
| `PARTIAL_EVIDENCE` | 20-50% of all expected markers expressed — abstain pending review |
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
5. **Review-ready drafts** — figures, evidence table, and draft methodology are included; human sign-off remains required
6. **Plugin, not Agent** — keep orchestration deterministic and task-bounded; do not add autonomous planning, self-directed analysis, or multi-agent behavior
6. **Fail closed** — insufficient/conflicting evidence writes `Unknown` and preserves a separate candidate

## Project layout

```
celltypepilot/
├── .codex-plugin/
│   └── plugin.json              # ← Codex plugin manifest (you are here)
├── .claude-plugin/
│   └── plugin.json              # Claude Code plugin manifest
├── AGENTS.md                    # Codex agent instructions (this file)
├── skills/celltypepilot/
│   ├── SKILL.md                 # Shared skill instructions
│   ├── agents/openai.yaml       # Codex agent interface config
│   └── reference/               # Reference docs
├── commands/                    # Claude Code slash commands (/annotate, /critic, etc.)
├── hooks/                       # Claude Code lifecycle hooks
├── rules/                       # Claude Code behavior rules
├── .mcp.json                    # MCP servers (PubMed, bioRxiv)
├── src/celltypepilot/           # Python package (shared backend)
│   ├── cli.py                   # CLI entry point (thin layer, delegates to orchestrator)
│   ├── orchestrator.py          # Pipeline business logic (shared by CLI & Web Inspector)
│   ├── data_adapter.py          # h5ad loading, robust species/tissue detection
│   ├── seurat_adapter.py        # Seurat .rds → AnnData conversion
│   ├── constants.py             # Thresholds, species/tissue constants
│   ├── data/
│   │   ├── marker_atlas.json    # Built-in marker knowledge graph (80+ types)
│   │   └── premium_atlas.json   # Premium atlas (tumor, brain, immune)
│   ├── templates/               # Jinja2 templates (HTML report, web dashboard)
│   ├── marker_scorer.py         # DE + marker overlap scoring
│   ├── reference_scorer.py      # Reference embedding scoring (4 backends)
│   ├── ensemble_scorer.py       # Adaptive ensemble fusion
│   ├── critic.py                # Annotation Critic review
│   ├── visualizer.py            # UMAP, dotplot, confidence figures
│   ├── web_inspector.py         # Flask web review panel
│   ├── literature.py            # PubMed literature validation
│   ├── license_manager.py       # Tiered license system
│   ├── provenance.py            # manifest.json generation
│   ├── reporter.py              # HTML report + methodology text
│   └── doctor.py                # Environment check
└── tests/                       # Unit, contract, and scientific-boundary tests
```
