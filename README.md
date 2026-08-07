# CellTypePilot

> **Single-cell annotation trust layer that parasitizes on your existing coding agent.**
> For individual researchers and small labs — no standalone app, no heavy infrastructure.

**CellTypePilot** is a **plugin** for Claude Code / OpenAI Codex. It turns pre-clustered
single-cell data into trusted, publication-ready cell-type annotations — with marker evidence,
confidence levels, an independent critic review, and a draft methodology paragraph for your paper.

It is **not** another cell-type annotation algorithm. It is a **trust layer** that lives inside
the coding agent you already use, turning a conversation into a reviewed annotation workflow.

## Why it exists

For individual researchers without dedicated bioinformatics support, cell-type annotation is
painful — not because the algorithms don't exist, but because:

| Pain point | What CellTypePilot does |
|---|---|
| **No one reviews your annotations** | Annotation Critic independently checks evidence sufficiency, negative marker conflicts, doublet signals, and ontology consistency — *before* you trust a label |
| **Config barrier is too high** (MCP, pixi, conda...) | `git clone` + `pip install -e .` + run. Zero MCP required for the basic path. `doctor` tells you what you have *before* anything fails |
| **Workflow fragmentation** (scripts here, tool there) | Runs inside your Claude Code / Codex session — no context switch to a separate app |
| **Can't explain *why* a cluster got its label** | Every annotation ships with: supporting markers, expression stats, critic flags, and a draft methods paragraph |
| **Rare / transitional states forced into a label** | Critic flags doublet signals and low-confidence calls instead of silently assigning a wrong label |
| **Token cost spirals** | Deterministic marker scoring is free (Tier 0); expensive multi-model consensus only triggers for genuinely ambiguous clusters (Tier 1) |
| **Results unreproducible** | `manifest.json` records knowledge graph version, parameters, data hash, and output hashes for every run |

## What you get

```
output/
├── data.annotated.h5ad          # ctp_cell_type, ctp_cl_id, ctp_confidence in obs
├── evidence_table.csv           # Per-cluster: scores, markers, critic flags, confidence
├── report_draft.html            # Self-contained HTML report with all figures embedded
├── methodology_draft.txt        # "We annotated N clusters using CellTypePilot v0.1.0..."
├── manifest.json                # Provenance: versions, params, data hash, output hashes
└── figures/
    ├── umap_cluster.png         # UMAP by cluster (colorblind-friendly Wong palette)
    ├── umap_celltype.png        # UMAP by annotated cell type
    ├── umap_confidence.png      # UMAP by critic confidence level
    ├── marker_dotplot.png       # Cell type x marker gene expression
    └── confidence_distribution.png
```

## Quick start

```bash
# 1. Install
pip install -e .

# 2. Environment check — tells you what works and what's missing
celltypepilot doctor

# 3. Inspect your data — auto-detects species, tissue, clusters, embeddings
celltypepilot inspect --input data.h5ad

# 4. Annotate — full pipeline in one command
celltypepilot annotate --input data.h5ad --cluster-key leiden --tissue blood

# 5. Deep-review a flagged cluster
celltypepilot critic --input data.h5ad --cluster-key leiden --focus cluster_7

# 6. Validate markers against literature (optional, needs network)
celltypepilot literature --cell-type "T cells" --markers "CD3E,CD4,CD8A"
```

Or just tell your agent: *"annotate my clusters in data.h5ad"* — the plugin instructions
guide Claude Code / Codex through the full workflow automatically.

## Plugin architecture

CellTypePilot **parasitizes** on the coding agent you already use. It does not ask you to
switch to a new app, learn a new UI, or configure a new environment.

```
┌──────────────────────────────────────────────────────────────────┐
│  Host agent:  Claude Code (SKILL.md)  │  Codex CLI (AGENTS.md)  │
│                 └──────────┬──────────┘                          │
│                            ↓                                     │
│  CellTypePilot plugin (shared Python backend)                    │
│  ┌─────────────────────────────────────────────────────────────┐ │
│  │  doctor          Environment gate — what can run, what can't│ │
│  │  inspect         Data intelligence — species, tissue, keys  │ │
│  │  annotate        Full pipeline:                             │ │
│  │    ├─ Data Adapter         Load .h5ad, detect, validate     │ │
│  │    ├─ Marker Knowledge Graph  80+ cell types, 11 tissues    │ │
│  │    ├─ Marker Scorer        Wilcoxon DE + 5-dim scoring      │ │
│  │    ├─ Annotation Critic    Independent evidence review      │ │
│  │    ├─ Visualizer           UMAP, dotplot, confidence (Wong) │ │
│  │    ├─ Reporter             HTML report + methods paragraph  │ │
│  │    ├─ Literature           PubMed validation (optional MCP) │ │
│  │    └─ Provenance           manifest.json versioning         │ │
│  │  critic          Deep-review a specific cluster             │ │
│  └─────────────────────────────────────────────────────────────┘ │
└──────────────────────────────────────────────────────────────────┘
```

### Dual-platform, single backend

| Platform | Entry point | How it works |
|---|---|---|
| **Claude Code** | `skills/celltypepilot/SKILL.md` | YAML frontmatter + 4-stage orchestration |
| **Codex CLI** | `AGENTS.md` (repo root) | Pure markdown instructions, auto-discovered |

Both platforms call the same `celltypepilot` CLI. The only difference is how the agent
discovers and orchestrates the workflow.

### Installation

```bash
# Claude Code
git clone https://github.com/HERRY423/CellTypePilot ~/.claude/skills/celltypepilot
cd ~/.claude/skills/celltypepilot && pip install -e .

# Codex CLI
git clone https://github.com/HERRY423/CellTypePilot ~/celltypepilot
cd ~/celltypepilot && pip install -e .
# Codex reads AGENTS.md automatically.

# Or just use it as a CLI tool
git clone https://github.com/HERRY423/CellTypePilot && cd CellTypePilot
pip install -e .
```

## Built-in Marker Knowledge Graph

The Marker Knowledge Graph (`mkg-2026.08`) covers 80+ cell types across 11 tissues,
with positive/negative markers, Cell Ontology IDs, and synonyms. Human and mouse are
both supported with automatic gene symbol conversion.

| Tissue | Cell Types |
|---|---|
| Blood/PBMC | T cells (CD4/CD8/naive/memory/Treg/Th1/Th17), B cells, NK cells, monocytes, DCs, platelets |
| Lung | Alveolar macrophages, AT1/AT2, ciliated, club, goblet, basal cells |
| Liver | Hepatocytes, Kupffer cells, cholangiocytes, endothelial, stellate cells |
| Brain | Excitatory/inhibitory neurons, astrocytes, oligodendrocytes, OPCs, microglia |
| Kidney | Proximal tubule, loop of Henle, DCT, podocytes, collecting duct |
| Gut | Enterocytes, goblet, enteroendocrine, Paneth, stem cells |
| Skin | Keratinocytes, melanocytes, Langerhans cells, fibroblasts |
| Heart | Cardiomyocytes, fibroblasts, endothelial, smooth muscle |
| Pancreas | Alpha/beta/delta cells, ductal, acinar cells |
| Skeletal muscle | Myofibers, satellite cells, FAPs |
| General | Endothelial, pericytes, fibroblasts, macrophages, mast cells, epithelial |

## Annotation Critic — the soul of the plugin

The Critic doesn't just score — it *doubts*. Every annotation is independently reviewed:

| Check | What it catches |
|---|---|
| Evidence sufficiency | < 20% marker coverage → `LOW_EVIDENCE`; 20–50% → `PARTIAL_EVIDENCE` |
| Negative marker conflict | Negative markers expressed in > 15% of cells → `NEG_MARKER_CONFLICT` |
| Doublet signal | Two mutually exclusive lineage signatures co-expressed → `POSSIBLE_DOUBLET` |
| Ontology consistency | Is the label a valid Cell Ontology term in the right tissue context? |

Confidence levels: **high** / **medium** / **low** / **needs_review**. The Critic can only
downgrade, never upgrade. A flagged cluster is a *success* — it means the system caught
something worth your attention.

## Roadmap

- [x] **Phase 1 (MVP)** — h5ad adapter, marker knowledge graph, Wilcoxon DE scoring, Annotation Critic, doctor, figures, JSON output, HTML report, methodology draft, manifest provenance, literature validation (PubMed)
- [ ] **Phase 2** — Tiered consensus orchestrator (Tier 0 → Tier 1 adaptive upgrade), doublet detection enhancement, ontology validation enhancement
- [ ] **Phase 3** — docx/pptx submission package, Seurat/R support, reference mapping methods
- [ ] **Phase 4** — Extended atlas subscription, team sharing, rare cell type mining

## License

MIT
