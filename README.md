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

# 5. Review interactively in browser (optional)
celltypepilot inspect-web --output output/

# 6. Deep-review a flagged cluster
celltypepilot critic --input data.h5ad --cluster-key leiden --focus cluster_7

# 7. Validate markers against literature (optional, needs network)
celltypepilot literature --cell-type "T cells" --markers "CD3E,CD4,CD8A"

# 8. Convert Seurat .rds → .h5ad (optional, needs R)
celltypepilot convert-rds --input data.rds --output data.h5ad
```

Or just tell your agent: *"annotate my clusters in data.h5ad"* — the plugin instructions
guide Claude Code / Codex through the full workflow automatically.

## Plugin architecture

CellTypePilot **parasitizes** on the coding agent you already use. It does not ask you to
switch to a new app, learn a new UI, or configure a new environment.

> **Why "plugin" and not "skill"?** "Plugin" is the product concept — a self-contained
> intelligence layer that attaches to a host agent. Each platform has its own native
> plugin format: Claude Code uses `.claude-plugin/plugin.json`, Codex uses
> `.codex-plugin/plugin.json`. Both share the same `skills/` directory and Python backend.
> CellTypePilot ships **both** plugin manifests plus platform-specific components.

```
┌──────────────────────────────────────────────────────────────────────┐
│  Host agent:                                                        │
│  Claude Code                          │  Codex                      │
│  .claude-plugin/plugin.json           │  .codex-plugin/plugin.json  │
│  commands/  agents/  hooks/  rules/   │  skills/*/agents/openai.yaml│
│                 └──────────┬──────────┘                              │
│                            ↓                                        │
│  CellTypePilot plugin (shared Python backend)                       │
│  ┌─────────────────────────────────────────────────────────────────┐│
│  │  skills/       SKILL.md (4-stage orchestration, shared)         ││
│  │  .mcp.json     PubMed, bioRxiv (optional)                       ││
│  │  ─────────────────────────────────────────────────────────────── ││
│  │  Python backend (src/celltypepilot/)                             ││
│  │    doctor          Environment gate — what can run, what can't   ││
│  │    inspect         Data intelligence — species, tissue, keys     ││
│  │    annotate        Full pipeline:                                ││
│  │      ├─ Data Adapter         Load .h5ad/.rds, detect, validate   ││
│  │      ├─ Marker Knowledge Graph  80+ cell types, 11 tissues       ││
│  │      ├─ Marker Scorer        Wilcoxon DE + 5-dim scoring         ││
│  │      ├─ Annotation Critic    Independent evidence review         ││
│  │      ├─ Web Inspector        Flask interactive review panel      ││
│  │      ├─ Visualizer           UMAP, dotplot, confidence (Wong)    ││
│  │      ├─ Reporter             HTML report + methods paragraph     ││
│  │      ├─ Literature           PubMed validation (optional MCP)    ││
│  │      ├─ License Manager      Tiered: free/academic/commercial    ││
│  │      └─ Provenance           manifest.json versioning            ││
│  │    critic          Deep-review a specific cluster                ││
│  └─────────────────────────────────────────────────────────────────┘│
└──────────────────────────────────────────────────────────────────────┘
```

### Dual-platform, single backend

| Platform | Plugin manifest | Agent config | Discovery |
|---|---|---|---|
| **Claude Code** | `.claude-plugin/plugin.json` | `agents/*.md` + `commands/*.md` + `hooks/` + `rules/` | Plugin manager discovers via `.claude-plugin/` |
| **Codex** | `.codex-plugin/plugin.json` | `skills/*/agents/openai.yaml` + `AGENTS.md` | Plugin cache or repo-root auto-discovery |
| **Standalone** | N/A | N/A | `celltypepilot` CLI command |

All three modes call the same Python backend. The `skills/` directory is shared between
both platforms — each reads `SKILL.md` for workflow orchestration.

### Installation

```bash
# Claude Code — install as a plugin
git clone https://github.com/HERRY423/CellTypePilot ~/.claude/plugins/marketplaces/local/plugins/celltypepilot
cd ~/.claude/plugins/marketplaces/local/plugins/celltypepilot
pip install -e .
# Claude Code discovers .claude-plugin/plugin.json → skills/ + commands/ + agents/ + hooks/

# Codex — install as a plugin
git clone https://github.com/HERRY423/CellTypePilot ~/.codex/plugins/cache/local/celltypepilot
cd ~/.codex/plugins/cache/local/celltypepilot
pip install -e .
# Codex discovers .codex-plugin/plugin.json → skills/ + agents/openai.yaml

# Standalone CLI — no agent needed
git clone https://github.com/HERRY423/CellTypePilot && cd CellTypePilot
pip install -e .
celltypepilot doctor
```

Optional extras: `pip install -e ".[web]"` (Web Inspector), `"[seurat]"` (.rds support), `"[all]"` (everything).

### Plugin structure

```
CellTypePilot/
├── .claude-plugin/
│   └── plugin.json              ← Claude Code plugin manifest (v0.2.0)
├── .codex-plugin/
│   └── plugin.json              ← Codex plugin manifest (v0.2.0, with interface block)
├── skills/
│   └── celltypepilot/
│       ├── SKILL.md              ← Shared skill instructions (4-stage workflow)
│       ├── agents/
│       │   └── openai.yaml       ← Codex agent interface config
│       └── reference/            ← Reference docs
├── commands/                     ← Claude Code slash commands: /annotate, /critic, /inspect, /doctor
├── agents/                       ← Claude Code sub-agents: annotation-critic.md
├── hooks/                        ← Claude Code lifecycle hooks: session-start check
├── rules/                        ← Claude Code behavior rules: annotation-workflow
├── .mcp.json                     ← MCP servers (PubMed, bioRxiv)
├── AGENTS.md                     ← Codex agent instructions
├── src/celltypepilot/            ← Python backend (shared by all platforms)
└── tests/                        ← 18 smoke tests
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
- [x] **Phase 2** — Dual-platform plugin (Claude Code `.claude-plugin/` + Codex `.codex-plugin/`), slash commands, sub-agents, hooks, rules, MCP integration
- [x] **Phase 3** — Web Inspector (Flask interactive panel), Seurat .rds adapter, tiered license system (free/academic/commercial), premium atlas (tumor/brain/immune)
- [ ] **Phase 4** — Tiered consensus orchestrator (Tier 0 → Tier 1 adaptive upgrade), doublet detection enhancement, ontology validation enhancement
- [ ] **Phase 5** — docx/pptx submission package, reference mapping methods, extended atlas subscription, team sharing, rare cell type mining

## License

MIT
