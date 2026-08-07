# CellTypePilot

[![CI](https://github.com/HERRY423/CellTypePilot/actions/workflows/ci.yml/badge.svg)](https://github.com/HERRY423/CellTypePilot/actions/workflows/ci.yml)
![Tests](https://img.shields.io/badge/tests-266%20passed-brightgreen)
![Coverage](https://img.shields.io/badge/coverage-75%25-blue)
![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Single-cell annotation trust layer that parasitizes on your existing coding agent.**
> For individual researchers and small labs — no standalone app, no heavy infrastructure.

**CellTypePilot** is a **plugin** for Claude Code / OpenAI Codex. It turns pre-clustered
single-cell data into trusted, publication-ready cell-type annotations — with **dual-engine**
scoring (marker overlap + reference embedding), adaptive ensemble fusion, confidence levels,
an independent critic review, and a draft methodology paragraph for your paper.

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
| **Rare / transitional states forced into a label** | Reference embedding mapping captures continuous trajectories; critic flags doublet signals and transitional states instead of silently assigning a wrong label |
| **Token cost spirals** | Deterministic marker scoring is free (Tier 0); expensive multi-model consensus only triggers for genuinely ambiguous clusters (Tier 1) |
| **Results unreproducible** | `manifest.json` records knowledge graph version, parameters, data hash, and output hashes for every run |
| **Metadata is messy** (Ensembl IDs, mixed gene naming, non-standard tissue columns) | Robust auto-detection: Ensembl prefix voting across 9 species, gene-symbol convention fallback, case-insensitive tissue matching with synonym groups |

## What you get

```
output/
├── data.annotated.h5ad          # ctp_cell_type, ctp_cl_id, ctp_confidence in obs
├── evidence_table.csv           # Per-cluster: scores, markers, critic flags, confidence
├── ensemble_scores.csv          # Per-cell-type: marker + ref + ensemble scores
├── transitional_states.csv      # Clusters flagged as differentiation intermediates
├── disagreements.csv            # Marker vs reference disagreement analysis
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

# 3. Inspect your data — robustly auto-detects species, tissue, clusters, embeddings
celltypepilot inspect --input data.h5ad

# 4. Annotate — full pipeline in one command
celltypepilot annotate --input data.h5ad --cluster-key leiden --tissue blood

# 5. Annotate with reference embedding (resolves trajectories & rare states)
celltypepilot annotate-embedding --input data.h5ad --cluster-key leiden \
    --reference atlas.h5ad --tissue blood

# 6. Check available reference scoring backends
celltypepilot backends

# 7. Review interactively in browser (optional)
celltypepilot inspect-web --output output/

# 8. Deep-review a flagged cluster
celltypepilot critic --input data.h5ad --cluster-key leiden --focus cluster_7

# 9. Validate markers against literature (optional, needs network)
celltypepilot literature --cell-type "T cells" --markers "CD3E,CD4,CD8A"

# 10. Convert Seurat .rds → .h5ad (optional, needs R)
celltypepilot convert-rds --input data.rds --output data.h5ad
```

Or just tell your agent: *"annotate my clusters in data.h5ad"* — the plugin instructions
guide Claude Code / Codex through the full workflow automatically.

## Robust automatic detection

Real-world datasets are messy. CellTypePilot's inspection layer is built to handle it:

**Species detection** — multi-signal voting, not brittle pattern matching:

| Signal | How it works |
|---|---|
| Ensembl ID prefixes | Longest-prefix voting over 9 species: `ENS` (human), `ENSMUS` (mouse), `ENSRNO` (rat), `ENSDAR` (zebrafish), `ENSGAL` (chicken), `ENSSSC` (pig), `ENSBTA` (cow), `ENSMMU` (macaque), `ENSCFA` (dog) |
| Gene symbol conventions | ALL-CAPS (e.g. `CD3E`, `S100A8`) → human; capitalized (e.g. `Cd3e`) → mouse — robust to symbols containing digits |
| Mixed naming | Per-gene voting with clear majority rule; ambiguous data defaults safely and is reported, never silently mislabeled |
| Explicit override | `--species human\|mouse` always wins |

**Tissue detection** — case-insensitive matching over `obs` metadata with synonym groups
(e.g. `pbmc` / `peripheral blood` / `blood` all resolve to the blood atlas), so columns
named `Tissue`, `tissue_type`, or `organ` all work. `--tissue` overrides when needed.

**Orchestrator architecture** — pipeline business logic (loading, scoring coordination,
override application, progress reporting) lives in a dedicated `orchestrator` module,
shared by both the CLI and the Web Inspector. The CLI only parses arguments and renders
output — no duplicated pipeline code.

**Templates, not hardcoded HTML** — the HTML report and Web Inspector dashboard are
rendered from Jinja2 templates (`src/celltypepilot/templates/`), keeping markup, JS,
and CSS cleanly separated from Python logic.

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
│  │    orchestrator    Pipeline logic shared by CLI & Web Inspector  ││
│  │    annotate        Full pipeline:                                ││
│  │      ├─ Data Adapter         Load .h5ad/.rds, detect, validate   ││
│  │      ├─ Marker Knowledge Graph  80+ cell types, 11 tissues       ││
│  │      ├─ Marker Scorer        Wilcoxon DE + 5-dim scoring         ││
│  │      ├─ Reference Scorer     CellTypist/scANVI/KNN/Correlation   ││
│  │      ├─ Ensemble Scorer      Adaptive fusion + disagreement      ││
│  │      ├─ Annotation Critic    Independent evidence review         ││
│  │      ├─ Web Inspector        Flask review panel (Jinja2 templates)││
│  │      ├─ Visualizer           UMAP, dotplot, confidence (Wong)    ││
│  │      ├─ Reporter             Jinja2 HTML report + methods text   ││
│  │      ├─ Literature           PubMed validation (optional MCP)    ││
│  │      ├─ License Manager      RSA-2048 signed, machine-bound      ││
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

Optional extras: `pip install -e ".[web]"` (Web Inspector), `"[seurat]"` (.rds support),
`"[reference]"` (CellTypist), `"[embedding]"` (scVI/scANVI), `"[all]"` (everything).

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
│   ├── orchestrator.py           ← Pipeline business logic (CLI + Web Inspector)
│   ├── templates/                ← Jinja2 templates (HTML report, dashboard)
│   ├── data/                     ← Marker knowledge graph + premium atlas
│   └── ...                       ← Scorers, critic, visualizer, reporter, etc.
├── tests/                        ← 266 tests (all passing, ~75% coverage)
└── .github/workflows/            ← CI: ruff lint/format + test matrix (3.10–3.12)
```

### Continuous integration

Every push runs GitHub Actions CI: `ruff check` + `ruff format --check` for code quality,
then the full test suite across Python 3.10 / 3.11 / 3.12 with coverage reporting.

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

## Reference Embedding + Ensemble Fusion

For continuous differentiation trajectories (stem → progenitor → mature) and rare
transitional states, pure marker overlap scoring can fail. CellTypePilot addresses this
with a **dual-engine** architecture:

**Engine 1 — Marker Scorer** (deterministic):
Wilcoxon DE + 5-dim scoring against the built-in marker knowledge graph. Zero cost,
fully reproducible, works offline.

**Engine 2 — Reference Scorer** (deep learning):
Projects query cells into a reference embedding space and transfers labels. Four backends:

| Backend | Method | Best for | Dependencies |
|---|---|---|---|
| **CellTypist** | Pre-trained logistic regression on CellxGene Census | Standard human/mouse tissues | `celltypist` |
| **scANVI** | Semi-supervised VAE with custom reference atlas | Custom atlases, cross-species | `scvi-tools` |
| **KNN** | PCA + inverse-distance KNN label transfer | Quick mapping, no model needed | `sklearn` (always available) |
| **Correlation** | Pearson correlation with reference mean profiles | Lightweight fallback | None (always available) |

**Ensemble Fusion** — Adaptive weighting combines both engines:

| Marker confidence | Marker weight | Reference weight | Rationale |
|---|---|---|---|
| ≥ 0.6 (high) | 0.70 | 0.30 | Markers are reliable → trust them |
| 0.3–0.6 (medium) | 0.50 | 0.50 | Balanced fusion |
| ≤ 0.3 + ref ≥ 0.5 | 0.20 | 0.80 | Markers fail → reference override |
| Strong disagreement | 0.15 | 0.85 | Reference confident, markers wrong |

**Transitional state detection** — Three signals identify differentiation intermediates:
1. **Cross-ranking**: each method's top-1 appears in the other's top-N
2. **Distribution entropy**: broad probability = diffuse identity
3. **Confidence asymmetry**: one method confident, other uncertain

Disagreements are flagged with biological interpretation (novel subtype, transitional
state, marker database gap, or low-quality cluster).

```bash
# Use CellTypist pre-trained model
celltypepilot annotate-embedding -i data.h5ad -k leiden -m Immune_All_Low.pkl

# Use custom reference atlas with auto-selected backend
celltypepilot annotate-embedding -i data.h5ad -k leiden -r atlas.h5ad

# Force specific backend
celltypepilot annotate-embedding -i data.h5ad -k leiden -r atlas.h5ad -b knn
```

## Annotation Critic — the soul of the plugin

The Critic doesn't just score — it *doubts*. Every annotation is independently reviewed
across **6 checks** (including ensemble agreement when reference embedding is available):

| Check | What it catches |
|---|---|
| Evidence sufficiency | < 20% marker coverage → `LOW_EVIDENCE`; 20–50% → `PARTIAL_EVIDENCE` |
| Negative marker conflict | Negative markers expressed in > 15% of cells → `NEG_MARKER_CONFLICT` |
| Doublet signal | Two mutually exclusive lineage signatures co-expressed → `POSSIBLE_DOUBLET` |
| Ontology consistency | Is the label a valid Cell Ontology term in the right tissue context? |
| Ensemble agreement | Marker vs reference disagreement → `ENSEMBLE_DISAGREEMENT` / `ENSEMBLE_MILD_DISAGREEMENT` |
| Weak reference | Reference-only support with low confidence → `WEAK_REFERENCE_ONLY` |

Confidence levels: **high** / **medium** / **low** / **needs_review**. The Critic can only
downgrade, never upgrade. A flagged cluster is a *success* — it means the system caught
something worth your attention.

## Roadmap

- [x] **Phase 1 (MVP)** — h5ad adapter, marker knowledge graph, Wilcoxon DE scoring, Annotation Critic, doctor, figures, JSON output, HTML report, methodology draft, manifest provenance, literature validation (PubMed)
- [x] **Phase 2** — Dual-platform plugin (Claude Code `.claude-plugin/` + Codex `.codex-plugin/`), slash commands, sub-agents, hooks, rules, MCP integration
- [x] **Phase 3** — Web Inspector (Flask interactive panel), Seurat .rds adapter, tiered license system (free/academic/commercial), premium atlas (tumor/brain/immune)
- [x] **Phase 4** — Reference Embedding + Ensemble fusion (CellTypist / scANVI / KNN / Correlation backends), adaptive weighting, transitional state detection, ensemble-aware critic, RSA-2048 license security, sparse-preserving Seurat conversion, Web Inspector override API
- [x] **Architecture hardening** — Orchestrator layer (pipeline logic extracted from CLI), Jinja2 templates (HTML/JS out of Python), robust multi-species Ensembl detection (9 species), synonym-based tissue detection, 266 tests at ~75% coverage, GitHub Actions CI (ruff + Python 3.10–3.12 matrix)
- [ ] **Phase 5** — Tiered consensus orchestrator (Tier 0 → Tier 1 adaptive upgrade), docx/pptx submission package, extended atlas subscription, team sharing, rare cell type mining

## License

MIT
