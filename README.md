# CellTypePilot

> Single-cell annotation intelligence layer for coding agents

**CellTypePilot** is a Claude Code / Codex skill that turns pre-clustered single-cell data into
trusted, publication-ready cell-type annotations — with evidence, confidence levels, and an
independent critic review.

## What makes it different

CellTypePilot is **not** another cell-type annotation algorithm. It's a **trust layer** that
runs inside your existing coding agent:

1. **Zero-friction entry** — `git clone` + `pip install -e .` + run. No MCP servers, no pixi,
   no conda environment to configure for the basic path.
2. **Evidence over black-box** — Every annotation comes with marker evidence, expression
   statistics, negative marker checks, and an independent critic verdict.
3. **Annotation Critic** — A skeptical review module that checks evidence sufficiency,
   negative marker conflicts, doublet signals, and ontology consistency. It flags problems
   instead of silently giving wrong answers.
4. **Publication-ready output** — UMAP, dotplot, confidence figures, evidence table, and a
   draft methodology paragraph ready for your paper.

## Quick start

```bash
# Install
pip install -e .

# Check your environment
celltypepilot doctor

# Inspect your data
celltypepilot inspect --input data.h5ad

# Run annotation
celltypepilot annotate --input data.h5ad --cluster-key leiden --tissue blood

# Deep-review a flagged cluster
celltypepilot critic --input data.h5ad --cluster-key leiden --focus cluster_7

# List available markers
celltypepilot markers --tissue blood
```

## Output

```
output/
├── data.annotated.h5ad          # Annotated data (ctp_cell_type, ctp_cl_id, ctp_confidence in obs)
├── evidence_table.csv           # Per-cluster evidence and critic results
├── report_draft.html            # Comprehensive HTML report with figures
├── methodology_draft.txt        # Draft methods paragraph for papers
├── manifest.json                # Run provenance (versions, params, hashes)
└── figures/
    ├── umap_cluster.png         # UMAP by cluster
    ├── umap_celltype.png        # UMAP by cell type
    ├── umap_confidence.png      # UMAP by confidence level
    ├── marker_dotplot.png       # Marker gene dotplot
    └── confidence_distribution.png
```

## Supported tissues (built-in)

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

Both human and mouse are supported with automatic species detection and gene symbol conversion.

## As a coding agent skill

### Claude Code

```bash
# Install as a Claude Code skill
git clone https://github.com/your-org/celltypepilot ~/.claude/skills/celltypepilot
cd ~/.claude/skills/celltypepilot && pip install -e .

# Then in Claude Code, just say:
# "annotate my clusters in data.h5ad"
# "what cell types are these?"
# "run celltypepilot on my PBMC data"
```

### OpenAI Codex CLI

```bash
# Clone the repo
git clone https://github.com/your-org/celltypepilot ~/celltypepilot
cd ~/celltypepilot && pip install -e .

# Codex reads AGENTS.md automatically from the repo root.
# Then in Codex, just say:
# "annotate my clusters in data.h5ad"
# "run celltypepilot doctor"

# Or use the quick-start script:
bash scripts/codex_annotate.sh data.h5ad leiden ./output blood
# Windows:
.\scripts\codex_annotate.ps1 -Input data.h5ad -ClusterKey leiden -OutputDir .\output -Tissue blood
```

Both platforms share the same Python backend (`src/celltypepilot/`) and CLI commands.
The only difference is the agent instruction format:
- Claude Code → `skills/celltypepilot/SKILL.md` (YAML frontmatter + orchestration)
- Codex CLI → `AGENTS.md` (pure markdown instructions)

## Architecture

```
User → Claude Code (SKILL.md)  OR  Codex CLI (AGENTS.md)
         │                              │
         └──────────┬───────────────────┘
                    ↓
         celltypepilot CLI (shared Python backend)
         │
         ├─ celltypepilot doctor          # Environment check
         ├─ celltypepilot inspect         # Data inspection
         ├─ celltypepilot annotate        # Full pipeline:
         │    ├─ Data Adapter             #   Load, detect, validate
         │    ├─ Marker Knowledge Graph   #   Built-in curated atlas
         │    ├─ Marker Scorer            #   DE + overlap + specificity
         │    ├─ Annotation Critic        #   Evidence review + flags
         │    ├─ Visualizer               #   UMAP, dotplot, confidence
         │    ├─ Reporter                 #   HTML report + evidence table
         │    └─ Provenance               #   Manifest + versioning
         └─ celltypepilot critic          # Deep-review specific cluster
```

## Roadmap

- [x] **Phase 1 (MVP)**: h5ad adapter + marker atlas + scoring + critic + doctor + figures + JSON
- [ ] **Phase 2**: Tiered consensus orchestrator + doublet detection + ontology validation
- [ ] **Phase 3**: docx/pptx export + manifest polish + Seurat/R support
- [ ] **Phase 4**: Extended atlas subscription + team sharing + rare cell type mining

## License

MIT
