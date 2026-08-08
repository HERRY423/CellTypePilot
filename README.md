# CellTypePilot

[![CI](https://github.com/HERRY423/CellTypePilot/actions/workflows/ci.yml/badge.svg)](https://github.com/HERRY423/CellTypePilot/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Local-first single-cell annotation review plugin for your existing coding workspace.**
> For individual researchers and small labs — no standalone app, no heavy infrastructure.

**CellTypePilot** is a **plugin** for Claude Code / OpenAI Codex. It turns pre-clustered
single-cell data into auditable draft cell-type annotations — with a governed context interface,
**dual-engine** identity scoring (marker overlap + reference embedding), an independent cell-state
lens, conservative critic review, and a draft methodology paragraph for your paper.

It is **not an autonomous analysis agent**. It is a deterministic, artifact-producing plugin
that adds evidence review, conservative abstention, benchmark hooks, and provenance to the
coding workspace you already use. A qualified human owns the final biological decision.

## Current validation boundary

| Available in v0.3.0 | Not yet claimed |
|---|---|
| Direction-, log2FC-, FDR-, and expression-fraction-gated DE evidence | Biological superiority over CellTypist, SingleR, Azimuth, or popV |
| Complete expected-marker denominators, with missing genes separated from present-but-silent genes | A completed public multi-study/donor benchmark |
| Explicit `Unknown`/abstain decisions; low-evidence clusters cannot pass | Calibrated accuracy for every tissue, platform, disease state, or species |
| One annotation pipeline for scoring, critic review, write-back, report, and manifest | Primary-source verification of every bundled marker relationship |
| Locked study/donor holdout runner and comparator adapters | Clinical-grade or fully automated biological decisions |
| Governed Context Pack: structured custom markers share the normal evidence gates; free text is provenance-only | That free-text biological context is validated evidence |
| Separate identity and state outputs; state scoring cannot overwrite identity or rescue an abstained identity | Calibrated cell-state accuracy or comprehensive state coverage |

For `mkg-2026.08`, all 733 bundled marker relationships are currently classified as
`aggregate_source_only_not_edge_verified`. They remain usable under the exploratory
`database` evidence policy, but are excluded by stricter `edge_verified` and `primary`
policies. The curation backlog is published in
[`docs/atlas_curation_queue.csv`](docs/atlas_curation_queue.csv).

## Why it exists

For individual researchers without dedicated bioinformatics support, cell-type annotation is
painful — not because the algorithms don't exist, but because:

| Pain point | What CellTypePilot does |
|---|---|
| **No one reviews your annotations** | A rules-based Annotation Critic checks evidence sufficiency, negative marker conflicts, doublet signals, and ontology consistency — *before* you trust a label |
| **Config barrier is too high** (MCP, pixi, conda...) | `git clone` + `pip install -e .` + run. Zero MCP required for the basic path. `doctor` tells you what you have *before* anything fails |
| **Workflow fragmentation** (scripts here, tool there) | Runs inside your Claude Code / Codex session — no context switch to a separate app |
| **Can't explain *why* a cluster got its label** | Every annotation ships with: supporting markers, expression stats, critic flags, and a draft methods paragraph |
| **Rare / transitional states forced into a label** | Low or conflicting evidence produces an explicit `Unknown`/`abstain`; the best candidate is retained separately for review |
| **Disease context is either ignored or blindly trusted** | Free text is recorded but never counted as evidence; structured custom markers are hashed, scoped, and subjected to the same DE and critic gates |
| **Identity and state are conflated** | Canonical identity and exploratory state are written on independent axes; an `Unknown` identity can retain a supported state without becoming a cell type |
| **Token cost spirals** | Annotation, critic checks, calibration, and reporting are deterministic local code; no LLM call is required |
| **Results unreproducible** | `manifest.json` records knowledge graph version, parameters, data hash, and output hashes for every run |
| **Metadata is messy** (Ensembl IDs, mixed gene naming, non-standard tissue columns) | Robust auto-detection: Ensembl prefix voting across 9 species, gene-symbol convention fallback, case-insensitive tissue matching with synonym groups |

## What you get

```
output/
├── data.annotated.h5ad          # final label, candidate, decision, reason, confidence in obs
├── evidence_table.csv           # Per-cluster: scores, markers, critic flags, confidence
├── state_results.csv            # Independent state candidate, decision, score, evidence
├── context_pack.normalized.json # Optional: normalized, scoped, hashed user context
├── ensemble_scores.csv          # Per-cell-type: marker + ref + ensemble scores
├── transitional_states.csv      # Clusters flagged as differentiation intermediates
├── disagreements.csv            # Marker vs reference disagreement analysis
├── report_draft.html            # Self-contained HTML report with all figures embedded
├── methodology_draft.txt        # "We annotated N clusters using CellTypePilot v0.3.0..."
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
celltypepilot annotate --input data.h5ad --cluster-key leiden \
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

# 11. Lock study/donor holdouts; evaluate imported out-of-fold predictions
celltypepilot benchmark -i benchmark.h5ad --truth-key truth \
    --study-key study --donor-key donor --output benchmark/

# 12. Actually execute CellTypePilot/CellTypist on every isolated fold.
# SingleR/Azimuth/popV use explicit JSON argv adapters under the same protocol.
celltypepilot benchmark-run -i benchmark.h5ad --truth-key truth \
    --study-key study --donor-key donor --cluster-key leiden \
    --species human --tissue blood --methods celltypepilot,celltypist

# 13. Fit a downgrade-only abstention policy on a separate calibration dataset
celltypepilot calibrate -i calibration.h5ad --truth-key truth \
    --predictions calibration_predictions.csv -o abstention_policy.json
celltypepilot annotate -i query.h5ad -k leiden -t blood \
    --calibration-policy abstention_policy.json

# 14. Add disease context without turning prose into evidence. Structured markers
# are accepted only through a versioned JSON Context Pack and/or marker CSV.
celltypepilot annotate -i kidney_iri.h5ad -k leiden -s human -t kidney \
    --context "Post-ischemic reperfusion injury" \
    --context-file context.json --custom-markers custom_markers.csv

# 15. Disable exploratory state scoring while keeping the identity pipeline unchanged
celltypepilot annotate -i data.h5ad -k leiden -t blood --no-states
```

Or invoke the plugin from Claude Code / Codex with *"annotate my clusters in data.h5ad"*.
The host follows the bundled workflow instructions while the local Python backend produces
the annotation artifacts.

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

CellTypePilot integrates with the coding workspace you already use. It does not require a
separate annotation application or autonomous analysis service.

> **Why "plugin" and not "agent"?** "Plugin" is the product concept — a self-contained,
> user-invoked capability bundle that attaches to a host workspace. Each platform has its own native
> plugin format: Claude Code uses `.claude-plugin/plugin.json`, Codex uses
> `.codex-plugin/plugin.json`. Both share the same `skills/` directory and Python backend.
> CellTypePilot ships **both** plugin manifests plus platform-specific components.

```
┌──────────────────────────────────────────────────────────────────────┐
│  Host workspace:                                                    │
│  Claude Code                          │  Codex                      │
│  .claude-plugin/plugin.json           │  .codex-plugin/plugin.json  │
│  commands/  hooks/  rules/            │  skills/*/agents/openai.yaml│
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
│  │      ├─ Governed Context Pack  Scoped, hashed custom hypotheses  ││
│  │      ├─ Marker Scorer        Wilcoxon DE + 5-dim scoring         ││
│  │      ├─ Reference Scorer     CellTypist/scANVI/KNN/Correlation   ││
│  │      ├─ Ensemble Scorer      Adaptive fusion + disagreement      ││
│  │      ├─ Annotation Critic    Rules-based same-run review         ││
│  │      ├─ State Lens           Independent, identity-invariant axis ││
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

| Platform | Plugin manifest | Host integration | Discovery |
|---|---|---|---|
| **Claude Code** | `.claude-plugin/plugin.json` | `commands/*.md` + `hooks/` + `rules/` | Plugin manager discovers via `.claude-plugin/` |
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
# Claude Code discovers .claude-plugin/plugin.json → skills/ + commands/ + hooks/

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
│   └── plugin.json              ← Claude Code plugin manifest (v0.3.0)
├── .codex-plugin/
│   └── plugin.json              ← Codex plugin manifest (v0.3.0, with interface block)
├── skills/
│   └── celltypepilot/
│       ├── SKILL.md              ← Shared skill instructions (4-stage workflow)
│       ├── agents/
│       │   └── openai.yaml       ← Codex agent interface config
│       └── reference/            ← Reference docs
├── commands/                     ← Claude Code slash commands: /annotate, /critic, /inspect, /doctor
├── hooks/                        ← Claude Code lifecycle hooks: session-start check
├── rules/                        ← Claude Code behavior rules: annotation-workflow
├── .mcp.json                     ← MCP servers (PubMed, bioRxiv)
├── AGENTS.md                     ← Codex agent instructions
├── src/celltypepilot/            ← Python backend (shared by all platforms)
│   ├── orchestrator.py           ← Pipeline business logic (CLI + Web Inspector)
│   ├── context_pack.py           ← Governed prior-context validation and hashing
│   ├── state_scorer.py           ← Independent cell-state scoring and invariants
│   ├── templates/                ← Jinja2 templates (HTML report, dashboard)
│   ├── data/                     ← Identity marker atlases + state atlas
│   └── ...                       ← Scorers, critic, visualizer, reporter, etc.
├── tests/                        ← Unit, contract, and scientific-boundary tests
└── .github/workflows/            ← CI: ruff lint/format + test matrix (3.10–3.12)
```

### Continuous integration

Every push runs GitHub Actions CI: `ruff check` + `ruff format --check` for code quality,
then the configured test suite across Python 3.10 / 3.11 / 3.12 with coverage reporting.
The CI badge reports the current default-branch state; it is not biological validation.

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

## Governed Context Pack

CellTypePilot accepts disease, region, timepoint, and experiment-specific knowledge without
letting a prompt bypass the evidence boundary:

- `--context` records free text for interpretation and provenance only. It never creates a
  marker, adds a score, unlocks an unsupported tissue, or changes an acceptance threshold.
- `--context-file` accepts the versioned `celltypepilot.context.v1` JSON schema.
- `--custom-markers` accepts a row-wise CSV with `axis,label,gene,polarity` plus optional
  `cl_id,parent_cell_types,source,review_status` columns.
- Structured identity markers enter the ordinary missing/silent, expression-fraction, positive
  log2FC, BH-FDR, negative-marker, and critic gates. Draft context-only support must abstain;
  reviewed context-only support is capped at medium confidence.
- The normalized pack, source-file hashes, schema version, scope, and canonical content hash are
  recorded in run artifacts. See [`docs/context_pack.md`](docs/context_pack.md).

## Identity × State

Cell identity and cell state are deliberately separate outputs. `ctp_cell_type` and `ctp_cl_id`
remain the conservative canonical identity. State fields (`ctp_cell_state_candidate`,
`ctp_state_decision`, `ctp_cell_state`, and `ctp_state_evidence`) are scored independently from a
versioned state atlas or structured context hypotheses. The merge step asserts that identity,
candidate, decision, and abstention-reason columns are byte-for-byte unchanged.

This permits `Unknown · interferon_responsive` when state evidence is sufficient but lineage
identity is not. It never permits a state to rescue or overwrite an abstained identity. Bundled
state modules currently have aggregate source-level provenance and are exploratory, not a claim
of calibrated state accuracy. See [`docs/state_lens.md`](docs/state_lens.md).

## Reference Embedding + Ensemble Fusion

For continuous differentiation trajectories (stem → progenitor → mature) and rare
transitional states, pure marker overlap scoring can fail. CellTypePilot addresses this
with a **dual-engine** architecture:

**Engine 1 — Marker Scorer** (deterministic):
Wilcoxon DE against the built-in marker knowledge graph. A supporting positive marker
must pass direction, log2FC ≥ 0.5, BH-FDR ≤ 0.05, and expression in ≥25% of cluster cells.
Coverage is divided by the complete expected marker panel, including genes absent from
the matrix. The scorer is deterministic, reproducible, and offline.

**Engine 2 — Reference Scorer** (deep learning):
Projects query cells into a reference embedding space and transfers labels. Four backends:

| Backend | Method | Best for | Dependencies |
|---|---|---|---|
| **CellTypist** | Pre-trained or explicitly selected logistic-regression model | Model-matched tissues only | `celltypist` |
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
celltypepilot annotate -i data.h5ad -k leiden -m Immune_All_Low.pkl

# Use custom reference atlas with auto-selected backend
celltypepilot annotate -i data.h5ad -k leiden -r atlas.h5ad

# Force specific backend
celltypepilot annotate -i data.h5ad -k leiden -r atlas.h5ad -b knn
```

## Annotation Critic — the soul of the plugin

The Critic doesn't just score — it *doubts*. Every annotation is reviewed by explicit rules
across **6 checks** (including ensemble agreement when reference embedding is available):

| Check | What it catches |
|---|---|
| Evidence sufficiency | < 20% expected-marker coverage → `LOW_EVIDENCE`; 20–50% → `PARTIAL_EVIDENCE`; both abstain |
| Negative marker conflict | Negative markers expressed in > 20% of cells → `NEG_MARKER_CONFLICT` |
| Doublet signal | Two mutually exclusive lineage signatures co-expressed → `POSSIBLE_DOUBLET` |
| Ontology consistency | Does the label exactly match the CL identifier declared in the versioned atlas? Live ontology resolution is not claimed |
| Ensemble agreement | Marker vs reference disagreement → `ENSEMBLE_DISAGREEMENT` / `ENSEMBLE_MILD_DISAGREEMENT` |
| Weak reference | Reference-only support with low confidence → `WEAK_REFERENCE_ONLY` |

Confidence levels: **high** / **medium** / **low** / **needs_review**. The Critic can only
downgrade, never upgrade. A flagged cluster is a *success* — it means the system caught
something worth your attention.

The atlas v2 schema records gene, polarity, species, tissue, state, atlas version,
PMID/DOI/URL, and verification status for every bundled marker relationship. Existing
relationships are honestly marked `aggregate_source_only_not_edge_verified`: database-paper
provenance is present, but a marker-specific primary experiment has not been claimed.
Use `--marker-evidence-policy edge_verified` or `primary` to exclude relationships below the
requested evidence tier; the default `database` policy keeps them for exploratory draft labeling
and marks the critic result `AGGREGATE_PROVENANCE_ONLY`.

## Roadmap

- [x] **Phase 1 (MVP)** — h5ad adapter, marker knowledge graph, Wilcoxon DE scoring, Annotation Critic, doctor, figures, JSON output, HTML report, methodology draft, manifest provenance, literature validation (PubMed)
- [x] **Phase 2** — Dual-platform plugin packaging (Claude Code `.claude-plugin/` + Codex `.codex-plugin/`), commands, hooks, rules, and optional literature integration
- [x] **Phase 3** — Web Inspector (Flask interactive panel), Seurat .rds adapter, tiered license system (free/academic/commercial), premium atlas (tumor/brain/immune)
- [x] **Phase 4** — Reference Embedding + Ensemble fusion (CellTypist / scANVI / KNN / Correlation backends), adaptive weighting, transitional state detection, ensemble-aware critic, RSA-2048 license security, sparse-preserving Seurat conversion, Web Inspector override API
- [x] **Architecture hardening** — Orchestrator layer (pipeline logic extracted from CLI), Jinja2 templates, multi-species detection, synonym-based tissue detection, and Python 3.10–3.12 CI
- [x] **Phase 5** — Governed Context Pack, custom marker trust boundary, legal identity ontology IDs, and independent Identity × State outputs
- [ ] **Validation release** — Published multi-study/donor benchmark, calibrated abstention card, verified reference registry, and marker-edge curation coverage report

## License

MIT
