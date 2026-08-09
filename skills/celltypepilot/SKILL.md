---
name: celltypepilot
description: >-
  Single-cell cell-type annotation with evidence, conservative abstention, and reviewable drafts.
  Use this whenever the user wants to annotate, label, or identify cell types for pre-clustered
  single-cell or spatial transcriptomics data (.h5ad), or says things like "annotate my clusters",
  "what cell types are these?", "label my scRNA-seq", "figure out the cell types", or "run CellTypePilot".
  CellTypePilot runs through the host coding workspace — no separate app, no MCP required for basic use.
  It provides governed study context, marker-based identity scoring, an independent state lens,
  a rules-based Annotation Critic that flags doubtful calls, colorblind-friendly figures, and
  a draft methodology paragraph ready for papers.
license: MIT
---

# CellTypePilot

Turn pre-clustered single-cell data into **auditable draft cell-type annotations** —
with evidence, confidence levels, and a rules-based critic review. CellTypePilot is a plugin,
not an autonomous analysis agent: its deterministic backend produces candidates and durable
artifacts; a qualified human remains responsible for final biological adjudication.

## Operating principles

- **Zero-friction entry.** `git clone` + `pip install -e .` + run. No MCP servers, no pixi, no
  conda environment to configure for the basic path. The `doctor` command tells users exactly
  what they have and what's missing, *before* anything fails.
- **Evidence over black-box.** Every annotation comes with: which markers support it, what
  fraction of cells express them, whether negative markers conflict, and a rules-based
  critic verdict. No label without evidence.
- **Governed context.** Free text is retained for provenance only. Only explicit structured
  marker hypotheses enter scoring, under the same statistical and critic gates as atlas markers.
- **Identity and state are separate.** State output cannot overwrite or rescue the canonical
  identity decision; `Unknown` identity may retain a supported state for review.
- **Critic is the soul.** The Annotation Critic doesn't just score — it *doubts*. It checks
  evidence sufficiency, negative marker conflicts, doublet signals, and ontology consistency.
  A cluster flagged by the critic is a success, not a failure.
- **Cost-aware.** Deterministic marker scoring is free and reproducible. LLM reasoning is
  available as an optional second opinion, not a mandatory expense.
- **Artifact-centric.** All heavy compute writes durable files (CSV, PNG, JSON). You reason
  over the files. Runs are resumable and auditable.
- **Review-ready draft output.** UMAP, dotplot, confidence figures, evidence table, and a
  draft methodology paragraph are generated for qualified human review and revision.

## Agent-host golden workflow: exactly four steps

When the local MCP facade is available, ordinary Agent hosts must use only this
stateful sequence:

1. `prepare_annotation` — inspect, resolve blockers, preflight evidence
   addressability, and write a hashed `annotation_plan.json`.
2. `annotate_from_plan` — execute that plan without changing species, tissue,
   clusters, packs, evidence policy, or input hash.
3. `review_uncertain_clusters` — create a bounded, read-only review queue. The
   Agent may explain evidence but must not select or apply a biological label.
4. `finalize_reviewed_annotations` — apply only explicitly supplied human
   decisions, require a named human signer, reconcile evidence rows, regenerate
   derived artifacts, and write a fresh review signature.

Step 2 also writes `contrastive_evidence.csv` and `evidence_gaps.json`. Step 3
must use those artifacts to explain the existing top-two ranking and turn each
`Unknown` into bounded follow-up actions; it must never infer or apply a
replacement identity.

Do not assemble an ordinary annotation workflow from advanced primitives. The
CLI path below remains the compatibility/manual path when MCP is unavailable.
All four MCP results use `celltypepilot.agent-decision.v1`. Follow only
`allowed_next_actions`, honor `blockers`, `forbidden_claims`, and
`human_action_required`, and treat `artifact_paths` as the durable evidence
handoff. A score margin is never a probability or an override instruction.

## CLI compatibility workflow

### Stage 1 — Environment gate & data inspection

**1a. Run `celltypepilot doctor` FIRST.** This checks Python version, all dependencies,
and reports the capability tier (full / degraded). If core deps are missing, tell the user
to `pip install celltypepilot` before proceeding. Do NOT attempt annotation without passing
the doctor check.

**1b. Inspect the data.** Run `celltypepilot inspect --input <path>` (add `--json` for
structured output). This reports:
- Species (auto-detected from gene naming conventions)
- Tissue (from obs metadata if available)
- Cluster keys found (leiden, louvain, etc.)
- Embedding keys found (UMAP, tSNE, etc.)
- Layer info (counts vs log-normalized)
- Gene ID convention (human symbols, mouse symbols, Ensembl)
- Cluster sizes and any warnings

**1c. Confirm parameters with the user.** Use AskUserQuestion to batch 2-3 questions:
- Confirm or correct the auto-detected species and tissue
- Confirm the cluster key (if multiple found, let user choose)
- Confirm the embedding key for visualization
- Ask about tissue context if not auto-detected (this determines which marker atlas to use)
- Ask whether the user has disease/region/timepoint context or a custom marker panel. Use
  `--context-file`/`--custom-markers` for evidence-bearing hypotheses; prose alone is not evidence.

### Stage 2 — Annotate

**2a. Run the annotation pipeline.** Execute:
```
celltypepilot annotate \
  --input <path> \
  --cluster-key <key> \
  --output <output_dir> \
  --species <human|mouse> \
  --tissue <tissue> \
  --embedding-key <key> \
  [--context <text>] \
  [--context-file <context.json>] \
  [--custom-markers <markers.csv>]
```

This runs the full pipeline:
1. Load data and compute hash (for provenance)
2. Auto-detect or confirm species/tissue
3. Load the built-in Marker Knowledge Graph (MKG) for the species/tissue
4. Compute marker scores: DE analysis (Wilcoxon) + marker overlap + specificity + negative marker check
5. Run the Annotation Critic: evidence sufficiency, negative marker conflict, doublet heuristic, ontology check
6. Score cell states on an independent, identity-invariant output axis
7. Generate figures: UMAP (clusters, cell types, confidence), dotplot, confidence distribution
8. Save outputs: annotated .h5ad, evidence/state tables, HTML report, manifest, methodology draft

**2b. Review the results.** Read the evidence table and HTML report. Present the key findings:
- How many clusters were annotated
- Confidence distribution (how many high/medium/low/needs_review)
- Which clusters were flagged by the critic and why
- Show the UMAP figures

**2c. Let the user review.** Point the user at:
- The HTML report (interactive, with all figures embedded)
- Flagged clusters that may need manual review
- The methodology draft paragraph they can adapt for their paper

### Stage 3 — Deep review (optional)

If the user wants to investigate a specific flagged cluster:
```
celltypepilot critic --input <path> --cluster-key <key> --focus <cluster_id>
```

This shows:
- The top-5 candidate cell types with scores
- Detailed critic evidence for the focus cluster
- Which markers support/conflict with each candidate

### Stage 4 — Literature validation (optional, requires MCP or network)

If the user wants literature-backed validation of annotations, or if the critic flagged clusters:

**4a. Check MCP/literature availability:**
```
celltypepilot doctor
```
Look for the "MCP / Literature Integration" section. If `pubmed_direct` is available, you can use the literature command.

**4b. Validate a specific annotation:**
```
celltypepilot literature --cell-type "T cells" --markers "CD3E,CD4,CD8A" --json
```
This searches PubMed for literature supporting each marker-cell_type association.

**4c. Use MCP tools for richer search.** If MCP servers are configured (PubMed, bioRxiv, etc.):
- Use the MCP search tools with queries like: `"CD3E" "T cells" marker expression single-cell`
- Cross-reference results with the built-in MKG
- Add literature citations to the methodology paragraph

**4d. Present literature evidence.** Show the user:
- Which markers have literature support
- Total references found per marker
- Any conflicting evidence
- Suggested citations for the paper

## Command reference

| Command | Purpose |
|---|---|
| `celltypepilot doctor` | Check environment, dependencies, capability tier, MCP status (`--json` for hosts) |
| `celltypepilot inspect -i <path>` | Inspect h5ad: species, tissue, clusters, embeddings, layers |
| `celltypepilot annotate -i <path> -k <key>` | Full annotation pipeline |
| `celltypepilot benchmark ... --evaluation-unit cell\|cluster\|both` | Keep cell and cluster benchmark endpoints separate |
| `celltypepilot critic -i <path> -k <key> -f <cluster>` | Deep-review a specific cluster |
| `celltypepilot markers -t <tissue>` | List available cell types and markers |
| `celltypepilot literature -c <type> -m <markers>` | Literature validation via PubMed |
| `celltypepilot observe -o <dir>` | Read-only run lifecycle + host metrics (`--json`) |
| `celltypepilot host-acceptance` | Codex/Claude/MCP discovery + lifecycle discrimination harness |
| `celltypepilot evidence propose-promotion/review-promotion/apply-promotion` | Two-reviewer, versioned marker-edge evidence promotion |
| `celltypepilot pack install/list/validate/remove` | Manage data-only evidence packs without executable code |

### Agent lifecycle discrimination (host acceptance)

When orchestrating `doctor` / `inspect` / `benchmark-run`, **do not** treat a zero exit
code as scientific success. Discriminate durable states from checkpoints and release
readiness:

| State | Meaning |
|-------|---------|
| `running` | Fold/method still executing — no final metrics |
| `completed` | Unit finished — not automatically a public claim |
| `failed` | Started then errored/timeout — keep as negative result |
| `unavailable` | Dependency/adapter missing — never impute predictions |
| `claim_ready` | Only when release gate is green |
| `incomplete_not_claim_ready` | Protocol blocked — no public robustness claim |

Run `celltypepilot host-acceptance --skip-worktree` (or
`python scripts/run_host_acceptance_worktree.py`) before trusting a new host
integration. See `docs/host_acceptance.md`.

All commands support `--json` for structured output that the host integration can parse and present.

## Output files

| File | Description |
|---|---|
| `data.annotated.h5ad` | AnnData with separate identity/state candidates, decisions, evidence, CL ID, and confidence in obs |
| `evidence_table.csv` | Per-cluster evidence: scores, markers, critic flags, confidence |
| `contrastive_evidence.csv` | Existing top-two ranking with shared and candidate-specific support, missing/silent markers, conflicts, and provenance; no reranking |
| `evidence_gaps.json` | Per-Unknown observed evidence gaps with bounded next actions, prohibited actions, and human-review requirement |
| `state_results.csv` | Independent state candidate, decision, missing/silent markers, score, and evidence |
| `identity_contract.json` | Gene-symbol, active-scope, alias/CL, and pack identity audit |
| `context_pack.normalized.json` | Normalized governed context and hashes, when context is supplied |
| `figures/umap_cluster.png` | UMAP colored by cluster |
| `figures/umap_celltype.png` | UMAP colored by annotated cell type |
| `figures/umap_confidence.png` | UMAP colored by critic confidence level |
| `figures/marker_dotplot.png` | Marker gene expression dotplot |
| `figures/confidence_distribution.png` | Bar chart of confidence levels |
| `report_draft.html` | Comprehensive HTML report with all figures and evidence |
| `methodology_draft.txt` | Draft methods paragraph for papers |
| `manifest.json` | Run manifest: versions, parameters, data hash, output hashes |

## Marker Knowledge Graph (MKG)

The built-in atlas (`mkg-2026.08`) covers:
- **Blood/PBMC**: T cells (CD4/CD8/naive/memory/Treg/Th1/Th17), B cells (naive/memory/plasma), NK cells, monocytes (classical/non-classical), dendritic cells (cDC1/cDC2/pDC), platelets
- **Lung**: alveolar macrophages, AT1/AT2, ciliated, club, goblet, basal cells
- **Liver**: hepatocytes, Kupffer cells, cholangiocytes, endothelial, stellate cells
- **Brain**: excitatory/inhibitory neurons, astrocytes, oligodendrocytes, OPCs, microglia
- **Kidney**: proximal tubule, loop of Henle, DCT, podocytes, collecting duct
- **Gut**: enterocytes, goblet, enteroendocrine, Paneth, intestinal stem cells
- **Skin**: keratinocytes, melanocytes, Langerhans cells, fibroblasts
- **Heart**: cardiomyocytes, fibroblasts, endothelial, smooth muscle
- **Pancreas**: alpha/beta/delta cells, ductal, acinar cells
- **Skeletal muscle**: myofibers, satellite cells, FAPs
- **General**: endothelial, pericytes, fibroblasts, macrophages, mast cells, epithelial

Each entry includes: Cell Ontology ID, positive markers, negative markers, synonyms, and subtypes.
Mouse gene symbols are auto-converted from human conventions.

## Critic flags explained

| Flag | Meaning | Action |
|---|---|---|
| `PASS` | All checks passed | Accept annotation |
| `LOW_EVIDENCE` | <20% of all expected markers expressed | Abstain as Unknown; manual review needed |
| `PARTIAL_EVIDENCE` | 20-50% of all expected markers expressed | Abstain as Unknown pending review |
| `NEG_MARKER_CONFLICT` | Negative markers unexpectedly expressed | Likely misannotation or doublet |
| `POSSIBLE_DOUBLET` | Two lineage signatures co-expressed | Sub-cluster or mark as doublet |
| `NO_CL_ID` | No legal Cell Ontology ID assigned | Abstain until a versioned mapping is supplied |
| `UNREVIEWED_CONTEXT_ONLY` | Identity depends only on a draft custom panel | Keep Unknown; review the panel |
| `REVIEWED_CONTEXT_SUPPORT` | Identity depends only on a reviewed custom panel | Accepted only if other gates pass; confidence capped at medium |

## Confidence levels

| Level | Meaning | When |
|---|---|---|
| `high` | Strong evidence, no conflicts | ≥70% markers detected, no neg conflicts, high specificity |
| `medium` | Good evidence, minor concerns | 50-70% markers, mild neg conflicts |
| `low` | Weak or ambiguous evidence | 30-50% markers, some concerns |
| `needs_review` | Insufficient or conflicting evidence | <30% markers, major conflicts, possible doublet |

## Reference map

| File | When to read |
|---|---|
| `reference/methodology.md` | The algorithm details: scoring formula, critic design, confidence calibration |
| `reference/outputs.md` | Detailed output file schemas and how to use them |
| `reference/extending.md` | How to add custom markers, new tissues, or extend the atlas |

## Installation

```bash
# As a Claude Code plugin (recommended)
git clone https://github.com/HERRY423/CellTypePilot ~/.claude/plugins/marketplaces/local/plugins/celltypepilot
cd ~/.claude/plugins/marketplaces/local/plugins/celltypepilot
pip install -e .

# Verify installation
celltypepilot doctor
```

## Design philosophy

CellTypePilot is **not** another cell-type annotation algorithm. It is a **trust layer** that
sits on top of your existing coding workspace. The key insight: for individual researchers
and small labs, the bottleneck is not "getting labels" — it is **knowing whether to trust them**.
CellTypePilot's Annotation Critic is designed to be the skeptical colleague who asks "but why
do you think that cluster is a T cell? Show me the evidence."
