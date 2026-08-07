# CellTypePilot — Codex subdirectory instructions

> These instructions apply when working within the CellTypePilot project.

## Development commands

```bash
# Install in development mode
pip install -e ".[dev]"

# Run all tests
python -m pytest tests/ -v

# Run specific test class
python -m pytest tests/test_smoke.py::TestCritic -v

# Check code quality
python -m pytest tests/ --tb=short
```

## Architecture notes

CellTypePilot is a dual-platform skill:
- **Claude Code**: Entry via `skills/celltypepilot/SKILL.md` (YAML frontmatter + orchestration)
- **Codex CLI**: Entry via root `AGENTS.md` (pure markdown instructions)
- **Shared backend**: `src/celltypepilot/` Python package (CLI + compute modules)

Both platforms call the same `celltypepilot` CLI commands. The difference is only in how
the agent discovers and orchestrates the workflow.

## Module responsibilities

| Module | Responsibility | Key functions |
|---|---|---|
| `data_adapter.py` | Load h5ad, detect species/tissue, find keys | `load_h5ad()`, `inspect_adata()`, `load_marker_atlas()` |
| `marker_scorer.py` | DE analysis + marker overlap scoring | `compute_marker_scores()`, `generate_annotation_summary()` |
| `critic.py` | Independent evidence review | `run_critic()`, `generate_critic_summary()` |
| `visualizer.py` | UMAP, dotplot, confidence figures | `generate_all_figures()` |
| `reporter.py` | HTML report + methodology text | `generate_html_report()`, `generate_methodology_text()` |
| `provenance.py` | manifest.json with versions/hashes | `create_manifest()`, `save_manifest()` |
| `doctor.py` | Environment dependency check | `run_doctor()`, `print_doctor()` |
| `cli.py` | Typer CLI with 5 commands | `doctor`, `inspect`, `annotate`, `critic`, `markers` |

## Data flow

```
h5ad → data_adapter.inspect() → parameters confirmed
     → data_adapter.load_marker_atlas(species, tissue) → marker knowledge graph
     → marker_scorer.compute_marker_scores() → candidate annotations
     → critic.run_critic() → reviewed annotations with flags
     → visualizer.generate_all_figures() → PNG figures
     → reporter.generate_html_report() → HTML report
     → provenance.save_manifest() → manifest.json
```

## Adding new tissues/markers

Edit `src/celltypepilot/data/marker_atlas.json`. Each tissue entry needs:
- `name`, `organ_system`, `cell_types` dict
- Each cell type: `cl_id`, `synonyms`, `positive_markers`, `negative_markers`, optional `subtypes`
- Mouse gene symbols are auto-converted from human conventions

## Testing

Tests use a synthetic PBMC dataset (`make_synthetic_pbmc()` in `tests/test_smoke.py`)
that generates 500 cells across 5 cell types with known markers. No external data needed.
