# Identity × State Lens

CellTypePilot treats lineage identity and transient cell state as independent axes. The State Lens
adds exploratory state evidence after the identity Critic has made its decision.

## Safety invariant

Attaching state results must not change these identity fields:

- `cell_type`
- `cl_id`
- `candidate_cell_type`
- `decision`
- `abstain_reason`

The implementation checks this invariant and fails the pipeline if it is violated. A supported
state may therefore produce a display label such as `Unknown · interferon_responsive`, but the
canonical identity remains `Unknown` and the identity decision remains `abstain`.

## Decisions

| State decision | Meaning | Canonical state output |
|---|---|---|
| `supported` | Sufficient expression and positive DE evidence, with no negative-marker conflict | State label |
| `hypothesis` | Some evidence, or a draft context-only state | `Unknown` |
| `abstain` | Insufficient, conflicting, or unavailable evidence | `Unknown` |

State scoring reports the full expected-marker denominator and separately records markers that are
missing from the matrix and markers that are present but silent. It applies the same positive
direction, log2FC, BH-FDR, and expression-fraction concepts used by identity marker scoring.

## AnnData fields

| `.obs` field | Description |
|---|---|
| `ctp_cell_state_candidate` | Best state hypothesis, including non-supported candidates |
| `ctp_state_decision` | `supported`, `hypothesis`, or `abstain` |
| `ctp_cell_state` | Supported state, otherwise `Unknown` |
| `ctp_state_score` | Heuristic state evidence score |
| `ctp_state_confidence` | Heuristic confidence category |
| `ctp_state_evidence` | Compact evidence summary |
| `ctp_display_label` | Identity plus supported state for visualization only |

Per-cluster state rows are also written to `state_results.csv`. Use `--no-states` to disable this
axis; identity scoring and critic behavior remain unchanged.

## Evidence boundary

The bundled `state-atlas-2026.08` contains cell-cycle, interferon-response, hypoxia, inflammatory,
and exhaustion modules with source identifiers and PMID/DOI metadata. Its relations are presently
marked `aggregate_source_only_not_edge_verified`. These modules support exploratory review; they do
not establish calibrated accuracy, comprehensive biological coverage, or primary verification of
each marker-state edge.
