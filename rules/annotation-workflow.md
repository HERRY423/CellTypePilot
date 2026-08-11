# Annotation Workflow Rules

## When the user mentions single-cell annotation

Automatically route to the appropriate CellTypePilot command:

| User says | Action |
|---|---|
| "annotate my clusters" / "label my data" | `/annotate` |
| "what cell types are these?" | `/annotate` (after `/ctp-inspect`) |
| "check this cluster" / "why is this a T cell?" | `/critic` |
| "what's in my h5ad?" / "inspect my data" | `/ctp-inspect` |
| "is everything installed?" / "doctor" | `/ctp-doctor` |
| "review flagged clusters" | `/critic` on each flagged cluster |

## Mandatory workflow order

1. **ALWAYS** run `/ctp-doctor` first in a new session
2. **ALWAYS** run `/ctp-inspect` before `/annotate`
3. **ALWAYS** confirm species/tissue/cluster-key with the user before annotating
4. **NEVER** skip the critic review for `needs_review` clusters

## Evidence-first principle

- No cell-type label without evidence
- Every annotation must have a confidence level
- Flagged clusters are successes, not failures
- Show the user the evidence table, not just the labels

## Cost control

- Marker scoring is deterministic and free — always run it first
- LLM-based reasoning is optional — only use when the user asks or when the critic is uncertain
- Prefer reading files over re-running computations
