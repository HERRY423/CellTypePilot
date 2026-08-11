# Backend-neutral hierarchical selective decisions

CellTypePilot separates three responsibilities:

1. Candidate backends propose identity labels.
2. The hierarchical selector chooses a leaf, a governed ancestor, or abstention.
3. Marker evidence and the Annotation Critic can downgrade that choice, never rescue it.

This is a draft-annotation decision protocol. It is not a calibrated probability model and does
not establish biological accuracy.

## Candidate artifact contract

CSV and JSON artifacts use `celltypepilot.backend-candidates.v1`. A cluster-level artifact needs:

| Field | Requirement |
|---|---|
| `cluster` or `cluster_id` | Required |
| `backend` or `method` | Required: `celltypist`, `popv`, `singler`, `scanvi`, `custom_reference`, `knn`, `correlation`, or `llm` |
| `cell_type`, `candidate_cell_type`, `predicted_label`, `prediction`, or `label` | Required |
| `cl_id` / `ontology_id` | Strongly recommended; unresolved or ambiguous labels cannot be accepted |
| `score` and `score_semantics` | Optional; retained without cross-backend score fusion |
| `rank` | Optional; inferred deterministically when absent |

A cell-level artifact may provide `cell_id`, `barcode`, or `obs_name` instead of `cluster`. It must
contain exactly one top-1 prediction per cell/backend. CellTypePilot joins query `obs_names` to the
locked cluster assignments, emits within-backend cluster vote fractions, and labels those fractions
as descriptive rather than probabilistic.

The golden Agent workflow hashes every candidate artifact and optional decision policy during
`prepare_annotation`; `annotate_from_plan` rejects changed files.

## Backend roles and independence

| Backend | Default role | Independence group |
|---|---|---|
| CellTypist | decision candidate | reference linear |
| popV | decision candidate | reference ensemble |
| SingleR | decision candidate | reference correlation |
| scANVI | decision candidate | reference latent |
| custom reference | decision candidate | custom reference |
| KNN / correlation over the same custom-reference family | decision candidate | custom reference; counted once |
| optional LLM | hypothesis only | does not vote by default |
| CellTypePilot marker scorer | evidence only | never votes |

Running multiple wrappers from the same independence group cannot manufacture consensus.

## Default decision policy

The default `hierarchical-selective-default-v1` policy requires two independent backend groups.

- Exact canonical agreement accepts the leaf candidate for critic review.
- Sibling disagreement may accept the closest common Atlas ancestor within two parent hops.
- Cross-lineage disagreement, an unresolved label, or insufficient independent support produces
  `Unknown` while retaining the candidate set and explicit abstention reason.
- The downstream marker critic can downgrade an accepted leaf/ancestor; it cannot promote an
  abstention.

`backend_agreement_fraction` is top-1 agreement among independent groups. It is not confidence,
accuracy, posterior probability, or a selective-risk guarantee. A separately locked held-out
calibration artifact remains necessary for risk/coverage claims.

## Three depth-validation domains

The machine-readable registry is `src/celltypepilot/data/validation_domains.json`. Claim-building is
concentrated on:

- lung homeostasis and injury;
- gut and inflammatory bowel disease;
- tumor microenvironment.

Each domain requires locked multi-study and donor holdouts, at least two platform families,
candidate-backend comparisons, hierarchy-aware error, selective risk/coverage, abstention audits,
separate calibration, and expert adjudication. All three currently fail closed as
`evidence_required`. Other Atlas tissues remain exploratory scopes, not validated product domains.
