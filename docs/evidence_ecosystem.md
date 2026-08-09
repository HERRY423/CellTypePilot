# CellTypePilot biological evidence layer

CellTypePilot treats biological evidence as a versioned plugin ecosystem, not
as a larger hard-coded marker dictionary. The host may discover and route
bounded operations, but it cannot promote evidence, invent identity mappings,
or turn a review signal into a biological claim.

## Priority 1 — gene and cell identity contract

Before scoring, CellTypePilot now audits marker addressability. It may use a
recognized gene-symbol column (for example `feature_name`) only when that
strictly improves overlap with the declared marker universe. Duplicate symbols
remain unresolved on their original feature identifiers. Original `var_names`
are restored before an annotated `.h5ad` is written or returned.

Runtime identity scope is composed from the requested tissue, `general`, and
only those extra scopes declared by a selected data-only pack. Labels resolve
through canonical names, synonyms, CL IDs, pack aliases, and explicit safe
parent fallbacks. Unknown or ambiguous labels remain unresolved.

Every annotation run writes `identity_contract.json` and records the same
contract in `manifest.json`. These artifacts establish representation and
scope; they do not establish biological accuracy.

## Priority 2 — first-party lung evidence pack

`lung_evidence_v0_1` is a first-party, data-only, hypothesis-trust pack. It
contains:

- a lung marker-atlas extension;
- CL-aware aliases and explicit safe-parent fallbacks;
- a reference manifest that declares no bundled reference and requires source
  disjointness for later benchmark evidence.

Its marker edges remain
`aggregate_source_only_not_edge_verified`. The pack is therefore useful for
coverage discovery and governed review, but is not a claim-ready source of
edge-level validation.

```text
celltypepilot pack validate lung_evidence_v0_1
celltypepilot annotate ... --pack lung_evidence_v0_1
```

## Priority 3 — human-gated evidence promotion

Automated search may draft a promotion proposal but cannot modify the atlas.
Database-record or primary-source promotion requires two distinct human
approvers; the requester cannot self-approve. The proposal locks the source
atlas version and marker-edge hash. Application fails if the edge changed and
always creates a new atlas version.

```text
celltypepilot evidence propose-promotion ...
celltypepilot evidence review-promotion --decision approve ...
celltypepilot evidence review-promotion --decision approve ...
celltypepilot evidence apply-promotion --new-version ...
```

## Priority 4 — separate cell and cluster benchmark tracks

CellTypePilot is a cluster annotator, so a cell-level comparator call is not a
fair substitute for a cluster endpoint. `benchmark` and `benchmark-run` now
support `--evaluation-unit cell|cluster|both`. The cluster track aggregates
truth and every method prediction to the same predeclared
`fold_id + cluster` unit, records truth purity and prediction support, and
abstains on ties. Cell and cluster results are written to different tables and
must not be blended into one performance claim.

Low-purity clusters remain visible in `cluster_track_diagnostics.csv`. A truth
tie is excluded with an explicit reason rather than silently tie-broken.

## Priority 5 — bounded Agent-host tools

The default local MCP facade exposes exactly four product operations:

- `prepare_annotation`;
- `annotate_from_plan`;
- `review_uncertain_clusters`;
- `finalize_reviewed_annotations`.

Each operation returns the same `celltypepilot.agent-decision.v1` envelope so
the host can distinguish blockers, evidence, allowed next actions, forbidden
claims, artifacts, and required human action without parsing prose. Annotation
writes two Agent-facing evidence products: `evidence_gaps.json` converts every
fail-closed `Unknown` into an observed, actionable gap without selecting a
replacement label; `contrastive_evidence.csv` explains the existing top-two
ranking through shared and candidate-specific marker support, missing/silent
markers, conflicts, and provenance. Neither artifact recalibrates or optimizes
the underlying score.

The diagnostic primitives below remain available only when maintainers opt in
with `CELLTYPEPILOT_MCP_SURFACE=advanced`:

- `tool_evidence_coverage`: runtime-eligible and candidate marker
  addressability, without annotation;
- `tool_evidence_trace`: CL/alias resolution and marker-edge provenance;
- `tool_resolve_evidence_packs`: compatible installed packs, with no automatic
  installation or trust upgrade;
- `tool_evidence_gap_queue`: read-only curation candidates from run artifacts;
- `tool_benchmark_card`: separate cell/cluster endpoints and claim readiness.

These tools are routing and artifact interfaces. They do not autonomously plan
experiments, install evidence, approve promotions, rename Novelty/OOD
candidates, or sign off biological conclusions.

## Claim boundary

Identifier reachability is not marker specificity. Aggregate database
provenance is not edge verification. Cluster purity is not annotation
accuracy. Benchmark tables are not a release claim unless the locked datasets,
comparators, label maps, calibration policy, and release manifest are complete.
Human review remains responsible for biological acceptance.
