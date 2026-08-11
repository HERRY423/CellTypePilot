# Changelog

All notable changes to CellTypePilot are documented here. The project follows
[Semantic Versioning](https://semver.org/). Release claims remain bounded by the
validation scope recorded in the repository and generated manifests.

## [0.4.0] - 2026-08-11

### Added

- Backend-neutral candidate artifacts and an ontology-aware selective decision layer for
  CellTypist, popV, SingleR, scANVI, custom references, and optional hypothesis-only LLM review.
- Native, resumable runners for those backends in ordinary `annotate`, with truth stripping and
  fold-local reuse in public benchmarks.
- Three depth-domain Atlas contracts for lung, gut/IBD, and tumor microenvironment, including
  governed scope packs and explicit evidence/claim boundaries.
- Outcome-blind donor-role locking for downgrade-only calibration and truth-free multi-lineage
  coverage audits.
- A hash-verified governance freeze spanning decision code, Atlas content, state/novelty policy,
  calibration policy, and domain contracts.

### Changed

- The marker scorer is evidence-only and no longer serves as the primary identity classifier when
  backend-neutral candidates are configured.
- Domain validation uses fold-isolated native backends, separate cell/cluster endpoints, atomic
  checkpoints, and retained unavailable/failed method statuses.

### Validation boundary

- The software release and lung multi-lineage addressability audit do not establish annotation
  accuracy, calibrated selective risk, independent-cohort calibration, or domain validation.
- Three-domain and five-cohort evidence claims remain fail-closed until their locked minimum cohort,
  runtime, calibration, expert-adjudication, and robustness requirements are satisfied.

## [0.3.1] - 2026-08-10

### Fixed

- Generate release checksums with a tested, cross-platform Python builder instead of shell
  redirections whose paths depended on the parent shell working directory.
- Fail closed unless the current version has exactly one wheel, one source distribution, and one
  Agent plugin bundle before `SHA256SUMS` is written.
- Rename the Claude Code `/doctor` slash command to `/ctp-doctor` so it no longer silently
  shadows Claude Code's built-in `/doctor` health checkup.
- Rename the Claude Code `/inspect` slash command to `/ctp-inspect` to avoid a future collision
  with a possible built-in `/inspect` and to stay distinct from the `claude inspect` CLI.

### Release note

- The immutable `v0.3.0` tag did not create a GitHub Release or publish a PyPI distribution because
  its verification job stopped before either publication stage.
- This patch changes release infrastructure and version metadata only; it adds no biological
  validation or annotation-accuracy claim.

## [0.3.0] - 2026-08-10

### Added

- Four-step Agent golden workflow with the `celltypepilot.agent-decision.v1` protocol.
- Actionable evidence gaps for `Unknown` decisions and contrastive top-two evidence.
- Independent Identity, State, and Novelty/OOD review axes with human sign-off boundaries.
- Governed context packs, data-only extension packs, atlas governance, and evidence promotion.
- Donor-aware benchmark infrastructure, comparator adapters, resumable execution, and QC diagnostics.
- Native MCP facade, Web Review audit artifacts, and dual Codex/Claude plugin manifests.
- Reproducible GitHub plugin bundle alongside the Python backend distribution.

### Changed

- The historical `premium` first-party atlas is explicitly MIT-licensed and available to all users.
- Python packaging now represents the deterministic backend; the complete Agent plugin is distributed
  as a separate GitHub Release bundle.
- Release automation verifies version alignment, tests, package metadata, and installed-wheel behavior
  before creating the GitHub Release and publishing to PyPI through Trusted Publishing.

### Validation boundary

- This is a technical preview release of an auditable annotation-review plugin.
- It does not establish biological superiority over CellTypist, SingleR, Azimuth, popV, or expert review.
- A qualified human remains responsible for final annotations and biological claims.

[0.4.0]: https://github.com/HERRY423/CellTypePilot/releases/tag/v0.4.0
[0.3.1]: https://github.com/HERRY423/CellTypePilot/releases/tag/v0.3.1
[0.3.0]: https://github.com/HERRY423/CellTypePilot/releases/tag/v0.3.0
