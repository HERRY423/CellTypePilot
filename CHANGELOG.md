# Changelog

All notable changes to CellTypePilot are documented here. The project follows
[Semantic Versioning](https://semver.org/). Release claims remain bounded by the
validation scope recorded in the repository and generated manifests.

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

[0.3.0]: https://github.com/HERRY423/CellTypePilot/releases/tag/v0.3.0
