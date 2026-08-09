# Immutable asset repository

This directory is the **catalog and local object-cache root** for benchmark assets that
must not be mixed with live fold workspaces.

| Path | Role |
|------|------|
| `catalog.json` | Immutable asset index (CELLxGENE, Azimuth refs, label maps, Docker images) |
| `storage_policy.json` | Object store / CDN key templates and immutability rules |
| `objects/` | Content-addressed local cache (`v1/{kind}/{asset_id}/{version}/sha256/{sha256}`) — gitignored |

## Hard rules

1. **Never write under `benchmarks/**/runs/**`.** Active fold checkpoints, predictions, and
   logs are not asset-store targets.
2. **Content-addressed keys only.** Object keys include the SHA-256. Overwrite is denied.
3. **Every asset records:** `url`, `version`, `sha256`, `license`, `species`, `tissue`,
   `training_study_provenance`, and `availability`.
4. **Availability is explicit.** Missing CDN bytes are `pending_upload` or
   `blocked_overlap_audit`, never silently treated as ready.
5. **Azimuth references** stay `blocked_overlap_audit` until `ref.Rds` + `idx.annoy` digests
   and a training-study overlap audit are frozen.

## Object store / CDN layout

```text
s3://celltypepilot-public-assets/v1/{kind}/{asset_id}/{version}/sha256/{sha256}
https://cdn.celltypepilot.example/assets/v1/{kind}/{asset_id}/{version}/sha256/{sha256}
```

Templates live in `storage_policy.json`. Replace the example CDN host when a real bucket
is provisioned; do not rewrite historical catalog rows—publish a new catalog version.

## Asset kinds

| kind | Purpose |
|------|---------|
| `cellxgene_dataset` | Immutable CELLxGENE H5AD dataset versions |
| `azimuth_reference` | Frozen `ref.Rds` + `idx.annoy` bundles with provenance audit |
| `label_map` | Predeclared evaluation label maps (SHA-locked) |
| `docker_image` | Dockerfile and/or registry image digests for comparator runtimes |

## Commands

```bash
# Schema + local verification (read-only w.r.t. runs/)
python scripts/sync_asset_catalog.py

# Materialize file: sources into objects/ (label maps, Dockerfile recipe)
python scripts/sync_asset_catalog.py --materialize

# Machine-readable summary
python scripts/sync_asset_catalog.py --json

# CLI
celltypepilot assets list
celltypepilot assets verify
celltypepilot assets materialize
```

## Relationship to `public_v1/registry.json`

The public cohort registry remains the **analysis-plan and cohort lock**. This catalog is
the **byte-level asset index and CDN contract**. They link via `related_cohort_ids` and
shared SHA-256 digests; neither module rewrites the other's files during a running fold.
