# Public multi-cohort benchmark release

This workflow is the evidence boundary between “the plugin executes” and a public
robustness claim. A release is claim-ready only when every immutable cohort, locked split,
predeclared label map, required comparator, and diagnostic artifact is present.

## Immutable asset catalog (object store / CDN)

Byte-level assets (CELLxGENE H5AD versions, Azimuth reference slots, label maps, Docker
images) are indexed under `benchmarks/assets/` with URL, version, SHA-256, license,
training-study provenance, tissue/species, and availability. See
`benchmarks/assets/README.md`, `catalog.json`, and `storage_policy.json`.

The asset catalog is parallel to `registry.json`. Sync and materialize tools write only
under `benchmarks/assets/objects/` and **never** under `benchmarks/**/runs/**` (active fold
workspaces).

```bash
python scripts/sync_asset_catalog.py
python scripts/sync_asset_catalog.py --materialize
celltypepilot assets list
celltypepilot assets verify
```

## Frozen v1 cohort set

The machine-readable registry is `benchmarks/public_v1/registry.json`. It locks exact
CELLxGENE dataset-version URLs and byte sizes for:

1. Garrido-Trigo IBD/healthy colon mucosa — inflamed-tissue stress test;
2. Pelka colorectal-cancer stromal cells — tumor and stromal-state stress test;
3. Travaglini lung 10x and Smart-seq2 assets — shared-study, cross-platform strata;
4. the Smillie portal subset — retained only as healthy gut because the downloaded asset
   contains 12 normal donors despite the UC collection title.

The last item is a deliberate negative result. A collection title is not outcome metadata.
The registry excludes this asset from inflamed-tissue claims.

## Reproduce without truth leakage

Install the benchmark dependencies and download the immutable source assets. Source H5AD
files are never rewritten. The registry supplies the versioned download URLs and expected
byte sizes; the release inventory additionally verifies SHA-256 hashes.

```bash
pip install -e ".[benchmark]"
python scripts/download_public_benchmark.py \
  --registry benchmarks/public_v1/registry.json
python scripts/prepare_public_benchmark.py \
  --registry benchmarks/public_v1/registry.json
python scripts/lock_public_holdouts.py \
  --registry benchmarks/public_v1/registry.json
python scripts/build_public_label_maps.py \
  --registry benchmarks/public_v1/registry.json
python scripts/preflight_public_comparators.py \
  --registry benchmarks/public_v1/registry.json
```

Preparation creates one cluster map per cohort. It freezes the 1,000 genes with highest
cohort-wide expression variance without reading labels, performs randomized TruncatedSVD,
and applies a deterministic KMeans rule separately inside each donor. `cell_type`, author
annotations, disease, condition, and other outcome fields are not used. The clustering audit
records the source hash, parameters, seed, and per-donor cluster counts. The global,
truth-blind feature-selection step is transductive and is declared as such; cluster boundaries
never pool donors. This preprocessing is a reproducible stress-test convention, not evidence
that KMeans is the preferred biological clustering method.

Run each cohort with leave-one-donor-out folds. For example:

```bash
celltypepilot benchmark-run \
  --input benchmarks/public_v1/data/garrido_trigo_ibd.h5ad \
  --truth-key cell_type --study-key ctp_study_id \
  --constant-study-id garrido_trigo_2023 --donor-key donor_id \
  --cluster-key ctp_cluster \
  --cluster-map benchmarks/public_v1/runs/garrido_trigo_ibd_2023/cluster_map.csv \
  --species human --tissue gut --strategy donor \
  --methods celltypepilot,celltypist,singler,azimuth,popv \
  --comparator-config benchmarks/public_v1/adapters/singler.json \
  --comparator-config benchmarks/public_v1/adapters/popv.json \
  --label-map benchmarks/public_v1/label_maps/garrido_trigo_ibd_2023.csv \
  --output benchmarks/public_v1/runs/garrido_trigo_ibd_2023
```

For immutable sources without a study column, the CLI injects the registry's constant study
ID in memory. Do not reuse donor IDs across studies unless they identify the same person.

## Comparator policy

- CellTypePilot and CellTypist are trained or configured separately inside every fold.
- SingleR and popV have shell-free fold-train-only adapter contracts. They count as completed
  only after the local execution status records every fold; adapter code alone is not evidence.
- Azimuth is unavailable until a compatible frozen reference (`ref.Rds`, `idx.annoy`) has a
  documented training-study overlap audit. No generic result is fabricated.
- Missing or failed tools are negative results. Predictions are never imputed.
- Confidence values retain method-specific semantics and are not called probabilities.

### GPU popV track (separate frozen batch)

CPU Docker popV (`celltypepilot-popv:0.6.1`) and any partial three-fold CPU run under
`runs/travaglini_lung_smartseq2_2020` are **not** mixed with GPU results. GPU re-execution
is a new batch:

- Pins: `benchmarks/public_v1/environment/gpu/stack_pins.json`
- Batch: `benchmarks/public_v1/batches/gpu_popv_retrain_v1/`
- Adapter: `adapters/popv_gpu.json` (requires NVIDIA Container Toolkit + `--gpus all`)
- Workers write only atomic `run/checkpoints/*`; aggregators are read-only over checkpoints
  and write only under `aggregate/`
- All three Travaglini Smart-seq2 donor folds must be re-run for claim-ready GPU popV

```bash
sudo bash scripts/bootstrap_gpu_linux_runner.sh
bash scripts/build_gpu_popv_image.sh
python scripts/run_gpu_fold_worker.py --batch-root benchmarks/public_v1/batches/gpu_popv_retrain_v1 --all-folds --write-aggregate
# or multi-node: --fold-id 'donor=travaglini_2020::1' (no --write-aggregate)
python scripts/aggregate_gpu_batch_checkpoints.py --batch-root benchmarks/public_v1/batches/gpu_popv_retrain_v1 --require-complete
```

## Build the release

```bash
celltypepilot benchmark-release \
  --registry benchmarks/public_v1/registry.json \
  --output benchmarks/public_v1/release \
  --n-boot 2000 --seed 42
```

The release contains cohort- and donor-level results, equal-donor bootstrap summaries,
paired donor comparisons with BH correction, batch/platform/condition sensitivity,
sample-enrichment diagnostics, predeclared low-quality/doublet/ambient strata,
`negative_results.csv`, a readable report, and a hashed `release_manifest.json`.

## Statistical interpretation

Cells from one donor are not independent replicates. Primary estimates therefore average
donor-level macro-F1 equally. Confidence intervals resample studies and then donors (or donors
only in a single-study track). Method comparisons use only donors with paired results.

Batch/platform/condition ranges and sample enrichment are descriptive diagnostics. A large
range identifies sensitivity; it does not identify a causal batch effect. Missing doublet or
ambient-RNA metadata is `not_assessed`, never “no artifact detected.” Expert-curated
`cell_type` values are reference labels rather than infallible biological truth. Label maps
are frozen before prediction scoring. Negative, non-significant, and failed results remain in
the release.
