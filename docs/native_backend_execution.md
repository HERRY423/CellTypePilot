# Native backend execution and depth-domain evidence pipeline

`annotate --native-backends` executes candidate generators inside the ordinary annotation
workflow. The configuration and each local reference/model are hashed. Every backend writes a
checkpoint, normalized candidates, version metadata, and a structured status. An unavailable or
failed backend remains visible and contributes no vote; it never causes marker scoring to become a
classifier.

## Configuration

```json
{
  "schema_version": "celltypepilot.native-backends.v1",
  "continue_on_failure": true,
  "resume": true,
  "backends": [
    {
      "backend": "celltypist",
      "mode": "retrain",
      "reference_path": "references/domain_reference.h5ad",
      "label_key": "cell_type"
    },
    {
      "backend": "popv",
      "mode": "retrain",
      "reference_path": "references/domain_reference.h5ad",
      "label_key": "cell_type",
      "query_counts_layer": "counts",
      "ref_counts_layer": "counts"
    },
    {
      "backend": "singler",
      "reference_path": "references/domain_reference.h5ad",
      "label_key": "cell_type",
      "rscript": "Rscript"
    },
    {
      "backend": "scanvi",
      "reference_path": "references/domain_reference.h5ad",
      "label_key": "cell_type",
      "query_counts_layer": "counts",
      "ref_counts_layer": "counts"
    },
    {
      "backend": "custom_reference",
      "method": "correlation",
      "reference_path": "references/domain_reference.h5ad",
      "label_key": "cell_type"
    },
    {
      "backend": "llm",
      "provider": "openai",
      "model": "YOUR_STRUCTURED_OUTPUT_MODEL",
      "allow_network": true,
      "api_key_env": "OPENAI_API_KEY"
    }
  ]
}
```

```bash
celltypepilot annotate -i query.h5ad -k leiden -s human -t lung \
  --native-backends native_backends.json
```

popV and scANVI require raw integer counts from the configured layer or `AnnData.raw`. SingleR
runs through a packaged `Rscript` adapter with shell-free argv. CellTypist can use a governed
pretrained model or train only on the configured reference. Query observation metadata is stripped
to the cluster key and explicitly needed batch keys before external execution.

The LLM runner requires both a config entry and `allow_network=true`. It receives only a bounded
per-cluster packet made from the top marker-evidence candidates and allowed labels, uses strict
JSON-schema output, and is always normalized as `hypothesis_only`. It cannot satisfy the
hierarchical selector's independent-backend requirement.

## Artifacts and resume contract

```text
output/native_backends/
  native_backend_status.csv
  <backend>/checkpoint.json
  <backend>/candidates.normalized.csv
  <backend>/run.json
  <backend>/stdout.log or stderr.log
```

The ordinary `manifest.json` records the config hash, dependency hashes, backend status, versions,
and bounded metadata. A completed checkpoint is reused only when the input, cluster key, backend
entry, and every local dependency hash still match.

## Three depth-validation domains

The evidence workflow is limited to lung, gut/IBD, and tumor microenvironment. It uses the same
native fold runners behind the leakage-resistant benchmark materializer:

```bash
celltypepilot domain-validation-plan \
  --registry benchmarks/public_v1/registry.json \
  --output benchmarks/domain_depth_v1

celltypepilot domain-validation-run \
  --plan benchmarks/domain_depth_v1/domain_validation_plan.json \
  --domain lung
```

Planning audits local assets, donor/study/platform depth, and locked label-map coverage. For
scANVI and custom-reference fold training, the emitted label universe is exactly the fold-training
truth label universe; the planner copies the already locked `__truth__` mapping into a new hashed
label-map artifact before execution. It does not infer new biological equivalences.

Execution uses atomic method/fold checkpoints, strips test truth before every backend call, and
keeps cell- and cluster-level endpoints separate. popV, SingleR, scANVI, custom reference, and the
hierarchical CellTypePilot product share the package-native execution code. LLM is excluded from
the accuracy comparator family. When the product and individual comparator are requested in the
same fold, the comparator reuses the product's hashed raw candidate artifact instead of retraining
the same backend twice.

`completed` means only that every requested method/fold artifact exists. Domain `claim_ready`
remains false until minimum independent cohorts/studies/donors/platforms, a separate calibration
cohort, external holdout, expert adjudication, and domain-specific stress tests are complete.
