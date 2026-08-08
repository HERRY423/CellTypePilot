# Locked comparator benchmark protocol

CellTypePilot is evaluated as a plugin, not as an autonomous agent. The benchmark runner
executes deterministic annotation tools behind one leakage-resistant file contract.

## Independence rules

1. Lock study/donor assignments before executing a comparator.
2. Each fold writes `train_reference.h5ad` and `test_query.h5ad`.
3. Test truth and label-like `obs` columns are removed from `test_query.h5ad`.
4. A fold-trained method may read labels only from `train_reference.h5ad`.
5. Every output must contain exactly one row for every test cell:
   `cell_id,predicted_label,confidence`.
6. Confidence is constrained to `[0,1]`; its semantics must be recorded by the adapter.
7. Label harmonization uses a CSV locked before result inspection:
   `method,raw_label,canonical_label`.
8. Missing or failed tools are reported as unavailable; their predictions are never imputed.

`celltypepilot benchmark-run` directly executes fold-trained CellTypePilot and CellTypist.
SingleR, Azimuth, and popV are invoked through shell-free JSON argv adapters because their
runtime/reference requirements are ecosystem-specific. Example:

```json
{
  "method": "singler",
  "reference_policy": "fold_train_only",
  "confidence_semantics": "adapter-specific value in [0,1]; describe the transformation",
  "argv": [
    "Rscript",
    "run_singler.R",
    "{train_h5ad}",
    "{test_h5ad}",
    "{output_csv}"
  ],
  "version_command": ["Rscript", "-e", "packageVersion('SingleR')"],
  "timeout_seconds": 3600
}
```

Placeholders must be separate argv entries. Commands run without a shell. Azimuth adapters
must construct an Azimuth reference from fold-train data or declare an external frozen
reference; a hosted reference must not be described as study-independent unless its training
studies have been audited against the test set.

## Reported metrics

- accuracy, macro-F1, balanced accuracy;
- abstain rate, coverage, selective accuracy;
- top-label correctness Brier score, ECE;
- risk-coverage curve and AURC;
- per-fold results and comparator execution status.

An OOF test benchmark estimates performance. It must not also be used to fit the released
abstention threshold. Thresholds are fitted with `celltypepilot calibrate` on a separately
designated calibration dataset.
