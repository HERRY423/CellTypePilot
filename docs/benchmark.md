# Independent holdout benchmark

CellTypePilot's benchmark path separates split design from prediction generation. It does
not treat a technical run on the same donors as evidence of generalization.

## 1. Lock the split plan

```bash
celltypepilot benchmark \
  --input benchmark.h5ad \
  --truth-key ground_truth \
  --study-key study_id \
  --donor-key donor_id \
  --strategy study \
  --output benchmark/
```

`study` creates leave-one-study-out test folds. `donor` creates leave-one-donor-unit-out
folds. Donor IDs must be globally unique: an ID appearing in more than one study is rejected
rather than silently treated as independent. Missing metadata, duplicate cell IDs, and fold
mismatches also fail closed. The assignment CSV hash is locked in the benchmark manifest.

## 2. Generate predictions without outcome leakage

For each fold, fit or configure CellTypePilot, CellTypist, SingleR, Azimuth, and popV using
training studies only. Do not inspect held-out labels while selecting markers, thresholds,
references, or label mappings. Store a long CSV:

```text
cell_id,fold_id,method,predicted_label
cell-1,study=study_A,celltypepilot,T cell
cell-1,study=study_A,celltypist,T cell
```

The harness intentionally does not pretend that an unavailable R/Python comparator ran.
Methods without supplied out-of-fold predictions are reported as `not_provided`.

## 3. Evaluate

```bash
celltypepilot benchmark \
  --input benchmark.h5ad \
  --truth-key ground_truth \
  --study-key study_id \
  --donor-key donor_id \
  --strategy study \
  --predictions predictions.csv \
  --output benchmark/
```

Outputs include aggregate and per-fold accuracy, macro-F1, balanced accuracy, coverage,
abstain rate, and selective accuracy. Report both coverage and selective accuracy so an
abstaining system cannot appear superior merely by declining difficult cells.
