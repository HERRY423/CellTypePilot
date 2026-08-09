# Batch: `gpu_popv_retrain_v1`

Frozen **GPU** re-execution track for popV on Travaglini lung Smart-seq2.

## Why a new batch?

The existing tree `runs/travaglini_lung_smartseq2_2020` holds **CPU** three-fold artifacts.
GPU results must **not** be mixed into that tree. This batch re-runs **all three** donor
folds for popV under:

```text
batches/gpu_popv_retrain_v1/run/checkpoints/   # workers write atomic files only
batches/gpu_popv_retrain_v1/aggregate/         # read-only merge of checkpoints
```

Holdout plan / cluster map / label map are **read** from the locked public assets; GPU
workers never rewrite the CPU run's predictions.

## Required folds

1. `donor=travaglini_2020::1`
2. `donor=travaglini_2020::2`
3. `donor=travaglini_2020::3`

## Workflow

1. Bootstrap Linux + NVIDIA CTK (`scripts/bootstrap_gpu_linux_runner.sh`)
2. Build & pin GPU image (`scripts/build_gpu_popv_image.sh`)
3. Run workers (one fold per node, or all folds on one node)
4. Aggregate checkpoints read-only
5. Evaluate / release as a **separate** track from CPU
