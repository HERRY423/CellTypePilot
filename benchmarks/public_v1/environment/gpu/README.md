# GPU Linux runner (NVIDIA Container Toolkit)

This directory freezes the **GPU** popV stack for a **new** public-benchmark batch.
It is intentionally separate from the CPU image `celltypepilot-popv:0.6.1` and from
`runs/travaglini_lung_smartseq2_2020` (CPU three-fold artifacts).

## Policy

| Rule | Detail |
|------|--------|
| New frozen batch | `batches/gpu_popv_retrain_v1/` |
| Re-run scope | **All three** Travaglini Smart-seq2 donor folds for popV |
| No mixing | Do not copy or merge GPU predictions into the CPU checkpoint tree |
| Worker writes | Atomic `checkpoints/{fold}__popv.{status.json,csv}` only when distributed |
| Aggregator | Read-only over checkpoints; writes only under the GPU batch `aggregate/` |

## Pin set

See `stack_pins.json`:

- NVIDIA driver ≥ 550.54.15 (recommended 550.90.07)
- CUDA 12.4 container base
- NVIDIA Container Toolkit ≥ 1.15
- Image tag `celltypepilot-popv-gpu:0.6.1-cu124`
- Digests filled by `scripts/freeze_gpu_image_identity.py` after build on the Linux node

## Bootstrap a Linux node

```bash
# On Ubuntu 22.04 x86_64 with root/sudo
sudo bash scripts/bootstrap_gpu_linux_runner.sh

# Build GPU image and record digests
bash scripts/build_gpu_popv_image.sh

# Smoke GPU visibility
docker run --rm --gpus all celltypepilot-popv-gpu:0.6.1-cu124 \
  python -c "import torch; assert torch.cuda.is_available(); print(torch.cuda.get_device_name(0))"
```

## Distributed donor-fold execution

```bash
# Node A — fold 1 only (no global OOF rewrite)
python scripts/run_gpu_fold_worker.py \
  --batch-root benchmarks/public_v1/batches/gpu_popv_retrain_v1 \
  --fold-id 'donor=travaglini_2020::1' \
  --worker-id gpu-node-a

# Node B / C similarly for folds 2 and 3

# Aggregator (read-only checkpoints → batch aggregate tables)
python scripts/aggregate_gpu_batch_checkpoints.py \
  --batch-root benchmarks/public_v1/batches/gpu_popv_retrain_v1
```

## Single-node re-run of all three folds

```bash
python scripts/run_gpu_fold_worker.py \
  --batch-root benchmarks/public_v1/batches/gpu_popv_retrain_v1 \
  --all-folds \
  --worker-id gpu-node-1 \
  --write-aggregate
```
