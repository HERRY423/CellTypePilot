"""Execute fold-local popV inside the frozen GPU container (NVIDIA Container Toolkit).

Requires --gpus and the GPU image pin. Must not be pointed at the CPU three-fold
run tree; use batches/gpu_popv_retrain_v1 instead.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

IMAGE = os.environ.get("CTP_POPV_GPU_IMAGE", "celltypepilot-popv-gpu:0.6.1-cu124")


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: run_popv_docker_gpu.py TRAIN_H5AD TEST_H5AD OUTPUT_CSV")
    train, test, output = (Path(value).resolve() for value in sys.argv[1:])
    if train.parent != test.parent or train.parent != output.parent:
        raise ValueError("popV GPU Docker adapter requires fold-local input and output paths")
    # Hard guard: refuse writing into the known CPU Smart-seq2 run root.
    parts = {part.lower() for part in train.parts}
    if "travaglini_lung_smartseq2_2020" in parts and "batches" not in parts:
        raise SystemExit(
            "Refusing GPU popV under the CPU run tree travaglini_lung_smartseq2_2020. "
            "Use benchmarks/public_v1/batches/gpu_popv_retrain_v1 instead."
        )
    adapter_dir = Path(__file__).resolve().parent
    command = [
        "docker",
        "run",
        "--rm",
        "--gpus",
        "all",
        "--mount",
        f"type=bind,src={train.parent},dst=/work",
        "--mount",
        f"type=bind,src={adapter_dir},dst=/adapter,readonly",
        IMAGE,
        "python",
        "/adapter/run_popv.py",
        f"/work/{train.name}",
        f"/work/{test.name}",
        f"/work/{output.name}",
    ]
    completed = subprocess.run(command, check=False)
    raise SystemExit(completed.returncode)


if __name__ == "__main__":
    main()
