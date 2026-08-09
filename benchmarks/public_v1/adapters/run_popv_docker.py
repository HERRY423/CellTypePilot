"""Execute the fold-local popV adapter in the frozen Linux container."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


IMAGE = "celltypepilot-popv:0.6.1"


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit("usage: run_popv_docker.py TRAIN_H5AD TEST_H5AD OUTPUT_CSV")
    train, test, output = (Path(value).resolve() for value in sys.argv[1:])
    if train.parent != test.parent or train.parent != output.parent:
        raise ValueError("popV Docker adapter requires fold-local input and output paths")
    adapter_dir = Path(__file__).resolve().parent
    command = [
        "docker",
        "run",
        "--rm",
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
