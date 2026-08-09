"""Freeze benchmark runtime versions and image identity into tracked artifacts."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path


def run(command: list[str], *, environment: dict[str, str] | None = None) -> str:
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=environment,
        check=False,
        timeout=300,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {command}\n{completed.stderr}"
        )
    return completed.stdout.strip() + "\n"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-root", required=True)
    args = parser.parse_args()
    root = Path(args.benchmark_root).resolve()
    output = root / "environment"
    output.mkdir(parents=True, exist_ok=True)

    python_executable = root / "envs" / "python" / "Scripts" / "python.exe"
    python_freeze = output / "python-celltypist.freeze.txt"
    atomic_write_text(
        python_freeze, run([str(python_executable), "-m", "pip", "freeze"])
    )

    rscript = Path("C:/R1/R-4.5.3/bin/Rscript.exe")
    r_library = root / "envs" / "R" / "library"
    r_environment = {**os.environ, "R_LIBS_USER": str(r_library)}
    r_packages = output / "r-singler.packages.csv"
    expression = (
        "p <- as.data.frame(installed.packages()[,c('Package','Version','LibPath')]); "
        "write.csv(p[order(p$Package),], stdout(), row.names=FALSE, quote=TRUE)"
    )
    atomic_write_text(
        r_packages,
        run([str(rscript), "-e", expression], environment=r_environment),
    )
    r_session = output / "r-session-info.txt"
    atomic_write_text(
        r_session,
        run([str(rscript), "-e", "sessionInfo()"], environment=r_environment),
    )

    image_record = output / "popv-image.json"
    popv_freeze = output / "python-popv.freeze.txt"
    try:
        inspect = run(
            [
                "docker",
                "image",
                "inspect",
                "celltypepilot-popv:0.6.1",
                "--format",
                "{{json .}}",
            ]
        )
        image_payload = json.loads(inspect)
        atomic_write_text(
            popv_freeze,
            run(
                [
                    "docker",
                    "run",
                    "--rm",
                    "--entrypoint",
                    "python",
                    "celltypepilot-popv:0.6.1",
                    "-m",
                    "pip",
                    "freeze",
                ]
            ),
        )
    except (RuntimeError, json.JSONDecodeError) as exc:
        image_payload = {"status": "unavailable", "error": str(exc)}
        atomic_write_text(popv_freeze, f"UNAVAILABLE: {exc}\n")
    atomic_write_text(image_record, json.dumps(image_payload, indent=2) + "\n")

    installation_attempts = output / "installation_attempts.jsonl"
    environment_readme = output / "README.md"
    artifacts = [
        python_freeze,
        r_packages,
        r_session,
        image_record,
        popv_freeze,
        installation_attempts,
        environment_readme,
    ]
    manifest = {
        "schema_version": "celltypepilot.benchmark-environment-freeze.v1",
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "artifacts": {path.name: digest(path) for path in artifacts},
    }
    atomic_write_text(
        output / "environment_manifest.json", json.dumps(manifest, indent=2) + "\n"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
