"""Execute or resume every frozen cohort using the isolated benchmark runtimes."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", required=True)
    parser.add_argument(
        "--methods", default="celltypepilot,celltypist,singler,popv"
    )
    parser.add_argument("--cohort", action="append", default=[])
    parser.add_argument("--continue-on-unavailable", action="store_true")
    args = parser.parse_args()

    registry_path = Path(args.registry).resolve()
    root = registry_path.parent
    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    executable = root / "envs" / "python" / "Scripts" / "celltypepilot.exe"
    methods = tuple(part.strip().lower() for part in args.methods.split(",") if part)
    selected = set(args.cohort)
    comparator_configs = {
        "singler": root / "adapters" / "singler.json",
        "popv": root / "adapters" / "popv.json",
    }
    logs = root / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    environment = {
        **os.environ,
        "R_LIBS_USER": str(root / "envs" / "R" / "library"),
    }

    failures: list[str] = []
    for cohort in registry["cohorts"]:
        cohort_id = cohort["cohort_id"]
        if selected and cohort_id not in selected:
            continue
        metadata = cohort["metadata"]
        command = [
            str(executable),
            "benchmark-run",
            "--input",
            str((root / cohort["local_path"]).resolve()),
            "--truth-key",
            metadata["truth_key"],
            "--study-key",
            "__ctp_study__",
            "--constant-study-id",
            cohort["constant_study_id"],
            "--donor-key",
            metadata["donor_key"],
            "--cluster-key",
            metadata["cluster_key"],
            "--cluster-map",
            str((root / cohort["cluster_map_path"]).resolve()),
            "--species",
            cohort["species"],
            "--tissue",
            cohort["tissue"],
            "--output",
            str((root / Path(cohort["predictions_path"]).parent).resolve()),
            "--strategy",
            "donor",
            "--methods",
            ",".join(methods),
            "--label-map",
            str((root / cohort["label_map_path"]).resolve()),
        ]
        for method in methods:
            if method in comparator_configs:
                command.extend(
                    ["--comparator-config", str(comparator_configs[method].resolve())]
                )
        if args.continue_on_unavailable:
            command.append("--continue-on-unavailable")
        completed = subprocess.run(
            command,
            cwd=root.parent.parent,
            env=environment,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
        log_path = logs / f"{cohort_id}.log"
        log_path.write_text(
            f"exit_code={completed.returncode}\nSTDOUT\n{completed.stdout}\nSTDERR\n{completed.stderr}",
            encoding="utf-8",
        )
        print(f"{cohort_id}: exit={completed.returncode} log={log_path}", flush=True)
        if completed.returncode:
            failures.append(cohort_id)
            if not args.continue_on_unavailable:
                break
    if failures:
        raise SystemExit(f"failed cohorts: {', '.join(failures)}")


if __name__ == "__main__":
    main()
