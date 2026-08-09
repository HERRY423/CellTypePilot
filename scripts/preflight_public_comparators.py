"""Record comparator availability without turning absence into silent omission."""

from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd


def python_status(package: str) -> tuple[bool, str]:
    if importlib.util.find_spec(package) is None:
        return False, "not_installed"
    try:
        return True, importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return True, "version_unknown"


def r_packages_status(packages: list[str]) -> tuple[bool, str]:
    rscript = shutil.which("Rscript")
    if not rscript:
        return False, "Rscript_not_found"
    expression = ";".join(
        f"cat('{package}=', requireNamespace('{package}', quietly=TRUE), '\\n')"
        for package in packages
    )
    completed = subprocess.run(
        [rscript, "-e", expression],
        check=False,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if completed.returncode != 0:
        return False, f"R_preflight_failed:{completed.stderr.strip()[:200]}"
    states = {
        line.split("=", 1)[0]: line.split("=", 1)[1].strip()
        for line in completed.stdout.splitlines()
        if "=" in line
    }
    missing = [package for package in packages if states.get(package) != "TRUE"]
    if missing:
        return False, "missing_R_packages:" + ",".join(missing)
    return True, "available"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", required=True)
    args = parser.parse_args()
    registry_path = Path(args.registry).resolve()
    registry = json.loads(registry_path.read_text(encoding="utf-8"))

    ctp_ok, ctp_detail = python_status("celltypepilot")
    celltypist_ok, celltypist_detail = python_status("celltypist")
    popv_ok, popv_detail = python_status("popv")
    singler_ok, singler_detail = r_packages_status(["SingleR", "zellkonverter"])
    azimuth_ok, azimuth_detail = r_packages_status(["Seurat", "Azimuth"])
    shared = {
        "celltypepilot": (
            "ready_not_run" if ctp_ok else "dependency_unavailable",
            ctp_detail,
        ),
        "celltypist": (
            "ready_not_run" if celltypist_ok else "dependency_unavailable",
            celltypist_detail,
        ),
        "singler": (
            "ready_not_run" if singler_ok else "dependency_unavailable",
            singler_detail,
        ),
        "popv": (
            "ready_not_run" if popv_ok else "dependency_unavailable",
            popv_detail,
        ),
        "azimuth": (
            "reference_not_configured" if azimuth_ok else "dependency_unavailable",
            (
                "A compatible fold-training reference and reference/idx.annoy pair "
                "must be audited per tissue; " + azimuth_detail
            ),
        ),
    }

    summaries = []
    for cohort in registry["cohorts"]:
        preflight_rows = [
            {
                "method": method,
                "status": shared[method][0],
                "detail": shared[method][1],
                "fold_id": "preflight",
            }
            for method in registry["required_methods"]
        ]
        output = (registry_path.parent / cohort["comparator_status_path"]).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        existing_rows: list[dict[str, object]] = []
        if output.exists():
            existing = pd.read_csv(output)
            existing = existing[existing["fold_id"].astype(str) != "preflight"]
            existing_rows = existing.to_dict(orient="records")
        executed_methods = {str(row["method"]) for row in existing_rows}
        rows = existing_rows + [
            row for row in preflight_rows if row["method"] not in executed_methods
        ]
        pd.DataFrame(rows).to_csv(output, index=False)
        summaries.append(
            {
                "cohort_id": cohort["cohort_id"],
                "path": str(output),
                "statuses": {row["method"]: row["status"] for row in rows},
            }
        )
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
