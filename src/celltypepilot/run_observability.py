"""Read-only run observability for plugin / Web Inspector dashboards.

Surfaces checkpoint status files, fold ETA, host CPU/GPU, failure reasons,
product hashes, and artifact stale state. Never writes into fold workspaces
or prediction tables.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

OBSERVABILITY_SCHEMA = "celltypepilot.run-observability.v1"
READ_ONLY_CONTRACT = (
    "This snapshot is read-only. It does not mutate checkpoints, predictions, "
    "or fold workspaces. Manual annotation overrides must use the append-only "
    "audit log and apply-overrides path — never the observability dashboard."
)

# Typical product files under a benchmark-run or annotation output directory.
PRODUCT_CANDIDATES = (
    "out_of_fold_predictions.csv",
    "comparator_status.csv",
    "holdout_assignments.csv",
    "holdout_plan.json",
    "benchmark_manifest.json",
    "benchmark_results.csv",
    "benchmark_results_by_fold.csv",
    "evidence_table.csv",
    "manifest.json",
    "data.annotated.h5ad",
    "artifact_status.json",
    "annotation_audit_log.jsonl",
)

_STATUS_STEM = re.compile(r"^(?P<stem>.+)__(?P<method>[A-Za-z0-9_.+-]+)\.status\.json$")


class ObservabilityError(ValueError):
    """Raised for invalid observability inputs (still never writes fold files)."""


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _parse_utc(value: Any) -> datetime | None:
    if value is None or value == "":
        return None
    text = str(value).strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError:
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _seconds_between(start: datetime | None, end: datetime | None) -> float | None:
    if start is None or end is None:
        return None
    return max(0.0, (end - start).total_seconds())


def file_sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_read_json(path: Path) -> dict[str, Any] | None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def discover_run_root(path: str | Path) -> Path:
    """Resolve a directory that may contain checkpoints/ directly or one level up."""
    root = Path(path).resolve()
    if not root.exists():
        raise ObservabilityError(f"Path does not exist: {root}")
    if not root.is_dir():
        raise ObservabilityError(f"Path is not a directory: {root}")
    if (root / "checkpoints").is_dir():
        return root
    # Allow pointing at checkpoints/ itself.
    if root.name == "checkpoints" and root.is_dir():
        return root.parent
    return root


def list_checkpoint_status_files(run_root: Path) -> list[Path]:
    checkpoint_dir = run_root / "checkpoints"
    if not checkpoint_dir.is_dir():
        return []
    return sorted(checkpoint_dir.glob("*.status.json"))


def _infer_method_fold_from_name(path: Path, payload: dict[str, Any]) -> tuple[str, str]:
    method = str(payload.get("method") or "").strip()
    fold_id = str(payload.get("fold_id") or "").strip()
    match = _STATUS_STEM.match(path.name)
    if match:
        if not method:
            method = match.group("method")
        if not fold_id:
            fold_id = match.group("stem")
    return method or "unknown", fold_id or path.stem


def load_checkpoint_records(run_root: Path) -> list[dict[str, Any]]:
    """Load all checkpoints/*.status.json as normalized read-only records."""
    now = _utc_now()
    records: list[dict[str, Any]] = []
    for path in list_checkpoint_status_files(run_root):
        payload = _safe_read_json(path)
        if payload is None:
            records.append(
                {
                    "status_path": str(path),
                    "status": "unreadable_checkpoint",
                    "method": "unknown",
                    "fold_id": path.stem,
                    "failure_reason": "status_json_unreadable",
                    "prediction_csv": None,
                    "prediction_csv_present": False,
                    "prediction_sha256": None,
                    "duration_seconds": None,
                    "elapsed_seconds": None,
                    "raw": None,
                }
            )
            continue

        method, fold_id = _infer_method_fold_from_name(path, payload)
        status = str(payload.get("status") or "unknown")
        started = _parse_utc(payload.get("started_at_utc"))
        completed = _parse_utc(payload.get("completed_at_utc"))
        failed_at = _parse_utc(payload.get("failed_at_utc"))
        end = completed or failed_at
        duration = _seconds_between(started, end)
        elapsed = _seconds_between(started, now) if status == "running" else duration

        # Sibling prediction product (read-only hash if present).
        pred_path = path.with_suffix("").with_suffix(".csv")
        # status path is *.status.json → stem ends with .status; sibling is *.csv
        pred_path = path.parent / (path.name.replace(".status.json", ".csv"))
        pred_present = pred_path.is_file()
        pred_sha = None
        if pred_present:
            try:
                pred_sha = file_sha256(pred_path)
            except OSError:
                pred_sha = None

        failure_reason = None
        if status in {"failed_or_unavailable", "failed", "error"}:
            failure_reason = (
                str(payload.get("error") or payload.get("detail") or status).strip()
                or status
            )
        elif status == "running" and not pred_present and started is not None:
            # Informational only — still running is not a failure.
            failure_reason = None
        elif status == "completed" and not pred_present:
            failure_reason = "completed_status_but_prediction_csv_missing"

        records.append(
            {
                "status_path": str(path),
                "status_filename": path.name,
                "method": method,
                "fold_id": fold_id,
                "status": status,
                "started_at_utc": payload.get("started_at_utc"),
                "completed_at_utc": payload.get("completed_at_utc"),
                "failed_at_utc": payload.get("failed_at_utc"),
                "previous_status": payload.get("previous_status"),
                "provenance": payload.get("provenance"),
                "failure_reason": failure_reason,
                "error": payload.get("error"),
                "duration_seconds": duration,
                "elapsed_seconds": elapsed,
                "prediction_csv": str(pred_path) if pred_present else str(pred_path),
                "prediction_csv_present": pred_present,
                "prediction_sha256": pred_sha,
                "raw": payload,
            }
        )
    return records


def estimate_fold_eta(records: list[dict[str, Any]]) -> dict[str, Any]:
    """Estimate remaining wall time from completed fold durations (read-only)."""
    completed_durations = [
        float(row["duration_seconds"])
        for row in records
        if row.get("status") == "completed" and row.get("duration_seconds") is not None
    ]
    by_method: dict[str, list[float]] = {}
    for row in records:
        if row.get("status") == "completed" and row.get("duration_seconds") is not None:
            by_method.setdefault(str(row["method"]), []).append(float(row["duration_seconds"]))

    mean_overall = (
        sum(completed_durations) / len(completed_durations) if completed_durations else None
    )
    mean_by_method = {
        method: sum(values) / len(values) for method, values in by_method.items() if values
    }

    running = [row for row in records if row.get("status") == "running"]
    pending_like = [
        row
        for row in records
        if row.get("status")
        not in {"completed", "failed_or_unavailable", "failed", "error", "unreadable_checkpoint"}
    ]
    # Counts for progress
    n_total = len(records)
    n_completed = sum(1 for row in records if row.get("status") == "completed")
    n_failed = sum(
        1
        for row in records
        if row.get("status") in {"failed_or_unavailable", "failed", "error"}
    )
    n_running = len(running)

    remaining_seconds = 0.0
    eta_components: list[dict[str, Any]] = []
    for row in running:
        method = str(row["method"])
        mean = mean_by_method.get(method, mean_overall)
        elapsed = float(row.get("elapsed_seconds") or 0.0)
        if mean is None:
            remaining_for_row = None
        else:
            remaining_for_row = max(0.0, mean - elapsed)
            remaining_seconds += remaining_for_row
        eta_components.append(
            {
                "method": method,
                "fold_id": row["fold_id"],
                "elapsed_seconds": elapsed,
                "expected_mean_seconds": mean,
                "eta_remaining_seconds": remaining_for_row,
            }
        )

    # If there are methods/folds not yet started we cannot know the plan from
    # status files alone; ETA covers observed running units only.
    return {
        "n_checkpoints": n_total,
        "n_completed": n_completed,
        "n_failed": n_failed,
        "n_running": n_running,
        "n_other": n_total - n_completed - n_failed - n_running,
        "progress_fraction": (n_completed / n_total) if n_total else None,
        "mean_completed_duration_seconds": mean_overall,
        "mean_duration_by_method_seconds": mean_by_method,
        "running_eta": eta_components,
        "estimated_remaining_seconds": remaining_seconds if running and mean_overall is not None else (
            remaining_seconds if any(c.get("eta_remaining_seconds") is not None for c in eta_components) else None
        ),
        "eta_basis": (
            "mean_completed_duration_minus_elapsed_for_running"
            if completed_durations
            else "insufficient_completed_history"
        ),
        "pending_observed_units": len(pending_like),
    }


def collect_host_resources() -> dict[str, Any]:
    """Best-effort CPU/GPU snapshot. Never raises for missing optional deps."""
    cpu: dict[str, Any] = {
        "platform": platform.platform(),
        "processor": platform.processor() or None,
        "python_version": platform.python_version(),
        "cpu_count_logical": os.cpu_count(),
    }
    try:
        load = os.getloadavg()  # type: ignore[attr-defined]
        cpu["load_avg_1_5_15"] = [float(load[0]), float(load[1]), float(load[2])]
    except (AttributeError, OSError):
        cpu["load_avg_1_5_15"] = None

    try:
        import psutil  # type: ignore

        cpu["cpu_percent"] = float(psutil.cpu_percent(interval=0.05))
        mem = psutil.virtual_memory()
        cpu["memory"] = {
            "total_bytes": int(mem.total),
            "available_bytes": int(mem.available),
            "percent": float(mem.percent),
        }
    except Exception:
        cpu["cpu_percent"] = None
        cpu["memory"] = None

    gpu: dict[str, Any] = {"available": False, "devices": [], "source": None}
    # Prefer nvidia-smi when present.
    nvsmi = shutil.which("nvidia-smi")
    if nvsmi:
        try:
            completed = subprocess.run(
                [
                    nvsmi,
                    "--query-gpu=index,name,utilization.gpu,memory.used,memory.total,temperature.gpu",
                    "--format=csv,noheader,nounits",
                ],
                check=False,
                capture_output=True,
                text=True,
                timeout=3,
            )
            if completed.returncode == 0 and completed.stdout.strip():
                devices = []
                for line in completed.stdout.strip().splitlines():
                    parts = [p.strip() for p in line.split(",")]
                    if len(parts) < 6:
                        continue
                    devices.append(
                        {
                            "index": int(parts[0]) if parts[0].isdigit() else parts[0],
                            "name": parts[1],
                            "utilization_gpu_percent": _float_or_none(parts[2]),
                            "memory_used_mib": _float_or_none(parts[3]),
                            "memory_total_mib": _float_or_none(parts[4]),
                            "temperature_c": _float_or_none(parts[5]),
                        }
                    )
                gpu = {"available": bool(devices), "devices": devices, "source": "nvidia-smi"}
        except (OSError, subprocess.SubprocessError, ValueError):
            pass

    if not gpu["available"]:
        try:
            import torch  # type: ignore

            if torch.cuda.is_available():
                devices = []
                for index in range(torch.cuda.device_count()):
                    props = torch.cuda.get_device_properties(index)
                    devices.append(
                        {
                            "index": index,
                            "name": props.name,
                            "total_memory_bytes": int(props.total_memory),
                        }
                    )
                gpu = {"available": True, "devices": devices, "source": "torch.cuda"}
        except Exception:
            pass

    return {"cpu": cpu, "gpu": gpu, "pid": os.getpid()}


def _float_or_none(value: str) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def hash_product_files(run_root: Path) -> list[dict[str, Any]]:
    """SHA-256 product inventory for known artifacts (read-only)."""
    rows: list[dict[str, Any]] = []
    for name in PRODUCT_CANDIDATES:
        path = run_root / name
        row: dict[str, Any] = {
            "name": name,
            "path": str(path),
            "present": path.is_file(),
            "byte_size": None,
            "sha256": None,
        }
        if path.is_file():
            try:
                row["byte_size"] = path.stat().st_size
                # Skip hashing multi-GB h5ad by default for dashboard latency;
                # still report size. Hash smaller products.
                if path.suffix.lower() in {".h5ad"} and path.stat().st_size > 50_000_000:
                    row["sha256"] = None
                    row["hash_skipped"] = "large_h5ad_skipped_for_latency"
                else:
                    row["sha256"] = file_sha256(path)
            except OSError as exc:
                row["error"] = str(exc)
        rows.append(row)

    # Also hash completed checkpoint prediction CSVs (bounded).
    checkpoint_dir = run_root / "checkpoints"
    if checkpoint_dir.is_dir():
        csvs = sorted(checkpoint_dir.glob("*.csv"))[:64]
        for path in csvs:
            try:
                rows.append(
                    {
                        "name": f"checkpoints/{path.name}",
                        "path": str(path),
                        "present": True,
                        "byte_size": path.stat().st_size,
                        "sha256": file_sha256(path),
                    }
                )
            except OSError as exc:
                rows.append(
                    {
                        "name": f"checkpoints/{path.name}",
                        "path": str(path),
                        "present": True,
                        "error": str(exc),
                    }
                )
    return rows


def load_stale_status(run_root: Path) -> dict[str, Any]:
    """Read artifact_status.json if present; otherwise derive a conservative view."""
    path = run_root / "artifact_status.json"
    if path.is_file():
        payload = _safe_read_json(path) or {}
        return {
            "source": "artifact_status.json",
            "path": str(path),
            "review_state": payload.get("review_state"),
            "derived_artifacts_stale": bool(
                payload.get("derived_artifacts_stale")
                or str(payload.get("review_state", "")).endswith("stale")
            ),
            "stale_reason": payload.get("stale_reason") or payload.get("message"),
            "stale_artifacts": payload.get("stale_artifacts")
            or payload.get("affected_artifacts"),
            "raw": payload,
        }

    # Benchmark-run heuristic: predictions present but no release assembly.
    has_predictions = (run_root / "out_of_fold_predictions.csv").is_file()
    has_results = (run_root / "benchmark_results.csv").is_file()
    incomplete = has_predictions and not has_results
    return {
        "source": "derived_heuristic",
        "path": None,
        "review_state": "not_applicable_no_artifact_status"
        if not (run_root / "data.annotated.h5ad").is_file()
        else "fresh_or_unknown",
        "derived_artifacts_stale": False,
        "stale_reason": (
            "No artifact_status.json; benchmark results incomplete relative to predictions"
            if incomplete
            else "No artifact_status.json present"
        ),
        "stale_artifacts": [],
        "raw": None,
    }


def summarize_failures(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    failures = []
    for row in records:
        reason = row.get("failure_reason") or row.get("error")
        if row.get("status") in {
            "failed_or_unavailable",
            "failed",
            "error",
            "unreadable_checkpoint",
        } or (row.get("status") == "completed" and not row.get("prediction_csv_present")):
            failures.append(
                {
                    "method": row.get("method"),
                    "fold_id": row.get("fold_id"),
                    "status": row.get("status"),
                    "failure_reason": reason or row.get("status"),
                    "status_path": row.get("status_path"),
                }
            )
    return failures


def build_observability_snapshot(
    path: str | Path,
    *,
    include_host: bool = True,
    include_product_hashes: bool = True,
) -> dict[str, Any]:
    """Assemble a full read-only observability snapshot for a run/output directory."""
    run_root = discover_run_root(path)
    records = load_checkpoint_records(run_root)
    eta = estimate_fold_eta(records)
    stale = load_stale_status(run_root)
    failures = summarize_failures(records)
    products = hash_product_files(run_root) if include_product_hashes else []
    host = collect_host_resources() if include_host else {"cpu": None, "gpu": None}

    by_status: dict[str, int] = {}
    for row in records:
        by_status[str(row["status"])] = by_status.get(str(row["status"]), 0) + 1

    return {
        "schema_version": OBSERVABILITY_SCHEMA,
        "generated_at_utc": _utc_now().isoformat(),
        "run_root": str(run_root),
        "read_only": True,
        "contract": READ_ONLY_CONTRACT,
        "prediction_mutation_allowed": False,
        "override_policy": "append_only_audit_log_via_web_review_or_cli_apply_overrides",
        "checkpoints": {
            "directory": str(run_root / "checkpoints"),
            "present": (run_root / "checkpoints").is_dir(),
            "n_status_files": len(records),
            "by_status": by_status,
            "records": records,
        },
        "fold_eta": eta,
        "host": host,
        "failures": failures,
        "products": products,
        "stale": stale,
    }
