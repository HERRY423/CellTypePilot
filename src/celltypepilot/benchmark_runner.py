"""Execute comparators behind a leakage-resistant plugin benchmark protocol."""

from __future__ import annotations

import importlib.metadata
import json
import os
import subprocess
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from .benchmark import BenchmarkValidationError, validate_out_of_fold_predictions


def configure_benchmark_runtime(output_dir: str | Path) -> dict[str, str]:
    """Use benchmark-owned writable temp and Numba cache directories.

    Scanpy imports Numba-backed kernels. On restricted Windows hosts, the
    process can otherwise spend minutes probing an unwritable default temp
    directory before any scientific computation starts.
    """
    root = Path(output_dir).resolve() / "_runtime_cache"
    temp_dir = root / "tmp"
    numba_dir = root / "numba"
    temp_dir.mkdir(parents=True, exist_ok=True)
    numba_dir.mkdir(parents=True, exist_ok=True)
    os.environ["TMP"] = str(temp_dir)
    os.environ["TEMP"] = str(temp_dir)
    os.environ["NUMBA_CACHE_DIR"] = str(numba_dir)
    tempfile.tempdir = str(temp_dir)
    return {
        "temp_dir": str(temp_dir),
        "numba_cache_dir": str(numba_dir),
    }


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def _atomic_write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False, lineterminator="\n")
    os.replace(temporary, path)


@dataclass(frozen=True)
class CommandComparator:
    """A shell-free argv contract for an external comparator adapter."""

    method: str
    argv: tuple[str, ...]
    timeout_seconds: int = 3600
    version_command: tuple[str, ...] = ()
    reference_policy: str = "fold_train_only"
    confidence_semantics: str = ""
    environment: tuple[tuple[str, str], ...] = ()
    config_dir: Path = Path(".")

    @classmethod
    def from_json(cls, path: str | Path) -> CommandComparator:
        config_path = Path(path).resolve()
        payload = json.loads(config_path.read_text(encoding="utf-8"))
        required = {"method", "argv", "reference_policy", "confidence_semantics"}
        missing = required - set(payload)
        if missing:
            raise BenchmarkValidationError(f"Comparator config missing: {sorted(missing)}")
        if payload["method"] not in {"singler", "azimuth", "popv"}:
            raise BenchmarkValidationError(
                "External comparator method must be singler, azimuth, or popv"
            )
        if not isinstance(payload["argv"], list) or not payload["argv"]:
            raise BenchmarkValidationError("Comparator argv must be a non-empty JSON list")
        environment = payload.get("environment", {})
        if not isinstance(environment, dict) or not all(
            isinstance(key, str) and isinstance(value, str) for key, value in environment.items()
        ):
            raise BenchmarkValidationError(
                "Comparator environment must be a JSON object of string values"
            )
        if payload["reference_policy"] != "fold_train_only":
            raise BenchmarkValidationError(
                "External comparator configs must declare reference_policy='fold_train_only'"
            )
        placeholders = set(payload["argv"])
        required_placeholders = {"{train_h5ad}", "{test_h5ad}", "{output_csv}"}
        if not required_placeholders <= placeholders:
            raise BenchmarkValidationError(
                "Comparator argv must include separate {train_h5ad}, {test_h5ad}, "
                "and {output_csv} entries"
            )
        return cls(
            method=payload["method"],
            argv=tuple(str(value) for value in payload["argv"]),
            timeout_seconds=int(payload.get("timeout_seconds", 3600)),
            version_command=tuple(str(value) for value in payload.get("version_command", [])),
            reference_policy=payload["reference_policy"],
            confidence_semantics=str(payload["confidence_semantics"]),
            environment=tuple(environment.items()),
            config_dir=config_path.parent,
        )


def _safe_query_obs(adata: ad.AnnData, cluster_key: str, truth_key: str) -> ad.AnnData:
    """Remove truth and label-like metadata before a comparator sees test data."""
    query = adata.copy()
    blocked_tokens = ("truth", "cell_type", "celltype", "annotation", "ground_truth")
    drop = [
        column
        for column in query.obs.columns
        if column != cluster_key
        and (column == truth_key or any(token in column.lower() for token in blocked_tokens))
    ]
    if drop:
        query.obs = query.obs.drop(columns=drop)
    return query


def materialize_fold(
    adata: ad.AnnData,
    assignments: pd.DataFrame,
    fold_id: str,
    truth_key: str,
    cluster_key: str,
    output_dir: str | Path,
    species: str,
    tissue: str,
) -> dict[str, Path]:
    """Write fold-train reference and truth-stripped test query artifacts."""
    if truth_key not in adata.obs or cluster_key not in adata.obs:
        raise BenchmarkValidationError("truth_key and cluster_key must be present in obs")
    test_ids = assignments.loc[assignments["fold_id"] == fold_id, "cell_id"].astype(str)
    if test_ids.empty:
        raise BenchmarkValidationError(f"Fold {fold_id!r} has no test cells")
    all_ids = pd.Index(adata.obs_names.astype(str))
    test_mask = all_ids.isin(test_ids)
    train_mask = ~test_mask
    if not np.any(train_mask):
        raise BenchmarkValidationError(f"Fold {fold_id!r} has no training cells")

    train_ids = set(all_ids[train_mask].astype(str))
    train_rows = assignments[assignments["cell_id"].astype(str).isin(train_ids)]
    fold_dir = Path(output_dir) / _safe_fold_name(fold_id)
    fold_dir.mkdir(parents=True, exist_ok=True)
    train_path = fold_dir / "train_reference.h5ad"
    test_path = fold_dir / "test_query.h5ad"

    if train_path.exists() and test_path.exists():
        train_backed = ad.read_h5ad(train_path, backed="r")
        test_backed = ad.read_h5ad(test_path, backed="r")
        try:
            observed_train = set(train_backed.obs_names.astype(str))
            observed_test = set(test_backed.obs_names.astype(str))
            blocked_tokens = ("truth", "cell_type", "celltype", "annotation", "ground_truth")
            leaked = [
                column
                for column in test_backed.obs.columns
                if column != cluster_key
                and (
                    column == truth_key or any(token in column.lower() for token in blocked_tokens)
                )
            ]
            if (
                observed_train == train_ids
                and observed_test == set(test_ids)
                and "cell_type" in train_backed.obs
                and not leaked
            ):
                return {"fold_dir": fold_dir, "train": train_path, "test": test_path}
        finally:
            if train_backed.file is not None:
                train_backed.file.close()
            if test_backed.file is not None:
                test_backed.file.close()

    train = adata[train_mask].copy()
    test = _safe_query_obs(adata[test_mask], cluster_key, truth_key)
    train.obs["cell_type"] = train.obs[truth_key].astype(str)
    train.uns["celltypepilot_reference"] = {
        "species": species,
        "tissues": [tissue],
        "source": "benchmark_fold_training_only",
        "version": "ephemeral-fold-v1",
        "label_ontology": "benchmark_locked_canonical_labels",
        "training_studies": sorted(set(train_rows["held_out_study"].astype(str))),
        "held_out_fold": fold_id,
    }
    train.write_h5ad(train_path)
    test.write_h5ad(test_path)
    return {"fold_dir": fold_dir, "train": train_path, "test": test_path}


def _safe_fold_name(fold_id: str) -> str:
    return "".join(
        character if character.isalnum() or character in "-_." else "_" for character in fold_id
    )


def _confidence_from_cluster_rows(
    test: ad.AnnData,
    cluster_key: str,
    rows: pd.DataFrame,
) -> pd.DataFrame:
    cluster = test.obs[cluster_key].astype(str)
    label_map = rows.set_index(rows["cluster"].astype(str))["cell_type"].to_dict()
    score_map = pd.to_numeric(rows["combined_score"], errors="coerce").fillna(0.0)
    score_map.index = rows["cluster"].astype(str)
    return pd.DataFrame(
        {
            "cell_id": test.obs_names.astype(str),
            "predicted_label": cluster.map(label_map).fillna("Unknown").to_numpy(),
            "confidence": cluster.map(score_map.to_dict()).fillna(0.0).clip(0, 1).to_numpy(),
        }
    )


def run_celltypepilot_fold(
    paths: dict[str, Path],
    cluster_key: str,
    species: str,
    tissue: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run the plugin with a fold-train reference; test truth is absent."""
    from .orchestrator import run_annotation_pipeline

    output = paths["fold_dir"] / "celltypepilot"
    result = run_annotation_pipeline(
        paths["test"],
        cluster_key,
        output,
        species=species,
        tissue=tissue,
        reference_path=paths["train"],
        ref_label_key="cell_type",
        reference_backend="correlation",
        no_figures=True,
    )
    test = ad.read_h5ad(paths["test"])
    predictions = _confidence_from_cluster_rows(test, cluster_key, result["critic_results"])
    try:
        version = importlib.metadata.version("celltypepilot")
    except importlib.metadata.PackageNotFoundError:
        from . import __version__

        version = __version__
    return predictions, {
        "implementation": "celltypepilot.orchestrator",
        "version": version,
        "reference_policy": "fold_train_only",
        "confidence_semantics": "cluster_evidence_score_not_probability_calibrated",
    }


def _celltypist_ready(adata: ad.AnnData) -> ad.AnnData:
    import scanpy as sc

    working = adata.copy()
    if "log1p" not in working.uns:
        sc.pp.normalize_total(working, target_sum=1e4)
        sc.pp.log1p(working)
    return working


def run_celltypist_fold(paths: dict[str, Path]) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Train CellTypist on fold-train cells and predict truth-stripped test cells."""
    try:
        import celltypist
    except ImportError as exc:
        raise BenchmarkValidationError("CellTypist is not installed") from exc

    train = _celltypist_ready(ad.read_h5ad(paths["train"]))
    test = _celltypist_ready(ad.read_h5ad(paths["test"]))
    shared = train.var_names.intersection(test.var_names)
    if len(shared) < 50:
        raise BenchmarkValidationError("CellTypist fold has fewer than 50 shared genes")
    train = train[:, shared].copy()
    test = test[:, shared].copy()
    model = celltypist.train(
        train,
        labels=train.obs["cell_type"].astype(str),
        check_expression=False,
        n_jobs=1,
    )
    result = celltypist.annotate(test, model=model, majority_voting=False)
    labels_frame = result.predicted_labels
    label_column = (
        "predicted_labels" if "predicted_labels" in labels_frame else labels_frame.columns[0]
    )
    probabilities = result.probability_matrix
    confidence = probabilities.max(axis=1).to_numpy(dtype=float)
    predictions = pd.DataFrame(
        {
            "cell_id": test.obs_names.astype(str),
            "predicted_label": labels_frame[label_column].astype(str).to_numpy(),
            "confidence": confidence,
        }
    )
    return predictions, {
        "implementation": "celltypist.train+annotate",
        "version": importlib.metadata.version("celltypist"),
        "reference_policy": "fold_train_only",
        "confidence_semantics": "maximum_multiclass_probability",
    }


def run_command_fold(
    spec: CommandComparator,
    paths: dict[str, Path],
    truth_key: str,
    cluster_key: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run SingleR/Azimuth/popV adapters without a shell and validate output."""
    output_path = paths["fold_dir"] / f"{spec.method}_predictions.csv"
    replacements = {
        "{train_h5ad}": str(paths["train"].resolve()),
        "{test_h5ad}": str(paths["test"].resolve()),
        "{output_csv}": str(output_path.resolve()),
        "{truth_key}": truth_key,
        "{cluster_key}": cluster_key,
    }
    argv = [
        replacements.get(value, value).replace("{config_dir}", str(spec.config_dir))
        for value in spec.argv
    ]
    adapter_environment = {
        key: value.replace("{config_dir}", str(spec.config_dir)) for key, value in spec.environment
    }
    environment = {
        **os.environ,
        **adapter_environment,
        "CELLTYPEPILOT_BENCHMARK_MODE": "1",
    }
    completed = subprocess.run(
        argv,
        cwd=paths["fold_dir"],
        env=environment,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=spec.timeout_seconds,
        check=False,
    )
    log_path = paths["fold_dir"] / f"{spec.method}.log"
    log_path.write_text(
        f"exit_code={completed.returncode}\nSTDOUT\n{completed.stdout}\nSTDERR\n{completed.stderr}",
        encoding="utf-8",
    )
    if completed.returncode != 0:
        raise BenchmarkValidationError(
            f"{spec.method} adapter failed with exit code {completed.returncode}; see {log_path}"
        )
    if not output_path.exists():
        raise BenchmarkValidationError(f"{spec.method} did not create {output_path.name}")
    predictions = pd.read_csv(output_path, dtype={"cell_id": str})
    required = {"cell_id", "predicted_label", "confidence"}
    missing = required - set(predictions)
    if missing:
        raise BenchmarkValidationError(f"{spec.method} output missing columns: {sorted(missing)}")
    version = "not_reported"
    if spec.version_command:
        version_result = subprocess.run(
            list(spec.version_command),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            cwd=paths["fold_dir"],
            env=environment,
            timeout=60,
            check=False,
        )
        version = (version_result.stdout or version_result.stderr).strip()[:500]
    return predictions, {
        "implementation": "external_argv_adapter",
        "argv": argv,
        "version": version,
        "reference_policy": spec.reference_policy,
        "confidence_semantics": spec.confidence_semantics,
        "log": str(log_path),
    }


def _validate_fold_output(predictions: pd.DataFrame, expected_cells: pd.Index) -> None:
    observed = pd.Index(predictions["cell_id"].astype(str))
    if observed.has_duplicates:
        raise BenchmarkValidationError("Comparator emitted duplicate cell identifiers")
    missing = expected_cells.difference(observed)
    extra = observed.difference(expected_cells)
    if len(missing) or len(extra):
        raise BenchmarkValidationError(
            f"Comparator cell set mismatch: missing={len(missing)}, extra={len(extra)}"
        )
    confidence = pd.to_numeric(predictions["confidence"], errors="coerce")
    if confidence.isna().any() or ((confidence < 0) | (confidence > 1)).any():
        raise BenchmarkValidationError("Comparator confidence must be numeric within [0, 1]")


def apply_locked_label_map(
    predictions: pd.DataFrame,
    label_map: pd.DataFrame | None,
) -> pd.DataFrame:
    """Map raw labels through a predeclared method-specific canonical table."""
    if label_map is None:
        return predictions
    required = {"method", "raw_label", "canonical_label"}
    missing = required - set(label_map)
    if missing:
        raise BenchmarkValidationError(f"Label map missing columns: {sorted(missing)}")
    if label_map.duplicated(["method", "raw_label"]).any():
        raise BenchmarkValidationError("Label map contains ambiguous method/raw_label pairs")
    mapping = {
        (str(row.method), str(row.raw_label)): str(row.canonical_label)
        for row in label_map.itertuples(index=False)
    }
    # Formatting-only case variants must not break a locked benchmark. Build a
    # folded fallback only when every variant resolves to the same canonical
    # label; biologically ambiguous folded aliases still fail closed.
    folded_candidates: dict[tuple[str, str], set[str]] = {}
    for (method, raw_label), canonical in mapping.items():
        key = (method.strip().casefold(), raw_label.strip().casefold())
        folded_candidates.setdefault(key, set()).add(canonical)
    ambiguous_folded = {key: values for key, values in folded_candidates.items() if len(values) > 1}
    if ambiguous_folded:
        examples = sorted(ambiguous_folded)[:5]
        raise BenchmarkValidationError(
            f"Label map contains case-insensitive ambiguities (examples: {examples})"
        )
    folded_mapping = {key: next(iter(values)) for key, values in folded_candidates.items()}

    def resolve(method: object, raw_label: object) -> str | None:
        exact = mapping.get((str(method), str(raw_label)))
        if exact is not None:
            return exact
        return folded_mapping.get(
            (str(method).strip().casefold(), str(raw_label).strip().casefold())
        )

    output = predictions.copy()
    if "raw_predicted_label" in output:
        raw_labels = output["raw_predicted_label"].astype(str)
        expected = [
            resolve(method, label)
            for method, label in zip(output["method"], raw_labels, strict=True)
        ]
        missing = [
            (str(method), str(label))
            for method, label, canonical in zip(output["method"], raw_labels, expected, strict=True)
            if canonical is None
        ]
        if missing:
            raise BenchmarkValidationError(
                f"Label map is not exhaustive (examples: {sorted(set(missing))[:5]})"
            )
        observed = output["predicted_label"].astype(str).tolist()
        if observed != expected:
            raise BenchmarkValidationError(
                "Previously mapped predictions do not match the locked label map"
            )
        return output

    output["raw_predicted_label"] = output["predicted_label"].astype(str)
    keys = list(
        zip(
            output["method"].astype(str),
            output["raw_predicted_label"].astype(str),
            strict=True,
        )
    )
    resolved = [resolve(method, label) for method, label in keys]
    missing = sorted(
        {key for key, canonical in zip(keys, resolved, strict=True) if canonical is None}
    )
    if missing:
        raise BenchmarkValidationError(f"Label map is not exhaustive (examples: {missing[:5]})")
    output["predicted_label"] = resolved
    return output


def run_benchmark_comparators(
    adata: ad.AnnData,
    assignments: pd.DataFrame,
    truth_key: str,
    cluster_key: str,
    output_dir: str | Path,
    species: str,
    tissue: str,
    methods: tuple[str, ...] = ("celltypepilot", "celltypist"),
    command_specs: tuple[CommandComparator, ...] = (),
    label_map: pd.DataFrame | None = None,
    continue_on_unavailable: bool = False,
    fold_ids: tuple[str, ...] | None = None,
    write_aggregate_tables: bool = True,
    worker_id: str | None = None,
    batch_id: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Execute requested comparators with atomic per-method/fold checkpoints.

    Distributed GPU/CPU workers should pass ``fold_ids`` for their assigned folds and
    set ``write_aggregate_tables=False`` so they only write atomic
    ``checkpoints/{fold}__{method}.{status.json,csv}`` files. Aggregators merge
    checkpoints read-only and must not re-execute folds.
    """
    output = Path(output_dir).resolve()
    configure_benchmark_runtime(output)
    checkpoint_dir = output / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    spec_by_method = {spec.method: spec for spec in command_specs}
    prediction_frames: dict[tuple[str, str], pd.DataFrame] = {}
    status_rows: dict[tuple[str, str], dict[str, Any]] = {}

    all_fold_ids = list(assignments["fold_id"].drop_duplicates().astype(str))
    if fold_ids is not None:
        wanted = {str(value) for value in fold_ids}
        unknown = sorted(wanted - set(all_fold_ids))
        if unknown:
            raise BenchmarkValidationError(
                f"Unknown fold_id(s) for this assignment plan: {unknown}"
            )
        selected_fold_ids = [fold for fold in all_fold_ids if fold in wanted]
        if not selected_fold_ids:
            raise BenchmarkValidationError("fold_ids filter selected zero folds")
    else:
        selected_fold_ids = all_fold_ids

    def checkpoint_paths(method: str, fold_id: str) -> tuple[Path, Path]:
        stem = f"{_safe_fold_name(fold_id)}__{method}"
        return checkpoint_dir / f"{stem}.status.json", checkpoint_dir / f"{stem}.csv"

    def worker_metadata() -> dict[str, Any]:
        meta: dict[str, Any] = {}
        if worker_id:
            meta["worker_id"] = worker_id
        if batch_id:
            meta["batch_id"] = batch_id
        return meta

    def persist_tables() -> None:
        # Partial workers must not rewrite global OOF/status tables; that would
        # clobber sibling nodes' completed folds under a shared output tree.
        if not write_aggregate_tables:
            return
        status_frame = pd.DataFrame(status_rows.values())
        if not status_frame.empty:
            _atomic_write_csv(status_frame, output / "comparator_status.csv")
        if prediction_frames:
            predictions_frame = pd.concat(prediction_frames.values(), ignore_index=True)
            _atomic_write_csv(predictions_frame, output / "out_of_fold_predictions.csv")

    for fold_id in selected_fold_ids:
        expected_cells = pd.Index(
            assignments.loc[assignments["fold_id"] == fold_id, "cell_id"].astype(str)
        )
        pending_methods = []
        for method in methods:
            key = (method, fold_id)
            status_path, prediction_path = checkpoint_paths(method, fold_id)
            if status_path.exists() and prediction_path.exists():
                try:
                    checkpoint = json.loads(status_path.read_text(encoding="utf-8"))
                    if checkpoint.get("status") == "completed":
                        frame = pd.read_csv(prediction_path, dtype={"cell_id": str})
                        _validate_fold_output(frame, expected_cells)
                        prediction_frames[key] = frame
                        status_rows[key] = {
                            "method": method,
                            "fold_id": fold_id,
                            "status": "completed",
                            "resumed_from_checkpoint": True,
                            **checkpoint.get("provenance", {}),
                            **worker_metadata(),
                        }
                        continue
                except (OSError, ValueError, json.JSONDecodeError, BenchmarkValidationError):
                    pass
            pending_methods.append(method)

        if not pending_methods:
            persist_tables()
            continue

        paths = materialize_fold(
            adata,
            assignments,
            fold_id,
            truth_key,
            cluster_key,
            output_dir,
            species,
            tissue,
        )
        for method in pending_methods:
            key = (method, fold_id)
            status_path, prediction_path = checkpoint_paths(method, fold_id)
            previous_status = None
            if status_path.exists():
                try:
                    previous_status = json.loads(status_path.read_text(encoding="utf-8")).get(
                        "status"
                    )
                except (OSError, ValueError, json.JSONDecodeError):
                    previous_status = "unreadable_checkpoint"
            running = {
                "method": method,
                "fold_id": fold_id,
                "status": "running",
                "started_at_utc": _utc_now(),
                "previous_status": previous_status,
                **worker_metadata(),
            }
            status_rows[key] = running
            _atomic_write_text(status_path, json.dumps(running, indent=2))
            persist_tables()
            try:
                if method == "celltypepilot":
                    frame, provenance = run_celltypepilot_fold(paths, cluster_key, species, tissue)
                elif method == "celltypist":
                    frame, provenance = run_celltypist_fold(paths)
                elif method in spec_by_method:
                    frame, provenance = run_command_fold(
                        spec_by_method[method], paths, truth_key, cluster_key
                    )
                else:
                    raise BenchmarkValidationError(
                        f"No executable adapter configured for {method!r}"
                    )
                _validate_fold_output(frame, expected_cells)
                frame["method"] = method
                frame["fold_id"] = fold_id
                frame = apply_locked_label_map(frame, label_map)
                _atomic_write_csv(frame, prediction_path)
                prediction_frames[key] = frame
                completed = {
                    "method": method,
                    "fold_id": fold_id,
                    "status": "completed",
                    "started_at_utc": running["started_at_utc"],
                    "completed_at_utc": _utc_now(),
                    "previous_status": previous_status,
                    "provenance": {**provenance, **worker_metadata()},
                    **worker_metadata(),
                }
                _atomic_write_text(status_path, json.dumps(completed, indent=2))
                status_rows[key] = {
                    "method": method,
                    "fold_id": fold_id,
                    "status": "completed",
                    "resumed_after_status": previous_status,
                    **provenance,
                    **worker_metadata(),
                }
                persist_tables()
            except Exception as exc:
                failed = {
                    "method": method,
                    "fold_id": fold_id,
                    "status": "failed_or_unavailable",
                    "started_at_utc": running["started_at_utc"],
                    "failed_at_utc": _utc_now(),
                    "previous_status": previous_status,
                    "error": str(exc),
                    **worker_metadata(),
                }
                _atomic_write_text(status_path, json.dumps(failed, indent=2))
                status_rows[key] = failed
                persist_tables()
                if not continue_on_unavailable:
                    raise

    predictions = (
        pd.concat(prediction_frames.values(), ignore_index=True)
        if prediction_frames
        else pd.DataFrame()
    )
    if not predictions.empty:
        predictions = apply_locked_label_map(predictions, label_map)
        validate_out_of_fold_predictions(assignments, predictions)
    status = pd.DataFrame(status_rows.values())
    persist_tables()
    return predictions, status
