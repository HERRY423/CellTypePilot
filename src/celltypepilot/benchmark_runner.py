"""Execute comparators behind a leakage-resistant plugin benchmark protocol."""

from __future__ import annotations

import importlib.metadata
import json
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from .benchmark import BenchmarkValidationError, validate_out_of_fold_predictions


@dataclass(frozen=True)
class CommandComparator:
    """A shell-free argv contract for an external comparator adapter."""

    method: str
    argv: tuple[str, ...]
    timeout_seconds: int = 3600
    version_command: tuple[str, ...] = ()
    reference_policy: str = "fold_train_only"
    confidence_semantics: str = ""

    @classmethod
    def from_json(cls, path: str | Path) -> CommandComparator:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
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

    train = adata[train_mask].copy()
    test = _safe_query_obs(adata[test_mask], cluster_key, truth_key)
    train.obs["cell_type"] = train.obs[truth_key].astype(str)
    train_ids = set(train.obs_names.astype(str))
    train_rows = assignments[assignments["cell_id"].astype(str).isin(train_ids)]
    train.uns["celltypepilot_reference"] = {
        "species": species,
        "tissues": [tissue],
        "source": "benchmark_fold_training_only",
        "version": "ephemeral-fold-v1",
        "label_ontology": "benchmark_locked_canonical_labels",
        "training_studies": sorted(set(train_rows["held_out_study"].astype(str))),
        "held_out_fold": fold_id,
    }

    fold_dir = Path(output_dir) / _safe_fold_name(fold_id)
    fold_dir.mkdir(parents=True, exist_ok=True)
    train_path = fold_dir / "train_reference.h5ad"
    test_path = fold_dir / "test_query.h5ad"
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
        "confidence_semantics": "cluster_combined_score_not_probability_calibrated",
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
        "{train_h5ad}": str(paths["train"]),
        "{test_h5ad}": str(paths["test"]),
        "{output_csv}": str(output_path),
        "{truth_key}": truth_key,
        "{cluster_key}": cluster_key,
    }
    argv = [replacements.get(value, value) for value in spec.argv]
    completed = subprocess.run(
        argv,
        cwd=paths["fold_dir"],
        env={**os.environ, "CELLTYPEPILOT_BENCHMARK_MODE": "1"},
        capture_output=True,
        text=True,
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
    output = predictions.copy()
    output["raw_predicted_label"] = output["predicted_label"].astype(str)
    output["predicted_label"] = [
        mapping.get((str(method), str(label)), str(label))
        for method, label in zip(output["method"], output["raw_predicted_label"], strict=True)
    ]
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
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Execute requested comparators on every locked fold and record failures."""
    spec_by_method = {spec.method: spec for spec in command_specs}
    prediction_frames = []
    status_rows = []
    for fold_id in assignments["fold_id"].drop_duplicates().astype(str):
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
        expected_cells = pd.Index(
            assignments.loc[assignments["fold_id"] == fold_id, "cell_id"].astype(str)
        )
        for method in methods:
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
                prediction_frames.append(frame)
                status_rows.append(
                    {"method": method, "fold_id": fold_id, "status": "completed", **provenance}
                )
            except Exception as exc:
                status_rows.append(
                    {
                        "method": method,
                        "fold_id": fold_id,
                        "status": "failed_or_unavailable",
                        "error": str(exc),
                    }
                )
                if not continue_on_unavailable:
                    raise

    predictions = (
        pd.concat(prediction_frames, ignore_index=True) if prediction_frames else pd.DataFrame()
    )
    if not predictions.empty:
        predictions = apply_locked_label_map(predictions, label_map)
        validate_out_of_fold_predictions(assignments, predictions)
    return predictions, pd.DataFrame(status_rows)
