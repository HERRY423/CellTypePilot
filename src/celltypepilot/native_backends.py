"""Package-native candidate backend execution for ordinary annotation runs.

Every backend is isolated behind the same fail-closed contract.  Backends may
generate candidates; only the downstream hierarchical selector can publish a
draft identity.  LLM output is permanently hypothesis-only.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import os
import shutil
import subprocess
import time
from importlib.resources import files
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .candidate_backends import (
    OUTPUT_COLUMNS,
    aggregate_cell_candidates,
    concatenate_candidates,
    normalize_candidate_table,
)
from .native_backend_config import hash_native_backend_dependencies

NATIVE_BACKEND_RUN_SCHEMA = "celltypepilot.native-backend-run.v1"
NATIVE_BACKEND_STATUS_SCHEMA = "celltypepilot.native-backend-status.v1"
STATUS_COLUMNS = [
    "schema_version",
    "backend",
    "status",
    "mode",
    "backend_version",
    "signature_sha256",
    "elapsed_seconds",
    "n_candidates",
    "candidate_artifact",
    "log_artifact",
    "error_type",
    "error_detail",
    "claim_boundary",
]


class NativeBackendRunError(RuntimeError):
    """Raised when a strict native-backend run cannot complete."""


class NativeBackendUnavailableError(NativeBackendRunError):
    """Raised when a requested optional runtime is unavailable."""


def check_native_backend_runtimes() -> dict[str, dict[str, Any]]:
    """Report import/executable availability without training or network calls."""
    rscript = shutil.which("Rscript")
    return {
        "celltypist": {
            "available": importlib.util.find_spec("celltypist") is not None,
            "check": "python_import_only",
        },
        "popv": {
            "available": importlib.util.find_spec("popv") is not None,
            "check": "python_import_only",
        },
        "singler": {
            "available": rscript is not None,
            "check": "Rscript_executable_only_R_packages_unverified",
            "executable": rscript,
        },
        "scanvi": {
            "available": importlib.util.find_spec("scvi") is not None,
            "check": "python_import_only",
        },
        "custom_reference": {"available": True, "check": "core_runtime"},
        "llm": {
            "available": importlib.util.find_spec("openai") is not None,
            "check": "sdk_import_only_credentials_and_network_unverified",
        },
    }


def _canonical_hash(payload: Any) -> str:
    content = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False, default=str
    )
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unavailable"


def _safe_error(exc: Exception) -> str:
    return " ".join(str(exc).replace("\r", " ").replace("\n", " ").split())[:1000]


def _write_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def _validate_counts(matrix, backend: str) -> None:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if values.size == 0:
        raise NativeBackendRunError(f"{backend} requires a non-empty raw count matrix")
    if not np.isfinite(values).all() or np.nanmin(values) < 0:
        raise NativeBackendRunError(f"{backend} requires finite non-negative raw counts")
    if not np.allclose(values, np.rint(values), atol=1e-6):
        raise NativeBackendRunError(
            f"{backend} requires integer raw counts; configure counts_layer or AnnData.raw"
        )


def _counts_view(adata: ad.AnnData, layer: str | None, backend: str) -> ad.AnnData:
    if layer:
        if layer not in adata.layers:
            raise NativeBackendRunError(f"{backend} counts_layer not found: {layer}")
        work = ad.AnnData(
            X=adata.layers[layer].copy(),
            obs=adata.obs.copy(),
            var=adata.var.copy(),
        )
    elif adata.raw is not None:
        work = adata.raw.to_adata()
        work.obs = adata.obs.copy()
    else:
        work = adata.copy()
    _validate_counts(work.X, backend)
    return work


def _sanitize_query(
    query: ad.AnnData,
    cluster_key: str,
    *,
    extra_obs_keys: tuple[str, ...] = (),
) -> ad.AnnData:
    keep = [cluster_key, *extra_obs_keys]
    keep = [key for key in dict.fromkeys(keep) if key in query.obs]
    work = query.copy()
    work.obs = work.obs.loc[:, keep].copy()
    work.obs["__ctp_cell_id__"] = work.obs_names.astype(str)
    work.uns = {}
    return work


def _sanitize_reference(reference: ad.AnnData, label_key: str, *obs_keys: str) -> ad.AnnData:
    if label_key not in reference.obs:
        raise NativeBackendRunError(f"Reference lacks label key: {label_key}")
    keep = [label_key, *(key for key in obs_keys if key and key in reference.obs)]
    work = reference.copy()
    work.obs = work.obs.loc[:, list(dict.fromkeys(keep))].copy()
    reference_contract = work.uns.get("celltypepilot_reference")
    work.uns = {}
    if reference_contract is not None:
        work.uns["celltypepilot_reference"] = reference_contract
    return work


def _shared_features(
    reference: ad.AnnData, query: ad.AnnData, backend: str, minimum: int = 50
) -> tuple[ad.AnnData, ad.AnnData]:
    shared = reference.var_names.intersection(query.var_names)
    if len(shared) < minimum:
        raise NativeBackendRunError(
            f"{backend} found {len(shared)} shared genes; at least {minimum} are required"
        )
    return reference[:, shared].copy(), query[:, shared].copy()


def _load_reference(entry: dict) -> ad.AnnData:
    path = Path(entry["reference_path"])
    if not path.exists():
        raise NativeBackendRunError(f"Reference not found: {path}")
    return ad.read_h5ad(path)


def _run_reference_backend(
    query: ad.AnnData,
    cluster_key: str,
    entry: dict,
    species: str,
    tissue: str,
) -> tuple[pd.DataFrame, dict]:
    backend = entry["backend"]
    method = entry.get("method", backend)
    reference = _load_reference(entry) if entry.get("reference_path") else None
    if backend == "scanvi":
        query = _sanitize_query(
            _counts_view(query, entry.get("query_counts_layer"), "scanvi"),
            cluster_key,
        )
        reference = _sanitize_reference(
            _counts_view(reference, entry.get("ref_counts_layer"), "scanvi"),
            entry.get("label_key", "cell_type"),
        )
    from .reference_scorer import score_by_reference

    result = score_by_reference(
        query,
        cluster_key,
        reference=reference,
        ref_label_key=entry.get("label_key", "cell_type"),
        model_path=entry.get("model_path"),
        backend=method,
        n_neighbors=int(entry.get("n_neighbors", 15)),
        species=species,
        tissue=tissue,
        allow_unverified_reference=bool(entry.get("allow_unverified_reference", False)),
    )
    raw = result.rename(columns={"ref_score": "score", "ref_rank": "rank"}).copy()
    raw["backend"] = "custom_reference" if backend == "custom_reference" else backend
    raw["score_semantics"] = "backend_cluster_support_not_cross_backend_probability"
    package = "scvi-tools" if backend == "scanvi" else "celltypist"
    return raw, {
        "backend_version": _package_version(package),
        "reference_contract": dict(result.attrs.get("reference_contract", {})),
    }


def _run_celltypist_retrain(
    query: ad.AnnData, cluster_key: str, entry: dict, run_dir: Path
) -> tuple[pd.DataFrame, dict]:
    try:
        import celltypist
        import scanpy as sc
    except ImportError as exc:
        raise NativeBackendUnavailableError(
            "CellTypist is unavailable; install celltypepilot[reference]"
        ) from exc

    reference = _sanitize_reference(_load_reference(entry), entry["label_key"])
    query_work = _sanitize_query(query, cluster_key)
    reference, query_work = _shared_features(reference, query_work, "celltypist")
    sc.pp.normalize_total(reference, target_sum=1e4)
    sc.pp.log1p(reference)
    sc.pp.normalize_total(query_work, target_sum=1e4)
    sc.pp.log1p(query_work)
    model = celltypist.train(
        reference,
        labels=entry["label_key"],
        n_jobs=int(entry.get("n_jobs", 1)),
        max_iter=int(entry.get("max_iter", 100)),
        use_SGD=bool(entry.get("use_sgd", False)),
        feature_selection=bool(entry.get("feature_selection", True)),
    )
    model_path = run_dir / "celltypist_model.pkl"
    model.write(str(model_path))
    predictions = celltypist.annotate(query_work, model=model, majority_voting=False)
    labels = predictions.predicted_labels
    if isinstance(labels, pd.DataFrame):
        labels = labels.iloc[:, 0]
    probabilities = predictions.probability_matrix
    confidence = probabilities.max(axis=1) if isinstance(probabilities, pd.DataFrame) else np.nan
    raw = pd.DataFrame(
        {
            "cell_id": query_work.obs["__ctp_cell_id__"].astype(str).to_numpy(),
            "predicted_label": pd.Series(labels).astype(str).to_numpy(),
            "confidence": np.asarray(confidence, dtype=float),
            "backend": "celltypist",
            "score_semantics": "maximum_model_probability_within_celltypist",
        }
    )
    return raw, {
        "backend_version": _package_version("celltypist"),
        "model_artifact": str(model_path),
    }


def _run_popv(
    query: ad.AnnData, cluster_key: str, entry: dict, run_dir: Path
) -> tuple[pd.DataFrame, dict]:
    try:
        import popv
    except ImportError as exc:
        raise NativeBackendUnavailableError(
            "popV is unavailable; install the optional popV runtime"
        ) from exc

    label_key = entry["label_key"]
    query_batch_key = str(entry.get("query_batch_key", "__ctp_query_batch__"))
    ref_batch_key = str(entry.get("ref_batch_key", "__ctp_ref_batch__"))
    query_counts = _counts_view(query, entry.get("query_counts_layer"), "popv")
    reference_counts = _counts_view(_load_reference(entry), entry.get("ref_counts_layer"), "popv")
    query_work = _sanitize_query(
        query_counts,
        cluster_key,
        extra_obs_keys=(query_batch_key,),
    )
    reference = _sanitize_reference(reference_counts, label_key, ref_batch_key)
    reference, query_work = _shared_features(reference, query_work, "popv")
    if query_batch_key not in query_work.obs:
        query_work.obs[query_batch_key] = "query"
    if ref_batch_key not in reference.obs:
        reference.obs[ref_batch_key] = "reference"
    model_dir = run_dir / "models"
    model_dir.mkdir(exist_ok=True)
    process = popv.preprocessing.Process_Query(
        query_work,
        reference,
        query_batch_key=query_batch_key,
        ref_labels_key=label_key,
        ref_batch_key=ref_batch_key,
        cl_obo_folder=False,
        unknown_celltype_label="Unknown",
        save_path_trained_models=str(model_dir),
        prediction_mode=entry.get("mode", "retrain"),
        hvg=min(int(entry.get("hvg", 4000)), reference.n_vars),
    )
    processed = process.adata
    output_dir = run_dir / "popv_output"
    output_dir.mkdir(exist_ok=True)
    popv.annotation.annotate_data(processed, save_path=str(output_dir))
    if "_dataset" not in processed.obs:
        raise NativeBackendRunError("popV output lacks the _dataset query/reference boundary")
    result = processed.obs.loc[processed.obs["_dataset"].astype(str).eq("query")]
    required = {"popv_prediction", "popv_prediction_score"}
    if not required.issubset(result.columns):
        raise NativeBackendRunError("popV did not emit prediction and agreement score columns")
    ids = (
        result["__ctp_cell_id__"].astype(str).to_numpy()
        if "__ctp_cell_id__" in result
        else result.index.astype(str).to_numpy()
    )
    score = pd.to_numeric(result["popv_prediction_score"], errors="raise").to_numpy(float)
    prediction_keys = processed.uns.get(
        "prediction_keys_seen", processed.uns.get("prediction_keys", [])
    )
    if score.size and np.nanmax(score) > 1:
        score = score / max(1, len(prediction_keys))
    raw = pd.DataFrame(
        {
            "cell_id": ids,
            "predicted_label": result["popv_prediction"].astype(str).to_numpy(),
            "confidence": np.clip(score, 0, 1),
            "backend": "popv",
            "score_semantics": "popv_expert_agreement_fraction_not_calibrated_probability",
        }
    )
    return raw, {
        "backend_version": _package_version("popv"),
        "prediction_keys": [str(key) for key in prediction_keys],
        "model_artifact": str(model_dir),
    }


def _run_singler(
    query: ad.AnnData, cluster_key: str, entry: dict, run_dir: Path
) -> tuple[pd.DataFrame, dict]:
    reference = _sanitize_reference(_load_reference(entry), entry["label_key"])
    query_work = _sanitize_query(query, cluster_key)
    reference, query_work = _shared_features(reference, query_work, "singler")
    reference_path = run_dir / "reference.sanitized.h5ad"
    query_path = run_dir / "query.truth_stripped.h5ad"
    output_path = run_dir / "raw_candidates.csv"
    reference.write_h5ad(reference_path)
    query_work.write_h5ad(query_path)
    adapter_path = entry.get("adapter_path") or str(
        files("celltypepilot.adapters") / "run_singler.R"
    )
    rscript = str(entry.get("rscript", "Rscript"))
    command = [
        rscript,
        str(adapter_path),
        str(reference_path),
        str(query_path),
        str(output_path),
        entry["label_key"],
    ]
    environment = {**os.environ, **entry.get("environment", {})}
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=entry["timeout_seconds"],
            shell=False,
            env=environment,
        )
    except FileNotFoundError as exc:
        raise NativeBackendUnavailableError(
            "Rscript is unavailable; SingleR requires R, SingleR, zellkonverter, and SummarizedExperiment"
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise NativeBackendRunError(
            f"SingleR timed out after {entry['timeout_seconds']} seconds"
        ) from exc
    (run_dir / "stdout.log").write_text(completed.stdout or "", encoding="utf-8")
    (run_dir / "stderr.log").write_text(completed.stderr or "", encoding="utf-8")
    if completed.returncode != 0:
        raise NativeBackendRunError(
            f"SingleR exited with code {completed.returncode}: {_safe_error(RuntimeError(completed.stderr))}"
        )
    if not output_path.exists():
        raise NativeBackendRunError("SingleR completed without writing candidates")
    raw = pd.read_csv(output_path)
    raw["backend"] = "singler"
    try:
        version = subprocess.run(
            [rscript, "-e", "cat(as.character(packageVersion('SingleR')))"],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
            shell=False,
            env=environment,
        ).stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        version = "unknown"
    return raw, {
        "backend_version": version or "unknown",
        "command": command[:2],
        "adapter_sha256": _sha256_file(Path(adapter_path)),
    }


def _llm_payload(scores: pd.DataFrame, resolver: dict) -> list[dict]:
    payload = []
    for cluster, group in scores.groupby(scores["cluster"].astype(str), sort=True):
        top = group.sort_values("combined_score", ascending=False).head(5)
        candidates = []
        for _, row in top.iterrows():
            label = str(row.get("cell_type", ""))
            candidates.append(
                {
                    "label": label,
                    "cl_id": resolver.get("cl_by_name", {}).get(label, ""),
                    "marker_evidence_score": float(row.get("combined_score", 0.0)),
                    "positive_markers_detected": str(row.get("positive_markers_detected", ""))[
                        :500
                    ],
                    "negative_markers_detected": str(row.get("negative_markers_detected", ""))[
                        :500
                    ],
                }
            )
        payload.append({"cluster": str(cluster), "allowed_candidates": candidates})
    return payload


def _run_llm(scores: pd.DataFrame, resolver: dict, entry: dict) -> tuple[pd.DataFrame, dict]:
    if entry.get("allow_network") is not True:
        raise NativeBackendRunError("LLM backend was not explicitly authorized for network use")
    api_key = os.environ.get(entry["api_key_env"])
    if not api_key:
        raise NativeBackendUnavailableError(
            f"LLM API key environment variable is unset: {entry['api_key_env']}"
        )
    try:
        from openai import OpenAI
    except ImportError as exc:
        raise NativeBackendUnavailableError(
            "OpenAI Python SDK is unavailable; install celltypepilot[llm]"
        ) from exc
    packets = _llm_payload(scores, resolver)
    request_hash = _canonical_hash(packets)
    schema = {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "decisions": {
                "type": "array",
                "items": {
                    "type": "object",
                    "additionalProperties": False,
                    "properties": {
                        "cluster": {"type": "string"},
                        "selected_label": {"type": "string"},
                        "abstain": {"type": "boolean"},
                        "rationale": {"type": "string", "maxLength": 500},
                    },
                    "required": ["cluster", "selected_label", "abstain", "rationale"],
                },
            }
        },
        "required": ["decisions"],
    }
    client = OpenAI(api_key=api_key)
    response = client.responses.create(
        model=entry["model"],
        input=[
            {
                "role": "system",
                "content": (
                    "You are a bounded cell-type candidate reviewer. For each cluster, either "
                    "select exactly one label from allowed_candidates or abstain. Marker scores "
                    "are evidence-ranking signals, not probabilities. Do not invent labels."
                ),
            },
            {"role": "user", "content": json.dumps(packets, ensure_ascii=False)},
        ],
        text={
            "format": {
                "type": "json_schema",
                "name": "celltypepilot_llm_candidates",
                "schema": schema,
                "strict": True,
            }
        },
    )
    if not getattr(response, "output_text", ""):
        raise NativeBackendRunError("LLM returned no structured output")
    try:
        decoded = json.loads(response.output_text)
    except json.JSONDecodeError as exc:
        raise NativeBackendRunError("LLM output was not valid JSON") from exc
    allowed = {
        item["cluster"]: {candidate["label"] for candidate in item["allowed_candidates"]}
        for item in packets
    }
    rows = []
    seen: set[str] = set()
    for decision in decoded.get("decisions", []):
        cluster = str(decision.get("cluster", ""))
        if cluster not in allowed or cluster in seen:
            raise NativeBackendRunError("LLM returned an unknown or duplicate cluster")
        seen.add(cluster)
        if decision.get("abstain"):
            continue
        label = str(decision.get("selected_label", ""))
        if label not in allowed[cluster]:
            raise NativeBackendRunError("LLM selected a label outside allowed_candidates")
        rows.append(
            {
                "cluster": cluster,
                "cell_type": label,
                "backend": "llm",
                "rank": 1,
                "score_semantics": "llm_hypothesis_without_numeric_probability",
            }
        )
    usage = getattr(response, "usage", None)
    usage_payload = usage.model_dump() if hasattr(usage, "model_dump") else str(usage or "")
    return pd.DataFrame(rows), {
        "backend_version": _package_version("openai"),
        "provider": "openai",
        "model": entry["model"],
        "request_sha256": request_hash,
        "response_id": str(getattr(response, "id", "")),
        "usage": usage_payload,
    }


def _normalize_result(
    raw: pd.DataFrame,
    query: ad.AnnData,
    cluster_key: str,
    resolver: dict,
    backend: str,
    source_artifact: Path,
    version: str,
) -> pd.DataFrame:
    if raw.empty:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)
    work = raw
    if "cluster" not in work and "cluster_id" not in work:
        work = aggregate_cell_candidates(work, query.obs[cluster_key])
    return normalize_candidate_table(
        work,
        resolver,
        source_artifact=str(source_artifact),
        default_backend=backend,
        source_version=version,
    )


def run_native_backends(
    query: ad.AnnData,
    cluster_key: str,
    marker_scores: pd.DataFrame,
    resolver: dict,
    output_dir: str | Path,
    config: dict,
    *,
    species: str,
    tissue: str,
    input_sha256: str,
) -> dict:
    """Run configured backends with checkpoints and structured failure artifacts."""
    root = Path(output_dir) / "native_backends"
    root.mkdir(parents=True, exist_ok=True)
    dependency_hashes = hash_native_backend_dependencies(config)
    candidate_tables: list[pd.DataFrame] = []
    statuses: list[dict] = []
    backend_metadata: dict[str, dict] = {}

    for entry in config.get("backends", []):
        if not entry.get("enabled", True):
            continue
        backend = entry["backend"]
        run_dir = root / backend
        run_dir.mkdir(parents=True, exist_ok=True)
        raw_candidate_path = run_dir / "raw_candidates.csv"
        candidate_path = run_dir / "candidates.normalized.csv"
        checkpoint_path = run_dir / "checkpoint.json"
        log_path = run_dir / "run.json"
        signature = _canonical_hash(
            {
                "schema_version": NATIVE_BACKEND_RUN_SCHEMA,
                "backend": entry,
                "input_sha256": input_sha256,
                "dependencies": dependency_hashes,
                "cluster_key": cluster_key,
            }
        )
        started = time.monotonic()
        if config.get("resume", True) and checkpoint_path.exists() and candidate_path.exists():
            try:
                checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                checkpoint = {}
            if (
                checkpoint.get("signature_sha256") == signature
                and checkpoint.get("status") == "completed"
            ):
                cached = pd.read_csv(candidate_path)
                candidate_tables.append(cached)
                metadata = checkpoint.get("metadata", {})
                backend_metadata[backend] = metadata
                statuses.append(
                    {
                        "schema_version": NATIVE_BACKEND_STATUS_SCHEMA,
                        "backend": backend,
                        "status": "completed_from_checkpoint",
                        "mode": entry.get("mode", entry.get("method", "native")),
                        "backend_version": metadata.get("backend_version", "unknown"),
                        "signature_sha256": signature,
                        "elapsed_seconds": 0.0,
                        "n_candidates": len(cached),
                        "candidate_artifact": str(candidate_path),
                        "log_artifact": str(log_path),
                        "error_type": "",
                        "error_detail": "",
                        "claim_boundary": "candidate_generation_only_not_final_identity",
                    }
                )
                continue

        try:
            if backend == "celltypist" and entry.get("mode") == "retrain":
                raw, metadata = _run_celltypist_retrain(query, cluster_key, entry, run_dir)
            elif backend in {"celltypist", "scanvi", "custom_reference"}:
                raw, metadata = _run_reference_backend(query, cluster_key, entry, species, tissue)
            elif backend == "popv":
                raw, metadata = _run_popv(query, cluster_key, entry, run_dir)
            elif backend == "singler":
                raw, metadata = _run_singler(query, cluster_key, entry, run_dir)
            elif backend == "llm":
                raw, metadata = _run_llm(marker_scores, resolver, entry)
            else:  # validated config makes this unreachable
                raise NativeBackendRunError(f"Unsupported backend: {backend}")
            normalized = _normalize_result(
                raw,
                query,
                cluster_key,
                resolver,
                backend,
                candidate_path,
                metadata.get("backend_version", "unknown"),
            )
            raw.to_csv(raw_candidate_path, index=False)
            normalized.to_csv(candidate_path, index=False)
            candidate_tables.append(normalized)
            backend_metadata[backend] = metadata
            checkpoint = {
                "schema_version": NATIVE_BACKEND_RUN_SCHEMA,
                "status": "completed",
                "signature_sha256": signature,
                "metadata": metadata,
                "n_candidates": len(normalized),
            }
            _write_json(checkpoint_path, checkpoint)
            _write_json(log_path, checkpoint)
            status = "completed"
            error_type = ""
            error_detail = ""
            version = metadata.get("backend_version", "unknown")
            n_candidates = len(normalized)
        except Exception as exc:
            unavailable = isinstance(exc, NativeBackendUnavailableError)
            status = "unavailable" if unavailable else "failed"
            error_type = type(exc).__name__
            error_detail = _safe_error(exc)
            version = "unavailable" if unavailable else "unknown"
            n_candidates = 0
            failure = {
                "schema_version": NATIVE_BACKEND_RUN_SCHEMA,
                "status": status,
                "signature_sha256": signature,
                "error_type": error_type,
                "error_detail": error_detail,
            }
            _write_json(checkpoint_path, failure)
            _write_json(log_path, failure)
            if not config.get("continue_on_failure", True):
                raise NativeBackendRunError(
                    f"Native backend {backend} {status}: {error_detail}"
                ) from exc
        statuses.append(
            {
                "schema_version": NATIVE_BACKEND_STATUS_SCHEMA,
                "backend": backend,
                "status": status,
                "mode": entry.get("mode", entry.get("method", "native")),
                "backend_version": version,
                "signature_sha256": signature,
                "elapsed_seconds": round(time.monotonic() - started, 3),
                "n_candidates": n_candidates,
                "candidate_artifact": str(candidate_path) if candidate_path.exists() else "",
                "log_artifact": str(log_path),
                "error_type": error_type,
                "error_detail": error_detail,
                "claim_boundary": "candidate_generation_only_not_final_identity",
            }
        )

    status_frame = pd.DataFrame(statuses, columns=STATUS_COLUMNS)
    status_path = root / "native_backend_status.csv"
    status_frame.to_csv(status_path, index=False)
    return {
        "candidates": concatenate_candidates(candidate_tables),
        "status": status_frame,
        "status_path": status_path,
        "root": root,
        "config_sha256": config["config_sha256"],
        "dependency_sha256": dependency_hashes,
        "backend_metadata": backend_metadata,
    }


def run_fold_native_backend(
    train_path: str | Path,
    query_path: str | Path,
    cluster_key: str,
    backend: str,
    run_dir: str | Path,
    *,
    species: str,
    tissue: str,
    entry_overrides: dict | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Execute one fold-train-only backend using the ordinary native runners.

    This bridge keeps public-domain validation on the same implementation as
    ``annotate --native-backends`` while preserving truth stripping performed
    by :func:`benchmark_runner.materialize_fold`.
    """
    if backend not in {"celltypist", "popv", "singler", "scanvi", "custom_reference"}:
        raise NativeBackendRunError(f"Unsupported fold-native backend: {backend}")
    train_path = Path(train_path).resolve()
    query_path = Path(query_path).resolve()
    work_dir = Path(run_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    overrides = dict(entry_overrides or {})
    reuse_raw_path = overrides.pop("reuse_raw_path", None)
    entry = {
        "backend": backend,
        "reference_path": str(train_path),
        "label_key": "cell_type",
        "mode": "retrain",
        "timeout_seconds": 14400,
        "allow_unverified_reference": False,
        **overrides,
    }
    if backend == "custom_reference":
        entry.setdefault("method", "correlation")
    query = ad.read_h5ad(query_path)
    reuse_path = Path(reuse_raw_path) if reuse_raw_path else None
    if reuse_path is not None and reuse_path.is_file():
        raw = pd.read_csv(reuse_path, dtype={"cell_id": str})
        checkpoint_path = reuse_path.parent / "checkpoint.json"
        try:
            checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            checkpoint = {}
        metadata = {
            **checkpoint.get("metadata", {}),
            "reused_product_native_artifact": str(reuse_path),
        }
    elif backend == "celltypist":
        raw, metadata = _run_celltypist_retrain(query, cluster_key, entry, work_dir)
    elif backend == "popv":
        raw, metadata = _run_popv(query, cluster_key, entry, work_dir)
    elif backend == "singler":
        raw, metadata = _run_singler(query, cluster_key, entry, work_dir)
    else:
        raw, metadata = _run_reference_backend(query, cluster_key, entry, species, tissue)

    if "cell_id" in raw:
        label_column = next(
            name for name in ("predicted_label", "cell_type", "prediction", "label") if name in raw
        )
        if "confidence" in raw:
            confidence = pd.to_numeric(raw["confidence"], errors="coerce")
        elif "similarity" in raw:
            similarity = pd.to_numeric(raw["similarity"], errors="coerce")
            confidence = (similarity + 1.0) / 2.0
        else:
            confidence = pd.Series(0.0, index=raw.index)
        predictions = pd.DataFrame(
            {
                "cell_id": raw["cell_id"].astype(str),
                "predicted_label": raw[label_column].astype(str),
                "confidence": confidence.clip(0, 1).fillna(0.0),
            }
        )
    else:
        ranks = (
            pd.to_numeric(raw["rank"], errors="coerce").fillna(1)
            if "rank" in raw
            else pd.Series(1, index=raw.index)
        )
        top = raw[ranks.eq(1)].copy()
        if top.duplicated("cluster").any():
            raise NativeBackendRunError(f"{backend} emitted duplicate cluster top predictions")
        labels = top.set_index(top["cluster"].astype(str))["cell_type"].to_dict()
        score_column = "score" if "score" in top else "raw_score"
        scores = (
            pd.to_numeric(top[score_column], errors="coerce").fillna(0.0)
            if score_column in top
            else pd.Series(0.0, index=top.index)
        )
        scores.index = top["cluster"].astype(str)
        clusters = query.obs[cluster_key].astype(str)
        predictions = pd.DataFrame(
            {
                "cell_id": query.obs_names.astype(str),
                "predicted_label": clusters.map(labels).fillna("Unknown").to_numpy(),
                "confidence": clusters.map(scores.to_dict()).fillna(0.0).clip(0, 1).to_numpy(),
            }
        )
    return predictions, {
        **metadata,
        "implementation": "celltypepilot.native_backends.run_fold_native_backend",
        "reference_policy": "fold_train_only",
        "confidence_semantics": (
            "backend_specific_rank_or_similarity_transformed_for_selective_evaluation; "
            "not_cross_backend_probability"
        ),
    }
