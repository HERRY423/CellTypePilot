"""Governed configuration for package-native annotation backends."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

NATIVE_BACKEND_CONFIG_SCHEMA = "celltypepilot.native-backends.v1"
SUPPORTED_NATIVE_BACKENDS = {
    "celltypist",
    "popv",
    "singler",
    "scanvi",
    "custom_reference",
    "llm",
}
PATH_FIELDS = {"reference_path", "model_path", "adapter_path"}


class NativeBackendConfigError(ValueError):
    """Raised when a native-backend configuration is unsafe or ambiguous."""


def _canonical_json(payload: Any) -> bytes:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode(
        "utf-8"
    )


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve_local_path(value: Any, base: Path, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise NativeBackendConfigError(f"{field} must be a non-empty path string")
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return str(path.resolve())


def _validate_backend(entry: Any, base: Path, index: int) -> dict:
    if not isinstance(entry, dict):
        raise NativeBackendConfigError(f"backends[{index}] must be an object")
    item = dict(entry)
    backend = str(item.get("backend", "")).strip().casefold().replace("-", "_")
    if backend not in SUPPORTED_NATIVE_BACKENDS:
        allowed = ", ".join(sorted(SUPPORTED_NATIVE_BACKENDS))
        raise NativeBackendConfigError(f"backends[{index}].backend must be one of: {allowed}")
    item["backend"] = backend
    item["enabled"] = item.get("enabled", True)
    if not isinstance(item["enabled"], bool):
        raise NativeBackendConfigError(f"backends[{index}].enabled must be boolean")

    for field in PATH_FIELDS:
        if item.get(field) is not None:
            item[field] = _resolve_local_path(item[field], base, field)

    timeout = item.get("timeout_seconds", 7200 if backend == "singler" else 14400)
    if not isinstance(timeout, int) or not 30 <= timeout <= 86400:
        raise NativeBackendConfigError(
            f"backends[{index}].timeout_seconds must be an integer from 30 to 86400"
        )
    item["timeout_seconds"] = timeout

    environment = item.get("environment", {})
    if not isinstance(environment, dict) or not all(
        isinstance(key, str) and key.replace("_", "").isalnum() and isinstance(value, str)
        for key, value in environment.items()
    ):
        raise NativeBackendConfigError(
            f"backends[{index}].environment must map environment-variable names to strings"
        )
    item["environment"] = environment

    if backend in {"popv", "singler", "scanvi", "custom_reference"}:
        if not item.get("reference_path"):
            raise NativeBackendConfigError(f"{backend} requires reference_path")
        item["label_key"] = str(item.get("label_key", "cell_type")).strip()
        if not item["label_key"]:
            raise NativeBackendConfigError(f"{backend}.label_key cannot be empty")

    if backend == "celltypist":
        mode = str(item.get("mode", "pretrained")).casefold()
        if mode not in {"pretrained", "retrain"}:
            raise NativeBackendConfigError("celltypist.mode must be pretrained or retrain")
        if mode == "retrain" and not item.get("reference_path"):
            raise NativeBackendConfigError("celltypist retrain mode requires reference_path")
        item["mode"] = mode
        item["label_key"] = str(item.get("label_key", "cell_type")).strip()

    if backend == "popv":
        mode = str(item.get("mode", "retrain")).casefold()
        if mode not in {"retrain", "inference", "fast"}:
            raise NativeBackendConfigError("popv.mode must be retrain, inference, or fast")
        item["mode"] = mode

    if backend == "custom_reference":
        method = str(item.get("method", "correlation")).casefold()
        if method not in {"correlation", "knn"}:
            raise NativeBackendConfigError("custom_reference.method must be correlation or knn")
        item["method"] = method

    if backend == "llm":
        provider = str(item.get("provider", "openai")).casefold()
        if provider != "openai":
            raise NativeBackendConfigError("llm.provider currently supports only openai")
        if not item.get("model"):
            raise NativeBackendConfigError("llm.model is required")
        if item.get("allow_network") is not True:
            raise NativeBackendConfigError(
                "llm requires explicit allow_network=true; otherwise omit or disable it"
            )
        item["provider"] = provider
        item["api_key_env"] = str(item.get("api_key_env", "OPENAI_API_KEY"))
        if not item["api_key_env"].replace("_", "").isalnum():
            raise NativeBackendConfigError("llm.api_key_env must be an environment variable name")

    return item


def load_native_backend_config(path: str | Path) -> dict:
    """Load, resolve, and validate a native-backend JSON configuration."""
    config_path = Path(path).expanduser().resolve()
    if not config_path.exists():
        raise NativeBackendConfigError(f"Native backend config not found: {config_path}")
    try:
        payload = json.loads(config_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise NativeBackendConfigError(f"Invalid native backend JSON: {exc}") from exc
    if not isinstance(payload, dict):
        raise NativeBackendConfigError("Native backend config must be a JSON object")
    if payload.get("schema_version") != NATIVE_BACKEND_CONFIG_SCHEMA:
        raise NativeBackendConfigError(f"schema_version must be {NATIVE_BACKEND_CONFIG_SCHEMA!r}")
    entries = payload.get("backends")
    if not isinstance(entries, list) or not entries:
        raise NativeBackendConfigError("backends must be a non-empty list")
    backends = [_validate_backend(item, config_path.parent, i) for i, item in enumerate(entries)]
    names = [item["backend"] for item in backends if item["enabled"]]
    if len(names) != len(set(names)):
        raise NativeBackendConfigError("Enabled backend names must be unique")
    continue_on_failure = payload.get("continue_on_failure", True)
    if not isinstance(continue_on_failure, bool):
        raise NativeBackendConfigError("continue_on_failure must be boolean")
    resume = payload.get("resume", True)
    if not isinstance(resume, bool):
        raise NativeBackendConfigError("resume must be boolean")
    normalized = {
        "schema_version": NATIVE_BACKEND_CONFIG_SCHEMA,
        "continue_on_failure": continue_on_failure,
        "resume": resume,
        "backends": backends,
        "config_path": str(config_path),
    }
    normalized["config_sha256"] = hashlib.sha256(_canonical_json(normalized)).hexdigest()
    return normalized


def hash_native_backend_dependencies(config: dict) -> dict[str, str]:
    """Hash every local reference/model/adapter dependency declared by the config."""
    hashes: dict[str, str] = {}
    for item in config.get("backends", []):
        if not item.get("enabled", True):
            continue
        for field in PATH_FIELDS:
            value = item.get(field)
            if value is None:
                continue
            path = Path(value)
            if not path.exists():
                raise NativeBackendConfigError(f"{item['backend']}.{field} not found: {path}")
            if str(path) in hashes:
                continue
            if path.is_dir():
                digest = hashlib.sha256()
                for child in sorted(
                    candidate for candidate in path.rglob("*") if candidate.is_file()
                ):
                    digest.update(str(child.relative_to(path)).replace("\\", "/").encode("utf-8"))
                    digest.update(_sha256_file(child).encode("ascii"))
                hashes[str(path)] = digest.hexdigest()
            else:
                hashes[str(path)] = _sha256_file(path)
    return hashes
