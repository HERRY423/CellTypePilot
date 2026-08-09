"""Immutable benchmark asset catalog for object-storage / CDN mirrors.

The catalog records CELLxGENE datasets, Azimuth references, label maps, and Docker
images with URL, version, SHA-256, license, training-study provenance, tissue/species,
and availability status.

Hard boundaries:
- Catalog metadata is versioned and fail-closed. Records are not mutated in place;
  new versions get new ``asset_id``/``version`` pairs.
- Local object cache lives under the catalog root (e.g. ``objects/``) and must never
  write into active fold workspaces under ``benchmarks/**/runs/``.
- Availability is an explicit enum. Missing bytes are ``pending_upload`` or
  ``declared_unavailable``, never silently treated as ready.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import urlparse

CATALOG_SCHEMA = "celltypepilot.immutable-asset-catalog.v1"
STORAGE_POLICY_SCHEMA = "celltypepilot.asset-storage-policy.v1"

ASSET_KINDS = (
    "cellxgene_dataset",
    "azimuth_reference",
    "label_map",
    "docker_image",
)

AVAILABILITY_STATUSES = (
    "available",  # URL resolvable and SHA-256 verified or verified local cache
    "source_available",  # original publisher URL is the source of truth; CDN mirror optional
    "pending_upload",  # declared for CDN/object store but bytes not yet published
    "declared_unavailable",  # intentionally not provided for this release track
    "blocked_overlap_audit",  # present metadata but barred until training-study audit passes
    "local_build_only",  # Dockerfile/build recipe frozen; image digest not yet published to registry
)

_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_OBJECT_STORE_SAFE = re.compile(r"^[a-z0-9][a-z0-9._-]{0,127}$")


class AssetCatalogError(ValueError):
    """Raised when the asset catalog or storage policy is invalid."""


def file_sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_mapping(payload: Any, label: str) -> dict[str, Any]:
    if not isinstance(payload, dict):
        raise AssetCatalogError(f"{label} must be a JSON object")
    return payload


def _require_str(record: dict[str, Any], field: str, asset_id: str) -> str:
    value = record.get(field)
    if not isinstance(value, str) or not value.strip():
        raise AssetCatalogError(f"Asset {asset_id!r} missing non-empty string field {field!r}")
    return value.strip()


def _require_sha256(value: str, field: str, asset_id: str) -> str:
    digest = value.strip().lower()
    if not _SHA256_RE.match(digest):
        raise AssetCatalogError(
            f"Asset {asset_id!r} field {field!r} must be lowercase hex SHA-256"
        )
    return digest


def _validate_url(url: str, asset_id: str, field: str = "url") -> None:
    parsed = urlparse(url)
    if parsed.scheme not in {"http", "https", "s3", "gs", "oci", "file"}:
        raise AssetCatalogError(
            f"Asset {asset_id!r} {field} must use http(s)/s3/gs/oci/file, got {parsed.scheme!r}"
        )
    if parsed.scheme != "file" and not parsed.netloc and not parsed.path:
        raise AssetCatalogError(f"Asset {asset_id!r} {field} is not a usable URL")


def _validate_training_study_provenance(
    provenance: Any, asset_id: str, kind: str
) -> dict[str, Any]:
    mapping = _require_mapping(provenance, f"Asset {asset_id!r} training_study_provenance")
    required = {
        "summary",
        "source_studies",
        "overlap_audit_status",
        "eligible_for_primary_holdout_track",
    }
    missing = required - set(mapping)
    if missing:
        raise AssetCatalogError(
            f"Asset {asset_id!r} training_study_provenance missing: {sorted(missing)}"
        )
    if not isinstance(mapping["summary"], str) or not mapping["summary"].strip():
        raise AssetCatalogError(f"Asset {asset_id!r} training_study_provenance.summary required")
    studies = mapping["source_studies"]
    if not isinstance(studies, list):
        raise AssetCatalogError(
            f"Asset {asset_id!r} training_study_provenance.source_studies must be a list"
        )
    for index, study in enumerate(studies):
        if not isinstance(study, dict):
            raise AssetCatalogError(
                f"Asset {asset_id!r} source_studies[{index}] must be an object"
            )
        for key in ("study_id", "role"):
            if not str(study.get(key, "")).strip():
                raise AssetCatalogError(
                    f"Asset {asset_id!r} source_studies[{index}] needs {key}"
                )
    audit = str(mapping["overlap_audit_status"]).strip()
    allowed_audit = {
        "not_applicable",
        "pending",
        "passed_no_eval_overlap",
        "failed_eval_overlap",
        "not_audited",
    }
    if audit not in allowed_audit:
        raise AssetCatalogError(
            f"Asset {asset_id!r} overlap_audit_status must be one of {sorted(allowed_audit)}"
        )
    eligible = mapping["eligible_for_primary_holdout_track"]
    if not isinstance(eligible, bool):
        raise AssetCatalogError(
            f"Asset {asset_id!r} eligible_for_primary_holdout_track must be bool"
        )
    # Reference assets that can leak evaluation studies must stay explicit.
    if kind in {"azimuth_reference", "docker_image"} and eligible and audit not in {
        "passed_no_eval_overlap",
        "not_applicable",
    }:
        raise AssetCatalogError(
            f"Asset {asset_id!r} cannot be primary-holdout-eligible without a passed "
            "or not_applicable overlap audit"
        )
    return mapping


def validate_asset_record(record: dict[str, Any]) -> dict[str, Any]:
    """Validate one asset record and return a normalized copy."""
    asset_id = _require_str(record, "asset_id", record.get("asset_id", "<unknown>"))
    if not _OBJECT_STORE_SAFE.match(asset_id):
        raise AssetCatalogError(
            f"asset_id {asset_id!r} must be lowercase object-store safe "
            "(start with alnum; only ._- thereafter)"
        )
    kind = _require_str(record, "kind", asset_id)
    if kind not in ASSET_KINDS:
        raise AssetCatalogError(f"Asset {asset_id!r} kind must be one of {ASSET_KINDS}")
    version = _require_str(record, "version", asset_id)
    url = _require_str(record, "url", asset_id)
    _validate_url(url, asset_id, "url")
    sha256 = _require_sha256(_require_str(record, "sha256", asset_id), "sha256", asset_id)
    license_id = _require_str(record, "license", asset_id)
    species = _require_str(record, "species", asset_id)
    tissue = _require_str(record, "tissue", asset_id)
    availability = _require_str(record, "availability", asset_id)
    if availability not in AVAILABILITY_STATUSES:
        raise AssetCatalogError(
            f"Asset {asset_id!r} availability must be one of {AVAILABILITY_STATUSES}"
        )
    provenance = _validate_training_study_provenance(
        record.get("training_study_provenance"), asset_id, kind
    )

    normalized = {
        "asset_id": asset_id,
        "kind": kind,
        "version": version,
        "url": url,
        "sha256": sha256,
        "license": license_id,
        "species": species,
        "tissue": tissue,
        "availability": availability,
        "training_study_provenance": provenance,
    }

    # Optional but validated when present.
    for optional_url_field in ("cdn_url", "source_url", "mirror_url"):
        value = record.get(optional_url_field)
        if value is None or value == "":
            continue
        if not isinstance(value, str):
            raise AssetCatalogError(
                f"Asset {asset_id!r} {optional_url_field} must be a string when set"
            )
        _validate_url(value, asset_id, optional_url_field)
        normalized[optional_url_field] = value

    if "byte_size" in record and record["byte_size"] is not None:
        try:
            byte_size = int(record["byte_size"])
        except (TypeError, ValueError) as exc:
            raise AssetCatalogError(
                f"Asset {asset_id!r} byte_size must be an integer"
            ) from exc
        if byte_size < 0:
            raise AssetCatalogError(f"Asset {asset_id!r} byte_size must be >= 0")
        normalized["byte_size"] = byte_size

    for optional_str in (
        "title",
        "media_type",
        "object_key",
        "local_cache_relpath",
        "notes",
        "image_tag",
        "image_digest",
        "dockerfile_sha256",
    ):
        value = record.get(optional_str)
        if value is None or value == "":
            continue
        if not isinstance(value, str):
            raise AssetCatalogError(f"Asset {asset_id!r} {optional_str} must be a string")
        if optional_str.endswith("sha256") or optional_str == "dockerfile_sha256":
            value = _require_sha256(value, optional_str, asset_id)
        if optional_str == "image_digest" and value.startswith("sha256:"):
            digest_hex = value.split(":", 1)[1]
            _require_sha256(digest_hex, optional_str, asset_id)
        normalized[optional_str] = value

    if "related_cohort_ids" in record and record["related_cohort_ids"] is not None:
        related = record["related_cohort_ids"]
        if not isinstance(related, list) or not all(isinstance(x, str) for x in related):
            raise AssetCatalogError(
                f"Asset {asset_id!r} related_cohort_ids must be a list of strings"
            )
        normalized["related_cohort_ids"] = list(related)

    if "components" in record and record["components"] is not None:
        components = record["components"]
        if not isinstance(components, list):
            raise AssetCatalogError(f"Asset {asset_id!r} components must be a list")
        normalized_components = []
        for index, component in enumerate(components):
            if not isinstance(component, dict):
                raise AssetCatalogError(
                    f"Asset {asset_id!r} components[{index}] must be an object"
                )
            name = str(component.get("name", "")).strip()
            if not name:
                raise AssetCatalogError(
                    f"Asset {asset_id!r} components[{index}] needs name"
                )
            component_sha = _require_sha256(
                str(component.get("sha256", "")),
                f"components[{index}].sha256",
                asset_id,
            )
            entry = {"name": name, "sha256": component_sha}
            if component.get("url"):
                _validate_url(str(component["url"]), asset_id, f"components[{index}].url")
                entry["url"] = str(component["url"])
            if component.get("byte_size") is not None:
                entry["byte_size"] = int(component["byte_size"])
            normalized_components.append(entry)
        normalized["components"] = normalized_components

    # Preserve non-conflicting extras for forward compatibility after required checks.
    reserved = set(normalized)
    for key, value in record.items():
        if key not in reserved:
            normalized[key] = value
    return normalized


def validate_asset_catalog(payload: dict[str, Any]) -> dict[str, Any]:
    """Validate a full catalog document."""
    payload = _require_mapping(payload, "catalog")
    if payload.get("schema_version") != CATALOG_SCHEMA:
        raise AssetCatalogError(f"Catalog schema_version must be {CATALOG_SCHEMA!r}")
    for field in ("catalog_id", "frozen_at_utc", "assets"):
        if field not in payload:
            raise AssetCatalogError(f"Catalog missing field {field!r}")
    assets = payload["assets"]
    if not isinstance(assets, list) or not assets:
        raise AssetCatalogError("Catalog assets must be a non-empty list")

    seen: set[tuple[str, str]] = set()
    normalized_assets = []
    for record in assets:
        if not isinstance(record, dict):
            raise AssetCatalogError("Each asset record must be an object")
        normalized = validate_asset_record(record)
        key = (normalized["asset_id"], normalized["version"])
        if key in seen:
            raise AssetCatalogError(
                f"Duplicate asset_id/version pair: {key[0]}@{key[1]}"
            )
        seen.add(key)
        normalized_assets.append(normalized)

    out = dict(payload)
    out["assets"] = normalized_assets
    return out


def validate_storage_policy(payload: dict[str, Any]) -> dict[str, Any]:
    """Validate object-storage / CDN path policy."""
    payload = _require_mapping(payload, "storage_policy")
    if payload.get("schema_version") != STORAGE_POLICY_SCHEMA:
        raise AssetCatalogError(
            f"Storage policy schema_version must be {STORAGE_POLICY_SCHEMA!r}"
        )
    for field in (
        "policy_id",
        "object_store_uri_template",
        "cdn_url_template",
        "object_key_template",
        "immutability",
        "forbidden_write_globs",
    ):
        if field not in payload:
            raise AssetCatalogError(f"Storage policy missing field {field!r}")
    immutability = _require_mapping(payload["immutability"], "immutability")
    for field in ("content_addressed", "overwrite_policy", "versioning"):
        if field not in immutability:
            raise AssetCatalogError(f"immutability missing field {field!r}")
    if immutability.get("overwrite_policy") != "deny":
        raise AssetCatalogError("immutability.overwrite_policy must be 'deny'")
    if not immutability.get("content_addressed") is True:
        raise AssetCatalogError("immutability.content_addressed must be true")
    globs = payload["forbidden_write_globs"]
    if not isinstance(globs, list) or not globs:
        raise AssetCatalogError("forbidden_write_globs must be a non-empty list")
    # Always protect active fold workspaces.
    required_forbidden = "benchmarks/**/runs/**"
    if required_forbidden not in globs:
        raise AssetCatalogError(
            f"forbidden_write_globs must include {required_forbidden!r}"
        )
    return payload


def load_asset_catalog(path: str | Path) -> tuple[dict[str, Any], Path]:
    catalog_path = Path(path).resolve()
    payload = json.loads(catalog_path.read_text(encoding="utf-8"))
    return validate_asset_catalog(payload), catalog_path


def load_storage_policy(path: str | Path) -> tuple[dict[str, Any], Path]:
    policy_path = Path(path).resolve()
    payload = json.loads(policy_path.read_text(encoding="utf-8"))
    return validate_storage_policy(payload), policy_path


def object_key_for_asset(asset: dict[str, Any], policy: dict[str, Any]) -> str:
    """Render the content-addressed object key for an asset."""
    if asset.get("object_key"):
        return str(asset["object_key"])
    template = str(policy["object_key_template"])
    return template.format(
        kind=asset["kind"],
        asset_id=asset["asset_id"],
        version=asset["version"],
        sha256=asset["sha256"],
        sha256_prefix=asset["sha256"][:12],
    )


def cdn_url_for_asset(asset: dict[str, Any], policy: dict[str, Any]) -> str:
    if asset.get("cdn_url"):
        return str(asset["cdn_url"])
    template = str(policy["cdn_url_template"])
    return template.format(
        kind=asset["kind"],
        asset_id=asset["asset_id"],
        version=asset["version"],
        sha256=asset["sha256"],
        sha256_prefix=asset["sha256"][:12],
        object_key=object_key_for_asset(asset, policy),
    )


def object_store_uri_for_asset(asset: dict[str, Any], policy: dict[str, Any]) -> str:
    template = str(policy["object_store_uri_template"])
    return template.format(
        kind=asset["kind"],
        asset_id=asset["asset_id"],
        version=asset["version"],
        sha256=asset["sha256"],
        sha256_prefix=asset["sha256"][:12],
        object_key=object_key_for_asset(asset, policy),
    )


def resolve_file_url(url: str, catalog_root: Path) -> Path | None:
    """Resolve a file: URL relative to the catalog root, or None if not a file URL."""
    parsed = urlparse(url)
    if parsed.scheme != "file":
        return None
    # Support file:relative/path and file:///absolute/path forms.
    raw = parsed.path if parsed.path else url[len("file:") :]
    if raw.startswith("/") and len(raw) > 2 and raw[2] == ":":
        # Windows absolute path leaked as /C:/...
        raw = raw[1:]
    path = Path(raw)
    if not path.is_absolute():
        path = (catalog_root / path).resolve()
    else:
        path = path.resolve()
    return path


def local_object_path(
    asset: dict[str, Any],
    catalog_root: Path,
    policy: dict[str, Any],
) -> Path:
    """Path under the catalog root for a cached immutable object (never under runs/)."""
    rel = asset.get("local_cache_relpath")
    if rel:
        path = (catalog_root / rel).resolve()
    else:
        key = object_key_for_asset(asset, policy)
        path = (catalog_root / "objects" / key).resolve()
    # Fail closed if the resolved path escapes the catalog root or hits runs/.
    root = catalog_root.resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise AssetCatalogError(
            f"Local object path for {asset['asset_id']} escapes catalog root"
        ) from exc
    assert_path_not_fold_workspace(path)
    return path


def assert_path_not_fold_workspace(path: Path) -> None:
    """Refuse any write/sync target that sits under an active fold workspace tree."""
    resolved = path.resolve()
    parts = list(resolved.parts)
    for index, part in enumerate(parts):
        if part == "runs" and index + 1 < len(parts):
            # benchmarks/.../runs/<cohort>/... is forbidden
            raise AssetCatalogError(
                f"Refusing path under fold workspace runs/: {resolved}"
            )


def materialize_source_to_object_cache(
    asset: dict[str, Any],
    catalog_root: Path,
    policy: dict[str, Any],
    *,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Copy a local file: source into the content-addressed object cache.

    Never reads or writes fold workspaces under ``runs/``. Remote HTTP sources are
    not downloaded here (use a separate CDN publish pipeline).
    """
    source_url = str(asset.get("source_url") or asset.get("url") or "")
    source = resolve_file_url(source_url, catalog_root)
    target = local_object_path(asset, catalog_root, policy)
    assert_path_not_fold_workspace(target)
    result: dict[str, Any] = {
        "asset_id": asset["asset_id"],
        "action": "materialize",
        "source": str(source) if source is not None else None,
        "target": str(target),
        "dry_run": dry_run,
        "status": "skipped_not_file_source",
    }
    if source is None:
        return result
    assert_path_not_fold_workspace(source)
    if not source.is_file():
        result["status"] = "source_missing"
        return result
    observed = file_sha256(source)
    if observed != asset["sha256"]:
        result["status"] = "source_sha256_mismatch"
        result["observed_sha256"] = observed
        result["expected_sha256"] = asset["sha256"]
        return result
    if asset.get("byte_size") is not None and source.stat().st_size != int(asset["byte_size"]):
        result["status"] = "source_byte_size_mismatch"
        return result
    if target.is_file():
        if file_sha256(target) == asset["sha256"]:
            result["status"] = "already_cached_verified"
            return result
        result["status"] = "target_exists_with_mismatch"
        return result
    if dry_run:
        result["status"] = "would_materialize"
        return result
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + f".partial")
    temporary.write_bytes(source.read_bytes())
    if file_sha256(temporary) != asset["sha256"]:
        temporary.unlink(missing_ok=True)
        result["status"] = "write_verify_failed"
        return result
    temporary.replace(target)
    result["status"] = "materialized"
    return result


def verify_local_asset(
    asset: dict[str, Any],
    catalog_root: Path,
    policy: dict[str, Any],
) -> dict[str, Any]:
    """Verify a local cache object (or file: source) against the catalog record."""
    target = local_object_path(asset, catalog_root, policy)
    assert_path_not_fold_workspace(target)
    candidates: list[Path] = [target]
    for field in ("source_url", "url"):
        resolved = resolve_file_url(str(asset.get(field) or ""), catalog_root)
        if resolved is not None:
            assert_path_not_fold_workspace(resolved)
            if resolved not in candidates:
                candidates.append(resolved)

    result: dict[str, Any] = {
        "asset_id": asset["asset_id"],
        "version": asset["version"],
        "kind": asset["kind"],
        "availability": asset["availability"],
        "expected_sha256": asset["sha256"],
        "local_path": str(target),
        "local_present": target.is_file(),
        "status": "missing_local",
    }

    verified_path: Path | None = None
    for candidate in candidates:
        if not candidate.is_file():
            continue
        size = candidate.stat().st_size
        if asset.get("byte_size") is not None and size != int(asset["byte_size"]):
            result["status"] = "byte_size_mismatch"
            result["checked_path"] = str(candidate)
            result["local_bytes"] = size
            result["expected_bytes"] = int(asset["byte_size"])
            return result
        observed = file_sha256(candidate)
        if observed != asset["sha256"]:
            result["status"] = "sha256_mismatch"
            result["checked_path"] = str(candidate)
            result["observed_sha256"] = observed
            return result
        verified_path = candidate
        break

    if verified_path is not None:
        result["status"] = "verified"
        result["checked_path"] = str(verified_path)
        result["local_bytes"] = verified_path.stat().st_size
        result["local_present"] = True
        return result

    if asset["availability"] in {
        "source_available",
        "pending_upload",
        "declared_unavailable",
        "blocked_overlap_audit",
        "local_build_only",
    }:
        result["status"] = f"ok_remote_or_declared:{asset['availability']}"
        return result
    return result


def summarize_catalog(
    catalog: dict[str, Any],
    policy: dict[str, Any] | None = None,
    catalog_root: Path | None = None,
    verify_local: bool = False,
) -> dict[str, Any]:
    """Summarize availability by kind; optionally verify local caches."""
    assets = catalog["assets"]
    by_kind: dict[str, dict[str, int]] = {kind: {} for kind in ASSET_KINDS}
    by_status: dict[str, int] = {}
    verifications: list[dict[str, Any]] = []

    for asset in assets:
        kind = asset["kind"]
        status = asset["availability"]
        by_kind.setdefault(kind, {})
        by_kind[kind][status] = by_kind[kind].get(status, 0) + 1
        by_status[status] = by_status.get(status, 0) + 1
        if verify_local:
            if policy is None or catalog_root is None:
                raise AssetCatalogError("verify_local requires policy and catalog_root")
            verifications.append(verify_local_asset(asset, catalog_root, policy))

    summary: dict[str, Any] = {
        "schema_version": catalog.get("schema_version"),
        "catalog_id": catalog.get("catalog_id"),
        "n_assets": len(assets),
        "by_kind": by_kind,
        "by_availability": by_status,
        "assets": [
            {
                "asset_id": asset["asset_id"],
                "kind": asset["kind"],
                "version": asset["version"],
                "url": asset["url"],
                "sha256": asset["sha256"],
                "license": asset["license"],
                "species": asset["species"],
                "tissue": asset["tissue"],
                "availability": asset["availability"],
                "object_key": (
                    object_key_for_asset(asset, policy) if policy is not None else None
                ),
                "cdn_url": cdn_url_for_asset(asset, policy) if policy is not None else asset.get("cdn_url"),
                "training_study_provenance": asset["training_study_provenance"],
            }
            for asset in assets
        ],
    }
    if verify_local:
        summary["local_verification"] = verifications
        summary["local_verified_count"] = sum(
            1 for row in verifications if row["status"] == "verified"
        )
        summary["local_failure_count"] = sum(
            1
            for row in verifications
            if row["status"] in {"sha256_mismatch", "byte_size_mismatch"}
        )
    return summary


def assets_by_kind(catalog: dict[str, Any], kind: str) -> list[dict[str, Any]]:
    if kind not in ASSET_KINDS:
        raise AssetCatalogError(f"Unknown kind {kind!r}")
    return [asset for asset in catalog["assets"] if asset["kind"] == kind]


def filter_assets(
    catalog: dict[str, Any],
    *,
    kinds: Iterable[str] | None = None,
    availability: Iterable[str] | None = None,
    species: str | None = None,
    tissue: str | None = None,
) -> list[dict[str, Any]]:
    kind_set = set(kinds) if kinds is not None else None
    avail_set = set(availability) if availability is not None else None
    rows = []
    for asset in catalog["assets"]:
        if kind_set is not None and asset["kind"] not in kind_set:
            continue
        if avail_set is not None and asset["availability"] not in avail_set:
            continue
        if species is not None and asset["species"] != species:
            continue
        if tissue is not None and asset["tissue"] != tissue:
            continue
        rows.append(asset)
    return rows
