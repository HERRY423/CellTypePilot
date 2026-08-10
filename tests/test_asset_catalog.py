"""Tests for the immutable asset catalog (no fold workspace writes)."""

from __future__ import annotations

from pathlib import Path

import pytest

from celltypepilot.asset_catalog import (
    AssetCatalogError,
    assert_path_not_fold_workspace,
    load_asset_catalog,
    load_storage_policy,
    materialize_source_to_object_cache,
    summarize_catalog,
    validate_asset_catalog,
    validate_asset_record,
    validate_storage_policy,
)

REPO = Path(__file__).resolve().parents[1]
CATALOG = REPO / "benchmarks" / "assets" / "catalog.json"
POLICY = REPO / "benchmarks" / "assets" / "storage_policy.json"


def test_repo_catalog_and_policy_validate():
    catalog, catalog_path = load_asset_catalog(CATALOG)
    policy, _ = load_storage_policy(POLICY)
    assert catalog["schema_version"].startswith("celltypepilot.immutable-asset-catalog")
    assert policy["immutability"]["overwrite_policy"] == "deny"
    assert "benchmarks/**/runs/**" in policy["forbidden_write_globs"]
    kinds = {asset["kind"] for asset in catalog["assets"]}
    assert kinds == {
        "cellxgene_dataset",
        "azimuth_reference",
        "label_map",
        "docker_image",
    }
    for asset in catalog["assets"]:
        for field in (
            "url",
            "version",
            "sha256",
            "license",
            "species",
            "tissue",
            "availability",
            "training_study_provenance",
        ):
            assert field in asset
    assert catalog_path.parent.name == "assets"


def test_label_maps_verify_against_file_sources():
    catalog, catalog_path = load_asset_catalog(CATALOG)
    policy, _ = load_storage_policy(POLICY)
    summary = summarize_catalog(
        catalog,
        policy=policy,
        catalog_root=catalog_path.parent,
        verify_local=True,
    )
    label_rows = [row for row in summary["local_verification"] if row["kind"] == "label_map"]
    assert label_rows
    assert all(row["status"] == "verified" for row in label_rows)


def test_azimuth_slots_are_blocked_until_audit():
    catalog, _ = load_asset_catalog(CATALOG)
    azimuth = [a for a in catalog["assets"] if a["kind"] == "azimuth_reference"]
    assert azimuth
    for asset in azimuth:
        assert asset["availability"] == "blocked_overlap_audit"
        assert asset["training_study_provenance"]["eligible_for_primary_holdout_track"] is False


def test_refuse_fold_workspace_path():
    # Avoid system temp dirs (permission issues on some Windows agents).
    workspace = REPO / "scratch" / "pytest_temp" / "asset_catalog_refuse_runs"
    bad = workspace / "benchmarks" / "public_v1" / "runs" / "cohort" / "file.h5ad"
    bad.parent.mkdir(parents=True, exist_ok=True)
    bad.write_text("x", encoding="utf-8")
    with pytest.raises(AssetCatalogError, match="runs/"):
        assert_path_not_fold_workspace(bad)


def test_materialize_does_not_target_runs():
    import hashlib
    import shutil

    catalog_root = REPO / "scratch" / "pytest_temp" / "asset_catalog_materialize"
    if catalog_root.exists():
        shutil.rmtree(catalog_root)
    catalog_root.mkdir(parents=True)
    # Source relative to catalog root, never under runs/
    external = catalog_root / "incoming" / "src_label.csv"
    external.parent.mkdir(parents=True)
    external.write_text("method,raw_label,canonical_label\n__truth__,A,a\n", encoding="utf-8")
    digest = hashlib.sha256(external.read_bytes()).hexdigest()

    asset = validate_asset_record(
        {
            "asset_id": "labelmap-test",
            "kind": "label_map",
            "version": "v1",
            "url": "file:incoming/src_label.csv",
            "source_url": "file:incoming/src_label.csv",
            "sha256": digest,
            "byte_size": external.stat().st_size,
            "license": "MIT",
            "species": "human",
            "tissue": "blood",
            "availability": "available",
            "training_study_provenance": {
                "summary": "unit test",
                "source_studies": [{"study_id": "t", "role": "evaluation_label_map"}],
                "overlap_audit_status": "not_applicable",
                "eligible_for_primary_holdout_track": True,
            },
        }
    )
    policy = validate_storage_policy(
        {
            "schema_version": "celltypepilot.asset-storage-policy.v1",
            "policy_id": "test",
            "object_store_uri_template": "s3://bucket/{object_key}",
            "cdn_url_template": "https://cdn.example/{object_key}",
            "object_key_template": "v1/{kind}/{asset_id}/{version}/sha256/{sha256}",
            "immutability": {
                "content_addressed": True,
                "overwrite_policy": "deny",
                "versioning": "test",
            },
            "forbidden_write_globs": ["benchmarks/**/runs/**"],
        }
    )
    row = materialize_source_to_object_cache(asset, catalog_root, policy)
    assert row["status"] == "materialized"
    target = Path(row["target"])
    assert target.is_file()
    assert "runs" not in target.parts
    assert target.read_bytes() == external.read_bytes()


def test_eligible_reference_requires_audit():
    with pytest.raises(AssetCatalogError, match="primary-holdout-eligible"):
        validate_asset_record(
            {
                "asset_id": "azimuth-bad",
                "kind": "azimuth_reference",
                "version": "v1",
                "url": "https://example.com/ref",
                "sha256": "a" * 64,
                "license": "x",
                "species": "human",
                "tissue": "lung",
                "availability": "blocked_overlap_audit",
                "training_study_provenance": {
                    "summary": "bad",
                    "source_studies": [{"study_id": "s", "role": "ref"}],
                    "overlap_audit_status": "not_audited",
                    "eligible_for_primary_holdout_track": True,
                },
            }
        )


def test_duplicate_asset_version_rejected():
    base = {
        "asset_id": "cellxgene-x",
        "kind": "cellxgene_dataset",
        "version": "v1",
        "url": "https://example.com/x.h5ad",
        "sha256": "b" * 64,
        "license": "CC-BY-4.0",
        "species": "human",
        "tissue": "lung",
        "availability": "source_available",
        "training_study_provenance": {
            "summary": "query",
            "source_studies": [{"study_id": "s", "role": "evaluation_query_asset"}],
            "overlap_audit_status": "not_applicable",
            "eligible_for_primary_holdout_track": True,
        },
    }
    with pytest.raises(AssetCatalogError, match="Duplicate"):
        validate_asset_catalog(
            {
                "schema_version": "celltypepilot.immutable-asset-catalog.v1",
                "catalog_id": "t",
                "frozen_at_utc": "2026-08-09T00:00:00Z",
                "assets": [base, dict(base)],
            }
        )
