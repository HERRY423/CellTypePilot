"""Tests for the data-only extension pack mechanism (pack spec v1).

Covers: manifest validation, fail-closed install, trust tiers, license
gating, content-hash tamper detection, atlas merging precedence, state
atlas trust mapping, and orchestrator/CLI integration.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from celltypepilot import pack_manager
from celltypepilot.pack_manager import (
    PACK_SCHEMA_VERSION,
    PackError,
    collect_pack_state_definitions_input,
    install_pack,
    list_installed_packs,
    merge_marker_atlas,
    pack_manifest_parameters,
    remove_pack,
    resolve_extension_packs,
    validate_pack,
    validate_pack_manifest,
)

PACK_VERSION = "test-pack-1.0.0"


# ──────────────────────────────────────────────
# Fixtures: data-only pack builders
# ──────────────────────────────────────────────
def _marker_evidence_record(gene: str, polarity: str, tissue: str) -> dict:
    return {
        "gene": gene,
        "polarity": polarity,
        "species": "human",
        "tissue": tissue,
        "state": "identity",
        "atlas_version": PACK_VERSION,
        "sources": [
            {
                "source_id": "test_source",
                "name": "Test Source",
                "pmid": "12345678",
                "doi": "10.1000/test",
                "url": "https://example.org/test",
            }
        ],
        "evidence_scope": "cell_type",
        "verification_status": "aggregate_source_only_not_edge_verified",
    }


def _pack_marker_atlas(tissue: str = "test_tissue") -> dict:
    return {
        "version": PACK_VERSION,
        "tissues": {
            tissue: {
                "cell_types": {
                    "test_cell": {
                        "cl_id": "CL:0000000",
                        "positive_markers": ["TESTG1", "TESTG2"],
                        "negative_markers": ["TESTN1"],
                        "marker_evidence": [
                            _marker_evidence_record("TESTG1", "positive", tissue),
                            _marker_evidence_record("TESTG2", "positive", tissue),
                            _marker_evidence_record("TESTN1", "negative", tissue),
                        ],
                    }
                }
            }
        },
    }


def _provenanceless_atlas(tissue: str = "draft_tissue") -> dict:
    """An atlas that fails validate_atlas_provenance (no marker_evidence)."""
    return {
        "version": PACK_VERSION,
        "tissues": {
            tissue: {
                "cell_types": {
                    "draft_cell": {
                        "cl_id": "CL:0000001",
                        "positive_markers": ["DRAFTG1"],
                        "negative_markers": [],
                        "marker_evidence": [],
                    }
                }
            }
        },
    }


def _pack_state_atlas() -> dict:
    return {
        "schema_version": "celltypepilot.state-atlas.v1",
        "version": PACK_VERSION,
        "sources": {
            "test_source": {
                "name": "Test Source",
                "pmid": "12345678",
                "doi": "10.1000/test",
                "url": "https://example.org/test",
            }
        },
        "states": {
            "test_state": {
                "species": ["human"],
                "tissues": ["general"],
                "parent_cell_types": ["test_cell"],
                "positive_markers": ["TESTS1"],
                "negative_markers": [],
                "marker_evidence": [
                    {
                        "gene": "TESTS1",
                        "polarity": "positive",
                        "state": "test_state",
                        "species": ["human"],
                        "tissues": ["general"],
                        "atlas_version": PACK_VERSION,
                        "source_ids": ["test_source"],
                        "verification_status": "aggregate_source_only_not_edge_verified",
                    }
                ],
            }
        },
    }


def _write_pack(
    root: Path,
    name: str = "testpack",
    tissues: list[str] | None = None,
    species: list[str] | None = None,
    license_tier: str = "community",
    marker_atlas: dict | None = "default",
    state_atlas: dict | None = None,
    extra_files: list[str] | None = None,
) -> Path:
    pack_dir = root / name
    pack_dir.mkdir(parents=True, exist_ok=True)
    files = []
    if marker_atlas == "default":
        marker_atlas = _pack_marker_atlas()
    if marker_atlas is not None:
        (pack_dir / "marker_atlas.json").write_text(
            json.dumps(marker_atlas, indent=2), encoding="utf-8"
        )
        files.append("marker_atlas.json")
    if state_atlas is not None:
        (pack_dir / "state_atlas.json").write_text(
            json.dumps(state_atlas, indent=2), encoding="utf-8"
        )
        files.append("state_atlas.json")
    for filename in extra_files or []:
        (pack_dir / filename).write_text("{}", encoding="utf-8")
    manifest = {
        "schema_version": PACK_SCHEMA_VERSION,
        "name": name,
        "version": PACK_VERSION,
        "description": "Test pack",
        "species": species or ["human"],
        "tissues": tissues or ["test_tissue"],
        "license_tier": license_tier,
        "files": files or ["marker_atlas.json"],
    }
    (pack_dir / "pack.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return pack_dir


@pytest.fixture
def pack_env(tmp_path, monkeypatch):
    """Redirect the user pack directory into the test workspace."""
    install_root = tmp_path / "installed_packs"
    monkeypatch.setenv(pack_manager.PACKS_ENV_VAR, str(install_root))
    source_root = tmp_path / "pack_sources"
    source_root.mkdir()
    return {"install_root": install_root, "source_root": source_root}


# ──────────────────────────────────────────────
# Manifest and pack validation
# ──────────────────────────────────────────────
def test_validate_pack_manifest_rejects_schema_name_and_files():
    issues = validate_pack_manifest(
        {
            "schema_version": "wrong.v0",
            "name": "Bad Name",
            "version": "",
            "species": ["fruitfly"],
            "license_tier": "platinum",
            "files": [],
        }
    )
    joined = "; ".join(issues)
    assert "schema_version" in joined
    assert "name" in joined
    assert "version" in joined
    assert "species" in joined
    assert "license_tier" in joined
    assert "files" in joined


def test_validate_pack_manifest_accepts_valid_manifest():
    manifest = {
        "schema_version": PACK_SCHEMA_VERSION,
        "name": "tme-pack",
        "version": "1.0.0",
        "species": ["human", "mouse"],
        "tissues": ["tumor_microenvironment"],
        "license_tier": "community",
        "files": ["marker_atlas.json"],
    }
    assert validate_pack_manifest(manifest) == []


def test_validate_pack_missing_manifest(tmp_path):
    issues = validate_pack(tmp_path / "nonexistent")
    assert issues and "pack.json not found" in issues[0]


def test_validate_pack_provenance_gate(pack_env):
    source = _write_pack(pack_env["source_root"], marker_atlas=_provenanceless_atlas())
    issues = validate_pack(source)
    assert any("evidence records differ" in issue for issue in issues)

    clean = _write_pack(pack_env["source_root"], name="cleanpack")
    assert validate_pack(clean) == []


# ──────────────────────────────────────────────
# Install: fail closed, trust tiers, license
# ──────────────────────────────────────────────
def test_install_fails_closed_without_provenance(pack_env):
    source = _write_pack(pack_env["source_root"], marker_atlas=_provenanceless_atlas())
    with pytest.raises(PackError, match="fail closed"):
        install_pack(str(source))
    assert not (pack_env["install_root"] / "testpack").exists()


def test_install_hypothesis_trust_downgrades_pack(pack_env):
    source = _write_pack(pack_env["source_root"], marker_atlas=_provenanceless_atlas())
    summary = install_pack(str(source), trust="hypothesis")
    assert summary["trust"] == "hypothesis"
    metadata = json.loads(
        (pack_env["install_root"] / "testpack" / "_installed.json").read_text(encoding="utf-8")
    )
    assert metadata["trust"] == "hypothesis"
    assert metadata["validation_issues_at_install"]
    assert metadata["content_sha256"]["marker_atlas.json"]

    records, warnings = resolve_extension_packs(["testpack"], "human")
    assert warnings == []
    assert records[0]["trust"] == "hypothesis"
    merged, _ = merge_marker_atlas({"version": "base", "tissues": {}}, records, "human")
    cell = merged["tissues"]["draft_tissue"]["cell_types"]["draft_cell"]
    assert cell["context_origin"] is True
    assert cell["context_review_status"] == "draft"
    assert cell["atlas_positive_markers"] == []
    assert cell["context_positive_markers"] == ["DRAFTG1"]


def test_install_rejects_invalid_trust_and_missing_source(pack_env):
    source = _write_pack(pack_env["source_root"])
    with pytest.raises(PackError, match="trust must be"):
        install_pack(str(source), trust="superuser")
    with pytest.raises(PackError, match="not found"):
        install_pack(str(pack_env["source_root"] / "missing"))


def test_install_requires_force_to_reinstall(pack_env):
    source = _write_pack(pack_env["source_root"])
    install_pack(str(source))
    with pytest.raises(PackError, match="already installed"):
        install_pack(str(source))
    install_pack(str(source), force=True)


def test_install_reserved_first_party_name_is_rejected(pack_env):
    source = _write_pack(pack_env["source_root"], name="premium")
    with pytest.raises(PackError, match="reserved"):
        install_pack(str(source))


def test_license_gate_blocks_academic_pack_on_free_tier(pack_env):
    source = _write_pack(pack_env["source_root"], license_tier="academic")
    with pytest.raises(PackError, match="academic license"):
        install_pack(str(source))

    # With an academic license the same pack installs and resolves.
    from celltypepilot import license_manager
    from celltypepilot.license_manager import LicenseTier

    class _AcademicLicense:
        tier = LicenseTier.ACADEMIC

    original = license_manager.load_license
    license_manager.load_license = lambda *args, **kwargs: _AcademicLicense()
    try:
        summary = install_pack(str(source))
        assert summary["name"] == "testpack"
        records, _ = resolve_extension_packs(["testpack"], "human")
        assert records[0]["license_tier"] == "academic"
    finally:
        license_manager.load_license = original


# ──────────────────────────────────────────────
# Resolve: unknown names, species scope, tampering
# ──────────────────────────────────────────────
def test_resolve_unknown_pack_fails_closed(pack_env):
    _write_pack(pack_env["source_root"])
    install_pack(str(pack_env["source_root"] / "testpack"))
    with pytest.raises(PackError, match="not installed"):
        resolve_extension_packs(["ghost"], "human")


def test_resolve_species_mismatch_skips_with_warning(pack_env):
    source = _write_pack(pack_env["source_root"], species=["mouse"])
    install_pack(str(source))
    records, warnings = resolve_extension_packs(["testpack"], "human")
    assert records == []
    assert any("skipped" in warning for warning in warnings)
    records, warnings = resolve_extension_packs(["testpack"], "mouse")
    assert len(records) == 1 and warnings == []


def test_resolve_detects_tampered_pack_file(pack_env):
    source = _write_pack(pack_env["source_root"])
    install_pack(str(source))
    target = pack_env["install_root"] / "testpack" / "marker_atlas.json"
    payload = json.loads(target.read_text(encoding="utf-8"))
    payload["tissues"]["test_tissue"]["cell_types"]["test_cell"]["positive_markers"].append(
        "SMUGGLED"
    )
    target.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(PackError, match="modified after install"):
        resolve_extension_packs(["testpack"], "human")


# ──────────────────────────────────────────────
# List / remove / first-party discovery
# ──────────────────────────────────────────────
def test_list_includes_first_party_premium_pack(pack_env):
    entries = list_installed_packs()
    premium = [entry for entry in entries if entry["name"] == "premium"]
    assert len(premium) == 1
    assert premium[0]["origin"] == "first_party"
    assert premium[0]["license_tier"] == "academic"
    assert "tumor_microenvironment" in premium[0]["tissues"]


def test_list_and_remove_user_pack(pack_env):
    source = _write_pack(pack_env["source_root"])
    install_pack(str(source))
    names = {entry["name"] for entry in list_installed_packs()}
    assert "testpack" in names

    remove_pack("testpack")
    assert "testpack" not in {entry["name"] for entry in list_installed_packs()}
    with pytest.raises(PackError, match="not installed"):
        remove_pack("testpack")


def test_remove_first_party_pack_is_rejected(pack_env):
    with pytest.raises(PackError, match="first-party"):
        remove_pack("premium")


# ──────────────────────────────────────────────
# Atlas merge semantics
# ──────────────────────────────────────────────
def _installed_record(pack_env, **pack_kwargs) -> dict:
    source = _write_pack(pack_env["source_root"], **pack_kwargs)
    install_pack(str(source), force=True)
    records, _ = resolve_extension_packs(["testpack"], "human")
    return records[0]


def test_merge_adds_new_tissue_and_shadows_collision(pack_env):
    record = _installed_record(pack_env)
    base = {
        "version": "builtin",
        "tissues": {
            "blood": {
                "cell_types": {"T cell": {"cl_id": "CL:0000084", "positive_markers": ["CD3E"]}}
            }
        },
    }
    merged, warnings = merge_marker_atlas(base, [record], "human")
    assert "test_tissue" in merged["tissues"]
    assert "test_cell" in merged["tissues"]["test_tissue"]["cell_types"]
    assert warnings == []
    # Base atlas must not be mutated in place.
    assert "test_tissue" not in base["tissues"]

    # Collision: pack tries to add a cell type the base already owns.
    colliding = _pack_marker_atlas(tissue="blood")
    colliding["tissues"]["blood"]["cell_types"] = {
        "T cell": {
            "cl_id": "CL:0009999",
            "positive_markers": ["FAKE"],
            "negative_markers": [],
            "marker_evidence": [
                _marker_evidence_record("FAKE", "positive", "blood"),
            ],
        }
    }
    record["marker_atlas"] = colliding
    merged, warnings = merge_marker_atlas(base, [record], "human")
    assert merged["tissues"]["blood"]["cell_types"]["T cell"]["cl_id"] == "CL:0000084"
    assert any("shadowed" in warning for warning in warnings)


def test_merge_converts_genes_for_mouse(pack_env):
    record = _installed_record(pack_env)
    # Conversion applies to purely alphabetic symbols (built-in convention).
    alphabetic = _pack_marker_atlas()
    cell = alphabetic["tissues"]["test_tissue"]["cell_types"]["test_cell"]
    cell["positive_markers"] = ["TESTGA", "TESTGB"]
    cell["marker_evidence"] = [
        _marker_evidence_record("TESTGA", "positive", "test_tissue"),
        _marker_evidence_record("TESTGB", "positive", "test_tissue"),
        _marker_evidence_record("TESTN1", "negative", "test_tissue"),
    ]
    record["marker_atlas"] = alphabetic
    merged, _ = merge_marker_atlas({"version": "builtin", "tissues": {}}, [record], "mouse")
    converted = merged["tissues"]["test_tissue"]["cell_types"]["test_cell"]
    assert converted["positive_markers"] == ["Testga", "Testgb"]
    assert converted["negative_markers"] == ["TESTN1"]  # digit kept as-is


def test_state_pack_trust_maps_to_review_status(pack_env):
    source = _write_pack(
        pack_env["source_root"], state_atlas=_pack_state_atlas(), marker_atlas=None
    )
    install_pack(str(source))
    records, _ = resolve_extension_packs(["testpack"], "human")
    entries = collect_pack_state_definitions_input(records)
    assert len(entries) == 1

    from celltypepilot.state_scorer import load_state_definitions

    definitions = load_state_definitions("human", "blood", pack_states=entries)
    pack_defs = [item for item in definitions if item["source"].startswith("pack:testpack")]
    assert len(pack_defs) == 1
    assert pack_defs[0]["review_status"] == "reviewed"

    entries[0]["trust"] = "hypothesis"
    definitions = load_state_definitions("human", "blood", pack_states=entries)
    pack_defs = [item for item in definitions if item["source"].startswith("pack:testpack")]
    assert pack_defs[0]["review_status"] == "draft"


def test_invalid_pack_state_atlas_fails_closed(pack_env):
    from celltypepilot.state_scorer import StateScoringError, load_state_definitions

    broken = _pack_state_atlas()
    broken["states"]["test_state"]["marker_evidence"] = []  # evidence/relationship mismatch
    with pytest.raises(StateScoringError, match="validation failed"):
        load_state_definitions(
            "human",
            "blood",
            pack_states=[
                {
                    "pack_name": "testpack",
                    "pack_version": PACK_VERSION,
                    "trust": "atlas",
                    "atlas": broken,
                }
            ],
        )


def test_pack_manifest_parameters_record_hashes(pack_env):
    record = _installed_record(pack_env)
    params = pack_manifest_parameters([record])
    assert params[0]["name"] == "testpack"
    assert params[0]["trust"] == "atlas"
    assert len(params[0]["files"]["marker_atlas.json"]) == 64


# ──────────────────────────────────────────────
# Orchestrator + CLI integration
# ──────────────────────────────────────────────
def test_pipeline_unknown_pack_fails_closed(h5ad_path, tmp_output_dir, pack_env):
    from celltypepilot.orchestrator import PipelineError, run_annotation_pipeline

    with pytest.raises(PipelineError, match="Extension pack safety gate failed"):
        run_annotation_pipeline(
            input_path=str(h5ad_path),
            cluster_key="leiden",
            output_dir=str(tmp_output_dir / "out"),
            species="human",
            tissue="blood",
            no_figures=True,
            packs=["ghost"],
        )


def test_pipeline_records_pack_in_manifest(h5ad_path, tmp_output_dir, pack_env):
    from celltypepilot.orchestrator import run_annotation_pipeline

    source = _write_pack(pack_env["source_root"])
    install_pack(str(source))
    result = run_annotation_pipeline(
        input_path=str(h5ad_path),
        cluster_key="leiden",
        output_dir=str(tmp_output_dir / "out"),
        species="human",
        tissue="blood",
        no_figures=True,
        packs=["testpack"],
    )
    pack_entries = result["manifest"]["parameters"]["extension_packs"]
    assert pack_entries[0]["name"] == "testpack"
    assert pack_entries[0]["trust"] == "atlas"

    written = json.loads((tmp_output_dir / "out" / "manifest.json").read_text(encoding="utf-8"))
    assert written["parameters"]["extension_packs"][0]["name"] == "testpack"


def test_cli_pack_install_list_validate_remove(pack_env):
    from celltypepilot.cli import app

    runner = CliRunner()
    source = _write_pack(pack_env["source_root"])

    result = runner.invoke(app, ["pack", "install", str(source)])
    assert result.exit_code == 0, result.output
    assert "testpack" in result.output

    result = runner.invoke(app, ["pack", "list"])
    assert result.exit_code == 0
    assert "testpack" in result.output

    result = runner.invoke(app, ["pack", "validate", "testpack"])
    assert result.exit_code == 0, result.output

    result = runner.invoke(app, ["pack", "remove", "testpack"])
    assert result.exit_code == 0, result.output

    result = runner.invoke(app, ["pack", "validate", "testpack"])
    assert result.exit_code != 0
