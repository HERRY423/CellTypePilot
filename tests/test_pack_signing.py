"""Data-only signed pack ecosystem tests."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from celltypepilot.pack_manager import PackError, install_pack, validate_pack
from celltypepilot.pack_signing import (
    scan_pack_for_code,
    scaffold_pack,
    sign_pack,
    validate_pack_ecosystem,
    verify_pack_signature,
)

REPO = Path(__file__).resolve().parents[1]
SCRATCH = REPO / "scratch" / "pytest_temp" / "pack_signing"


def _evidence_edge(gene: str, polarity: str, tissue: str, version: str) -> dict:
    return {
        "gene": gene,
        "polarity": polarity,
        "species": ["human"],
        "tissue": tissue,
        "state": "baseline_or_unspecified",
        "atlas_version": version,
        "sources": [
            {
                "source_id": "test",
                "name": "Test",
                "pmid": "12345678",
                "doi": "10.0/test",
                "url": "https://example.org",
            }
        ],
        "evidence_scope": "database_level_source",
        "verification_status": "aggregate_source_only_not_edge_verified",
    }


def _fill_min_atlas(pack_dir: Path, version: str = "1.0.0") -> None:
    atlas = {
        "version": version,
        "schema_version": "celltypepilot.marker-atlas.v2",
        "tissues": {
            "gut": {
                "cell_types": {
                    "test_epithelial": {
                        "cl_id": "CL:0000066",
                        "positive_markers": ["EPCAM"],
                        "negative_markers": ["PTPRC"],
                        "marker_evidence": [
                            _evidence_edge("EPCAM", "positive", "gut", version),
                            _evidence_edge("PTPRC", "negative", "gut", version),
                        ],
                    }
                }
            }
        },
    }
    (pack_dir / "marker_atlas.json").write_text(
        json.dumps(atlas, indent=2) + "\n", encoding="utf-8"
    )
    manifest = json.loads((pack_dir / "pack.json").read_text(encoding="utf-8"))
    manifest["version"] = version
    manifest["tissues"] = ["gut"]
    manifest["diseases"] = ["IBD"]
    (pack_dir / "pack.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def test_scaffold_sign_verify_roundtrip():
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    pack_dir = SCRATCH / "ibd-pack"
    scaffold_pack(
        pack_dir,
        name="ibd-pack",
        version="1.0.0",
        tissues=["gut"],
        diseases=["IBD"],
        pack_kind="evidence",
    )
    _fill_min_atlas(pack_dir)
    assert not scan_pack_for_code(pack_dir)
    sig = sign_pack(pack_dir, signer="tester", hmac_secret="unit-test-secret")
    assert sig["algorithm"] == "hmac-sha256"
    verify = verify_pack_signature(pack_dir, hmac_secret="unit-test-secret")
    assert verify["valid"] is True
    issues = validate_pack_ecosystem(
        pack_dir, require_signature=True, hmac_secret="unit-test-secret"
    )
    assert issues == [], issues


def test_code_file_rejected():
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    pack_dir = SCRATCH / "bad-pack"
    scaffold_pack(pack_dir, name="bad-pack", tissues=["blood"])
    (pack_dir / "evil.py").write_text("print('nope')\n", encoding="utf-8")
    issues = scan_pack_for_code(pack_dir)
    assert any("forbidden" in i for i in issues)
    with pytest.raises(PackError):
        sign_pack(pack_dir, hmac_secret="x")


def test_tamper_breaks_signature():
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    pack_dir = SCRATCH / "tamper-pack"
    scaffold_pack(pack_dir, name="tamper-pack", tissues=["blood"])
    _fill_min_atlas(pack_dir, version="2.0.0")
    sign_pack(pack_dir, hmac_secret="secret")
    atlas = json.loads((pack_dir / "marker_atlas.json").read_text(encoding="utf-8"))
    atlas["version"] = "2.0.1-tampered"
    (pack_dir / "marker_atlas.json").write_text(json.dumps(atlas) + "\n", encoding="utf-8")
    verify = verify_pack_signature(pack_dir, hmac_secret="secret")
    assert verify["valid"] is False
    assert "fingerprint" in verify.get("reason", "").lower() or verify.get("status") == "tampered"


def test_install_rejects_code_pack(monkeypatch):
    if SCRATCH.exists():
        shutil.rmtree(SCRATCH)
    pack_dir = SCRATCH / "codey"
    install_root = SCRATCH / "installed"
    monkeypatch.setenv("CELLTYPEPILOT_PACKS_DIR", str(install_root))
    scaffold_pack(pack_dir, name="codey", tissues=["blood"])
    _fill_min_atlas(pack_dir)
    (pack_dir / "hack.sh").write_text("echo x\n", encoding="utf-8")
    with pytest.raises(PackError):
        install_pack(str(pack_dir), trust="atlas")
