"""Signed, versioned, data-only extension packs for the plugin ecosystem.

Packs are community/tissue/disease contributions of marker and state knowledge.
They may never contain executable code. Signatures bind license, provenance,
ontology identifiers, and content hashes. Runtime still routes every pack marker
through ordinary DE evidence, critic, abstention, and conflict gates.
"""

from __future__ import annotations

import base64
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .pack_manager import (
    ATLAS_FILE,
    PACK_SCHEMA_VERSION,
    STATE_FILE,
    PackError,
    _file_sha256,
    read_pack_manifest,
    validate_pack,
)

# Data-only allowlist: no executables, no scripts.
ALLOWED_DATA_NAMES = frozenset(
    {
        "pack.json",
        ATLAS_FILE,
        STATE_FILE,
        "LICENSE",
        "LICENSE.txt",
        "README.md",
        "ontology_map.json",
        "reference_manifest.json",
        "pack.sig.json",
    }
)
FORBIDDEN_SUFFIXES = frozenset(
    {
        ".py",
        ".pyc",
        ".pyo",
        ".so",
        ".dll",
        ".dylib",
        ".exe",
        ".bat",
        ".cmd",
        ".ps1",
        ".sh",
        ".js",
        ".mjs",
        ".wasm",
        ".r",
        ".R",
        ".ipynb",
        ".jar",
        ".class",
    }
)

PACK_KINDS = ("evidence", "reference", "mixed")
SIGNATURE_SCHEMA = "celltypepilot.pack-signature.v1"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def scan_pack_for_code(pack_dir: str | Path) -> list[str]:
    """Reject any non-data / executable artifacts (fail closed)."""
    root = Path(pack_dir)
    issues: list[str] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(root).as_posix()
        if rel.startswith(".git/"):
            continue
        if path.name == "_installed.json":
            continue
        suffix = path.suffix.lower()
        if suffix in FORBIDDEN_SUFFIXES or path.name.endswith(".py"):
            issues.append(f"forbidden executable/code artifact: {rel}")
            continue
        # Only allow known data basenames at pack root (no nested code trees).
        if path.parent != root:
            issues.append(f"nested files are not allowed in data-only packs: {rel}")
            continue
        if path.name not in ALLOWED_DATA_NAMES:
            issues.append(f"unknown data file {rel!r}; allowed: {sorted(ALLOWED_DATA_NAMES)}")
    return issues


def compute_pack_content_hashes(pack_dir: str | Path) -> dict[str, str]:
    root = Path(pack_dir)
    hashes: dict[str, str] = {}
    for name in sorted(ALLOWED_DATA_NAMES):
        path = root / name
        if path.is_file() and name != "pack.sig.json":
            hashes[name] = _file_sha256(path)
    return hashes


def compute_pack_fingerprint(pack_dir: str | Path) -> dict[str, Any]:
    """Canonical fingerprint used for signing and install verification."""
    root = Path(pack_dir)
    manifest = read_pack_manifest(root)
    content_sha256 = compute_pack_content_hashes(root)
    # Fingerprint excludes the signature file itself.
    body = {
        "schema_version": PACK_SCHEMA_VERSION,
        "name": manifest.get("name"),
        "version": manifest.get("version"),
        "pack_kind": manifest.get("pack_kind", "evidence"),
        "license": manifest.get("license", manifest.get("license_spdx", "")),
        "license_tier": manifest.get("license_tier", "community"),
        "species": manifest.get("species", ["human"]),
        "tissues": manifest.get("tissues", []),
        "diseases": manifest.get("diseases", []),
        "ontology": manifest.get("ontology", {}),
        "provenance": manifest.get("provenance", {}),
        "files": manifest.get("files", []),
        "content_sha256": content_sha256,
    }
    canonical = json.dumps(body, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    return {
        "fingerprint_sha256": digest,
        "canonical_payload": body,
        "canonical_json": canonical,
    }


def enhance_manifest_fields(manifest: dict) -> list[str]:
    """Validate ecosystem fields for community packs (optional but recommended)."""
    issues: list[str] = []
    kind = manifest.get("pack_kind", "evidence")
    if kind not in PACK_KINDS:
        issues.append(f"pack_kind must be one of {list(PACK_KINDS)}, got {kind!r}")
    if "license" not in manifest and "license_spdx" not in manifest:
        issues.append("license or license_spdx is required for ecosystem packs")
    diseases = manifest.get("diseases", [])
    if diseases is not None and not isinstance(diseases, list):
        issues.append("diseases must be a list of disease/context labels when present")
    ontology = manifest.get("ontology", {})
    if ontology is not None and not isinstance(ontology, dict):
        issues.append("ontology must be an object (e.g. {cell_ontology: version})")
    provenance = manifest.get("provenance", {})
    if provenance is not None and not isinstance(provenance, dict):
        issues.append("provenance must be an object with source/curator/date fields")
    return issues


def validate_pack_ecosystem(
    pack_dir: str | Path,
    *,
    require_signature: bool = False,
    hmac_secret: str | None = None,
    public_key_pem: str | None = None,
) -> list[str]:
    """Full ecosystem validation: data-only + schema + provenance + optional signature."""
    pack_dir = Path(pack_dir)
    issues = scan_pack_for_code(pack_dir)
    issues.extend(validate_pack(pack_dir))
    try:
        manifest = read_pack_manifest(pack_dir)
        issues.extend(enhance_manifest_fields(manifest))
    except PackError as exc:
        issues.append(str(exc))
        return issues

    sig_path = pack_dir / "pack.sig.json"
    if require_signature and not sig_path.is_file():
        issues.append("pack.sig.json missing (require_signature=True)")
    if sig_path.is_file():
        verify = verify_pack_signature(
            pack_dir, hmac_secret=hmac_secret, public_key_pem=public_key_pem
        )
        if not verify.get("valid"):
            issues.append(f"signature invalid: {verify.get('reason')}")
    return issues


def _sign_bytes_hmac(message: bytes, secret: bytes) -> str:
    import hmac

    return hmac.new(secret, message, hashlib.sha256).hexdigest()


def _try_rsa_sign(message: bytes, private_key_pem: str) -> str | None:
    try:
        from cryptography.hazmat.primitives import hashes, serialization
        from cryptography.hazmat.primitives.asymmetric import padding
    except ImportError:
        return None
    try:
        key = serialization.load_pem_private_key(private_key_pem.encode("utf-8"), password=None)
        signature = key.sign(
            message,
            padding.PKCS1v15(),
            hashes.SHA256(),
        )
        return base64.b64encode(signature).decode("ascii")
    except Exception:
        return None


def _try_rsa_verify(message: bytes, signature_b64: str, public_key_pem: str) -> bool:
    try:
        from cryptography.hazmat.primitives import hashes, serialization
        from cryptography.hazmat.primitives.asymmetric import padding
    except ImportError:
        return False
    try:
        key = serialization.load_pem_public_key(public_key_pem.encode("utf-8"))
        key.verify(
            base64.b64decode(signature_b64.encode("ascii")),
            message,
            padding.PKCS1v15(),
            hashes.SHA256(),
        )
        return True
    except Exception:
        return False


def sign_pack(
    pack_dir: str | Path,
    *,
    private_key_pem: str | None = None,
    hmac_secret: str | None = None,
    signer: str = "local-curator",
) -> dict[str, Any]:
    """Write pack.sig.json binding content hashes + license/provenance/ontology.

    Prefer RSA private key PEM when cryptography is available; otherwise HMAC
    secret (community/dev). Either way packs remain data-only.
    """
    pack_dir = Path(pack_dir)
    issues = scan_pack_for_code(pack_dir)
    if issues:
        raise PackError("Cannot sign pack with code/unknown files: " + "; ".join(issues[:5]))
    base_issues = validate_pack(pack_dir)
    if base_issues:
        raise PackError("Cannot sign invalid pack: " + "; ".join(base_issues[:5]))

    fingerprint = compute_pack_fingerprint(pack_dir)
    message = fingerprint["canonical_json"].encode("utf-8")
    algorithm = None
    signature = None
    if private_key_pem:
        signature = _try_rsa_sign(message, private_key_pem)
        if signature:
            algorithm = "rsa-sha256"
    if signature is None:
        secret = (hmac_secret or "celltypepilot-dev-pack-hmac").encode("utf-8")
        signature = _sign_bytes_hmac(message, secret)
        algorithm = "hmac-sha256"

    payload = {
        "schema_version": SIGNATURE_SCHEMA,
        "signed_at_utc": _utc_now(),
        "signer": signer,
        "algorithm": algorithm,
        "fingerprint_sha256": fingerprint["fingerprint_sha256"],
        "content_sha256": fingerprint["canonical_payload"]["content_sha256"],
        "signature": signature,
        "manifest_name": fingerprint["canonical_payload"]["name"],
        "manifest_version": fingerprint["canonical_payload"]["version"],
        "license": fingerprint["canonical_payload"].get("license"),
        "ontology": fingerprint["canonical_payload"].get("ontology"),
        "provenance": fingerprint["canonical_payload"].get("provenance"),
        "note": (
            "Signature binds data files only. Runtime still applies marker evidence, "
            "critic, abstention, and conflict gates to every pack relationship."
        ),
    }
    path = pack_dir / "pack.sig.json"
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return payload


def verify_pack_signature(
    pack_dir: str | Path,
    *,
    public_key_pem: str | None = None,
    hmac_secret: str | None = None,
) -> dict[str, Any]:
    pack_dir = Path(pack_dir)
    sig_path = pack_dir / "pack.sig.json"
    if not sig_path.is_file():
        return {"valid": False, "reason": "pack.sig.json missing", "status": "unsigned"}
    try:
        sig = json.loads(sig_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        return {"valid": False, "reason": f"signature unreadable: {exc}", "status": "invalid"}

    fingerprint = compute_pack_fingerprint(pack_dir)
    if sig.get("fingerprint_sha256") != fingerprint["fingerprint_sha256"]:
        return {
            "valid": False,
            "reason": "fingerprint mismatch (content or manifest changed after signing)",
            "status": "tampered",
            "expected": fingerprint["fingerprint_sha256"],
            "recorded": sig.get("fingerprint_sha256"),
        }
    message = fingerprint["canonical_json"].encode("utf-8")
    algorithm = str(sig.get("algorithm", ""))
    signature = str(sig.get("signature", ""))
    ok = False
    if algorithm == "rsa-sha256" and public_key_pem:
        ok = _try_rsa_verify(message, signature, public_key_pem)
    elif algorithm == "hmac-sha256":
        secret = (hmac_secret or "celltypepilot-dev-pack-hmac").encode("utf-8")
        ok = _sign_bytes_hmac(message, secret) == signature
    elif algorithm == "rsa-sha256" and not public_key_pem:
        # Without public key, accept fingerprint match only as "fingerprint_ok".
        return {
            "valid": True,
            "status": "fingerprint_ok_rsa_unverified",
            "reason": "fingerprint matches; RSA public key not provided for crypto verify",
            "algorithm": algorithm,
            "signer": sig.get("signer"),
        }
    else:
        return {
            "valid": False,
            "reason": f"unsupported algorithm {algorithm!r}",
            "status": "invalid",
        }

    return {
        "valid": ok,
        "status": "valid" if ok else "invalid_signature",
        "reason": "ok" if ok else "cryptographic signature verification failed",
        "algorithm": algorithm,
        "signer": sig.get("signer"),
        "signed_at_utc": sig.get("signed_at_utc"),
    }


def scaffold_pack(
    output_dir: str | Path,
    *,
    name: str,
    version: str = "0.1.0",
    tissues: list[str] | None = None,
    diseases: list[str] | None = None,
    pack_kind: str = "evidence",
    license_spdx: str = "CC-BY-4.0",
) -> Path:
    """Create a data-only pack skeleton for community contribution."""
    if not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", name):
        raise PackError(f"name must match [a-z0-9][a-z0-9_-]*, got {name!r}")
    if pack_kind not in PACK_KINDS:
        raise PackError(f"pack_kind must be one of {list(PACK_KINDS)}")
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    tissues = tissues or ["general"]
    diseases = diseases or []
    manifest = {
        "schema_version": PACK_SCHEMA_VERSION,
        "name": name,
        "version": version,
        "pack_kind": pack_kind,
        "description": f"Data-only {pack_kind} pack for {', '.join(tissues)}",
        "species": ["human"],
        "tissues": tissues,
        "diseases": diseases,
        "license": license_spdx,
        "license_spdx": license_spdx,
        "license_tier": "community",
        "files": ["marker_atlas.json"],
        "ontology": {
            "cell_ontology": "CL",
            "note": "All cl_id values must be valid Cell Ontology identifiers",
        },
        "provenance": {
            "curator": "community",
            "created_at_utc": _utc_now(),
            "sources": [],
            "note": "Every marker edge needs sources + verification_status",
        },
        "runtime_gates": [
            "marker_evidence",
            "critic",
            "abstention",
            "conflict_detection",
        ],
        "code_policy": "data_only_no_executables",
    }
    (root / "pack.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (root / "LICENSE").write_text(f"{license_spdx}\n", encoding="utf-8")
    (root / "README.md").write_text(
        f"# {name}\n\nData-only CellTypePilot pack. No code allowed.\n"
        f"Gates: marker evidence, critic, abstention, conflicts.\n",
        encoding="utf-8",
    )
    # Minimal empty atlas shell (user fills markers with provenance).
    atlas = {
        "version": version,
        "schema_version": "celltypepilot.marker-atlas.v2",
        "description": f"{name} marker atlas",
        "tissues": {tissues[0]: {"cell_types": {}}},
        "provenance_policy": "pack_edges_require_sources",
    }
    (root / ATLAS_FILE).write_text(
        json.dumps(atlas, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return root
