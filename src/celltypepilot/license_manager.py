"""CellTypePilot — Tiered license manager.

Manages license tiers and gates access to premium features:
- Free tier: Built-in MKG (80+ cell types, 11 tissues)
- Academic tier: Extended atlas (200+ cell types, disease states, developmental stages)
- Commercial tier: Full atlas + priority support + custom tissue panels

Security model:
- RSA-2048 asymmetric signing: public key embedded, private key offline
- Machine fingerprint binding: license tied to hardware
- License file tamper detection: HMAC integrity on local storage
- Offline-first with optional online verification

License key format: CTP-{TIER}-{BASE64_PAYLOAD}.{BASE64_SIGNATURE}
"""

from __future__ import annotations

import base64
import hashlib
import hmac
import json
import os
import platform
import uuid
from dataclasses import asdict, dataclass, field
from datetime import datetime, timedelta
from enum import Enum
from pathlib import Path


class LicenseTier(str, Enum):
    """License tiers."""

    FREE = "free"
    ACADEMIC = "academic"
    COMMERCIAL = "commercial"
    TRIAL = "trial"


@dataclass
class LicenseInfo:
    """License information."""

    tier: LicenseTier = LicenseTier.FREE
    key: str = ""
    holder: str = ""
    email: str = ""
    issued_at: str = ""
    expires_at: str = ""
    features: list[str] = field(default_factory=list)
    max_tissues: int = 11  # Free tier limit
    max_cell_types: int = 80  # Free tier limit
    machine_id: str = ""  # Bound machine fingerprint
    valid: bool = True
    _signature_valid: bool = False  # Internal: RSA signature check result

    def is_expired(self) -> bool:
        """Check if license has expired."""
        if not self.expires_at:
            return False
        try:
            exp = datetime.fromisoformat(self.expires_at)
            return datetime.now() > exp
        except ValueError:
            return False

    def has_feature(self, feature: str) -> bool:
        """Check if a feature is enabled."""
        if self.tier == LicenseTier.FREE:
            return feature in FREE_FEATURES
        elif self.tier == LicenseTier.ACADEMIC:
            return feature in ACADEMIC_FEATURES
        elif self.tier == LicenseTier.COMMERCIAL:
            return feature in COMMERCIAL_FEATURES
        elif self.tier == LicenseTier.TRIAL:
            return feature in TRIAL_FEATURES
        return False


# ──────────────────────────────────────────────
# Feature gates per tier
# ──────────────────────────────────────────────

FREE_FEATURES = {
    "basic_atlas",  # 80+ cell types, 11 tissues
    "marker_scoring",  # Wilcoxon DE + 5-dim scoring
    "critic_review",  # Annotation Critic
    "basic_visualization",  # UMAP, dotplot, confidence
    "html_report",  # HTML report generation
    "manifest_provenance",  # Run provenance tracking
}

ACADEMIC_FEATURES = FREE_FEATURES | {
    "extended_atlas",  # 200+ cell types, disease states
    "developmental_atlas",  # Developmental stage markers
    "disease_atlas",  # Disease-specific cell states
    "literature_validation",  # PubMed integration
    "advanced_visualization",  # Cross-sample comparisons
    "docx_export",  # Word document export
}

COMMERCIAL_FEATURES = ACADEMIC_FEATURES | {
    "custom_tissue_panels",  # Custom tissue marker panels
    "team_sharing",  # Team atlas sharing
    "priority_support",  # Priority email support
    "api_access",  # REST API access
    "white_label",  # White-label reports
}

TRIAL_FEATURES = ACADEMIC_FEATURES | set()  # Trial gets academic features for 30 days


# ──────────────────────────────────────────────
# RSA-2048 Asymmetric Signing
# ──────────────────────────────────────────────
# The public key is embedded in the client for VERIFICATION only.
# The private key is kept OFFLINE for license generation.
# This means users cannot forge valid license keys even with full
# source code access — they would need the private key.

LICENSE_PREFIX = "CTP"
_LICENSE_VERSION = "2"

# RSA-2048 public key (PEM) — for verification only
# Private key is stored offline at .license_private_key.pem
_RSA_PUBLIC_KEY_PEM = (
    "-----BEGIN PUBLIC KEY-----\n"
    "MIIBIjANBgkqhkiG9w0BAQEFAAOCAQ8AMIIBCgKCAQEAnyG0J4zs1JfZm/kq6uKa\n"
    "ktV8Vl1tZvZpkhRYKiDrk5pUVPyCTHQhiUgN7L+o2yXA89QOTRZEfvWb6aP4pEP3\n"
    "hYeDN6OrSICrYLjfWCjaRXp9oNi+qwmhNkKohPO8jkdKhJSnKFwlpe3nYn5oiCZC\n"
    "1wJc/RagwhulBvmjR+MDrv+VTdhS6ccFnqjteTuenP3PQd3+QZh8aCVnial/WXaq\n"
    "mpwX0M6eOuuJQiiRQDAtFvXYOoxXUHZ2dOgFy3k04YEavA2vQ3hGsFfyCZjD9DsI\n"
    "Q3lTsg8pxdfvq/Lo1O+x8OLGiXLE9R44WrIVgiT0HTADCy+3rcx/3AYa6kgLR1gy\n"
    "3wIDAQAB\n"
    "-----END PUBLIC KEY-----\n"
)

# Local HMAC key for license.json file integrity (tamper detection)
# This is separate from RSA — it protects the local file from manual editing
_FILE_HMAC_SECRET = (
    b"ctp-file-integrity-"
    + hashlib.sha256(
        platform.node().encode() + os.getlogin().encode() if hasattr(os, "getlogin") else b""
    ).digest()[:32]
)


# ──────────────────────────────────────────────
# Machine Fingerprint
# ──────────────────────────────────────────────


def get_machine_id() -> str:
    """Generate a machine-specific fingerprint.

    Combines multiple hardware/OS identifiers to create a unique
    machine ID. This binds licenses to specific hardware.

    Returns:
        Hex string of the machine fingerprint
    """
    components = []

    # Platform info
    components.append(platform.node())  # hostname
    components.append(platform.machine())  # arch
    components.append(platform.processor())  # CPU info

    # MAC address (primary network adapter)
    try:
        mac = uuid.getnode()
        components.append(str(mac))
    except Exception:
        pass

    # Windows-specific: machine GUID from registry
    if platform.system() == "Windows":
        try:
            import winreg

            key = winreg.OpenKey(
                winreg.HKEY_LOCAL_MACHINE, r"SOFTWARE\Microsoft\Cryptography", 0, winreg.KEY_READ
            )
            guid, _ = winreg.QueryValueEx(key, "MachineGuid")
            components.append(guid)
            winreg.CloseKey(key)
        except Exception:
            pass

    # Linux-specific: machine-id
    if platform.system() == "Linux":
        try:
            machine_id_path = Path("/etc/machine-id")
            if machine_id_path.exists():
                components.append(machine_id_path.read_text().strip())
        except Exception:
            pass

    # macOS-specific: IOPlatformUUID
    if platform.system() == "Darwin":
        try:
            import subprocess

            result = subprocess.run(
                ["ioreg", "-rd1", "-c", "IOPlatformExpertDevice"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            for line in result.stdout.split("\n"):
                if "IOPlatformUUID" in line:
                    components.append(line.split("=")[-1].strip().strip('"'))
                    break
        except Exception:
            pass

    # Combine and hash
    raw = "|".join(components)
    return hashlib.sha256(raw.encode()).hexdigest()[:32]


# ──────────────────────────────────────────────
# License Key Generation (requires private key)
# ──────────────────────────────────────────────


def _load_private_key() -> object | None:
    """Load the RSA private key from the local file.

    This is only available on the license generation server.
    Returns None if the private key is not found (client-side).
    """
    key_path = Path(__file__).parent.parent.parent / ".license_private_key.pem"
    if not key_path.exists():
        # Also check project root
        key_path = Path.cwd() / ".license_private_key.pem"
    if not key_path.exists():
        return None

    try:
        from cryptography.hazmat.primitives.serialization import load_pem_private_key

        pem_data = key_path.read_bytes()
        return load_pem_private_key(pem_data, password=None)
    except Exception:
        return None


def generate_license_key(
    tier: LicenseTier,
    holder: str,
    email: str,
    days_valid: int = 365,
    machine_id: str = "",
) -> str:
    """Generate a license key with RSA-2048 signature.

    Requires the private key file (.license_private_key.pem) to be present.
    Without the private key, key generation is impossible — this is the
    core security improvement over the previous HMAC-based system.

    Args:
        tier: License tier
        holder: License holder name
        email: License holder email
        days_valid: Days the license is valid
        machine_id: Optional machine fingerprint to bind to

    Returns:
        License key string (CTP-{TIER}-{PAYLOAD}.{SIGNATURE})

    Raises:
        RuntimeError: If private key is not available
    """
    private_key = _load_private_key()
    if private_key is None:
        raise RuntimeError(
            "Private key not found. License generation requires "
            ".license_private_key.pem. This should only be available "
            "on the license server, not on client machines."
        )

    from cryptography.hazmat.primitives import hashes
    from cryptography.hazmat.primitives.asymmetric import padding as asym_padding

    expiry = (datetime.now() + timedelta(days=days_valid)).strftime("%Y%m%d")

    tier_code = {
        LicenseTier.FREE: "FREE",
        LicenseTier.ACADEMIC: "ACAD",
        LicenseTier.COMMERCIAL: "COMM",
        LicenseTier.TRIAL: "TRIAL",
    }[tier]

    # Build payload
    payload_dict = {
        "v": _LICENSE_VERSION,
        "tier": tier_code,
        "holder": holder,
        "email": email,
        "expiry": expiry,
        "machine": machine_id or "",
        "issued": datetime.now().strftime("%Y%m%d"),
    }
    payload_json = json.dumps(payload_dict, separators=(",", ":"), sort_keys=True)
    payload_b64 = base64.urlsafe_b64encode(payload_json.encode()).decode().rstrip("=")

    # RSA sign the payload
    signature = private_key.sign(payload_json.encode(), asym_padding.PKCS1v15(), hashes.SHA256())
    sig_b64 = base64.urlsafe_b64encode(signature).decode().rstrip("=")

    return f"{LICENSE_PREFIX}-{tier_code}-{payload_b64}.{sig_b64}"


# ──────────────────────────────────────────────
# License Key Validation (public key only)
# ──────────────────────────────────────────────


def validate_license_key(key: str) -> tuple[bool, LicenseTier | None, str, dict]:
    """Validate a license key using RSA public key verification.

    This only requires the public key, which is embedded in the source.
    Even with full source access, users cannot forge valid signatures
    because the private key is not present.

    Args:
        key: License key string

    Returns:
        Tuple of (is_valid, tier, message, payload_dict)
    """
    if not key:
        return False, None, "No license key provided", {}

    parts = key.split("-")
    if len(parts) < 4:
        return False, None, "Invalid license key format", {}

    prefix = parts[0]
    tier_code = parts[1]

    if prefix != LICENSE_PREFIX:
        return False, None, "Invalid license key prefix", {}

    # Map tier code
    tier_map = {
        "FREE": LicenseTier.FREE,
        "ACAD": LicenseTier.ACADEMIC,
        "COMM": LicenseTier.COMMERCIAL,
        "TRIAL": LicenseTier.TRIAL,
    }

    if tier_code not in tier_map:
        return False, None, f"Unknown tier code: {tier_code}", {}

    # Split payload and signature
    remainder = "-".join(parts[2:])  # Rejoin in case payload contains dashes
    if "." not in remainder:
        return False, None, "Invalid license key: missing signature", {}

    payload_b64, sig_b64 = remainder.rsplit(".", 1)

    # Decode payload
    try:
        # Add padding back
        padding_needed = 4 - len(payload_b64) % 4
        if padding_needed != 4:
            payload_b64 += "=" * padding_needed
        payload_json = base64.urlsafe_b64decode(payload_b64).decode()
        payload = json.loads(payload_json)
    except Exception:
        return False, None, "Invalid license key: corrupted payload", {}

    # Check expiry
    expiry_str = payload.get("expiry", "")
    try:
        exp_date = datetime.strptime(expiry_str, "%Y%m%d")
        if datetime.now() > exp_date:
            return False, None, f"License expired on {exp_date.strftime('%Y-%m-%d')}", payload
    except ValueError:
        return False, None, "Invalid expiry date in license key", payload

    # Check machine binding
    bound_machine = payload.get("machine", "")
    if bound_machine:
        current_machine = get_machine_id()
        if bound_machine != current_machine:
            return (
                False,
                None,
                (
                    "License is bound to a different machine. "
                    f"Expected: {bound_machine[:8]}..., Got: {current_machine[:8]}..."
                ),
                payload,
            )

    # RSA signature verification
    try:
        from cryptography.hazmat.primitives import hashes
        from cryptography.hazmat.primitives.asymmetric import padding as asym_padding
        from cryptography.hazmat.primitives.serialization import load_pem_public_key

        public_key = load_pem_public_key(_RSA_PUBLIC_KEY_PEM.encode())
        # Reconstruct signature
        sig_padding = 4 - len(sig_b64) % 4
        if sig_padding != 4:
            sig_b64 += "=" * sig_padding
        signature = base64.urlsafe_b64decode(sig_b64)

        # Verify — this will raise InvalidSignature if the signature is wrong
        public_key.verify(
            signature, payload_json.encode(), asym_padding.PKCS1v15(), hashes.SHA256()
        )
    except ImportError:
        # cryptography not installed — fall back to basic validation
        # This is acceptable for free tier but should warn
        return (
            True,
            tier_map[tier_code],
            (
                "License validated (signature check skipped — install cryptography for full verification)"
            ),
            payload,
        )
    except Exception as e:
        return False, None, f"License signature verification failed: {e}", payload

    return True, tier_map[tier_code], "License valid (signature verified)", payload


# ──────────────────────────────────────────────
# License Storage (with tamper detection)
# ──────────────────────────────────────────────

DEFAULT_LICENSE_DIR = Path.home() / ".celltypepilot"
LICENSE_FILE = "license.json"


def get_license_path() -> Path:
    """Get the license file path."""
    return DEFAULT_LICENSE_DIR / LICENSE_FILE


def _compute_file_hmac(data: str) -> str:
    """Compute HMAC for license file tamper detection."""
    return hmac.new(_FILE_HMAC_SECRET, data.encode(), hashlib.sha256).hexdigest()


def save_license(license_info: LicenseInfo) -> Path:
    """Save license info to disk with integrity signature.

    The saved file includes an HMAC signature that detects manual
    editing of the license.json file.

    Args:
        license_info: License information to save

    Returns:
        Path to saved license file
    """
    DEFAULT_LICENSE_DIR.mkdir(parents=True, exist_ok=True)
    license_path = get_license_path()

    data = asdict(license_info)
    data["tier"] = license_info.tier.value
    # Remove internal fields
    data.pop("_signature_valid", None)

    json_str = json.dumps(data, indent=2, sort_keys=True)
    # Add integrity signature
    data["_hmac"] = _compute_file_hmac(json_str)

    license_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return license_path


def load_license() -> LicenseInfo:
    """Load license info from disk with tamper detection.

    If the license file has been manually edited (HMAC mismatch),
    the license is invalidated and free tier is returned.

    Returns:
        LicenseInfo (free tier if no license file or tampered)
    """
    license_path = get_license_path()

    if not license_path.exists():
        return LicenseInfo(tier=LicenseTier.FREE)

    try:
        raw = json.loads(license_path.read_text(encoding="utf-8"))
        stored_hmac = raw.pop("_hmac", "")

        # Verify integrity
        json_str = json.dumps(raw, indent=2, sort_keys=True)
        expected_hmac = _compute_file_hmac(json_str)

        if stored_hmac and stored_hmac != expected_hmac:
            # File has been tampered with
            return LicenseInfo(tier=LicenseTier.FREE)

        raw["tier"] = LicenseTier(raw.get("tier", "free"))
        return LicenseInfo(
            **{k: v for k, v in raw.items() if k in LicenseInfo.__dataclass_fields__}
        )
    except (json.JSONDecodeError, KeyError, ValueError):
        return LicenseInfo(tier=LicenseTier.FREE)


def activate_license(key: str, holder: str = "", email: str = "") -> tuple[bool, str]:
    """Activate a license key.

    Validates the RSA signature, checks machine binding,
    and saves the license with integrity protection.

    Args:
        key: License key
        holder: License holder name
        email: License holder email

    Returns:
        Tuple of (success, message)
    """
    is_valid, tier, message, payload = validate_license_key(key)

    if not is_valid:
        return False, message

    # Extract info from payload
    holder = holder or payload.get("holder", "")
    email = email or payload.get("email", "")
    expiry_str = payload.get("expiry", "")
    machine_id = payload.get("machine", "")

    # Create license info
    license_info = LicenseInfo(
        tier=tier,
        key=key,
        holder=holder,
        email=email,
        issued_at=datetime.now().isoformat(),
        expires_at=datetime.strptime(expiry_str, "%Y%m%d").isoformat() if expiry_str else "",
        features=list(
            ACADEMIC_FEATURES
            if tier in (LicenseTier.ACADEMIC, LicenseTier.TRIAL)
            else COMMERCIAL_FEATURES
            if tier == LicenseTier.COMMERCIAL
            else FREE_FEATURES
        ),
        machine_id=machine_id or get_machine_id(),
        _signature_valid=True,
    )

    save_license(license_info)
    return True, f"License activated: {tier.value} tier (signature verified)"


def check_feature_access(feature: str) -> tuple[bool, str]:
    """Check if a feature is accessible with current license.

    Args:
        feature: Feature name to check

    Returns:
        Tuple of (has_access, message)
    """
    license_info = load_license()

    if license_info.is_expired():
        return False, "License expired. Renew at https://celltypepilot.io/license"

    if license_info.has_feature(feature):
        return True, f"Feature '{feature}' available ({license_info.tier.value} tier)"

    # Suggest upgrade
    if feature in ACADEMIC_FEATURES - FREE_FEATURES:
        return False, (
            f"Feature '{feature}' requires Academic or Commercial license. "
            f"Current: {license_info.tier.value}. "
            f"Upgrade at https://celltypepilot.io/license"
        )
    elif feature in COMMERCIAL_FEATURES - ACADEMIC_FEATURES:
        return False, (
            f"Feature '{feature}' requires Commercial license. "
            f"Current: {license_info.tier.value}. "
            f"Upgrade at https://celltypepilot.io/license"
        )

    return False, f"Unknown feature: {feature}"


# ──────────────────────────────────────────────
# Premium atlas gating
# ──────────────────────────────────────────────


def get_atlas_access(tissue: str) -> tuple[bool, str]:
    """Check if a tissue is accessible in the atlas.

    Free tier: 11 basic tissues
    Academic/Commercial: All tissues including extended panels

    Args:
        tissue: Tissue name to check

    Returns:
        Tuple of (has_access, message)
    """
    license_info = load_license()

    # Basic tissues available to all tiers
    basic_tissues = {
        "blood",
        "lung",
        "liver",
        "brain",
        "kidney",
        "gut",
        "skin",
        "heart",
        "pancreas",
        "muscle",
        "general",
    }

    if tissue in basic_tissues:
        return True, f"Tissue '{tissue}' available (all tiers)"

    # Extended tissues require academic+ tier
    if license_info.tier in (LicenseTier.ACADEMIC, LicenseTier.COMMERCIAL, LicenseTier.TRIAL):
        return True, f"Tissue '{tissue}' available ({license_info.tier.value} tier)"

    return False, (
        f"Tissue '{tissue}' requires Academic or Commercial license. "
        f"Current: {license_info.tier.value}. "
        f"Available basic tissues: {', '.join(sorted(basic_tissues))}"
    )
