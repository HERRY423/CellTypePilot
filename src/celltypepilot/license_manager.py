"""CellTypePilot — Tiered license manager.

Manages license tiers and gates access to premium features:
- Free tier: Built-in MKG (80+ cell types, 11 tissues)
- Academic tier: Extended atlas (200+ cell types, disease states, developmental stages)
- Commercial tier: Full atlas + priority support + custom tissue panels

License keys are validated locally (offline-first) with optional
online verification for team/commercial licenses.
"""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass, field, asdict
from datetime import datetime, timedelta
from enum import Enum
from pathlib import Path
from typing import Optional


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
    valid: bool = True

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
    "basic_atlas",           # 80+ cell types, 11 tissues
    "marker_scoring",        # Wilcoxon DE + 5-dim scoring
    "critic_review",         # Annotation Critic
    "basic_visualization",   # UMAP, dotplot, confidence
    "html_report",           # HTML report generation
    "manifest_provenance",   # Run provenance tracking
}

ACADEMIC_FEATURES = FREE_FEATURES | {
    "extended_atlas",        # 200+ cell types, disease states
    "developmental_atlas",   # Developmental stage markers
    "disease_atlas",         # Disease-specific cell states
    "literature_validation", # PubMed integration
    "advanced_visualization", # Cross-sample comparisons
    "docx_export",           # Word document export
}

COMMERCIAL_FEATURES = ACADEMIC_FEATURES | {
    "custom_tissue_panels",  # Custom tissue marker panels
    "team_sharing",          # Team atlas sharing
    "priority_support",      # Priority email support
    "api_access",            # REST API access
    "white_label",           # White-label reports
}

TRIAL_FEATURES = ACADEMIC_FEATURES | set()  # Trial gets academic features for 30 days


# ──────────────────────────────────────────────
# License key generation/validation
# ──────────────────────────────────────────────

# Simple license key format: CTP-{TIER}-{HASH}-{EXPIRY}
# Example: CTP-ACAD-A1B2C3D4-20271231

LICENSE_PREFIX = "CTP"
LICENSE_SECRET = "celltypepilot-license-v1"  # Not a real secret, just for checksum


def generate_license_key(
    tier: LicenseTier,
    holder: str,
    email: str,
    days_valid: int = 365,
) -> str:
    """Generate a license key.

    This is a simple offline license system. For production use,
    you'd want proper cryptographic signing.

    Args:
        tier: License tier
        holder: License holder name
        email: License holder email
        days_valid: Days the license is valid

    Returns:
        License key string
    """
    expiry = (datetime.now() + timedelta(days=days_valid)).strftime("%Y%m%d")

    # Create checksum
    payload = f"{tier.value}:{holder}:{email}:{expiry}:{LICENSE_SECRET}"
    checksum = hashlib.sha256(payload.encode()).hexdigest()[:8].upper()

    # Format key
    tier_code = {
        LicenseTier.FREE: "FREE",
        LicenseTier.ACADEMIC: "ACAD",
        LicenseTier.COMMERCIAL: "COMM",
        LicenseTier.TRIAL: "TRIAL",
    }[tier]

    return f"{LICENSE_PREFIX}-{tier_code}-{checksum}-{expiry}"


def validate_license_key(key: str) -> tuple[bool, Optional[LicenseTier], str]:
    """Validate a license key.

    Args:
        key: License key string

    Returns:
        Tuple of (is_valid, tier, message)
    """
    if not key:
        return False, None, "No license key provided"

    parts = key.split("-")
    if len(parts) != 4:
        return False, None, "Invalid license key format"

    prefix, tier_code, checksum, expiry = parts

    if prefix != LICENSE_PREFIX:
        return False, None, "Invalid license key prefix"

    # Check expiry
    try:
        exp_date = datetime.strptime(expiry, "%Y%m%d")
        if datetime.now() > exp_date:
            return False, None, f"License expired on {exp_date.strftime('%Y-%m-%d')}"
    except ValueError:
        return False, None, "Invalid expiry date in license key"

    # Map tier code
    tier_map = {
        "FREE": LicenseTier.FREE,
        "ACAD": LicenseTier.ACADEMIC,
        "COMM": LicenseTier.COMMERCIAL,
        "TRIAL": LicenseTier.TRIAL,
    }

    if tier_code not in tier_map:
        return False, None, f"Unknown tier code: {tier_code}"

    return True, tier_map[tier_code], "License valid"


# ──────────────────────────────────────────────
# License storage
# ──────────────────────────────────────────────

DEFAULT_LICENSE_DIR = Path.home() / ".celltypepilot"
LICENSE_FILE = "license.json"


def get_license_path() -> Path:
    """Get the license file path."""
    return DEFAULT_LICENSE_DIR / LICENSE_FILE


def save_license(license_info: LicenseInfo) -> Path:
    """Save license info to disk.

    Args:
        license_info: License information to save

    Returns:
        Path to saved license file
    """
    DEFAULT_LICENSE_DIR.mkdir(parents=True, exist_ok=True)
    license_path = get_license_path()

    data = asdict(license_info)
    data["tier"] = license_info.tier.value

    license_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
    return license_path


def load_license() -> LicenseInfo:
    """Load license info from disk.

    Returns:
        LicenseInfo (free tier if no license file exists)
    """
    license_path = get_license_path()

    if not license_path.exists():
        return LicenseInfo(tier=LicenseTier.FREE)

    try:
        data = json.loads(license_path.read_text(encoding="utf-8"))
        data["tier"] = LicenseTier(data.get("tier", "free"))
        return LicenseInfo(**data)
    except (json.JSONDecodeError, KeyError, ValueError):
        return LicenseInfo(tier=LicenseTier.FREE)


def activate_license(key: str, holder: str = "", email: str = "") -> tuple[bool, str]:
    """Activate a license key.

    Args:
        key: License key
        holder: License holder name
        email: License holder email

    Returns:
        Tuple of (success, message)
    """
    is_valid, tier, message = validate_license_key(key)

    if not is_valid:
        return False, message

    # Create license info
    license_info = LicenseInfo(
        tier=tier,
        key=key,
        holder=holder,
        email=email,
        issued_at=datetime.now().isoformat(),
        # Extract expiry from key
        expires_at=datetime.strptime(key.split("-")[-1], "%Y%m%d").isoformat(),
        features=list(ACADEMIC_FEATURES if tier in (LicenseTier.ACADEMIC, LicenseTier.TRIAL)
                      else COMMERCIAL_FEATURES if tier == LicenseTier.COMMERCIAL
                      else FREE_FEATURES),
    )

    save_license(license_info)
    return True, f"License activated: {tier.value} tier"


def check_feature_access(feature: str) -> tuple[bool, str]:
    """Check if a feature is accessible with current license.

    Args:
        feature: Feature name to check

    Returns:
        Tuple of (has_access, message)
    """
    license_info = load_license()

    if license_info.is_expired():
        return False, f"License expired. Renew at https://celltypepilot.io/license"

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
        "blood", "lung", "liver", "brain", "kidney",
        "gut", "skin", "heart", "pancreas", "muscle", "general"
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
