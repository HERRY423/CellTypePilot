"""Tests for the license manager module."""

import json
import tempfile
from datetime import datetime, timedelta
from pathlib import Path
from unittest.mock import patch

import pytest

from celltypepilot.license_manager import (
    LicenseTier, LicenseInfo,
    FREE_FEATURES, ACADEMIC_FEATURES, COMMERCIAL_FEATURES, TRIAL_FEATURES,
    get_machine_id, validate_license_key,
    save_license, load_license,
    check_feature_access, get_atlas_access,
    _compute_file_hmac,
    DEFAULT_LICENSE_DIR,
)


class TestLicenseTier:
    def test_tier_values(self):
        assert LicenseTier.FREE.value == "free"
        assert LicenseTier.ACADEMIC.value == "academic"
        assert LicenseTier.COMMERCIAL.value == "commercial"
        assert LicenseTier.TRIAL.value == "trial"

    def test_tier_is_string_enum(self):
        for tier in LicenseTier:
            assert isinstance(tier.value, str)


class TestLicenseInfo:
    def test_default_is_free(self):
        info = LicenseInfo()
        assert info.tier == LicenseTier.FREE
        assert info.valid is True

    def test_is_expired_no_expiry(self):
        info = LicenseInfo(expires_at="")
        assert not info.is_expired()

    def test_is_expired_future(self):
        future = (datetime.now() + timedelta(days=365)).isoformat()
        info = LicenseInfo(expires_at=future)
        assert not info.is_expired()

    def test_is_expired_past(self):
        past = (datetime.now() - timedelta(days=1)).isoformat()
        info = LicenseInfo(expires_at=past)
        assert info.is_expired()

    def test_is_expired_invalid_format(self):
        info = LicenseInfo(expires_at="not-a-date")
        assert not info.is_expired()

    def test_has_feature_free(self):
        info = LicenseInfo(tier=LicenseTier.FREE)
        assert info.has_feature("basic_atlas")
        assert info.has_feature("marker_scoring")
        assert not info.has_feature("extended_atlas")
        assert not info.has_feature("custom_tissue_panels")

    def test_has_feature_academic(self):
        info = LicenseInfo(tier=LicenseTier.ACADEMIC)
        assert info.has_feature("basic_atlas")
        assert info.has_feature("extended_atlas")
        assert info.has_feature("literature_validation")
        assert not info.has_feature("custom_tissue_panels")

    def test_has_feature_commercial(self):
        info = LicenseInfo(tier=LicenseTier.COMMERCIAL)
        assert info.has_feature("basic_atlas")
        assert info.has_feature("extended_atlas")
        assert info.has_feature("custom_tissue_panels")
        assert info.has_feature("white_label")

    def test_has_feature_trial(self):
        info = LicenseInfo(tier=LicenseTier.TRIAL)
        assert info.has_feature("extended_atlas")
        assert info.has_feature("literature_validation")

    def test_has_feature_unknown(self):
        info = LicenseInfo(tier=LicenseTier.FREE)
        assert not info.has_feature("nonexistent_feature")


class TestFeatureSets:
    def test_free_is_subset_of_academic(self):
        assert FREE_FEATURES.issubset(ACADEMIC_FEATURES)

    def test_academic_is_subset_of_commercial(self):
        assert ACADEMIC_FEATURES.issubset(COMMERCIAL_FEATURES)

    def test_trial_equals_academic(self):
        assert TRIAL_FEATURES == ACADEMIC_FEATURES

    def test_free_has_core_features(self):
        core = {"basic_atlas", "marker_scoring", "critic_review", "html_report"}
        assert core.issubset(FREE_FEATURES)

    def test_academic_adds_extended(self):
        extra = ACADEMIC_FEATURES - FREE_FEATURES
        assert "extended_atlas" in extra
        assert "literature_validation" in extra

    def test_commercial_adds_exclusive(self):
        extra = COMMERCIAL_FEATURES - ACADEMIC_FEATURES
        assert "custom_tissue_panels" in extra
        assert "white_label" in extra


class TestMachineId:
    def test_returns_string(self):
        mid = get_machine_id()
        assert isinstance(mid, str)
        assert len(mid) == 32  # SHA-256 hex truncated to 32 chars

    def test_deterministic(self):
        mid1 = get_machine_id()
        mid2 = get_machine_id()
        assert mid1 == mid2


class TestValidateLicenseKey:
    def test_empty_key(self):
        valid, tier, msg, payload = validate_license_key("")
        assert not valid
        assert tier is None
        assert "No license key" in msg

    def test_invalid_format(self):
        valid, tier, msg, _ = validate_license_key("INVALID")
        assert not valid

    def test_wrong_prefix(self):
        valid, tier, msg, _ = validate_license_key("XXX-FREE-payload.sig")
        assert not valid

    def test_unknown_tier_code(self):
        valid, tier, msg, _ = validate_license_key("CTP-UNKNOWN-payload.sig")
        assert not valid

    def test_missing_signature(self):
        valid, tier, msg, _ = validate_license_key("CTP-FREE-payload_no_dot")
        assert not valid


class TestLicenseSaveLoad:
    def test_save_and_load_roundtrip(self):
        info = LicenseInfo(
            tier=LicenseTier.FREE,
            holder="Test User",
            email="test@example.com",
        )
        with patch("celltypepilot.license_manager.DEFAULT_LICENSE_DIR",
                    new=Path(tempfile.mkdtemp())):
            saved_path = save_license(info)
            assert saved_path.exists()

            loaded = load_license()
            assert loaded.tier == LicenseTier.FREE
            assert loaded.holder == "Test User"
            assert loaded.email == "test@example.com"

    def test_load_nonexistent_returns_free(self):
        with patch("celltypepilot.license_manager.DEFAULT_LICENSE_DIR",
                    new=Path(tempfile.mkdtemp())):
            info = load_license()
            assert info.tier == LicenseTier.FREE

    def test_tamper_detection(self):
        with patch("celltypepilot.license_manager.DEFAULT_LICENSE_DIR",
                    new=Path(tempfile.mkdtemp())) as lic_dir:
            info = LicenseInfo(tier=LicenseTier.ACADEMIC, holder="Hacker")
            save_license(info)

            # Tamper with the file
            lic_path = lic_dir / "license.json"
            data = json.loads(lic_path.read_text(encoding="utf-8"))
            data["holder"] = "Modified"
            lic_path.write_text(json.dumps(data), encoding="utf-8")

            loaded = load_license()
            assert loaded.tier == LicenseTier.FREE  # Tampered → fallback


class TestFileHmac:
    def test_deterministic(self):
        h1 = _compute_file_hmac("test data")
        h2 = _compute_file_hmac("test data")
        assert h1 == h2

    def test_different_data_different_hmac(self):
        h1 = _compute_file_hmac("data A")
        h2 = _compute_file_hmac("data B")
        assert h1 != h2


class TestAtlasAccess:
    def test_basic_tissue_free_tier(self):
        with patch("celltypepilot.license_manager.load_license",
                    return_value=LicenseInfo(tier=LicenseTier.FREE)):
            ok, msg = get_atlas_access("blood")
            assert ok

    def test_extended_tissue_free_tier_denied(self):
        with patch("celltypepilot.license_manager.load_license",
                    return_value=LicenseInfo(tier=LicenseTier.FREE)):
            ok, msg = get_atlas_access("tumor_microenvironment")
            assert not ok

    def test_extended_tissue_academic_allowed(self):
        with patch("celltypepilot.license_manager.load_license",
                    return_value=LicenseInfo(tier=LicenseTier.ACADEMIC)):
            ok, msg = get_atlas_access("tumor_microenvironment")
            assert ok


class TestCheckFeatureAccess:
    def test_free_feature_available(self):
        with patch("celltypepilot.license_manager.load_license",
                    return_value=LicenseInfo(tier=LicenseTier.FREE)):
            ok, msg = check_feature_access("basic_atlas")
            assert ok

    def test_academic_feature_blocked_free(self):
        with patch("celltypepilot.license_manager.load_license",
                    return_value=LicenseInfo(tier=LicenseTier.FREE)):
            ok, msg = check_feature_access("extended_atlas")
            assert not ok
            assert "Academic" in msg

    def test_expired_license_blocked(self):
        past = (datetime.now() - timedelta(days=1)).isoformat()
        with patch("celltypepilot.license_manager.load_license",
                    return_value=LicenseInfo(tier=LicenseTier.ACADEMIC, expires_at=past)):
            ok, msg = check_feature_access("basic_atlas")
            assert not ok
            assert "expired" in msg.lower()
