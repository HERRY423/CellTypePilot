"""Three-domain validation focus and claim-readiness boundaries."""

from __future__ import annotations

import json
from importlib.resources import files
from typing import Any

from .constants import FIRST_PARTY_PACKS_DIR

DOMAIN_SCHEMA = "celltypepilot.validation-domains.v1"


class ValidationDomainError(ValueError):
    """Raised when the bundled validation-domain registry is invalid."""


def load_validation_domains() -> dict:
    path = files("celltypepilot.data").joinpath("validation_domains.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    validate_validation_domains(payload)
    return payload


def validate_validation_domains(payload: dict) -> None:
    if payload.get("schema_version") != DOMAIN_SCHEMA:
        raise ValidationDomainError(f"Expected schema_version {DOMAIN_SCHEMA}")
    domains = payload.get("domains")
    required = {"lung", "gut_ibd", "tumor_microenvironment"}
    if not isinstance(domains, dict) or set(domains) != required:
        raise ValidationDomainError(
            "Registry must define exactly lung, gut_ibd, and tumor_microenvironment"
        )
    for name, domain in domains.items():
        if domain.get("claim_ready") is not False:
            raise ValidationDomainError(f"Development registry domain {name} must fail closed")
        if len(domain.get("required_candidate_backends", [])) < 5:
            raise ValidationDomainError(f"Domain {name} lacks the required comparator depth")
        if not domain.get("required_validation_axes"):
            raise ValidationDomainError(f"Domain {name} lacks validation axes")
        atlas_contract = domain.get("atlas_contract")
        if not isinstance(atlas_contract, dict):
            raise ValidationDomainError(f"Domain {name} lacks atlas_contract")
        if atlas_contract.get("status") != "scope_complete_evidence_not_claim_ready":
            raise ValidationDomainError(f"Domain {name} atlas scope must remain non-claim-ready")
        lineages = atlas_contract.get("required_lineages")
        if not isinstance(lineages, list) or len(lineages) < 4:
            raise ValidationDomainError(f"Domain {name} requires at least four lineage scopes")
        required_packs = atlas_contract.get("required_packs")
        if not isinstance(required_packs, list) or not required_packs:
            raise ValidationDomainError(f"Domain {name} lacks required Atlas packs")
        missing_packs = [
            pack
            for pack in required_packs
            if not (FIRST_PARTY_PACKS_DIR / str(pack) / "pack.json").is_file()
        ]
        if missing_packs:
            raise ValidationDomainError(
                f"Domain {name} references missing Atlas packs: {missing_packs}"
            )


def assess_validation_domain(
    tissue: str | None, packs: list[str] | None = None, registry: dict | None = None
) -> dict[str, Any]:
    """Resolve the focused domain; out-of-focus tissues remain exploratory only."""
    payload = registry or load_validation_domains()
    tissue_key = str(tissue or "general").strip().casefold().replace("-", "_").replace(" ", "_")
    pack_keys = {str(value).casefold() for value in (packs or [])}
    if tissue_key == "lung":
        domain_id = "lung"
    elif tissue_key in {"gut", "ibd", "inflamed_gut", "intestine", "colon"}:
        domain_id = "gut_ibd"
    elif tissue_key in {"tumor_microenvironment", "tumor", "cancer", "tme"} or (
        "premium" in pack_keys and tissue_key == "tumor_microenvironment"
    ):
        domain_id = "tumor_microenvironment"
    else:
        return {
            "schema_version": DOMAIN_SCHEMA,
            "domain_id": None,
            "status": "out_of_focus_exploratory",
            "claim_ready": False,
            "blockers": ["OUTSIDE_THREE_DEPTH_VALIDATION_DOMAINS"],
            "claim_boundary": (
                "Atlas coverage outside the three focus domains is exploratory and not domain-validated"
            ),
        }
    domain = payload["domains"][domain_id]
    return {
        "schema_version": DOMAIN_SCHEMA,
        "domain_id": domain_id,
        "status": domain["status"],
        "claim_ready": False,
        "blockers": list(domain["current_blockers"]),
        "requirements": {
            "atlas_contract": domain["atlas_contract"],
            "required_candidate_backends": domain["required_candidate_backends"],
            "required_validation_axes": domain["required_validation_axes"],
            "minimum_evidence": domain["minimum_evidence"],
        },
        "claim_boundary": domain["claim_boundary"],
    }
