"""Biological identity contracts for marker and reference evidence.

This module resolves two representation boundaries that must be correct before
CellTypePilot can interpret biological evidence:

* feature identifiers (for example Ensembl IDs with symbols in ``var``), and
* cell-type identities (canonical atlas labels, synonyms, CL IDs, and explicit
  safe parent fallbacks supplied by a governed evidence pack).

The contract is deliberately conservative. Ambiguous gene symbols and aliases
are left unresolved, and a fine label is only collapsed to a parent when an
explicit pack mapping authorizes that loss of resolution.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import Counter
from collections.abc import Iterable
from copy import deepcopy
from typing import Any

import pandas as pd

GENE_SYMBOL_COLUMNS = (
    "gene_symbol",
    "gene_symbols",
    "feature_name",
    "symbol",
    "gene_name",
)
ONTOLOGY_ID_COLUMNS = (
    "cell_type_ontology_term_id",
    "cell_ontology_id",
    "cl_id",
)
GENE_TOKEN = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


class IdentityContractError(ValueError):
    """Raised when an identity contract is ambiguous or malformed."""


def _norm_label(value: Any) -> str:
    return " ".join(str(value or "").strip().casefold().replace("_", " ").split())


def _sha256_lines(values: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _iter_nodes(cell_types: dict, tissue: str, parents: tuple[str, ...] = ()):
    for name, info in cell_types.items():
        path = (*parents, name)
        yield tissue, path, info
        yield from _iter_nodes(info.get("subtypes", {}), tissue, path)


def collect_pack_identity_contract(records: list[dict] | None) -> dict:
    """Merge data-only ontology maps from selected packs, failing on conflicts."""
    merged: dict[str, Any] = {
        "include_tissues": [],
        "aliases": {},
        "safe_parent_fallbacks": {},
        "sources": [],
    }
    for record in records or []:
        mapping = record.get("ontology_map") or {}
        if not mapping:
            continue
        if not isinstance(mapping, dict):
            raise IdentityContractError(
                f"Pack {record.get('name')!r} ontology_map must be an object"
            )
        merged["sources"].append({"name": record.get("name"), "version": record.get("version")})
        for tissue in mapping.get("include_tissues", []):
            tissue = str(tissue)
            if tissue not in merged["include_tissues"]:
                merged["include_tissues"].append(tissue)
        for field in ("aliases", "safe_parent_fallbacks"):
            values = mapping.get(field, {})
            if not isinstance(values, dict):
                raise IdentityContractError(
                    f"Pack {record.get('name')!r} ontology_map.{field} must be an object"
                )
            for raw, canonical in values.items():
                key = _norm_label(raw)
                canonical = str(canonical).strip()
                prior = merged[field].get(key)
                if prior is not None and prior != canonical:
                    raise IdentityContractError(
                        f"Conflicting {field} mapping for {raw!r}: {prior!r} vs {canonical!r}"
                    )
                merged[field][key] = canonical
    return merged


def compose_marker_definitions(
    atlas: dict,
    tissue: str,
    evidence_policy: str = "database",
    pack_contract: dict | None = None,
    *,
    include_unverified_candidates: bool = False,
) -> tuple[dict[str, dict], dict]:
    """Compose tissue marker scope with general and explicitly included scopes.

    Tissue-specific definitions have precedence. ``general`` is included for
    every non-general tissue; additional scopes (for example blood in a lung
    pack) must be declared by a selected pack ontology map.
    """
    from .data_adapter import get_all_markers_for_tissue

    available = atlas.get("tissues", {})
    requested = [str(tissue)]
    if tissue != "general":
        requested.append("general")
    requested.extend((pack_contract or {}).get("include_tissues", []))
    scopes: list[str] = []
    for scope in requested:
        if scope in available and scope not in scopes:
            scopes.append(scope)

    merged: dict[str, dict] = {}
    collisions: list[dict] = []
    for scope in scopes:
        definitions = get_all_markers_for_tissue(atlas, scope, evidence_policy=evidence_policy)
        if include_unverified_candidates:
            for info in definitions.values():
                records = info.get("marker_evidence", [])
                info["positive_markers"] = [
                    str(record["gene"])
                    for record in records
                    if record.get("polarity") == "positive" and record.get("gene")
                ]
                info["negative_markers"] = [
                    str(record["gene"])
                    for record in records
                    if record.get("polarity") == "negative" and record.get("gene")
                ]
                info["evidence_policy"] = "candidate_inventory_not_runtime_eligible"
        for name, info in definitions.items():
            if name in merged:
                if merged[name].get("cl_id") != info.get("cl_id"):
                    collisions.append(
                        {
                            "cell_type": name,
                            "kept_scope": merged[name].get("evidence_tissue"),
                            "ignored_scope": scope,
                            "reason": "cl_id_conflict",
                        }
                    )
                continue
            entry = deepcopy(info)
            entry["evidence_tissue"] = scope
            merged[name] = entry
    return merged, {
        "schema_version": "celltypepilot.identity-scope.v1",
        "requested_tissue": tissue,
        "active_tissues": scopes,
        "n_cell_types": len(merged),
        "collisions": collisions,
        "pack_contract_sources": (pack_contract or {}).get("sources", []),
        "include_unverified_candidates": include_unverified_candidates,
    }


def build_identity_resolver(
    atlas: dict,
    active_tissues: Iterable[str],
    pack_contract: dict | None = None,
) -> dict:
    """Build an ambiguity-aware label/CL resolver for the active scope."""
    by_label: dict[str, str | None] = {}
    by_cl: dict[str, str | None] = {}
    parent_by_name: dict[str, str] = {}

    def add(index: dict[str, str | None], key: str, canonical: str) -> None:
        if not key:
            return
        prior = index.get(key)
        if prior is None and key in index:
            return
        if prior is not None and prior != canonical:
            index[key] = None
        else:
            index[key] = canonical

    for tissue in active_tissues:
        tissue_data = atlas.get("tissues", {}).get(tissue, {})
        for _tissue, path, info in _iter_nodes(tissue_data.get("cell_types", {}), tissue):
            canonical = path[-1]
            add(by_label, _norm_label(canonical), canonical)
            for synonym in info.get("synonyms", []):
                add(by_label, _norm_label(synonym), canonical)
            cl_id = str(info.get("cl_id", "")).strip()
            if cl_id:
                add(by_cl, cl_id, canonical)
            if len(path) > 1:
                parent_by_name[canonical] = path[-2]

    for raw, canonical in (pack_contract or {}).get("aliases", {}).items():
        add(by_label, raw, canonical)

    return {
        "schema_version": "celltypepilot.identity-resolver.v1",
        "by_label": by_label,
        "by_cl": by_cl,
        "parent_by_name": parent_by_name,
        "safe_parent_fallbacks": (pack_contract or {}).get("safe_parent_fallbacks", {}),
    }


def resolve_identity_label(label: Any, resolver: dict, cl_id: Any = "") -> dict:
    """Resolve one label without inventing a biological equivalence."""
    raw = str(label or "").strip()
    cl_value = str(cl_id or "").strip()
    if cl_value:
        canonical = resolver.get("by_cl", {}).get(cl_value)
        if canonical:
            return {
                "raw_label": raw,
                "canonical_label": canonical,
                "resolution": "cell_ontology_id",
                "cl_id": cl_value,
            }

    key = _norm_label(raw)
    canonical = resolver.get("by_label", {}).get(key)
    if canonical:
        return {
            "raw_label": raw,
            "canonical_label": canonical,
            "resolution": "canonical_or_alias",
            "cl_id": cl_value,
        }
    if key in resolver.get("by_label", {}) and canonical is None:
        return {
            "raw_label": raw,
            "canonical_label": raw,
            "resolution": "ambiguous_alias_unresolved",
            "cl_id": cl_value,
        }

    parent = resolver.get("safe_parent_fallbacks", {}).get(key)
    if parent:
        parent_resolved = resolver.get("by_label", {}).get(_norm_label(parent))
        if not parent_resolved:
            raise IdentityContractError(
                f"Safe parent fallback target {parent!r} is absent or ambiguous in active scope"
            )
        return {
            "raw_label": raw,
            "canonical_label": parent_resolved,
            "resolution": "explicit_safe_parent_fallback",
            "cl_id": cl_value,
        }
    return {
        "raw_label": raw,
        "canonical_label": raw,
        "resolution": "unresolved",
        "cl_id": cl_value,
    }


def canonicalize_reference_adata(reference, label_key: str, resolver: dict) -> dict:
    """Canonicalize reference labels in memory and retain their raw values."""
    if label_key not in reference.obs:
        raise IdentityContractError(f"Reference label key {label_key!r} is missing")
    raw = reference.obs[label_key].astype(str)
    ontology_column = next((c for c in ONTOLOGY_ID_COLUMNS if c in reference.obs), None)
    ontology = (
        reference.obs[ontology_column].astype(str)
        if ontology_column is not None
        else pd.Series("", index=reference.obs.index)
    )
    resolved = [
        resolve_identity_label(label, resolver, cl_id)
        for label, cl_id in zip(raw, ontology, strict=True)
    ]
    reference.obs[f"ctp_raw_{label_key}"] = raw.to_numpy()
    reference.obs[label_key] = [item["canonical_label"] for item in resolved]
    counts = Counter(item["resolution"] for item in resolved)
    unresolved = sorted(
        {item["raw_label"] for item in resolved if item["resolution"].endswith("unresolved")}
    )
    return {
        "schema_version": "celltypepilot.reference-identity-audit.v1",
        "label_key": label_key,
        "ontology_column": ontology_column,
        "resolution_counts": dict(counts),
        "unresolved_labels": unresolved,
    }


def canonicalize_reference_scores(scores: pd.DataFrame, resolver: dict) -> pd.DataFrame:
    """Canonicalize model-only reference output and recompute deterministic ranks."""
    if scores.empty or "cell_type" not in scores:
        return scores
    attrs = dict(scores.attrs)
    frame = scores.copy()
    resolutions = [resolve_identity_label(value, resolver) for value in frame["cell_type"]]
    frame["raw_cell_type"] = frame["cell_type"].astype(str)
    frame["cell_type"] = [item["canonical_label"] for item in resolutions]
    frame["identity_resolution"] = [item["resolution"] for item in resolutions]
    if frame.duplicated(["cluster", "cell_type"]).any():
        frame = (
            frame.sort_values("ref_score", ascending=False)
            .groupby(["cluster", "cell_type"], as_index=False, sort=False)
            .first()
        )
    frame["ref_rank"] = (
        frame.groupby("cluster")["ref_score"].rank(ascending=False, method="first").astype(int)
    )
    frame = frame.sort_values(["cluster", "ref_rank"])
    frame.attrs.update(attrs)
    frame.attrs["identity_resolution"] = "canonicalized_max_on_collapsed_labels"
    return frame


def _valid_symbol(value: Any) -> str:
    value = str(value or "").strip()
    if not value or value.casefold() in {"nan", "none"} or not GENE_TOKEN.fullmatch(value):
        return ""
    return value


def apply_gene_identity_contract(adata, target_markers: Iterable[str]) -> dict:
    """Expose unique gene symbols for evidence scoring, retaining original IDs.

    A symbol-bearing ``var`` column is selected only when it improves overlap
    with the declared marker universe. Duplicate symbols remain on their
    original identifiers so a marker can never silently bind an arbitrary
    duplicate feature.
    """
    target = {str(value) for value in target_markers if str(value)}
    original = [str(value) for value in adata.var_names]
    before = len(target.intersection(original))
    candidates: list[tuple[int, int, str, list[str]]] = []
    for column in GENE_SYMBOL_COLUMNS:
        if column not in adata.var:
            continue
        values = [_valid_symbol(value) for value in adata.var[column]]
        overlap = len(target.intersection(values))
        valid = sum(bool(value) for value in values)
        candidates.append((overlap, valid, column, values))

    chosen = max(candidates, default=None, key=lambda item: (item[0], item[1]))
    applied = bool(chosen and chosen[0] > before)
    duplicate_symbols: list[str] = []
    if applied and chosen is not None:
        _overlap, _valid, column, values = chosen
        counts = Counter(value for value in values if value)
        duplicate_symbols = sorted(value for value, count in counts.items() if count > 1)
        proposed = [
            value if value and counts[value] == 1 else original[index]
            for index, value in enumerate(values)
        ]
        final_counts = Counter(proposed)
        proposed = [
            value if final_counts[value] == 1 else original[index]
            for index, value in enumerate(proposed)
        ]
        adata.var["ctp_original_var_name"] = original
        adata.var_names = proposed
        source = column
    else:
        source = "var_names"

    after_names = [str(value) for value in adata.var_names]
    audit = {
        "schema_version": "celltypepilot.gene-identity-contract.v1",
        "source": source,
        "applied": applied,
        "n_features": len(original),
        "marker_universe_size": len(target),
        "marker_overlap_before": before,
        "marker_overlap_after": len(target.intersection(after_names)),
        "ambiguous_duplicate_symbol_count": len(duplicate_symbols),
        "ambiguous_duplicate_symbols": duplicate_symbols[:100],
        "original_var_names_sha256": _sha256_lines(original),
        "active_var_names_sha256": _sha256_lines(after_names),
        "claim_boundary": (
            "Identifier normalization makes existing expression evidence addressable; "
            "it does not validate marker specificity or biological accuracy."
        ),
    }
    adata.uns["celltypepilot_gene_identity"] = json.loads(json.dumps(audit))
    return audit


def restore_original_gene_identifiers(adata) -> bool:
    """Restore input feature identifiers before returning/writing annotated data."""
    column = "ctp_original_var_name"
    if column not in adata.var:
        return False
    original = adata.var[column].astype(str).to_numpy()
    adata.var_names = original
    return True
