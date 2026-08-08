"""Governed user context for candidate expansion without prompt-driven acceptance."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from copy import deepcopy
from pathlib import Path

CONTEXT_SCHEMA_VERSION = "celltypepilot.context.v1"
REVIEW_STATUSES = {"draft", "reviewed"}
POLARITIES = ("positive", "negative")


class ContextPackError(ValueError):
    """Raised when user-supplied context is malformed or outside the run scope."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _as_marker_records(values: list, default_source: str) -> list[dict]:
    records = []
    seen = set()
    for value in values or []:
        if isinstance(value, str):
            gene = value.strip()
            source = default_source
        elif isinstance(value, dict):
            gene = str(value.get("gene", "")).strip()
            source = str(value.get("source", default_source)).strip()
        else:
            raise ContextPackError("Markers must be gene strings or objects")
        if not gene:
            raise ContextPackError("Context marker gene cannot be empty")
        if gene in seen:
            continue
        seen.add(gene)
        records.append({"gene": gene, "source": source or "user_asserted"})
    return records


def _normalize_hypothesis(raw: dict, axis: str, default_review: str) -> dict:
    if not isinstance(raw, dict):
        raise ContextPackError(f"{axis} hypotheses must be objects")
    label_key = "cell_type" if axis == "identity" else "state"
    label = str(raw.get(label_key, "")).strip()
    if not label:
        raise ContextPackError(f"{axis} hypothesis requires {label_key!r}")
    review_status = str(raw.get("review_status", default_review)).lower()
    if review_status not in REVIEW_STATUSES:
        raise ContextPackError(f"Invalid review_status {review_status!r}")
    default_source = str(raw.get("source", "user_asserted")).strip() or "user_asserted"
    positive = _as_marker_records(raw.get("positive_markers", []), default_source)
    negative = _as_marker_records(raw.get("negative_markers", []), default_source)
    if not positive:
        raise ContextPackError(f"{axis} hypothesis {label!r} needs positive_markers")

    normalized = {
        label_key: label,
        "positive_markers": positive,
        "negative_markers": negative,
        "review_status": review_status,
        "source": default_source,
    }
    if axis == "identity":
        cl_id = str(raw.get("cl_id", "")).strip()
        if cl_id and re.fullmatch(r"CL:\d{7}", cl_id) is None:
            raise ContextPackError(f"Invalid Cell Ontology ID {cl_id!r} for {label!r}")
        normalized["cl_id"] = cl_id
    else:
        parents = raw.get("parent_cell_types", [])
        if isinstance(parents, str):
            parents = [parents]
        normalized["parent_cell_types"] = [
            str(item).strip() for item in parents if str(item).strip()
        ]
    return normalized


def _rows_to_hypotheses(path: Path) -> tuple[list[dict], list[dict]]:
    identity: dict[str, dict] = {}
    states: dict[str, dict] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"axis", "label", "gene", "polarity"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ContextPackError(f"Custom marker CSV requires columns {sorted(required)}")
        for row_number, row in enumerate(reader, start=2):
            axis = str(row.get("axis", "")).strip().lower()
            polarity = str(row.get("polarity", "")).strip().lower()
            label = str(row.get("label", "")).strip()
            gene = str(row.get("gene", "")).strip()
            if (
                axis not in {"identity", "state"}
                or polarity not in POLARITIES
                or not label
                or not gene
            ):
                raise ContextPackError(f"Invalid custom marker row {row_number}")
            target = identity if axis == "identity" else states
            item = target.setdefault(
                label,
                {
                    "cell_type" if axis == "identity" else "state": label,
                    "cl_id": str(row.get("cl_id", "")).strip(),
                    "parent_cell_types": [
                        value.strip()
                        for value in str(row.get("parent_cell_types", "")).split(";")
                        if value.strip()
                    ],
                    "positive_markers": [],
                    "negative_markers": [],
                    "review_status": str(row.get("review_status", "draft")).strip() or "draft",
                    "source": str(row.get("source", "user_asserted")).strip() or "user_asserted",
                },
            )
            incoming_cl_id = str(row.get("cl_id", "")).strip()
            if item.get("cl_id") and incoming_cl_id and item["cl_id"] != incoming_cl_id:
                raise ContextPackError(
                    f"Conflicting CL IDs for {label!r} in custom marker row {row_number}"
                )
            if incoming_cl_id:
                item["cl_id"] = incoming_cl_id
            incoming_review = str(row.get("review_status", "draft")).strip() or "draft"
            if incoming_review.lower() not in REVIEW_STATUSES:
                raise ContextPackError(f"Invalid review_status in custom marker row {row_number}")
            if item["review_status"].lower() == "draft" or incoming_review.lower() == "draft":
                item["review_status"] = "draft"
            else:
                item["review_status"] = "reviewed"
            incoming_parents = [
                value.strip()
                for value in str(row.get("parent_cell_types", "")).split(";")
                if value.strip()
            ]
            item["parent_cell_types"] = list(
                dict.fromkeys([*item["parent_cell_types"], *incoming_parents])
            )
            item[f"{polarity}_markers"].append(
                {"gene": gene, "source": str(row.get("source", item["source"])).strip()}
            )
    return list(identity.values()), list(states.values())


def load_context_pack(
    context_text: str | None = None,
    context_file: str | Path | None = None,
    custom_markers_file: str | Path | None = None,
    species: str | None = None,
    tissue: str | None = None,
) -> dict:
    """Load, scope-check, normalize, and hash a governed context pack.

    Free text is retained for interpretation and provenance only. It never creates
    marker evidence or changes acceptance thresholds.
    """
    raw: dict = {}
    source_hashes: dict[str, str] = {}
    if context_file is not None:
        path = Path(context_file)
        if not path.exists():
            raise ContextPackError(f"Context file not found: {path}")
        try:
            raw = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise ContextPackError(f"Context file must be valid JSON: {exc}") from exc
        if not isinstance(raw, dict):
            raise ContextPackError("Context file root must be a JSON object")
        source_hashes["context_file_sha256"] = _sha256(path)

    schema_version = raw.get("schema_version", CONTEXT_SCHEMA_VERSION)
    if schema_version != CONTEXT_SCHEMA_VERSION:
        raise ContextPackError(
            f"Unsupported context schema {schema_version!r}; expected {CONTEXT_SCHEMA_VERSION!r}"
        )
    pack_species = str(raw.get("species", species or "")).strip().lower()
    pack_tissue = str(raw.get("tissue", tissue or "")).strip().lower()
    if species and pack_species and pack_species != species.lower():
        raise ContextPackError(
            f"Context species {pack_species!r} does not match run species {species!r}"
        )
    if tissue and pack_tissue and pack_tissue != tissue.lower():
        raise ContextPackError(
            f"Context tissue {pack_tissue!r} does not match run tissue {tissue!r}"
        )

    review_status = str(raw.get("review_status", "draft")).lower()
    if review_status not in REVIEW_STATUSES:
        raise ContextPackError(f"Invalid context review_status {review_status!r}")
    identity_value = raw.get("identity_hypotheses", [])
    state_value = raw.get("state_hypotheses", [])
    if not isinstance(identity_value, list) or not isinstance(state_value, list):
        raise ContextPackError("identity_hypotheses and state_hypotheses must be arrays")
    identity_raw = list(identity_value)
    state_raw = list(state_value)
    if custom_markers_file is not None:
        marker_path = Path(custom_markers_file)
        if not marker_path.exists():
            raise ContextPackError(f"Custom marker file not found: {marker_path}")
        extra_identity, extra_states = _rows_to_hypotheses(marker_path)
        identity_raw.extend(extra_identity)
        state_raw.extend(extra_states)
        source_hashes["custom_markers_sha256"] = _sha256(marker_path)

    free_text_parts = [str(raw.get("free_text", "")).strip(), str(context_text or "").strip()]
    free_text = "\n".join(part for part in free_text_parts if part)
    normalized = {
        "schema_version": CONTEXT_SCHEMA_VERSION,
        "species": pack_species,
        "tissue": pack_tissue,
        "condition": str(raw.get("condition", "")).strip(),
        "timepoint": str(raw.get("timepoint", "")).strip(),
        "anatomical_region": str(raw.get("anatomical_region", "")).strip(),
        "free_text": free_text,
        "review_status": review_status,
        "identity_hypotheses": [
            _normalize_hypothesis(item, "identity", review_status) for item in identity_raw
        ],
        "state_hypotheses": [
            _normalize_hypothesis(item, "state", review_status) for item in state_raw
        ],
        "source_hashes": source_hashes,
    }
    canonical = json.dumps(normalized, sort_keys=True, ensure_ascii=False).encode("utf-8")
    normalized["canonical_sha256"] = hashlib.sha256(canonical).hexdigest()
    return normalized


def resolve_atlas_tissue(requested_tissue: str, atlas: dict, context_pack: dict) -> str:
    """Resolve atlas scope without letting free text bypass a missing tissue panel."""
    if requested_tissue in atlas.get("tissues", {}):
        return requested_tissue
    if context_pack.get("identity_hypotheses"):
        return "general"
    raise ContextPackError(
        f"Tissue {requested_tissue!r} is not present in the bundled atlas. "
        "Provide a structured identity hypothesis or choose a supported tissue; "
        "free text alone cannot unlock an unsupported tissue."
    )


def context_manifest_parameters(context_pack: dict, enabled: bool) -> dict:
    """Return the non-secret, reproducibility-critical context manifest fields."""
    return {
        "context_enabled": enabled,
        "context_schema_version": context_pack["schema_version"],
        "context_sha256": context_pack["canonical_sha256"] if enabled else None,
        "context_source_hashes": context_pack["source_hashes"] if enabled else {},
        "context_review_status": context_pack["review_status"] if enabled else None,
        "context_free_text_present": bool(context_pack.get("free_text")),
        "context_identity_hypotheses": len(context_pack.get("identity_hypotheses", [])),
        "context_state_hypotheses": len(context_pack.get("state_hypotheses", [])),
    }


def merge_identity_hypotheses(markers: dict[str, dict], context_pack: dict) -> dict[str, dict]:
    """Merge structured identity hypotheses while preserving their trust boundary."""
    merged = deepcopy(markers)
    for hypothesis in context_pack.get("identity_hypotheses", []):
        label = hypothesis["cell_type"]
        target = merged.setdefault(
            label,
            {
                "positive_markers": [],
                "negative_markers": [],
                "cl_id": hypothesis.get("cl_id", ""),
                "marker_evidence": [],
            },
        )
        target.setdefault("atlas_positive_markers", list(target.get("positive_markers", [])))
        if hypothesis.get("cl_id"):
            existing = str(target.get("cl_id", ""))
            if existing and existing != hypothesis["cl_id"]:
                raise ContextPackError(
                    f"Context CL ID {hypothesis['cl_id']} conflicts with atlas ID {existing} for {label!r}"
                )
            target["cl_id"] = hypothesis["cl_id"]
        context_genes = []
        for polarity in POLARITIES:
            key = f"{polarity}_markers"
            for record in hypothesis[key]:
                gene = record["gene"]
                if gene not in target.setdefault(key, []):
                    target[key].append(gene)
                context_genes.append(gene)
                target.setdefault("marker_evidence", []).append(
                    {
                        "gene": gene,
                        "polarity": polarity,
                        "sources": [
                            {
                                "source_id": "user_context",
                                "name": record["source"],
                                "pmid": "",
                                "doi": "",
                                "url": "",
                            }
                        ],
                        "verification_status": "user_asserted",
                    }
                )
        target["context_origin"] = True
        previous_review = target.get("context_review_status")
        target["context_review_status"] = (
            "draft" if "draft" in {previous_review, hypothesis["review_status"]} else "reviewed"
        )
        target["context_positive_markers"] = list(
            dict.fromkeys(
                [
                    *target.get("context_positive_markers", []),
                    *(item["gene"] for item in hypothesis["positive_markers"]),
                ]
            )
        )
        target["context_negative_markers"] = list(
            dict.fromkeys(
                [
                    *target.get("context_negative_markers", []),
                    *(item["gene"] for item in hypothesis["negative_markers"]),
                ]
            )
        )
        target["context_marker_count"] = len(
            set(target["context_positive_markers"]) | set(target["context_negative_markers"])
        )
    return merged
