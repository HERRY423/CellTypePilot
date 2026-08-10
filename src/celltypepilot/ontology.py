"""Cell Ontology (CL) service — cached live ontology resolution and checks.

The bundled atlases declare Cell Ontology identifiers (``cl_id``) for every
cell type. This module makes those identifiers checkable against the *live*
Cell Ontology instead of trusting static strings:

- ``celltypepilot ontology update`` downloads ``cl.obo`` from the OBO
  Foundry into a local cache (overridable for tests via
  ``CELLTYPEPILOT_ONTOLOGY_DIR``).
- ``celltypepilot ontology check`` validates every bundled and installed
  atlas against the cached ontology: unknown terms, obsolete terms (with
  ``replaced_by`` / ``consider`` hints), and label mismatches.

Design contracts (same spirit as the rest of the plugin):

- Offline-safe: without a cache every check degrades to an explicit
  "no ontology cache" report instead of crashing or silently passing.
- Honest scope: co-occurrence with ontology labels is checked, not
  biological marker validity. Mismatches are warnings, obsolete or unknown
  identifiers are errors.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import tempfile
import urllib.request
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from . import __version__

ONTOLOGY_ENV_VAR = "CELLTYPEPILOT_ONTOLOGY_DIR"
ONTOLOGY_FILENAME = "cl.obo"
METADATA_FILENAME = "ontology_meta.json"
CL_OBO_URL = "http://purl.obolibrary.org/obo/cl.obo"
USER_AGENT = f"CellTypePilot/{__version__} (https://github.com/HERRY423/CellTypePilot)"


class OntologyError(ValueError):
    """Raised when the ontology cache cannot be created, read, or parsed."""


def ontology_dir() -> Path:
    """Ontology cache directory (override via env var for tests)."""
    override = os.environ.get(ONTOLOGY_ENV_VAR)
    if override:
        return Path(override)
    return Path.home() / ".celltypepilot" / "ontology"


def ontology_cache_path() -> Path:
    return ontology_dir() / ONTOLOGY_FILENAME


def _metadata_path() -> Path:
    return ontology_dir() / METADATA_FILENAME


def download_ontology(force: bool = False, timeout: int = 120) -> dict:
    """Download cl.obo into the cache. Returns cache status metadata."""
    target = ontology_cache_path()
    if target.is_file() and not force:
        return ontology_cache_status()
    target.parent.mkdir(parents=True, exist_ok=True)
    temp_fd, temp_name = tempfile.mkstemp(suffix=".obo", dir=str(target.parent))
    os.close(temp_fd)
    temp_path = Path(temp_name)
    try:
        request = urllib.request.Request(CL_OBO_URL, headers={"User-Agent": USER_AGENT})
        with (
            urllib.request.urlopen(request, timeout=timeout) as response,
            temp_path.open("wb") as handle,
        ):
            shutil.copyfileobj(response, handle)
        if temp_path.stat().st_size < 100_000:
            raise OntologyError("Downloaded cl.obo is suspiciously small; aborting")
        shutil.move(str(temp_path), target)
    except OntologyError:
        temp_path.unlink(missing_ok=True)
        raise
    except (urllib.error.URLError, OSError) as exc:
        temp_path.unlink(missing_ok=True)
        raise OntologyError(f"Failed to download Cell Ontology: {exc}") from exc
    metadata = {
        "source": CL_OBO_URL,
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "size_bytes": target.stat().st_size,
    }
    _metadata_path().write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return ontology_cache_status()


def ontology_cache_status() -> dict:
    """Report cache presence, age, and provenance without network access."""
    target = ontology_cache_path()
    if not target.is_file():
        return {
            "cached": False,
            "path": str(target),
            "detail": "No ontology cache. Run: celltypepilot ontology update",
        }
    status: dict = {
        "cached": True,
        "path": str(target),
        "size_bytes": target.stat().st_size,
    }
    try:
        metadata = json.loads(_metadata_path().read_text(encoding="utf-8"))
        status.update(metadata)
    except (OSError, json.JSONDecodeError):
        status["detail"] = "Cache present but metadata missing"
    return status


# ──────────────────────────────────────────────
# Minimal OBO parsing (Term stanzas only)
# ──────────────────────────────────────────────


@dataclass
class OntologyTerm:
    cl_id: str
    name: str = ""
    is_obsolete: bool = False
    replaced_by: list[str] = field(default_factory=list)
    consider: list[str] = field(default_factory=list)
    parents: list[str] = field(default_factory=list)
    synonyms: list[str] = field(default_factory=list)


@dataclass
class OntologyService:
    terms: dict[str, OntologyTerm]
    source: str = "cache"

    def resolve(self, cl_id: str) -> OntologyTerm | None:
        return self.terms.get(str(cl_id).strip())

    def label_of(self, cl_id: str) -> str:
        term = self.resolve(cl_id)
        return term.name if term else ""


def parse_obo(path: str | Path) -> dict[str, OntologyTerm]:
    """Parse Term stanzas from an .obo file into a term index."""
    terms: dict[str, OntologyTerm] = {}
    current: OntologyTerm | None = None
    in_term = False
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line == "[Term]":
                in_term = True
                current = None
                continue
            if line.startswith("[") and line.endswith("]"):
                in_term = False
                current = None
                continue
            if not in_term or not line:
                continue
            if line.startswith("id: "):
                current = OntologyTerm(cl_id=line[4:].strip())
                terms[current.cl_id] = current
                continue
            if current is None:
                continue
            if line.startswith("name: "):
                current.name = line[6:].strip()
            elif line.startswith("is_obsolete: true"):
                current.is_obsolete = True
            elif line.startswith("replaced_by: "):
                current.replaced_by.append(line[13:].split("!")[0].strip())
            elif line.startswith("consider: "):
                current.consider.append(line[10:].split("!")[0].strip())
            elif line.startswith("is_a: "):
                current.parents.append(line[6:].split("!")[0].strip())
            elif line.startswith("synonym: "):
                match = re.match(r'synonym:\s*"(.+?)"', line)
                if match:
                    current.synonyms.append(match.group(1))
    if not terms:
        raise OntologyError(f"No [Term] stanzas parsed from {path}")
    return terms


def load_ontology() -> OntologyService:
    """Load the cached ontology. Fails closed when the cache is missing."""
    path = ontology_cache_path()
    if not path.is_file():
        raise OntologyError("Cell Ontology cache not found. Run: celltypepilot ontology update")
    try:
        terms = parse_obo(path)
    except (OSError, UnicodeDecodeError) as exc:
        raise OntologyError(f"Failed to read ontology cache: {exc}") from exc
    return OntologyService(terms=terms, source=str(path))


# ──────────────────────────────────────────────
# Atlas checks against the live ontology
# ──────────────────────────────────────────────


def _normalize_label(text: str) -> str:
    return re.sub(r"\s+", " ", str(text).replace("_", " ").strip().lower())


def _term_label_matches(term: OntologyTerm, cell_type_key: str) -> bool:
    expected = _normalize_label(cell_type_key)
    candidates = [_normalize_label(term.name)] + [_normalize_label(syn) for syn in term.synonyms]
    # Allow the atlas key to be a refinement of the ontology label
    # (e.g. "CD4+ T cell" vs "CD4-positive, alpha-beta T cell" is NOT a
    # lexical match; such cases stay visible as warnings for curators).
    return expected in candidates


def check_atlas_ontology(service: OntologyService, atlas: dict) -> list[dict]:
    """Check every cl_id in an atlas against the ontology.

    Returns a list of findings: ``{"severity", "path", "cl_id", "issue"}``.
    Severities: ``error`` (unknown/obsolete identifier) and ``warning``
    (lexical label mismatch — curators decide whether the refinement is
    intentional).
    """
    findings: list[dict] = []

    def walk(cell_types: dict, path: str) -> None:
        for name, info in cell_types.items():
            node_path = f"{path}/{name}"
            cl_id = str(info.get("cl_id", ""))
            if not re.fullmatch(r"CL:\d{7}", cl_id):
                findings.append(
                    {
                        "severity": "error",
                        "path": node_path,
                        "cl_id": cl_id,
                        "issue": f"malformed cl_id {cl_id!r}",
                    }
                )
                walk(info.get("subtypes", {}), node_path)
                continue
            term = service.resolve(cl_id)
            if term is None:
                findings.append(
                    {
                        "severity": "error",
                        "path": node_path,
                        "cl_id": cl_id,
                        "issue": "CL identifier not found in the current Cell Ontology",
                    }
                )
            elif term.is_obsolete:
                hints = []
                if term.replaced_by:
                    hints.append(f"replaced_by: {', '.join(term.replaced_by)}")
                if term.consider:
                    hints.append(f"consider: {', '.join(term.consider)}")
                findings.append(
                    {
                        "severity": "error",
                        "path": node_path,
                        "cl_id": cl_id,
                        "issue": f"obsolete CL term ({'; '.join(hints) or 'no replacement hint'})",
                    }
                )
            elif not _term_label_matches(term, name):
                findings.append(
                    {
                        "severity": "warning",
                        "path": node_path,
                        "cl_id": cl_id,
                        "issue": (
                            f"atlas key {name!r} does not lexically match ontology label "
                            f"{term.name!r} (may be an intentional refinement)"
                        ),
                    }
                )
            walk(info.get("subtypes", {}), node_path)

    for tissue, tissue_info in atlas.get("tissues", {}).items():
        walk(tissue_info.get("cell_types", {}), tissue)
    return findings


def summarize_findings(findings: list[dict], checked_nodes: int = 0) -> dict:
    errors = [item for item in findings if item["severity"] == "error"]
    warnings = [item for item in findings if item["severity"] == "warning"]
    return {
        "checked_nodes": checked_nodes or len(findings),
        "errors": len(errors),
        "warnings": len(warnings),
        "ok": not errors,
    }
