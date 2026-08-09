"""Regenerate derived review artifacts and re-sign the run after human edits.

After append-only overrides are applied, evidence/report/figures/manifest are
stale. This module regenerates them from the annotated h5ad + evidence table
and writes a fresh signature binding content hashes (fail closed if inputs missing).
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

RESIGN_SCHEMA = "celltypepilot.review-resign.v1"
ARTIFACT_STATUS = "artifact_status.json"
SIGN_FILE = "review_signature.json"

DERIVED = (
    "evidence_table.csv",
    "contrastive_evidence.csv",
    "evidence_gaps.json",
    "report_draft.html",
    "methodology_draft.txt",
    "manifest.json",
    "figures/",
)


class ReviewResignError(ValueError):
    """Raised when re-sign cannot complete fail-closed."""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def collect_content_hashes(output_dir: Path) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for name in (
        "data.annotated.h5ad",
        "evidence_table.csv",
        "contrastive_evidence.csv",
        "evidence_gaps.json",
        "state_results.csv",
        "novelty_results.csv",
        "manifest.json",
        "report_draft.html",
        "methodology_draft.txt",
        "qc_diagnostics.json",
    ):
        path = output_dir / name
        if path.is_file():
            hashes[name] = _sha256(path)
    figures = output_dir / "figures"
    if figures.is_dir():
        for path in sorted(figures.glob("*")):
            if path.is_file():
                hashes[f"figures/{path.name}"] = _sha256(path)
    return hashes


def mark_current_after_resign(output_dir: Path, signature: dict[str, Any]) -> dict[str, Any]:
    status = {
        "schema_version": "celltypepilot.artifact-status.v1",
        "updated_at": _utc_now(),
        "review_state": "current_after_resign",
        "stale_artifacts": [],
        "current_artifacts": list(signature.get("content_sha256", {}).keys()),
        "message": "Derived artifacts regenerated and re-signed after human review.",
        "signature_file": SIGN_FILE,
        "signature_sha256": signature.get("signature_sha256"),
    }
    path = output_dir / ARTIFACT_STATUS
    path.write_text(json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return status


def resign_review_outputs(
    output_dir: str | Path,
    *,
    signer: str = "human-reviewer",
    regenerate: bool = True,
) -> dict[str, Any]:
    """Re-hash (and optionally regenerate) derived artifacts; write review_signature.json.

    Regeneration uses existing annotated h5ad + evidence when full pipeline re-run
    is not requested: refreshes HTML report and methodology from current tables.
    """
    root = Path(output_dir)
    if not root.is_dir():
        raise ReviewResignError(f"Output directory not found: {root}")
    annotated = root / "data.annotated.h5ad"
    evidence = root / "evidence_table.csv"
    if not annotated.is_file():
        raise ReviewResignError("data.annotated.h5ad required before re-sign")
    if not evidence.is_file():
        raise ReviewResignError("evidence_table.csv required before re-sign")

    regenerated: list[str] = []
    if regenerate:
        try:
            import pandas as pd

            from .agent_evidence import build_actionable_evidence_gaps
            from .critic import generate_critic_summary
            from .provenance import (
                load_manifest,
                save_manifest,
                update_manifest_outputs,
            )
            from .reporter import generate_html_report, generate_methodology_text

            critic_results = pd.read_csv(evidence)
            evidence_gaps = build_actionable_evidence_gaps(critic_results)
            (root / "evidence_gaps.json").write_text(
                json.dumps(evidence_gaps, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            regenerated.append("evidence_gaps.json")
            critic_summary = generate_critic_summary(critic_results)
            try:
                manifest = load_manifest(root)
            except Exception:
                manifest = {
                    "schema_version": "celltypepilot.manifest.v1",
                    "parameters": {},
                    "outputs": {},
                }
            figure_paths = (
                sorted((root / "figures").glob("*")) if (root / "figures").is_dir() else []
            )
            generate_html_report(
                critic_results,
                critic_results,
                critic_summary,
                manifest,
                [str(p) for p in figure_paths],
                root,
            )
            regenerated.append("report_draft.html")
            method = generate_methodology_text(manifest, critic_summary, critic_results)
            (root / "methodology_draft.txt").write_text(method, encoding="utf-8")
            regenerated.append("methodology_draft.txt")
            # Rebuild the manifest after regenerating derived artifacts. Mutable
            # envelope files are intentionally excluded: the audit log is
            # append-only after signing, while signature/status describe the
            # bundle and cannot be members of their own hash set.
            manifest["outputs"] = {}
            update_manifest_outputs(
                manifest,
                root,
                exclude={
                    "annotation_audit_log.jsonl",
                    ARTIFACT_STATUS,
                    SIGN_FILE,
                    "annotation_workflow_state.json",
                },
            )
            save_manifest(manifest, root)
            regenerated.append("manifest.json")
        except Exception as exc:
            # Fail closed on regenerate errors only when regenerate requested.
            raise ReviewResignError(f"Regenerate failed: {exc}") from exc

    content_sha256 = collect_content_hashes(root)
    canonical = json.dumps(
        {"content_sha256": content_sha256, "output_dir": root.name},
        sort_keys=True,
        separators=(",", ":"),
    )
    signature_sha256 = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    payload = {
        "schema_version": RESIGN_SCHEMA,
        "signed_at_utc": _utc_now(),
        "signer": signer,
        "algorithm": "content-sha256-bundle",
        "signature_sha256": signature_sha256,
        "content_sha256": content_sha256,
        "regenerated": regenerated,
        "stale_policy": "cleared_after_successful_resign",
        "append_only_audit": True,
        "audit_binding": "excluded_mutable_append_only_log",
        "note": (
            "Human edits remain in the append-only annotation_audit_log.jsonl. "
            "The signature binds current regenerated artifacts and intentionally "
            "excludes that mutable log and its own envelope files."
        ),
    }
    (root / SIGN_FILE).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    status = mark_current_after_resign(root, payload)
    # Append audit event if log exists or always create
    audit = root / "annotation_audit_log.jsonl"
    event = {
        "schema_version": "celltypepilot.web-audit.v1",
        "timestamp": _utc_now(),
        "event_type": "review_resigned",
        "payload": {
            "signer": signer,
            "signature_sha256": signature_sha256,
            "regenerated": regenerated,
        },
    }
    with audit.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(event, ensure_ascii=False) + "\n")
    return {"signature": payload, "artifact_status": status}
