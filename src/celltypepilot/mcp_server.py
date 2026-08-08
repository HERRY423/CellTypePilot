"""Native MCP facade for CellTypePilot.

This server exposes deterministic CellTypePilot operations to Agent hosts.
It does not plan experiments, choose biological conclusions, or bypass the
same fail-closed gates used by the CLI.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class MCPServerError(RuntimeError):
    """Raised when the optional MCP runtime is unavailable."""


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "to_dict"):
        try:
            return _jsonable(value.to_dict(orient="records"))
        except TypeError:
            return _jsonable(value.to_dict())
    if hasattr(value, "tolist"):
        return _jsonable(value.tolist())
    if hasattr(value, "item"):
        try:
            return _jsonable(value.item())
        except (TypeError, ValueError):
            pass
    return value


def tool_inspect_h5ad(
    input_path: str,
    cluster_key: str | None = None,
    embedding_key: str | None = None,
) -> dict:
    """Inspect an h5ad file and report support boundaries."""
    from .data_adapter import inspect_adata

    return _jsonable(inspect_adata(input_path, cluster_key, embedding_key))


def tool_list_marker_scope(
    species: str = "human",
    tissue: str | None = None,
    packs: list[str] | None = None,
) -> dict:
    """List supported tissues or marker definitions for one tissue."""
    from .data_adapter import get_all_markers_for_tissue, load_marker_atlas

    atlas = load_marker_atlas(species)
    warnings: list[str] = []
    if packs:
        from .pack_manager import merge_marker_atlas, resolve_extension_packs

        records, pack_warnings = resolve_extension_packs(packs, species)
        atlas, merge_warnings = merge_marker_atlas(atlas, records, species)
        warnings.extend(pack_warnings)
        warnings.extend(merge_warnings)
    if tissue is None:
        return {
            "species": species,
            "available_tissues": sorted(atlas.get("tissues", {})),
            "warnings": warnings,
        }
    markers = get_all_markers_for_tissue(atlas, tissue)
    return {"species": species, "tissue": tissue, "markers": _jsonable(markers), "warnings": warnings}


def tool_atlas_governance(include_packs: bool = True) -> dict:
    """Return the offline atlas governance report."""
    from .atlas_governance import build_atlas_governance_report

    return _jsonable(build_atlas_governance_report(include_packs=include_packs))


def tool_uncertainty_language(
    calibration_policy_path: str | None = None,
    uses_reference: bool = False,
) -> dict:
    """Return CellTypePilot's score/confidence/probability claim contract."""
    from .uncertainty import build_uncertainty_language_manifest

    policy = None
    if calibration_policy_path:
        policy = json.loads(Path(calibration_policy_path).read_text(encoding="utf-8"))
    return _jsonable(
        build_uncertainty_language_manifest(
            calibration_policy=policy,
            uses_reference=uses_reference,
        )
    )


def tool_read_manifest(output_dir: str) -> dict:
    """Read a CellTypePilot output manifest and related audit state if present."""
    from .provenance import load_manifest

    root = Path(output_dir)
    manifest_path = root / "manifest.json"
    payload = {"manifest": load_manifest(manifest_path)}
    for name in ("artifact_status.json", "annotation_overrides.json"):
        path = root / name
        if path.is_file():
            payload[name.removesuffix(".json")] = json.loads(path.read_text(encoding="utf-8"))
    audit_path = root / "annotation_audit_log.jsonl"
    if audit_path.is_file():
        lines = audit_path.read_text(encoding="utf-8").splitlines()
        payload["audit_log_tail"] = [json.loads(line) for line in lines[-20:] if line.strip()]
    return _jsonable(payload)


def tool_critic_review_h5ad(
    input_path: str,
    cluster_key: str,
    focus_cluster: str,
    species: str | None = None,
    tissue: str | None = None,
) -> dict:
    """Run deterministic critic review for one cluster in a local h5ad file."""
    from .data_adapter import load_h5ad
    from .orchestrator import critic_review

    adata = load_h5ad(input_path)
    return _jsonable(
        critic_review(
            adata,
            cluster_key,
            focus_cluster,
            species=species,
            tissue=tissue,
        )
    )


def tool_annotate_clusters(
    input_path: str,
    cluster_key: str,
    output_dir: str,
    species: str | None = None,
    tissue: str | None = None,
    embedding_key: str | None = None,
    no_figures: bool = False,
) -> dict:
    """Run the bounded annotation pipeline and write reviewable artifacts."""
    from .orchestrator import run_annotation_pipeline

    result = run_annotation_pipeline(
        input_path=input_path,
        cluster_key=cluster_key,
        output_dir=output_dir,
        species=species,
        tissue=tissue,
        embedding_key=embedding_key,
        no_figures=no_figures,
    )
    return _jsonable(
        {
            "species": result["species"],
            "tissue": result["tissue"],
            "critic_summary": result["critic_summary"],
            "validation_scope": result.get("validation_scope"),
            "paths": result["paths"],
            "manifest": result["manifest"],
        }
    )


def build_mcp_server():
    """Create the FastMCP server, or raise an actionable dependency error."""
    try:
        from fastmcp import FastMCP
    except ImportError as exc:
        raise MCPServerError(
            "CellTypePilot MCP support requires the optional MCP runtime. "
            "Install with: pip install -e .[mcp]"
        ) from exc

    mcp = FastMCP("CellTypePilot")
    mcp.tool()(tool_inspect_h5ad)
    mcp.tool()(tool_list_marker_scope)
    mcp.tool()(tool_atlas_governance)
    mcp.tool()(tool_uncertainty_language)
    mcp.tool()(tool_read_manifest)
    mcp.tool()(tool_critic_review_h5ad)
    mcp.tool()(tool_annotate_clusters)
    return mcp


def main() -> None:
    """Run the CellTypePilot MCP server over stdio."""
    server = build_mcp_server()
    server.run()


if __name__ == "__main__":
    main()
