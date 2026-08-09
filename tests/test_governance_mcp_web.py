import builtins
import json
from pathlib import Path

import pytest

from celltypepilot.atlas_governance import build_atlas_governance_report
from celltypepilot.mcp_server import (
    MCPServerError,
    RESOURCE_NAMES,
    build_mcp_server,
    read_output_resource,
    tool_atlas_governance,
    tool_list_artifacts,
    tool_list_marker_scope,
    tool_read_artifact_status,
    tool_read_audit_log,
    tool_read_evidence,
    tool_read_manifest,
    tool_read_novelty_results,
    tool_review_add_override,
    tool_review_apply_overrides,
    tool_review_clear_overrides,
    tool_review_delete_override,
    tool_review_list_overrides,
    tool_uncertainty_language,
)


@pytest.fixture
def mcp_output_dir(tmp_path):
    """A minimal CellTypePilot output dir with manifest/evidence/novelty/h5ad."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    out = tmp_path / "output"
    out.mkdir()

    (out / "manifest.json").write_text(
        json.dumps(
            {
                "celltypepilot_version": "0.3.0",
                "mkg_version": "mkg-2026.08.1",
                "parameters": {
                    "cluster_key": "leiden",
                    "species": "human",
                    "tissue": "blood",
                },
                "outputs": {"evidence_table.csv": {"sha256": "abc", "size_bytes": 10}},
                "validation_scope": "draft_annotation",
            }
        ),
        encoding="utf-8",
    )

    pd.DataFrame(
        [
            {
                "cluster": "0",
                "cell_type": "T cell",
                "combined_score": 0.9,
                "critic_confidence": "high",
                "critic_flags": "PASS",
                "novelty_decision": "not_prioritized",
            },
            {
                "cluster": "1",
                "cell_type": "Unknown",
                "combined_score": 0.2,
                "critic_confidence": "needs_review",
                "critic_flags": "LOW_EVIDENCE",
                "novelty_decision": "ood_novel_candidate",
            },
        ]
    ).to_csv(out / "evidence_table.csv", index=False)

    pd.DataFrame(
        [
            {"cluster": "1", "novelty_decision": "ood_novel_candidate"},
            {"cluster": "0", "novelty_decision": "not_prioritized"},
        ]
    ).to_csv(out / "novelty_results.csv", index=False)

    adata = ad.AnnData(np.zeros((2, 2), dtype=np.float32))
    adata.var_names = ["A", "B"]
    adata.obs["leiden"] = ["0", "1"]
    adata.obs["ctp_cell_type"] = ["T cell", "Unknown"]
    adata.write(out / "data.annotated.h5ad")

    return out


def test_atlas_governance_report_is_machine_readable_and_bounded():
    report = build_atlas_governance_report(include_packs=False)
    assert report["schema_version"] == "celltypepilot.atlas-governance.v1"
    assert report["supported_annotation_species"] == ["human", "mouse"]
    assert report["aggregate"]["n_marker_relationships"] > 0
    assert "not biological validation" in report["claim_boundary"]
    builtin = report["atlas_assets"][0]
    assert builtin["provenance_validation"] == "passed"
    assert builtin["counts"]["n_cell_type_nodes"] > 0


def test_mcp_facade_tools_are_deterministic_and_fail_closed():
    governance = tool_atlas_governance(include_packs=False)
    assert governance["schema_version"] == "celltypepilot.atlas-governance.v1"

    uncertainty = tool_uncertainty_language()
    assert uncertainty["schema_version"] == "celltypepilot.uncertainty-language.v1"
    assert uncertainty["score_columns"]["evidence_score"].endswith("not_probability")

    novelty = tool_read_novelty_results("__missing_output_dir__")
    assert novelty["schema_version"] == "celltypepilot.novelty-ood.v1"
    assert novelty["status"] == "not_found"

    scope = tool_list_marker_scope(species="human")
    assert "blood" in scope["available_tissues"]

    try:
        tool_list_marker_scope(species="rat")
    except ValueError as exc:
        assert "supports scoring only human, mouse" in str(exc)
    else:
        raise AssertionError("unsupported species must fail closed through MCP facade")


def test_mcp_readers_fail_closed_without_artifacts(tmp_path):
    artifacts = tool_list_artifacts(str(tmp_path))
    assert artifacts["schema_version"] == "celltypepilot.artifacts.v1"
    assert artifacts["exists"] is True
    assert artifacts["available"]["manifest.json"] is False
    assert artifacts["available"]["evidence_table.csv"] is False

    manifest = tool_read_manifest(str(tmp_path))
    assert manifest["schema_version"] == "celltypepilot.manifest.v1"
    assert manifest["status"] == "not_found"

    evidence = tool_read_evidence(str(tmp_path))
    assert evidence["schema_version"] == "celltypepilot.evidence.v1"
    assert evidence["status"] == "not_found"

    novelty = tool_read_novelty_results(str(tmp_path))
    assert novelty["status"] == "not_found"


def test_tool_read_manifest_reads_audit_state(mcp_output_dir):
    from celltypepilot import mcp_server

    mcp_server._append_audit_event(mcp_output_dir, "override_saved", {"cluster": "0"})

    manifest = tool_read_manifest(str(mcp_output_dir))
    assert manifest["status"] == "available"
    assert manifest["manifest"]["validation_scope"] == "draft_annotation"
    assert manifest["manifest"]["parameters"]["species"] == "human"
    assert manifest["audit_log_tail"][-1]["event_type"] == "override_saved"


def test_tool_read_evidence_cluster_filter_and_claim_boundary(mcp_output_dir):
    all_rows = tool_read_evidence(str(mcp_output_dir))
    assert all_rows["status"] == "available"
    assert all_rows["n_rows"] == 2
    assert "not calibrated probabilities" in all_rows["claim_boundary"]

    single = tool_read_evidence(str(mcp_output_dir), cluster="1")
    assert single["status"] == "available"
    assert len(single["rows"]) == 1
    assert single["rows"][0]["cluster"] == "1"

    missing = tool_read_evidence(str(mcp_output_dir), cluster="99")
    assert missing["status"] == "cluster_not_found"


def test_mcp_resource_schema_exposes_read_only_artifacts(mcp_output_dir):
    assert set(RESOURCE_NAMES) == {
        "manifest",
        "evidence",
        "novelty",
        "artifact_status",
        "audit_log",
    }

    evidence = json.loads(read_output_resource(str(mcp_output_dir), "evidence"))
    assert evidence["schema_version"] == "celltypepilot.evidence.v1"
    assert evidence["status"] == "available"
    assert evidence["n_rows"] == 2

    manifest = json.loads(read_output_resource(str(mcp_output_dir), "manifest"))
    assert manifest["schema_version"] == "celltypepilot.manifest.v1"
    assert manifest["manifest"]["validation_scope"] == "draft_annotation"

    novelty = json.loads(read_output_resource(str(mcp_output_dir), "novelty"))
    assert novelty["decision_counts"]["ood_novel_candidate"] == 1
    assert "not validated new cell-type discovery" in novelty["claim_boundary"]

    status = json.loads(read_output_resource(str(mcp_output_dir), "artifact_status"))
    assert status["artifact_status"]["review_state"] == "current"

    unknown = json.loads(read_output_resource(str(mcp_output_dir), "bogus"))
    assert unknown["schema_version"] == "celltypepilot.resource.v1"
    assert unknown["status"] == "error"


def test_tool_read_audit_log_and_artifact_status(mcp_output_dir):
    from celltypepilot import mcp_server

    mcp_server._append_audit_event(mcp_output_dir, "override_saved", {"cluster": "0"})
    mcp_server._append_audit_event(mcp_output_dir, "override_saved", {"cluster": "1"})

    audit = tool_read_audit_log(str(mcp_output_dir), event_type="override_saved")
    assert audit["schema_version"] == "celltypepilot.audit-log.v1"
    assert audit["n_events"] == 2
    assert audit["events"][0]["schema_version"] == "celltypepilot.web-audit.v1"

    status = tool_read_artifact_status(str(mcp_output_dir))
    assert status["artifact_status"]["review_state"] == "current"
    assert "data.annotated.h5ad" in status["artifact_status"]["current_artifacts"]


def test_review_gates_agent_inferred_overrides(mcp_output_dir):
    denied = tool_review_add_override(
        str(mcp_output_dir),
        "0",
        "T cell",
        "agent guess",
        confirm_human_review=False,
    )
    assert denied["status"] == "error"
    assert "confirm_human_review" in denied["error"]
    assert tool_review_list_overrides(str(mcp_output_dir))["count"] == 0


def test_review_stage_delete_list_clear_flow(mcp_output_dir):
    saved = tool_review_add_override(
        str(mcp_output_dir),
        "0",
        "CD4 T cell",
        "reviewed markers",
        confirm_human_review=True,
    )
    assert saved["status"] == "saved"
    assert saved["total_overrides"] == 1

    listed = tool_review_list_overrides(str(mcp_output_dir))
    assert listed["count"] == 1
    assert listed["overrides"]["0"]["new_type"] == "CD4 T cell"

    deleted = tool_review_delete_override(str(mcp_output_dir), "0")
    assert deleted["status"] == "deleted"
    assert tool_review_list_overrides(str(mcp_output_dir))["count"] == 0

    cleared = tool_review_clear_overrides(str(mcp_output_dir))
    assert cleared["status"] == "cleared"


def test_review_add_requires_annotated_h5ad(tmp_path):
    result = tool_review_add_override(
        str(tmp_path),
        "0",
        "T cell",
        "human",
        confirm_human_review=True,
    )
    assert result["status"] == "error"
    assert "data.annotated.h5ad" in result["error"]


def test_review_apply_requires_human_confirmation_and_marks_stale(
    mcp_output_dir, monkeypatch
):
    tool_review_add_override(
        str(mcp_output_dir),
        "0",
        "CD4 T cell",
        "human adjudication",
        confirm_human_review=True,
    )

    denied = tool_review_apply_overrides(str(mcp_output_dir), confirm_human_review=False)
    assert denied["status"] == "error"
    assert "confirm_human_review" in denied["error"]

    monkeypatch.setattr(
        "celltypepilot.orchestrator.apply_overrides_to_h5ad",
        lambda h5ad_path, overrides: {
            "applied": 1,
            "skipped": 0,
            "total": 1,
            "backup": str(Path(h5ad_path).parent / "b.h5ad"),
            "details": [],
        },
    )
    applied = tool_review_apply_overrides(str(mcp_output_dir), confirm_human_review=True)
    assert applied["status"] == "applied"
    assert applied["result"]["applied"] == 1
    assert applied["artifact_status"]["review_state"] == "applied_overrides_artifacts_stale"
    assert "manifest.json" in applied["artifact_status"]["stale_artifacts"]
    assert applied["artifact_status"]["current_artifacts"] == ["data.annotated.h5ad"]

    # Applied overrides are no longer pending.
    assert tool_review_list_overrides(str(mcp_output_dir))["count"] == 0

    audit = tool_read_audit_log(str(mcp_output_dir))
    assert audit["events"][-1]["event_type"] == "overrides_applied"


def test_review_apply_fails_closed_without_pending(mcp_output_dir):
    result = tool_review_apply_overrides(str(mcp_output_dir), confirm_human_review=True)
    assert result["status"] == "error"
    assert "No pending overrides" in result["error"]


def test_build_mcp_server_raises_actionable_error_without_runtime(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "fastmcp":
            raise ImportError("blocked for test")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(MCPServerError, match=r"pip install -e \.\[mcp\]"):
        build_mcp_server()


def test_web_review_audit_log_and_stale_status(tmp_path, monkeypatch):
    from celltypepilot import web_inspector

    web_inspector._output_dir = tmp_path
    web_inspector._overrides = {}
    web_inspector._adata_cache = None
    web_inspector._evidence_cache = None
    monkeypatch.setattr(
        web_inspector,
        "_apply_overrides_to_h5ad",
        lambda: {"applied": 1, "skipped": 0, "total": 1, "backup": str(tmp_path / "b.h5ad")},
    )

    client = web_inspector.app.test_client()
    response = client.post(
        "/api/override",
        json={"cluster": "0", "new_type": "T cell", "reason": "reviewed markers"},
    )
    assert response.status_code == 200
    assert (tmp_path / "annotation_overrides.json").is_file()

    audit = client.get("/api/audit").get_json()
    assert audit["events"][-1]["event_type"] == "override_saved"

    apply_response = client.post("/api/overrides/apply")
    assert apply_response.status_code == 200
    payload = apply_response.get_json()
    assert payload["artifact_status"]["review_state"] == "applied_overrides_artifacts_stale"
    assert "manifest.json" in payload["artifact_status"]["stale_artifacts"]

    status = client.get("/api/artifact-status").get_json()["artifact_status"]
    assert status["review_state"] == "applied_overrides_artifacts_stale"
    overrides = json.loads((tmp_path / "annotation_overrides.json").read_text(encoding="utf-8"))
    assert overrides == {}