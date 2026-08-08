import json

from celltypepilot.atlas_governance import build_atlas_governance_report
from celltypepilot.mcp_server import (
    tool_atlas_governance,
    tool_list_marker_scope,
    tool_uncertainty_language,
)


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

    scope = tool_list_marker_scope(species="human")
    assert "blood" in scope["available_tissues"]

    try:
        tool_list_marker_scope(species="rat")
    except ValueError as exc:
        assert "supports scoring only human, mouse" in str(exc)
    else:
        raise AssertionError("unsupported species must fail closed through MCP facade")


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
