from __future__ import annotations

import asyncio
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from celltypepilot.agent_protocol import AGENT_DECISION_SCHEMA, validate_agent_decision
from celltypepilot.golden_workflow import (
    GoldenWorkflowError,
    _load_plan,
    annotate_from_plan,
    prepare_annotation,
    review_uncertain_clusters,
)
from celltypepilot.mcp_server import build_mcp_server
from celltypepilot.review_resign import collect_content_hashes, resign_review_outputs


def _write_lung_input(path: Path) -> None:
    obs = pd.DataFrame(
        {"cluster": pd.Categorical(["0", "0", "1", "1"])},
        index=["c1", "c2", "c3", "c4"],
    )
    var = pd.DataFrame(
        {"feature_name": ["CA4", "RGCC", "PECAM1", "EPCAM"]},
        index=["ENSG00000167434", "ENSG00000102760", "ENSG00000261371", "ENSG00000119888"],
    )
    ad.AnnData(X=np.ones((4, 4)), obs=obs, var=var).write_h5ad(path)


def test_prepare_writes_locked_executable_plan(tmp_path: Path):
    input_path = tmp_path / "lung.h5ad"
    output_dir = tmp_path / "run"
    _write_lung_input(input_path)

    result = prepare_annotation(
        str(input_path),
        str(output_dir),
        cluster_key="cluster",
        species="human",
        tissue="lung",
        marker_evidence_policy="database",
    )

    assert result["status"] == "ready"
    validate_agent_decision(result)
    assert result["schema_version"] == AGENT_DECISION_SCHEMA
    assert result["allowed_next_actions"] == ["annotate_from_plan"]
    plan = _load_plan(result["plan_path"])
    assert plan["cluster_key"] == "cluster"
    assert plan["input_sha256"] == hashlib.sha256(input_path.read_bytes()).hexdigest()
    assert plan["evidence_coverage"]["gene_identity"]["marker_overlap_after"] > 0

    raw = json.loads(Path(result["plan_path"]).read_text(encoding="utf-8"))
    raw["tissue"] = "brain"
    Path(result["plan_path"]).write_text(json.dumps(raw), encoding="utf-8")
    with pytest.raises(GoldenWorkflowError, match="hash mismatch"):
        _load_plan(result["plan_path"])


def test_prepare_locks_native_backend_config_and_reference(tmp_path: Path):
    input_path = tmp_path / "lung.h5ad"
    reference_path = tmp_path / "reference.h5ad"
    output_dir = tmp_path / "run"
    _write_lung_input(input_path)
    _write_lung_input(reference_path)
    reference = ad.read_h5ad(reference_path)
    reference.obs["cell_type"] = ["endothelial", "endothelial", "epithelial", "epithelial"]
    reference.write_h5ad(reference_path)
    native_config = tmp_path / "native.json"
    native_config.write_text(
        json.dumps(
            {
                "schema_version": "celltypepilot.native-backends.v1",
                "backends": [
                    {
                        "backend": "custom_reference",
                        "method": "correlation",
                        "reference_path": str(reference_path),
                        "label_key": "cell_type",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    result = prepare_annotation(
        str(input_path),
        str(output_dir),
        cluster_key="cluster",
        species="human",
        tissue="lung",
        native_backend_config_path=str(native_config),
    )
    plan = result["plan"]
    assert plan["native_backend_config_sha256"]
    assert str(reference_path.resolve()) in plan["native_backend_dependency_sha256"]

    with reference_path.open("ab") as handle:
        handle.write(b"changed")
    with pytest.raises(GoldenWorkflowError, match="dependency changed"):
        annotate_from_plan(result["plan_path"])


def test_review_queue_is_bounded_and_read_only(tmp_path: Path):
    pd.DataFrame(
        [
            {
                "cluster": "0",
                "decision": "accepted",
                "cell_type": "T cell",
                "candidate_cell_type": "T cell",
                "critic_confidence": "high",
                "critic_flags": "PASS",
            },
            {
                "cluster": "1",
                "decision": "abstain",
                "cell_type": "Unknown",
                "candidate_cell_type": "Macrophage",
                "critic_confidence": "needs_review",
                "critic_flags": "LOW_EVIDENCE",
                "abstain_reason": "insufficient marker support",
            },
        ]
    ).to_csv(tmp_path / "evidence_table.csv", index=False)

    result = review_uncertain_clusters(str(tmp_path))

    validate_agent_decision(result)
    assert result["status"] == "human_review_required"
    assert result["n_clusters_for_review"] == 1
    assert result["clusters"][0]["cluster"] == "1"
    assert result["clusters"][0]["allowed_actions"] == [
        "keep_unknown",
        "human_override_with_reason",
    ]
    assert result["human_action_required"] is True
    assert result["mutation_policy"] == "read_only_until_confirmed_human_finalize"


def test_default_mcp_surface_is_exactly_four_golden_tools(monkeypatch):
    monkeypatch.delenv("CELLTYPEPILOT_MCP_SURFACE", raising=False)
    server = build_mcp_server()
    names = {tool.name for tool in asyncio.run(server.list_tools())}
    assert names == {
        "prepare_annotation",
        "annotate_from_plan",
        "review_uncertain_clusters",
        "finalize_reviewed_annotations",
    }


def test_resign_rebuilds_manifest_and_excludes_mutable_envelope(tmp_path: Path):
    obs = pd.DataFrame({"ctp_cl_id": ["0"], "ctp_cell_type": ["T cell"]}, index=["c1"])
    ad.AnnData(X=np.ones((1, 1)), obs=obs, var=pd.DataFrame(index=["CD3E"])).write_h5ad(
        tmp_path / "data.annotated.h5ad"
    )
    pd.DataFrame(
        [
            {
                "cluster": "0",
                "cell_type": "T cell",
                "combined_score": 0.5,
                "critic_confidence": "medium",
                "critic_flags": "PASS",
            }
        ]
    ).to_csv(tmp_path / "evidence_table.csv", index=False)
    (tmp_path / "manifest.json").write_text(
        json.dumps({"parameters": {}, "outputs": {"old.txt": {"sha256": "stale"}}}),
        encoding="utf-8",
    )
    (tmp_path / "annotation_audit_log.jsonl").write_text("", encoding="utf-8")

    result = resign_review_outputs(tmp_path, signer="reviewer-1", regenerate=True)
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    gaps = json.loads((tmp_path / "evidence_gaps.json").read_text(encoding="utf-8"))
    signed = result["signature"]["content_sha256"]

    assert "old.txt" not in manifest["outputs"]
    assert "report_draft.html" in manifest["outputs"]
    assert "manifest.json" in signed
    assert "evidence_gaps.json" in signed
    assert gaps["n_unknown_clusters"] == 0
    assert "annotation_audit_log.jsonl" not in signed
    assert collect_content_hashes(tmp_path) == signed
