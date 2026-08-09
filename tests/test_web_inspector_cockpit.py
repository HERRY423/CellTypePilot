import json
import pandas as pd
import pytest
from pathlib import Path


def test_web_inspector_cockpit_full_flow(tmp_path, monkeypatch):
    """Test full Cockpit workflow: Checklist, Cluster History/Notes, Override Diff, Sign-off, and Packet Export."""
    from celltypepilot import web_inspector

    # Initialize mock test environment
    web_inspector._output_dir = tmp_path
    web_inspector._overrides = {}
    web_inspector._adata_cache = None
    web_inspector._evidence_cache = None

    # Create dummy evidence table and mock h5ad file
    evidence_csv = tmp_path / "evidence_table.csv"
    evidence_csv.write_text(
        "cluster,cell_type,cl_id,combined_score,critic_confidence,critic_flags,novelty_decision,novelty_score,n_cells\n"
        "0,CD4 T cell,CL:0000624,0.85,high,PASS,not_assessed,0.0,500\n"
        "1,B cell,CL:0000236,0.45,low,LOW_EVIDENCE,not_assessed,0.0,300\n",
        encoding="utf-8",
    )

    # Mock _load_data to return dummy DataFrame & mock AnnData with pandas obs
    class DummyAdata:
        obs = pd.DataFrame(
            {
                "ctp_cl_id": ["0", "1"],
                "ctp_confidence": ["high", "low"],
                "ctp_cell_type": ["CD4 T cell", "B cell"],
            }
        )

    evidence_df = pd.read_csv(evidence_csv)
    monkeypatch.setattr(web_inspector, "_load_data", lambda: (DummyAdata(), evidence_df))

    client = web_inspector.app.test_client()

    # 1. Test Dashboard Landing Route
    dash_resp = client.get("/")
    assert dash_resp.status_code == 200
    assert "Review Cockpit".encode("utf-8") in dash_resp.data or "判读驾驶舱".encode("utf-8") in dash_resp.data

    # 2. Test Checklist API (GET & POST)
    chk_resp = client.get("/api/checklist")
    assert chk_resp.status_code == 200
    chk_data = chk_resp.get_json()["checklist"]
    assert "readiness_pct" in chk_data
    assert chk_data["automated_items"]["critic_flags_reviewed"]["completed"] is False  # cluster 1 is LOW_EVIDENCE and unreviewed

    # Mark manual checklist item as completed
    update_chk = client.post(
        "/api/checklist",
        json={"item_key": "marker_literature_alignment", "completed": True},
    )
    assert update_chk.status_code == 200
    assert update_chk.get_json()["checklist"]["manual_items"]["marker_literature_alignment"]["completed"] is True

    # 3. Test Cluster History & Notes API
    hist_resp = client.get("/api/clusters/1/history")
    assert hist_resp.status_code == 200
    hist_data = hist_resp.get_json()
    assert hist_data["cluster"] == "1"
    assert hist_data["baseline"]["cell_type"] == "B cell"
    assert hist_data["review_status"] == "unreviewed"

    # Add Note to Cluster 1
    note_resp = client.post(
        "/api/clusters/1/note",
        json={"author": "Dr. Watson", "text": "qPCR confirmed MS4A1 expression"},
    )
    assert note_resp.status_code == 200
    notes = note_resp.get_json()["notes"]
    assert len(notes) == 1
    assert notes[0]["author"] == "Dr. Watson"

    # Update Cluster 1 Status to 'reviewed'
    st_resp = client.post("/api/clusters/1/status", json={"status": "reviewed"})
    assert st_resp.status_code == 200
    assert st_resp.get_json()["status"] == "reviewed"

    # Re-check Checklist: critic flags check should now pass because cluster 1 is reviewed!
    chk_resp2 = client.get("/api/checklist")
    chk_data2 = chk_resp2.get_json()["checklist"]
    assert chk_data2["automated_items"]["critic_flags_reviewed"]["completed"] is True

    # 4. Test Override Diff API
    # Add an override for cluster 1
    ov_resp = client.post(
        "/api/override",
        json={"cluster": "1", "new_type": "Memory B cell", "reason": "Expert review"},
    )
    assert ov_resp.status_code == 200

    diff_resp = client.get("/api/overrides/diff")
    assert diff_resp.status_code == 200
    diff = diff_resp.get_json()["diff"]
    assert diff["modified_clusters"] == 1
    assert diff["affected_cells"] == 300
    assert len(diff["diff_rows"]) == 2

    # 5. Test Sign-off Workflow API
    signoff_get = client.get("/api/signoff")
    assert signoff_get.status_code == 200
    assert signoff_get.get_json()["signoff"]["signed_off"] is False

    # Attempt sign-off with incomplete checklist without force flag
    signoff_fail = client.post(
        "/api/signoff",
        json={
            "reviewer_name": "Dr. Jane",
            "reviewer_role": "PI",
            "decision": "APPROVED",
            "notes": "Good results",
            "force": False,
        },
    )
    assert signoff_fail.status_code == 400
    assert "Checklist incomplete" in signoff_fail.get_json()["error"]

    # Submit forced sign-off
    signoff_pass = client.post(
        "/api/signoff",
        json={
            "reviewer_name": "Dr. Jane",
            "reviewer_role": "PI",
            "decision": "APPROVED",
            "notes": "Forced approval after inspection",
            "force": True,
        },
    )
    assert signoff_pass.status_code == 200
    so_cert = signoff_pass.get_json()["signoff"]
    assert so_cert["signed_off"] is True
    assert so_cert["decision"] == "APPROVED"
    assert (tmp_path / "review_signoff.json").is_file()

    # 6. Test Review Packet Export Endpoint
    packet_resp = client.get("/api/export/review-packet")
    assert packet_resp.status_code == 200
    assert packet_resp.headers["Content-Type"].startswith("application/json")
    packet = json.loads(packet_resp.data.decode("utf-8"))
    assert packet["schema_version"] == "celltypepilot.review-packet.v1"
    assert packet["signoff_certificate"]["signed_off"] is True
    assert len(packet["clusters"]) == 2
    assert packet["override_diff"]["modified_clusters"] == 1
