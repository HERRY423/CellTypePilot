from __future__ import annotations

import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from celltypepilot.calibration_split import build_calibration_split
from celltypepilot.governance_freeze import (
    GovernanceFreezeError,
    build_governance_freeze,
    governed_paths,
    verify_governance_freeze,
)
from celltypepilot.lineage_coverage import build_selector_lineage_audit
from celltypepilot.pack_manager import resolve_extension_packs, validate_pack
from celltypepilot.validation_domains import load_validation_domains

REPO_ROOT = Path(__file__).parents[1]
PACK_ROOT = REPO_ROOT / "src/celltypepilot/data/packs"


def test_three_depth_atlas_contracts_resolve_to_valid_first_party_packs():
    domains = load_validation_domains()["domains"]

    assert set(domains) == {"lung", "gut_ibd", "tumor_microenvironment"}
    for domain in domains.values():
        contract = domain["atlas_contract"]
        assert contract["status"] == "scope_complete_evidence_not_claim_ready"
        assert len(contract["required_lineages"]) >= 4
        records, warnings = resolve_extension_packs(contract["required_packs"], "human")
        assert warnings == []
        assert {record["name"] for record in records} == set(contract["required_packs"])
        for pack in contract["required_packs"]:
            assert validate_pack(PACK_ROOT / pack) == []


def test_calibration_split_is_truth_blind_and_donor_disjoint(tmp_path):
    data = tmp_path / "cohort.h5ad"
    obs = pd.DataFrame(
        {
            "donor": ["d1", "d1", "d2", "d2", "d3", "d3", "d4", "d4"],
            "truth_must_not_be_read": ["A", "A", "B", "B", "A", "A", "B", "B"],
        },
        index=[f"c{i}" for i in range(8)],
    )
    ad.AnnData(np.ones((8, 2)), obs=obs, var=pd.DataFrame(index=["G1", "G2"])).write_h5ad(data)
    registry = tmp_path / "registry.json"
    registry.write_text(
        json.dumps(
            {
                "cohorts": [
                    {
                        "cohort_id": "cohort",
                        "dataset_version_id": "version-1",
                        "constant_study_id": "study",
                        "donor_namespace": "study",
                        "local_path": data.name,
                        "metadata": {"donor_key": "donor"},
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    result = build_calibration_split(registry, tmp_path / "split", calibration_fraction=0.25)
    assignments = result["assignments"]

    assert set(assignments["role"]) == {"calibration", "evaluation"}
    assert result["manifest"]["truth_access"] == "prohibited_not_read"
    calibration = set(assignments.loc[assignments["role"] == "calibration", "donor_unit"])
    evaluation = set(assignments.loc[assignments["role"] == "evaluation", "donor_unit"])
    assert calibration.isdisjoint(evaluation)
    assert len(evaluation) >= 2


def test_lung_selector_audit_demonstrates_four_lineages_without_truth(tmp_path):
    labels = {
        "0": "pulmonary alveolar type 2 cell",
        "1": "capillary endothelial cell",
        "2": "fibroblast",
        "3": "natural killer cell",
    }
    cluster_rows = []
    prediction_rows = []
    for cluster, label in labels.items():
        for offset in range(2):
            cell_id = f"c{cluster}-{offset}"
            cluster_rows.append({"cell_id": cell_id, "cluster": cluster})
            for method in ("celltypist", "popv", "singler"):
                prediction_rows.append(
                    {
                        "cell_id": cell_id,
                        "method": method,
                        "predicted_label": label,
                        "confidence": 0.9,
                    }
                )
    predictions = tmp_path / "predictions.csv"
    clusters = tmp_path / "clusters.csv"
    pd.DataFrame(prediction_rows).to_csv(predictions, index=False)
    pd.DataFrame(cluster_rows).to_csv(clusters, index=False)

    result = build_selector_lineage_audit(
        predictions, clusters, tmp_path / "audit", domain_id="lung"
    )

    assert result["manifest"]["status"] == "passed"
    assert set(result["manifest"]["observed_lineages"]) == {
        "epithelial",
        "immune",
        "stromal",
        "vascular",
    }
    assert result["manifest"]["truth_access"] == "prohibited_not_read"
    assert set(result["decisions"]["selective_decision"]) == {"accepted_leaf"}


def test_governance_freeze_detects_tampering(tmp_path):
    source_root = REPO_ROOT
    copied_root = tmp_path / "repo"
    for source in governed_paths(source_root):
        destination = copied_root / source.relative_to(source_root)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
    freeze = copied_root / "freeze.json"
    build_governance_freeze(freeze, root=copied_root, release_id="test-release")

    verified = verify_governance_freeze(freeze, root=copied_root)
    assert verified["status"] == "verified"

    target = copied_root / "src/celltypepilot/data/validation_domains.json"
    target.write_text(target.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(GovernanceFreezeError, match="changed after freeze"):
        verify_governance_freeze(freeze, root=copied_root)


def test_governance_freeze_is_cross_platform_line_ending_stable(tmp_path):
    source_root = REPO_ROOT
    copied_root = tmp_path / "repo"
    for source in governed_paths(source_root):
        destination = copied_root / source.relative_to(source_root)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    target = copied_root / "src/celltypepilot/pack_manager.py"
    lf_content = target.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    target.write_bytes(lf_content)
    freeze = copied_root / "freeze.json"
    payload = build_governance_freeze(freeze, root=copied_root, release_id="cross-platform-test")

    target.write_bytes(lf_content.replace(b"\n", b"\r\n"))
    verified = verify_governance_freeze(freeze, root=copied_root)

    assert payload["hash_semantics"] == "utf8_text_lf_v1"
    assert verified["status"] == "verified"
