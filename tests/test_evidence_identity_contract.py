from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from celltypepilot.data_adapter import load_marker_atlas
from celltypepilot.identity_contract import (
    apply_gene_identity_contract,
    build_identity_resolver,
    collect_pack_identity_contract,
    compose_marker_definitions,
    resolve_identity_label,
    restore_original_gene_identifiers,
)
from celltypepilot.mcp_server import tool_evidence_coverage, tool_evidence_trace
from celltypepilot.pack_manager import merge_marker_atlas, resolve_extension_packs, validate_pack

PACK = Path(__file__).parents[1] / "src" / "celltypepilot" / "data" / "packs" / "lung_evidence_v0_1"


def test_gene_identity_uses_unique_symbol_column_and_restores_input_ids():
    var = pd.DataFrame(
        {"feature_name": ["PECAM1", "EPCAM", "DUP", "DUP"]},
        index=["ENSG1", "ENSG2", "ENSG3", "ENSG4"],
    )
    adata = ad.AnnData(X=np.ones((2, 4)), var=var)

    audit = apply_gene_identity_contract(adata, {"PECAM1", "EPCAM", "DUP"})

    assert audit["source"] == "feature_name"
    assert audit["marker_overlap_before"] == 0
    assert audit["marker_overlap_after"] == 2
    assert audit["ambiguous_duplicate_symbols"] == ["DUP"]
    assert list(adata.var_names) == ["PECAM1", "EPCAM", "ENSG3", "ENSG4"]
    assert restore_original_gene_identifiers(adata)
    assert list(adata.var_names) == ["ENSG1", "ENSG2", "ENSG3", "ENSG4"]


def test_lung_pack_composes_scope_and_only_uses_explicit_parent_fallbacks():
    assert validate_pack(PACK) == []
    records, warnings = resolve_extension_packs(["lung_evidence_v0_1"], "human")
    assert warnings == []
    record = records[0]
    contract = collect_pack_identity_contract([record])
    atlas = load_marker_atlas("human")
    atlas, merge_warnings = merge_marker_atlas(atlas, records, "human")
    assert isinstance(merge_warnings, list)
    markers, scope = compose_marker_definitions(atlas, "lung", pack_contract=contract)
    resolver = build_identity_resolver(atlas, scope["active_tissues"], contract)

    assert {"lung", "general", "blood"} <= set(scope["active_tissues"])
    assert "Capillary endothelial cell" in markers
    assert resolve_identity_label("natural killer cell", resolver)["canonical_label"] == "NK cell"
    parent = resolve_identity_label("endothelial cell of artery", resolver)
    assert parent["canonical_label"] == "Endothelial cell"
    assert parent["resolution"] == "explicit_safe_parent_fallback"
    unresolved = resolve_identity_label("invented ultra-fine cell", resolver)
    assert unresolved["resolution"] == "unresolved"
    assert unresolved["canonical_label"] == "invented ultra-fine cell"


def test_agent_preflight_reports_addressability_without_accuracy_claim(tmp_path):
    path = tmp_path / "symbols_in_var.h5ad"
    var = pd.DataFrame(
        {"feature_name": ["PECAM1", "EMCN", "CA4", "RGCC", "EPCAM", "PTPRC"]},
        index=[f"ENSG{i}" for i in range(6)],
    )
    ad.AnnData(X=np.ones((3, 6)), var=var).write_h5ad(path)

    result = tool_evidence_coverage(
        str(path), "human", "lung", packs=["lung_evidence_v0_1"], evidence_policy="database"
    )
    capillary = next(
        item for item in result["cell_types"] if item["cell_type"] == "Capillary endothelial cell"
    )
    assert result["gene_identity"]["marker_overlap_before"] == 0
    assert result["gene_identity"]["marker_overlap_after"] >= 6
    assert capillary["reachable_fraction"] == 1.0
    assert capillary["n_runtime_eligible_positive_markers"] == 4
    assert "not annotation accuracy" in result["claim_boundary"]

    trace = tool_evidence_trace(
        "capillary EC", species="human", tissue="lung", packs=["lung_evidence_v0_1"]
    )
    assert trace["status"] == "resolved"
    assert trace["resolution"]["canonical_label"] == "Capillary endothelial cell"
