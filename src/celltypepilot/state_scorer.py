"""Independent cell-state scoring that cannot change canonical identity decisions."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .constants import (
    MARKER_FC_THRESHOLD,
    MARKER_FDR_THRESHOLD,
    MARKER_PCT_THRESHOLD,
    STATE_ATLAS_PATH,
)

STATE_DECISIONS = {"supported", "hypothesis", "abstain"}


class StateScoringError(ValueError):
    """Raised when state definitions or identity invariants are invalid."""


def _marker_genes(values: list) -> list[str]:
    return [str(item.get("gene", "")) if isinstance(item, dict) else str(item) for item in values]


def validate_state_atlas(atlas: dict) -> list[str]:
    """Validate state-marker provenance and reject incomplete relationships."""
    issues = []
    sources = atlas.get("sources", {})
    for state_name, definition in atlas.get("states", {}).items():
        expected = {
            (gene, polarity)
            for polarity, key in (
                ("positive", "positive_markers"),
                ("negative", "negative_markers"),
            )
            for gene in _marker_genes(definition.get(key, []))
        }
        records = definition.get("marker_evidence", [])
        observed = {(item.get("gene"), item.get("polarity")) for item in records}
        if expected != observed:
            issues.append(f"{state_name}: marker relationships and evidence records differ")
        for index, record in enumerate(records):
            if record.get("state") != state_name:
                issues.append(f"{state_name}/marker_evidence/{index}: state mismatch")
            if record.get("species") != definition.get("species", []):
                issues.append(f"{state_name}/marker_evidence/{index}: species mismatch")
            if record.get("tissues") != definition.get("tissues", []):
                issues.append(f"{state_name}/marker_evidence/{index}: tissue mismatch")
            if record.get("atlas_version") != atlas.get("version"):
                issues.append(f"{state_name}/marker_evidence/{index}: atlas version mismatch")
            if not record.get("source_ids"):
                issues.append(f"{state_name}/marker_evidence/{index}: no source_ids")
            for source_id in record.get("source_ids", []):
                source = sources.get(source_id, {})
                if not source:
                    issues.append(
                        f"{state_name}/marker_evidence/{index}: unknown source {source_id}"
                    )
                elif not all(source.get(field) for field in ("name", "pmid", "doi", "url")):
                    issues.append(
                        f"{state_name}/marker_evidence/{index}: incomplete source {source_id}"
                    )
            if not record.get("verification_status"):
                issues.append(f"{state_name}/marker_evidence/{index}: missing verification_status")
    return issues


def load_state_definitions(
    species: str, tissue: str, context_pack: dict | None = None
) -> list[dict]:
    """Load built-in and context state definitions for the declared scope."""
    atlas = json.loads(Path(STATE_ATLAS_PATH).read_text(encoding="utf-8"))
    issues = validate_state_atlas(atlas)
    if issues:
        raise StateScoringError("State atlas validation failed: " + "; ".join(issues[:5]))
    definitions = []
    for state_name, raw in atlas.get("states", {}).items():
        species_scope = raw.get("species", [])
        tissue_scope = raw.get("tissues", ["general"])
        if species_scope and species not in species_scope:
            continue
        if "general" not in tissue_scope and tissue not in tissue_scope:
            continue
        positive = _marker_genes(raw.get("positive_markers", []))
        negative = _marker_genes(raw.get("negative_markers", []))
        if species == "mouse":
            positive = [gene[:1].upper() + gene[1:].lower() for gene in positive]
            negative = [gene[:1].upper() + gene[1:].lower() for gene in negative]
        definitions.append(
            {
                "state": state_name,
                "parent_cell_types": list(raw.get("parent_cell_types", [])),
                "positive_markers": positive,
                "negative_markers": negative,
                "source": f"builtin:{atlas.get('version', 'unknown')}",
                "review_status": "reviewed",
                "marker_evidence": list(raw.get("marker_evidence", [])),
            }
        )
    for raw in (context_pack or {}).get("state_hypotheses", []):
        definitions.append(
            {
                "state": raw["state"],
                "parent_cell_types": list(raw.get("parent_cell_types", [])),
                "positive_markers": _marker_genes(raw.get("positive_markers", [])),
                "negative_markers": _marker_genes(raw.get("negative_markers", [])),
                "source": f"context:{raw.get('source', 'user_asserted')}",
                "review_status": raw.get("review_status", "draft"),
                "marker_evidence": [],
            }
        )
    return definitions


def _expression_fractions(
    adata: ad.AnnData,
    cluster: str,
    cluster_key: str,
    genes: list[str],
    layer: str | None,
) -> dict[str, float]:
    matrix = adata.layers[layer] if layer else adata.X
    mask = adata.obs[cluster_key].astype(str).to_numpy() == str(cluster)
    gene_index = {str(gene): index for index, gene in enumerate(adata.var_names)}
    fractions = {}
    for gene in genes:
        if gene not in gene_index:
            continue
        values = matrix[mask, gene_index[gene]]
        values = values.toarray().ravel() if sparse.issparse(values) else np.asarray(values).ravel()
        fractions[gene] = float(np.mean(values > 0)) if values.size else 0.0
    return fractions


def _parent_matches(candidate: str, parents: list[str]) -> bool:
    if not parents:
        return True
    normalized = candidate.casefold().replace("+", " positive ")
    return any(parent.casefold().replace("+", " positive ") in normalized for parent in parents)


def score_cell_states(
    adata: ad.AnnData,
    cluster_key: str,
    identity_results: pd.DataFrame,
    definitions: list[dict],
    layer: str | None = None,
    de_results: dict[str, pd.DataFrame] | None = None,
) -> pd.DataFrame:
    """Score state modules independently from identity and return one decision per cluster."""
    if de_results is None:
        from .marker_scorer import _extract_de_results

        de_results = _extract_de_results(adata, cluster_key)
    candidates = identity_results.set_index("cluster").get(
        "candidate_cell_type", pd.Series(dtype=str)
    )
    rows = []
    for cluster in identity_results["cluster"].astype(str):
        candidate = str(candidates.get(cluster, "Unknown"))
        de = de_results.get(cluster, pd.DataFrame())
        if de.empty:
            significant = pd.DataFrame(columns=["gene", "logfoldchange"])
        else:
            significant = de[
                (pd.to_numeric(de["pval_adj"], errors="coerce") <= MARKER_FDR_THRESHOLD)
                & (pd.to_numeric(de["logfoldchange"], errors="coerce") >= MARKER_FC_THRESHOLD)
            ]
        de_fc = dict(
            zip(
                significant.get("gene", pd.Series(dtype=str)).astype(str),
                pd.to_numeric(
                    significant.get("logfoldchange", pd.Series(dtype=float)), errors="coerce"
                ),
                strict=True,
            )
        )
        scored = []
        for definition in definitions:
            if not _parent_matches(candidate, definition.get("parent_cell_types", [])):
                continue
            positive = list(dict.fromkeys(definition.get("positive_markers", [])))
            negative = list(dict.fromkeys(definition.get("negative_markers", [])))
            present = [gene for gene in positive if gene in adata.var_names]
            missing = [gene for gene in positive if gene not in adata.var_names]
            fractions = _expression_fractions(adata, cluster, cluster_key, present, layer)
            expressed = [
                gene for gene, fraction in fractions.items() if fraction >= MARKER_PCT_THRESHOLD
            ]
            silent = [gene for gene in present if gene not in expressed]
            supporting = [gene for gene in expressed if gene in de_fc]
            negative_fractions = _expression_fractions(adata, cluster, cluster_key, negative, layer)
            negative_expressed = [
                gene
                for gene, fraction in negative_fractions.items()
                if fraction >= MARKER_PCT_THRESHOLD
            ]
            denominator = max(len(positive), 1)
            coverage = len(supporting) / denominator
            expression_coverage = len(expressed) / denominator
            mean_fc = float(np.mean([de_fc[gene] for gene in supporting])) if supporting else 0.0
            negative_conflict = len(negative_expressed) / max(len(negative), 1) if negative else 0.0
            score = (
                0.45 * coverage
                + 0.25 * expression_coverage
                + 0.30 * max(0.0, min(mean_fc / 2.0, 1.0))
                - 0.30 * negative_conflict
            )
            scored.append(
                {
                    "cluster": cluster,
                    "cell_state_candidate": definition["state"],
                    "state_score": round(max(0.0, min(score, 1.0)), 4),
                    "state_marker_coverage": round(coverage, 4),
                    "state_expression_coverage": round(expression_coverage, 4),
                    "state_mean_log2fc": round(mean_fc, 4),
                    "state_negative_conflict": round(negative_conflict, 4),
                    "n_state_expected_markers": len(positive),
                    "n_state_present_markers": len(present),
                    "n_state_missing_markers": len(missing),
                    "n_state_silent_markers": len(silent),
                    "n_state_supporting_markers": len(supporting),
                    "state_present_markers": ";".join(present),
                    "state_missing_markers": ";".join(missing),
                    "state_silent_markers": ";".join(silent),
                    "state_supporting_markers": ";".join(supporting),
                    "state_negative_markers_expressed": ";".join(negative_expressed),
                    "state_source": definition.get("source", "unknown"),
                    "state_review_status": definition.get("review_status", "draft"),
                }
            )
        if not scored:
            rows.append(_abstain_state(cluster, "no_state_definition_matches_identity"))
            continue
        ranked = sorted(scored, key=lambda row: row["state_score"], reverse=True)
        top = ranked[0]
        required_support = min(2, max(top["n_state_expected_markers"], 1))
        supported = (
            top["n_state_supporting_markers"] >= required_support
            and top["state_marker_coverage"] >= 0.3
            and top["state_expression_coverage"] >= 0.4
            and top["state_score"] >= 0.4
            and top["state_negative_conflict"] < 0.2
        )
        any_signal = (
            top["n_state_supporting_markers"] >= 1 or top["state_expression_coverage"] >= 0.25
        )
        ambiguous = (
            len(ranked) > 1
            and ranked[1]["state_score"] >= 0.3
            and (top["state_score"] - ranked[1]["state_score"] < 0.1)
        )
        flags = []
        if ambiguous:
            flags.append("AMBIGUOUS_STATE")
        if top["state_review_status"] != "reviewed":
            flags.append("UNREVIEWED_CONTEXT_STATE")
        if top["n_state_missing_markers"]:
            flags.append("STATE_MARKERS_MISSING")
        if top["state_negative_conflict"] >= 0.2:
            flags.append("STATE_NEGATIVE_CONFLICT")
        if supported and not ambiguous and top["state_review_status"] == "reviewed":
            decision = "supported"
        elif any_signal:
            decision = "hypothesis"
        else:
            decision = "abstain"
            top["cell_state_candidate"] = "Unknown"
        top["state_decision"] = decision
        top["state_confidence"] = (
            "high"
            if decision == "supported" and top["state_score"] >= 0.7
            else "medium"
            if decision == "supported"
            else "low"
            if decision == "hypothesis"
            else "needs_review"
        )
        top["state_flags"] = "; ".join(flags) if flags else "PASS"
        top["state_evidence"] = (
            f"support={top['n_state_supporting_markers']}/{top['n_state_expected_markers']}; "
            f"expressed={top['state_expression_coverage']:.0%}; "
            f"missing={top['n_state_missing_markers']}; silent={top['n_state_silent_markers']}"
        )
        rows.append(top)
    return pd.DataFrame(rows)


def _abstain_state(cluster: str, reason: str) -> dict:
    return {
        "cluster": cluster,
        "cell_state_candidate": "Unknown",
        "state_decision": "abstain",
        "state_score": 0.0,
        "state_confidence": "needs_review",
        "state_flags": reason,
        "state_evidence": reason,
        "state_source": "none",
        "state_review_status": "not_applicable",
    }


def attach_state_results(
    identity_results: pd.DataFrame, state_results: pd.DataFrame
) -> pd.DataFrame:
    """Attach state columns and prove identity columns are unchanged."""
    identity_columns = [
        column
        for column in ("cell_type", "cl_id", "candidate_cell_type", "decision", "abstain_reason")
        if column in identity_results
    ]
    before = identity_results.set_index("cluster")[identity_columns].copy()
    merged = identity_results.merge(state_results, on="cluster", how="left", validate="one_to_one")
    after = merged.set_index("cluster")[identity_columns]
    if not before.equals(after):
        raise StateScoringError("State scoring changed canonical identity columns")
    merged["cell_state_candidate"] = merged["cell_state_candidate"].fillna("Unknown")
    merged["state_decision"] = merged["state_decision"].fillna("abstain")
    merged["state_evidence"] = merged["state_evidence"].fillna("state_not_scored")
    merged["display_label"] = merged.apply(
        lambda row: (
            f"{row['cell_type']} · {row['cell_state_candidate']}"
            if row["state_decision"] == "supported"
            else str(row["cell_type"])
        ),
        axis=1,
    )
    return merged
