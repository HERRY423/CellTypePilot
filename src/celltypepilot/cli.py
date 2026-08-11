"""CellTypePilot CLI — command-line interface."""

from __future__ import annotations

import json
from pathlib import Path

import typer
from rich.console import Console

from . import __version__

console = Console()
app = typer.Typer(
    name="celltypepilot",
    help="CellTypePilot — Local-first single-cell annotation review plugin",
    add_completion=False,
)


def _print_json(payload, **kwargs) -> None:
    """Write machine-readable JSON without Rich line wrapping."""
    typer.echo(json.dumps(payload, **kwargs))


def version_callback(value: bool):
    if value:
        console.print(f"CellTypePilot v{__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        "-v",
        callback=version_callback,
        is_eager=True,
        help="Show version and exit.",
    ),
):
    """CellTypePilot — Local-first single-cell annotation review plugin."""
    pass


# ──────────────────────────────────────────────
# doctor command
# ──────────────────────────────────────────────
@app.command()
def doctor(
    json_output: bool = typer.Option(
        False, "--json", help="Structured doctor report for Agent hosts"
    ),
):
    """Check environment: Python version, dependencies, capability level."""
    from .agent_lifecycle import doctor_report_to_dict
    from .doctor import print_doctor, run_doctor

    if json_output:
        report = run_doctor()
        print(json.dumps(doctor_report_to_dict(report), indent=2))
        if not report.python_ok or any(
            not dep.installed and dep.required for dep in report.dependencies
        ):
            raise typer.Exit(1)
        return
    print_doctor()


# ──────────────────────────────────────────────
# inspect command
# ──────────────────────────────────────────────
@app.command()
def inspect(
    input: str = typer.Option(..., "--input", "-i", help="Path to .h5ad file"),
    cluster_key: str | None = typer.Option(None, "--cluster-key", "-k", help="Cluster key in obs"),
    embedding_key: str | None = typer.Option(
        None, "--embedding-key", "-e", help="Embedding key in obsm"
    ),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """Inspect an h5ad file: detect species, tissue, clusters, embeddings, layers."""
    from .data_adapter import format_inspect_report, inspect_adata

    report = inspect_adata(input, cluster_key, embedding_key)

    if json_output:
        _print_json(report, indent=2)
    else:
        console.print(format_inspect_report(report))


# ──────────────────────────────────────────────
# annotate command (main pipeline)
# ──────────────────────────────────────────────
@app.command()
def annotate(
    input: str = typer.Option(..., "--input", "-i", help="Path to .h5ad file"),
    cluster_key: str = typer.Option(..., "--cluster-key", "-k", help="Cluster key in obs"),
    output_dir: str = typer.Option(".", "--output", "-o", help="Output directory"),
    species: str | None = typer.Option(
        None, "--species", "-s", help="Species: human/mouse (auto-detect if omitted)"
    ),
    tissue: str | None = typer.Option(
        None, "--tissue", "-t", help="Tissue context (e.g., blood, lung, brain)"
    ),
    embedding_key: str | None = typer.Option(
        None, "--embedding-key", "-e", help="Embedding key in obsm"
    ),
    layer: str | None = typer.Option(
        None, "--layer", help="Layer to use for expression (default: X)"
    ),
    json_output: bool = typer.Option(False, "--json", help="Output results as JSON"),
    no_figures: bool = typer.Option(False, "--no-figures", help="Skip figure generation"),
    reference: str | None = typer.Option(
        None, "--reference", "-r", help="Optional reference .h5ad with cell type labels"
    ),
    ref_label_key: str = typer.Option(
        "cell_type", "--ref-label", help="Cell type column in reference.obs"
    ),
    model_path: str | None = typer.Option(
        None, "--model", "-m", help="Optional CellTypist model path"
    ),
    backend: str = typer.Option(
        "auto", "--backend", "-b", help="Reference backend: auto/celltypist/scanvi/knn/correlation"
    ),
    marker_weight: float = typer.Option(
        0.5, "--marker-weight", help="Legacy diagnostic ensemble weight; never drives identity"
    ),
    no_ensemble: bool = typer.Option(
        False, "--no-ensemble", help="Skip the legacy diagnostic marker/reference ensemble"
    ),
    allow_unverified_reference: bool = typer.Option(
        False,
        "--allow-unverified-reference",
        help="Explicitly override the reference provenance safety gate; recorded in manifest",
    ),
    marker_evidence_policy: str = typer.Option(
        "database",
        "--marker-evidence-policy",
        help="database, literature, edge_verified, or primary; stricter policies exclude unverified edges",
    ),
    calibration_policy: str | None = typer.Option(
        None,
        "--calibration-policy",
        help="Optional abstention policy JSON fitted on a separate calibration dataset",
    ),
    context: str | None = typer.Option(
        None,
        "--context",
        help="Free biological context; recorded for interpretation but never treated as evidence",
    ),
    context_file: str | None = typer.Option(
        None,
        "--context-file",
        help="Governed celltypepilot.context.v1 JSON with structured identity/state hypotheses",
    ),
    custom_markers: str | None = typer.Option(
        None,
        "--custom-markers",
        help="CSV marker hypotheses with axis,label,gene,polarity columns",
    ),
    no_states: bool = typer.Option(
        False,
        "--no-states",
        help="Disable the independent State Lens without changing identity scoring",
    ),
    pack: list[str] = typer.Option(
        [],
        "--pack",
        help="Installed extension pack name to merge (repeatable), e.g. --pack premium",
    ),
    doublet_table: Path | None = typer.Option(
        None,
        "--doublet-table",
        help="External doublet tool CSV (diagnostic only; cannot rescue identity)",
    ),
    ambient_table: Path | None = typer.Option(
        None,
        "--ambient-table",
        help="External ambient-RNA tool CSV (diagnostic only; cannot rescue identity)",
    ),
    doublet_score_column: str = typer.Option(
        "doublet_score", "--doublet-score-column", help="Score column in --doublet-table"
    ),
    ambient_score_column: str = typer.Option(
        "ambient_score", "--ambient-score-column", help="Score column in --ambient-table"
    ),
    doublet_flag_column: str | None = typer.Option(
        None, "--doublet-flag-column", help="Optional boolean/class column in --doublet-table"
    ),
    ambient_flag_column: str | None = typer.Option(
        None, "--ambient-flag-column", help="Optional boolean/class column in --ambient-table"
    ),
    candidate_table: list[Path] = typer.Option(
        [],
        "--candidate-table",
        help=(
            "Repeatable CSV/JSON from CellTypist, popV, SingleR, scANVI, "
            "custom_reference, or optional LLM"
        ),
    ),
    native_backends: Path | None = typer.Option(
        None,
        "--native-backends",
        help=(
            "celltypepilot.native-backends.v1 JSON; executes configured CellTypist, "
            "popV, SingleR, scANVI, custom-reference, and opt-in LLM runners"
        ),
    ),
    decision_policy: Path | None = typer.Option(
        None,
        "--decision-policy",
        help="Optional hierarchical selective-decision policy JSON; not a calibration artifact",
    ),
):
    """Run candidate backends -> hierarchical selection -> critic -> report."""
    from .orchestrator import PipelineError, run_annotation_pipeline

    def _progress(step: int, total: int, message: str):
        console.print(f"[bold blue]Step {step}/{total}:[/bold blue] {message}")

    try:
        result = run_annotation_pipeline(
            input_path=input,
            cluster_key=cluster_key,
            output_dir=output_dir,
            species=species,
            tissue=tissue,
            embedding_key=embedding_key,
            layer=layer,
            no_figures=no_figures,
            reference_path=reference,
            ref_label_key=ref_label_key,
            model_path=model_path,
            reference_backend=backend,
            marker_weight=marker_weight,
            no_ensemble=no_ensemble,
            allow_unverified_reference=allow_unverified_reference,
            marker_evidence_policy=marker_evidence_policy,
            calibration_policy_path=calibration_policy,
            context_text=context,
            context_file_path=context_file,
            custom_markers_path=custom_markers,
            enable_states=not no_states,
            packs=pack or None,
            doublet_table_path=doublet_table,
            ambient_table_path=ambient_table,
            doublet_score_column=doublet_score_column,
            ambient_score_column=ambient_score_column,
            doublet_flag_column=doublet_flag_column,
            ambient_flag_column=ambient_flag_column,
            candidate_artifact_paths=[str(path) for path in candidate_table] or None,
            native_backend_config_path=native_backends,
            decision_policy_path=decision_policy,
            progress=_progress,
        )
    except PipelineError as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1) from e

    critic_summary = result["critic_summary"]
    console.print(f"  Detected species: [cyan]{result['species']}[/cyan]")
    console.print(f"  Tissue: [cyan]{result['tissue']}[/cyan]")
    if result["embedding_key"]:
        console.print(f"  Embedding: [cyan]{result['embedding_key']}[/cyan]")
    else:
        console.print("[yellow]  No embedding found. Figures were skipped.[/yellow]")
    console.print(f"  Annotated {len(result['critic_results'])} clusters")
    console.print(
        f"  Passed: [green]{critic_summary['pass']}[/green] | "
        f"Flagged: [red]{critic_summary['flagged']}[/red]"
    )
    qc = result.get("qc_diagnostics") or {}
    if qc:
        console.print(
            f"  QC diagnostics: [cyan]{qc.get('rollup_status')}[/cyan] "
            f"(identity rescue forbidden; missing axes = not_assessed)"
        )
    if critic_summary.get("narrative"):
        console.print(f"  [dim]{critic_summary['narrative']}[/dim]")
    if "state_decision" in result["critic_results"]:
        state_counts = result["critic_results"]["state_decision"].value_counts().to_dict()
        console.print(
            "  State Lens: "
            + " | ".join(f"{decision}={count}" for decision, count in state_counts.items())
        )
    console.print(f"  Generated {len(result['figure_paths'])} figures")
    for label, path in result["paths"].items():
        console.print(f"  {label.replace('_', ' ').title()}: {path}")

    # JSON output
    if json_output:
        output_json = {
            "annotations": result["critic_results"].to_dict(orient="records"),
            "critic_summary": critic_summary,
            "manifest": result["manifest"],
        }
        _print_json(output_json, indent=2, default=str)

    console.print("\n[bold green]Done![/bold green] CellTypePilot annotation complete.")
    console.print(f"Output directory: {result['output_path'].resolve()}")


# ──────────────────────────────────────────────
# critic command (re-review a specific cluster)
# ──────────────────────────────────────────────
@app.command()
def critic(
    input: str = typer.Option(..., "--input", "-i", help="Path to .h5ad file"),
    cluster_key: str = typer.Option(..., "--cluster-key", "-k", help="Cluster key in obs"),
    focus: str = typer.Option(..., "--focus", "-f", help="Cluster ID to deep-review"),
    species: str | None = typer.Option(None, "--species", "-s", help="Species"),
    tissue: str | None = typer.Option(None, "--tissue", "-t", help="Tissue context"),
):
    """Deep-review a specific cluster flagged by the critic."""
    from .critic import run_critic
    from .data_adapter import detect_species, detect_tissue, load_h5ad, load_marker_atlas
    from .marker_scorer import compute_marker_scores, generate_annotation_summary

    adata = load_h5ad(input)
    if species is None:
        species = detect_species(adata)
    if tissue is None:
        tissue = detect_tissue(adata) or "general"

    try:
        atlas = load_marker_atlas(species)
    except ValueError as exc:
        console.print(f"[red]Error: {exc}[/red]")
        raise typer.Exit(1) from exc
    from .data_adapter import get_all_markers_for_tissue

    markers = get_all_markers_for_tissue(atlas, tissue)

    scores = compute_marker_scores(adata, cluster_key, markers)
    summary = generate_annotation_summary(scores, cluster_key)

    # Filter to focus cluster
    focus_rows = summary[summary["cluster"] == focus]
    if focus_rows.empty:
        console.print(f"[red]Cluster '{focus}' not found in annotations.[/red]")
        raise typer.Exit(1)

    critic_results = run_critic(adata, cluster_key, focus_rows, atlas, tissue)

    console.print(f"\n[bold]Deep Review: Cluster {focus}[/bold]")
    console.print("=" * 50)
    for _, row in critic_results.iterrows():
        console.print(f"  Cell Type:   {row.get('cell_type', 'N/A')}")
        console.print(f"  CL ID:       {row.get('cl_id', 'N/A')}")
        console.print(f"  Evidence score: {row.get('combined_score', 0):.3f}")
        console.print(f"  Review level:   {row.get('critic_confidence', 'N/A')}")
        console.print(f"  Flags:       {row.get('critic_flags', 'PASS')}")
        if row.get("evidence_summary"):
            console.print(f"  [bold]Summary:[/bold]    {row.get('evidence_summary')}")
        console.print(f"  Evidence:    {row.get('critic_evidence', '')}")
        console.print(f"  Notes:       {row.get('critic_notes', '')}")

    # Show top-5 candidates
    console.print(f"\n[bold]Top 5 Candidates for Cluster {focus}:[/bold]")
    cluster_scores = scores[scores["cluster"] == focus].head(5)
    for _, row in cluster_scores.iterrows():
        console.print(
            f"  #{int(row['rank'])} {row['cell_type']} (evidence_score={row['combined_score']:.3f}, "
            f"overlap={row['pct_overlap']:.0%}, neg_conflict={row['neg_conflict']:.0%})"
        )


# ──────────────────────────────────────────────
# markers command (list available markers)
# ──────────────────────────────────────────────
@app.command()
def markers(
    tissue: str | None = typer.Option(None, "--tissue", "-t", help="Tissue to list markers for"),
    species: str = typer.Option("human", "--species", "-s", help="Species: human/mouse"),
    pack: list[str] = typer.Option(
        [],
        "--pack",
        help="Include cell types from this installed extension pack (repeatable)",
    ),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """List available cell types and markers in the knowledge graph."""
    from .data_adapter import get_all_markers_for_tissue, load_marker_atlas

    try:
        atlas = load_marker_atlas(species)
    except ValueError as exc:
        console.print(f"[red]Error: {exc}[/red]")
        raise typer.Exit(1) from exc
    if pack:
        from .pack_manager import (
            PackError,
            merge_marker_atlas,
            resolve_extension_packs,
        )

        try:
            records, warnings = resolve_extension_packs(list(pack), species)
            atlas, merge_warnings = merge_marker_atlas(atlas, records, species)
        except PackError as exc:
            console.print(f"[red]Extension pack error: {exc}[/red]")
            raise typer.Exit(1) from exc
        for warning in warnings + merge_warnings:
            console.print(f"[yellow]{warning}[/yellow]")
    available_tissues = list(atlas.get("tissues", {}).keys())

    if tissue is None:
        console.print("[bold]Available tissues:[/bold]")
        for t in available_tissues:
            n_types = len(atlas["tissues"][t].get("cell_types", {}))
            console.print(f"  {t}: {n_types} cell types")
        return

    markers_dict = get_all_markers_for_tissue(atlas, tissue)
    if not markers_dict:
        console.print(f"[yellow]No markers found for tissue '{tissue}'.[/yellow]")
        console.print(f"Available: {available_tissues}")
        return

    if json_output:
        print(json.dumps(markers_dict, indent=2))
    else:
        console.print(f"\n[bold]Cell types and markers for '{tissue}':[/bold]\n")
        for ct_name, info in markers_dict.items():
            pos = info.get("positive_markers", [])
            neg = info.get("negative_markers", [])
            cl_id = info.get("cl_id", "")
            console.print(f"  [cyan]{ct_name}[/cyan] ({cl_id})")
            console.print(f"    + markers: {', '.join(pos[:8])}{'...' if len(pos) > 8 else ''}")
            if neg:
                console.print(f"    - markers: {', '.join(neg[:5])}{'...' if len(neg) > 5 else ''}")
            console.print()


# ──────────────────────────────────────────────
# literature command (MCP-backed literature search)
# ──────────────────────────────────────────────
@app.command("atlas-governance")
def atlas_governance(
    output: str | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Optional JSON path for the atlas governance report",
    ),
    no_packs: bool = typer.Option(
        False,
        "--no-packs",
        help="Skip first-party/user extension pack inventory",
    ),
    diff_previous: str | None = typer.Option(
        None,
        "--diff-previous",
        help="Path to previous atlas JSON to compute structural diff",
    ),
    check_conflicts: bool = typer.Option(
        False,
        "--check-conflicts",
        help="Check for intra-atlas marker conflicts",
    ),
    json_output: bool = typer.Option(False, "--json", help="Print full report as JSON"),
):
    """Build an offline atlas governance report for release and Agent hosts."""
    from .atlas_governance import build_atlas_governance_report, write_atlas_governance_report

    report = build_atlas_governance_report(include_packs=not no_packs)

    if check_conflicts:
        from .atlas_conflict import detect_marker_conflicts
        from .constants import ATLAS_PATH
        from .data_adapter import load_marker_atlas

        atlas = load_marker_atlas(ATLAS_PATH)
        conflicts = detect_marker_conflicts(atlas)
        report["conflicts"] = [
            {
                "type": c.conflict_type,
                "gene": c.gene,
                "type_a": c.cell_type_a,
                "type_b": c.cell_type_b,
                "severity": c.severity,
            }
            for c in conflicts
        ]

    if diff_previous:
        from .atlas_diff import diff_atlases, format_diff_json
        from .data_adapter import load_marker_atlas

        prev_atlas = load_marker_atlas(diff_previous)
        from .constants import ATLAS_PATH

        curr_atlas = load_marker_atlas(ATLAS_PATH)
        diff = diff_atlases(prev_atlas, curr_atlas)
        report["diff_from_previous"] = format_diff_json(diff)

    if output:
        path = write_atlas_governance_report(output, include_packs=not no_packs)
        console.print(f"[green]Atlas governance report written:[/green] {path}")
    if json_output:
        print(json.dumps(report, indent=2, ensure_ascii=False))
        return
    aggregate = report["aggregate"]
    console.print("[bold]Atlas governance[/bold]")
    console.print(f"  Assets: {aggregate['n_assets']}")
    console.print(f"  Marker relationships: {aggregate['n_marker_relationships']}")
    console.print(f"  Needs edge curation: {aggregate['needs_edge_curation']}")
    if "governance_health_score" in report:
        console.print(f"  Governance health score: {report['governance_health_score']:.2f}")
    console.print(
        "  Supported annotation species: " + ", ".join(report["supported_annotation_species"])
    )
    ontology = report.get("ontology", {})
    cache = ontology.get("cache", {})
    console.print(f"  Ontology cache: {'present' if cache.get('cached') else 'missing'}")
    if "conflicts" in report:
        console.print(f"  Detected conflicts: {len(report['conflicts'])}")


@app.command("verify-novelty")
def verify_novelty(
    input_path: str = typer.Option(..., "--input", "-i", help="Path to input .h5ad file"),
    cluster_key: str = typer.Option(..., "--cluster-key", "-k", help="Cluster column in obs"),
    focus_cluster: str = typer.Option(..., "--focus", "-f", help="Cluster ID to evaluate"),
    tissue: str = typer.Option("general", "--tissue", "-t", help="Tissue context"),
    output_dir: str | None = typer.Option(None, "--output", "-o", help="Optional output directory"),
    json_output: bool = typer.Option(False, "--json", help="Print verification packet as JSON"),
):
    """Run 5-gate audit protocol on a focused OOD/novel cell-type candidate."""
    from .constants import ATLAS_PATH
    from .critic import run_critic
    from .data_adapter import get_all_markers_for_tissue, load_h5ad, load_marker_atlas
    from .marker_scorer import compute_marker_scores, generate_annotation_summary
    from .novelty_verification import verify_novelty_candidate

    adata = load_h5ad(input_path)
    atlas = load_marker_atlas(ATLAS_PATH)
    markers = get_all_markers_for_tissue(atlas, tissue)
    scores = compute_marker_scores(adata, cluster_key, markers)
    summary = generate_annotation_summary(scores, cluster_key)
    critic_df = run_critic(adata, cluster_key, summary, atlas, tissue)

    focus_row = critic_df[critic_df["cluster"].astype(str) == str(focus_cluster)]
    row_dict = focus_row.iloc[0].to_dict() if not focus_row.empty else {}

    packet = verify_novelty_candidate(adata, cluster_key, focus_cluster, row_dict, atlas, tissue)

    if output_dir:
        out_p = Path(output_dir)
        out_p.mkdir(parents=True, exist_ok=True)
        (out_p / f"novelty_verification_cluster_{focus_cluster}.json").write_text(
            json.dumps(packet, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        console.print(f"[green]Verification packet saved to {out_p}[/green]")

    if json_output:
        print(json.dumps(packet, indent=2, ensure_ascii=False))
        return

    console.print(f"[bold]Novelty Verification Packet for Cluster {focus_cluster}[/bold]")
    console.print(f"  Verification Passed: {packet['verification_passed']}")
    console.print(f"  Suggested Classification: {packet['suggested_classification']}")
    console.print(f"  Adjudication Status: {packet['adjudication_status']}")


@app.command("adjudicate-novelty")
def adjudicate_novelty(
    output_dir: str = typer.Option(
        ..., "--output", "-o", help="Output directory containing run artifacts"
    ),
    cluster: str = typer.Option(..., "--cluster", "-c", help="Cluster ID being adjudicated"),
    verdict: str = typer.Option(
        ...,
        "--verdict",
        "-v",
        help="Verdict: validated_novel_cell_type, novel_cell_state, atlas_gap_resolved, rejected_technical_artifact, rejected_mixed_cluster",
    ),
    reviewer: str = typer.Option(..., "--reviewer", "-r", help="Name/ID of expert reviewer"),
    notes: str | None = typer.Option(None, "--notes", help="Optional adjudication notes"),
    pmid: str | None = typer.Option(None, "--pmid", help="Optional PMID or evidence link"),
):
    """Log a signed human expert adjudication verdict for a novelty candidate."""
    from .novelty_verification import log_novelty_adjudication

    try:
        log_novelty_adjudication(output_dir, cluster, verdict, reviewer, notes=notes, pmid=pmid)
        console.print(f"[green]Adjudication logged for cluster {cluster}: {verdict}[/green]")
    except ValueError as exc:
        console.print(f"[red]Adjudication error: {exc}[/red]")
        raise typer.Exit(1) from exc


@app.command()
def literature(
    cell_type: str = typer.Option(..., "--cell-type", "-c", help="Cell type to search"),
    markers: str | None = typer.Option(
        None, "--markers", "-m", help="Comma-separated marker genes"
    ),
    max_refs: int = typer.Option(5, "--max-refs", help="Max references per query"),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """Search literature for marker validation (PubMed/bioRxiv)."""
    from .literature import (
        check_mcp_availability,
        generate_mcp_search_queries,
        validate_annotation_with_literature,
    )

    # Check MCP availability
    mcp_status = check_mcp_availability()
    if not mcp_status.get("pubmed_direct"):
        console.print(
            "[yellow]Warning: PubMed direct access not available. Check network.[/yellow]"
        )

    marker_list = [m.strip() for m in markers.split(",")] if markers else []

    if marker_list:
        # Validate specific markers for a cell type
        results = validate_annotation_with_literature(
            cell_type, marker_list, max_refs_per_marker=max_refs
        )

        if json_output:
            console.print(json.dumps(results, indent=2))
        else:
            console.print(f"\n[bold]Literature Validation for '{cell_type}':[/bold]\n")
            console.print(f"  Positive markers checked: {results['positive_markers_checked']}")
            console.print(f"  Markers supported: {results['positive_markers_supported']}")
            console.print(f"  Total refs found: {results['total_literature_refs']}")
            console.print(f"  Assessment: {results['overall_assessment']}\n")

            for ev in results.get("positive_evidence", []):
                status = (
                    "[green]OK[/green]" if ev["consensus"] == "supported" else "[yellow]?[/yellow]"
                )
                console.print(f"  {status} {ev['gene']}: {ev['total_refs']} refs")
                for hit in ev.get("top_hits", []):
                    console.print(
                        f"      - {hit['authors']} ({hit['year']}). {hit['title'][:60]}..."
                    )
    else:
        # Generate search queries for manual MCP use
        queries = generate_mcp_search_queries(cell_type, [])
        if json_output:
            console.print(json.dumps({"queries": queries, "mcp_status": mcp_status}, indent=2))
        else:
            console.print(f"\n[bold]Suggested search queries for '{cell_type}':[/bold]\n")
            for i, q in enumerate(queries, 1):
                console.print(f"  {i}. {q}")
            console.print("\n[bold]MCP Status:[/bold]")
            for tool, available in mcp_status.items():
                status = "[green]available[/green]" if available else "[red]not available[/red]"
                console.print(f"  {tool}: {status}")


# ──────────────────────────────────────────────
# inspect-web command (Web Inspector)
# ──────────────────────────────────────────────
@app.command()
def inspect_web(
    output_dir: Path = typer.Option(
        "./ctp_output",
        "--output",
        "-o",
        help="Path to CellTypePilot output directory (annotation review)",
    ),
    run_dir: Path | None = typer.Option(
        None,
        "--run-dir",
        help="Optional benchmark-run directory with checkpoints/ (read-only observability)",
    ),
    host: str = typer.Option("127.0.0.1", "--host", help="Host to bind to"),
    port: int = typer.Option(8765, "--port", "-p", help="Port to listen on"),
):
    """Launch the Web Inspector — interactive annotation review panel."""
    from .web_inspector import run_inspector

    if not output_dir.exists() and (run_dir is None or not run_dir.exists()):
        console.print(f"[red]Output directory not found: {output_dir}[/red]")
        console.print(
            "Run 'celltypepilot annotate' first, or pass --run-dir for observability-only."
        )
        raise typer.Exit(1)

    if run_dir is not None and not run_dir.exists():
        console.print(f"[red]Run directory not found: {run_dir}[/red]")
        raise typer.Exit(1)

    # Observability-only: allow pointing -o at a benchmark run that has checkpoints.
    effective_output = output_dir if output_dir.exists() else run_dir
    effective_run = run_dir if run_dir is not None else effective_output

    console.print("[bold]Launching Web Inspector...[/bold]")
    console.print(f"  Output dir: {effective_output}")
    console.print(f"  Observability run dir (read-only): {effective_run}")
    console.print(f"  URL: http://{host}:{port}")
    console.print("  Observability never mutates predictions; overrides use audit log.")
    console.print()

    try:
        run_inspector(effective_output, host=host, port=port, run_dir=effective_run)
    except KeyboardInterrupt:
        console.print("\n[yellow]Web Inspector stopped.[/yellow]")
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1) from e


@app.command("observe")
def observe_run(
    output_dir: Path = typer.Option(
        ...,
        "--output",
        "-o",
        help="Benchmark-run or annotation output directory (read-only)",
    ),
    json_output: bool = typer.Option(False, "--json", help="Machine-readable snapshot"),
    no_host: bool = typer.Option(False, "--no-host", help="Skip CPU/GPU probe"),
    no_hashes: bool = typer.Option(False, "--no-hashes", help="Skip product hashing"),
):
    """Read-only run observability: checkpoints, ETA, host, failures, hashes, stale.

    Never writes fold workspaces or prediction tables. Manual overrides remain on
    the append-only audit log + apply-overrides path.
    """
    from .agent_lifecycle import build_agent_status_report
    from .run_observability import ObservabilityError, build_observability_snapshot

    try:
        snapshot = build_observability_snapshot(
            output_dir,
            include_host=not no_host,
            include_product_hashes=not no_hashes,
        )
    except ObservabilityError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(1) from exc

    checkpoint_dir = Path(snapshot["run_root"]) / "checkpoints"
    snapshot["agent_lifecycle"] = build_agent_status_report(
        checkpoint_dir=checkpoint_dir if checkpoint_dir.is_dir() else None,
        release_manifest=None,
    )

    if json_output:
        print(json.dumps(snapshot, indent=2))
        return

    cp = snapshot.get("checkpoints") or {}
    eta = snapshot.get("fold_eta") or {}
    stale = snapshot.get("stale") or {}
    host = snapshot.get("host") or {}
    console.print(f"[bold]Run observability[/bold] (read-only) — {snapshot.get('run_root')}")
    console.print(
        f"  checkpoints: {cp.get('n_status_files', 0)}  "
        f"completed={eta.get('n_completed', 0)} running={eta.get('n_running', 0)} "
        f"failed={eta.get('n_failed', 0)}"
    )
    rem = eta.get("estimated_remaining_seconds")
    console.print(f"  ETA remaining (s): {rem if rem is not None else '—'}")
    console.print(
        f"  stale: {stale.get('derived_artifacts_stale')}  "
        f"({stale.get('review_state')}) — {stale.get('stale_reason')}"
    )
    cpu = host.get("cpu") or {}
    gpu = host.get("gpu") or {}
    console.print(
        f"  CPU%: {cpu.get('cpu_percent')}  cores={cpu.get('cpu_count_logical')}  "
        f"GPU: {'yes' if gpu.get('available') else 'no'}"
    )
    failures = snapshot.get("failures") or []
    if failures:
        console.print(f"  [red]failures: {len(failures)}[/red]")
        for row in failures[:10]:
            console.print(
                f"    - {row.get('method')} / {row.get('fold_id')}: {row.get('failure_reason')}"
            )
    console.print(
        "  prediction mutation: [green]denied[/green] "
        "(overrides → audit log + apply-overrides only)"
    )
    life = snapshot.get("agent_lifecycle") or {}
    ck = life.get("checkpoints") or {}
    if ck:
        console.print(
            f"  agent rollup: [bold]{ck.get('rollup_agent_state')}[/bold]  "
            f"counts={ck.get('counts_by_agent_state')}"
        )


@app.command("qc-diagnostics")
def qc_diagnostics_cmd(
    input: Path = typer.Option(..., "--input", "-i", help="Path to .h5ad file"),
    cluster_key: str | None = typer.Option(
        None, "--cluster-key", "-k", help="Cluster key for sample-enrichment axis"
    ),
    output_dir: Path = typer.Option(
        "qc_output", "--output", "-o", help="Directory for qc_diagnostics.json/csv"
    ),
    doublet_table: Path | None = typer.Option(
        None, "--doublet-table", help="External doublet tool CSV"
    ),
    ambient_table: Path | None = typer.Option(
        None, "--ambient-table", help="External ambient-RNA tool CSV"
    ),
    doublet_score_column: str = typer.Option("doublet_score", "--doublet-score-column"),
    ambient_score_column: str = typer.Option("ambient_score", "--ambient-score-column"),
    doublet_flag_column: str | None = typer.Option(None, "--doublet-flag-column"),
    ambient_flag_column: str | None = typer.Option(None, "--ambient-flag-column"),
    json_output: bool = typer.Option(False, "--json", help="Print full QC report JSON"),
):
    """Assemble composable QC diagnostic axes (never rewrites identity labels).

    Missing metadata yields not_assessed_*, never clean. External doublet/ambient
    tools are optional plugins on a diagnostic axis only.
    """
    from .data_adapter import load_h5ad
    from .qc_diagnostics import (
        QCDiagnosticError,
        assemble_qc_diagnostics,
        load_external_tool_table,
        write_qc_diagnostics,
    )

    try:
        adata = load_h5ad(str(input))
        doublet = (
            load_external_tool_table(
                doublet_table,
                axis="doublet",
                score_column=doublet_score_column,
                flag_column=doublet_flag_column,
            )
            if doublet_table is not None
            else None
        )
        ambient = (
            load_external_tool_table(
                ambient_table,
                axis="ambient_rna",
                score_column=ambient_score_column,
                flag_column=ambient_flag_column,
            )
            if ambient_table is not None
            else None
        )
        report = assemble_qc_diagnostics(
            adata,
            cluster_key=cluster_key,
            doublet_table=doublet,
            ambient_table=ambient,
        )
        paths = write_qc_diagnostics(report, output_dir)
    except (QCDiagnosticError, FileNotFoundError, OSError) as exc:
        console.print(f"[red]QC diagnostics failed: {exc}[/red]")
        raise typer.Exit(1) from exc

    if json_output:
        print(json.dumps(report, indent=2))
    else:
        console.print("[green]QC diagnostics written[/green] (can_rescue_identity=false)")
        console.print(f"  rollup: {report['rollup_status']} / {report['rollup_flag']}")
        for axis, payload in report["axes"].items():
            console.print(
                f"  {axis}: {payload['status']}  flag={payload['flag']}  "
                f"flagged={payload['n_cells_flagged']}"
            )
        for key, path in paths.items():
            console.print(f"  {key}: {path}")


@app.command("host-acceptance")
def host_acceptance(
    worktree: Path | None = typer.Option(
        None,
        "--worktree",
        help="Optional independent git worktree path (created if missing via script)",
    ),
    json_output: bool = typer.Option(False, "--json", help="Machine-readable report"),
    skip_worktree: bool = typer.Option(
        False,
        "--skip-worktree",
        help="Run harness in-repo under scratch/ (still isolated fixtures)",
    ),
):
    """Host acceptance: discovery + lifecycle discrimination for Agent hosts.

    Verifies Codex / Claude Code / MCP discovery surfaces and that an Agent can
    distinguish running / completed / failed / unavailable / claim_ready — not
    merely that CLI entry points exist.
    """
    from .host_acceptance import run_host_acceptance

    report = run_host_acceptance(
        worktree=worktree,
        skip_worktree=skip_worktree,
    )
    if json_output:
        print(json.dumps(report, indent=2))
    else:
        console.print(f"[bold]Host acceptance[/bold] status={report['overall_status']}")
        for check in report["checks"]:
            mark = "green" if check["passed"] else "red"
            console.print(f"  [{mark}]{check['id']}[/{mark}]: {check['detail']}")
        if report.get("agent_discrimination"):
            console.print(
                f"  discrimination: {report['agent_discrimination'].get('states_observed')}"
            )
    if report["overall_status"] != "passed":
        raise typer.Exit(1)


# ──────────────────────────────────────────────
# convert-rds command (Seurat support)
# ──────────────────────────────────────────────
@app.command()
def convert_rds(
    input_rds: Path = typer.Option(..., "--input", "-i", help="Path to Seurat .rds file"),
    output_h5ad: Path | None = typer.Option(None, "--output", "-o", help="Output .h5ad path"),
):
    """Convert Seurat .rds to .h5ad for CellTypePilot annotation."""
    from .seurat_adapter import check_seurat_support, seurat_to_h5ad

    # Check support
    support = check_seurat_support()
    if not support["seurat_rds_supported"]:
        console.print("[red]Seurat .rds support not available.[/red]")
        console.print("Install one of:")
        console.print("  pip install rpy2          # Option 1: rpy2 (requires R installed)")
        console.print("  Or install R + SeuratDisk  # Option 2: Rscript in PATH")
        raise typer.Exit(1)

    if output_h5ad is None:
        output_h5ad = input_rds.with_suffix(".h5ad")

    console.print("[bold]Converting Seurat .rds to .h5ad...[/bold]")
    console.print(f"  Input:  {input_rds}")
    console.print(f"  Output: {output_h5ad}")

    try:
        result_path = seurat_to_h5ad(input_rds, output_h5ad)
        console.print(f"[green]Conversion complete: {result_path}[/green]")
        console.print(f"\nNow run: celltypepilot annotate --input {result_path}")
    except Exception as e:
        console.print(f"[red]Conversion failed: {e}[/red]")
        raise typer.Exit(1) from e


# ──────────────────────────────────────────────
# apply-overrides command (write corrections to .h5ad)
# ──────────────────────────────────────────────
@app.command()
def apply_overrides(
    output_dir: Path = typer.Option(
        "./ctp_output",
        "--output",
        "-o",
        help="CellTypePilot output directory containing data.annotated.h5ad",
    ),
    overrides_file: Path = typer.Option(
        ...,
        "--overrides",
        "-f",
        help="Path to annotation_overrides.json",
    ),
    regenerate: bool = typer.Option(
        False,
        "--regenerate",
        "-r",
        help="Regenerate figures after applying overrides",
    ),
):
    """Apply annotation overrides from Web Inspector to .h5ad file.

    Reads the overrides JSON exported from the Web Inspector,
    writes corrected labels back to the annotated .h5ad, and
    creates a timestamped backup of the original file.

    Example:
        celltypepilot apply-overrides -o ./ctp_output -f annotation_overrides.json
    """
    from .orchestrator import (
        apply_overrides_to_h5ad,
        find_cluster_column,
        regenerate_figures_after_override,
    )

    h5ad_path = output_dir / "data.annotated.h5ad"
    if not h5ad_path.exists():
        console.print(f"[red]No annotated data found at {h5ad_path}[/red]")
        console.print("Run 'celltypepilot annotate' first.")
        raise typer.Exit(1)

    if not overrides_file.exists():
        console.print(f"[red]Overrides file not found: {overrides_file}[/red]")
        console.print("Export overrides from the Web Inspector first.")
        raise typer.Exit(1)

    # Load overrides
    try:
        overrides = json.loads(overrides_file.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, ValueError) as e:
        console.print(f"[red]Invalid overrides JSON: {e}[/red]")
        raise typer.Exit(1) from e

    if not overrides:
        console.print("[yellow]No overrides to apply.[/yellow]")
        return

    console.print(f"[bold]Applying {len(overrides)} override(s)...[/bold]")

    result = apply_overrides_to_h5ad(h5ad_path, overrides)
    console.print(f"  Backup: {result['backup']}")

    for detail in result["details"]:
        if detail["status"] == "applied":
            console.print(
                f"  [green]OK[/green] Cluster {detail['cluster']}: "
                f"{detail.get('old_type', 'Unknown')} → {detail['new_type']} "
                f"({detail['n_cells']} cells)"
            )
        else:
            console.print(
                f"  [yellow]Cluster {detail['cluster']}: "
                f"{detail.get('reason', detail['status'])}, skipped[/yellow]"
            )

    console.print(
        f"\n[bold green]Applied {result['applied']} override(s)[/bold green], "
        f"{result['skipped']} skipped"
    )

    # Regenerate figures if requested
    if regenerate:
        console.print("\n[bold blue]Regenerating figures...[/bold blue]")
        try:
            import scanpy as sc

            adata = sc.read_h5ad(h5ad_path)
            cluster_col = find_cluster_column(adata.obs)
            if cluster_col is None:
                console.print(
                    "  [yellow]No cluster column found, skipping figure regeneration[/yellow]"
                )
            else:
                figure_paths = regenerate_figures_after_override(output_dir, adata, cluster_col)
                if figure_paths:
                    console.print(f"  Regenerated {len(figure_paths)} figures")
                else:
                    console.print(
                        "  [yellow]No embedding found, skipping figure regeneration[/yellow]"
                    )
        except Exception as e:
            console.print(f"  [yellow]Figure regeneration failed: {e}[/yellow]")

    console.print(f"\nOutput: {h5ad_path.resolve()}")


# ──────────────────────────────────────────────
# annotate-embedding command (reference embedding + ensemble)
# ──────────────────────────────────────────────
@app.command()
def annotate_embedding(
    input: str = typer.Option(..., "--input", "-i", help="Path to .h5ad file"),
    cluster_key: str = typer.Option(..., "--cluster-key", "-k", help="Cluster key in obs"),
    reference: str | None = typer.Option(
        None, "--reference", "-r", help="Reference .h5ad with cell type labels"
    ),
    ref_label_key: str = typer.Option(
        "cell_type", "--ref-label", help="Cell type column in reference.obs"
    ),
    model_path: str | None = typer.Option(
        None, "--model", "-m", help="CellTypist model path (.pkl)"
    ),
    backend: str = typer.Option(
        "auto", "--backend", "-b", help="Backend: auto/celltypist/scanvi/knn/correlation"
    ),
    output_dir: str = typer.Option(".", "--output", "-o", help="Output directory"),
    species: str | None = typer.Option(None, "--species", "-s", help="Species: human/mouse"),
    tissue: str | None = typer.Option(None, "--tissue", "-t", help="Tissue context"),
    marker_weight: float = typer.Option(0.5, "--marker-weight", help="Base marker weight (0-1)"),
    no_ensemble: bool = typer.Option(
        False, "--no-ensemble", help="Skip ensemble, use reference only"
    ),
    allow_unverified_reference: bool = typer.Option(
        False,
        "--allow-unverified-reference",
        help="Explicitly override the reference provenance safety gate; recorded in manifest",
    ),
    marker_evidence_policy: str = typer.Option(
        "database",
        "--marker-evidence-policy",
        help="database, literature, edge_verified, or primary; stricter policies exclude unverified edges",
    ),
    calibration_policy: str | None = typer.Option(
        None,
        "--calibration-policy",
        help="Optional abstention policy JSON fitted on a separate calibration dataset",
    ),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """Compatibility wrapper for one in-process reference candidate backend.

    CellTypist/scANVI/KNN/correlation generates candidates; marker evidence and
    the critic may downgrade them. The default selector still requires an
    additional independent backend before publishing an accepted identity.

    Examples:
        # Auto-detect backend with reference atlas
        celltypepilot annotate-embedding -i data.h5ad -k leiden -r reference.h5ad

        # Use CellTypist pre-trained model
        celltypepilot annotate-embedding -i data.h5ad -k leiden -m models/Immune_All_Low.pkl

        # Force correlation backend (no extra deps)
        celltypepilot annotate-embedding -i data.h5ad -k leiden -r ref.h5ad -b correlation
    """
    # Compatibility wrapper: all scoring, critic review, writeback, report, and
    # provenance now run through the single annotation orchestrator.
    from .orchestrator import PipelineError, run_annotation_pipeline

    try:
        result = run_annotation_pipeline(
            input_path=input,
            cluster_key=cluster_key,
            output_dir=output_dir,
            species=species,
            tissue=tissue,
            reference_path=reference,
            ref_label_key=ref_label_key,
            model_path=model_path,
            reference_backend=backend,
            marker_weight=marker_weight,
            no_ensemble=no_ensemble,
            allow_unverified_reference=allow_unverified_reference,
            marker_evidence_policy=marker_evidence_policy,
            calibration_policy_path=calibration_policy,
        )
    except PipelineError as exc:
        console.print(f"[red]Error: {exc}[/red]")
        raise typer.Exit(1) from exc

    if json_output:
        console.print(
            json.dumps(
                {
                    "annotations": result["critic_results"].to_dict(orient="records"),
                    "critic_summary": result["critic_summary"],
                    "manifest": result["manifest"],
                },
                indent=2,
                default=str,
            )
        )
    console.print("[bold green]Done![/bold green] Reference candidates reviewed selectively.")
    console.print(f"Output directory: {result['output_path'].resolve()}")
    return

    from .data_adapter import (
        detect_species,
        detect_tissue,
        get_all_markers_for_tissue,
        load_h5ad,
        load_marker_atlas,
        require_supported_annotation_species,
    )
    from .ensemble_scorer import analyze_disagreements, ensemble_scores, generate_ensemble_summary
    from .marker_scorer import compute_marker_scores, generate_annotation_summary
    from .reference_scorer import (
        check_reference_backends,
        detect_transitional_states,
        score_by_reference,
    )

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Step 1: Load query
    console.print("[bold blue]Step 1/5:[/bold blue] Loading data...")
    adata = load_h5ad(input)

    if species is None:
        species = detect_species(adata)
    try:
        require_supported_annotation_species(species)
    except ValueError as exc:
        console.print(f"[red]Error: {exc}[/red]")
        raise typer.Exit(1) from exc
    if tissue is None:
        tissue = detect_tissue(adata) or "general"

    # Step 2: Marker scoring (existing pipeline)
    console.print("[bold blue]Step 2/5:[/bold blue] Marker-based scoring...")
    atlas = load_marker_atlas(species)
    markers = get_all_markers_for_tissue(atlas, tissue)
    marker_scores = compute_marker_scores(adata, cluster_key, markers)
    marker_summary = generate_annotation_summary(marker_scores, cluster_key)
    console.print(f"  Marker: {len(marker_summary)} clusters scored")

    # Step 3: Reference embedding scoring
    console.print("[bold blue]Step 3/5:[/bold blue] Reference embedding scoring...")

    # Check available backends
    backends = check_reference_backends()
    console.print(f"  Available backends: {', '.join(k for k, v in backends.items() if v)}")

    # Load reference if provided
    ref_adata = None
    if reference:
        console.print(f"  Loading reference: {reference}")
        ref_adata = load_h5ad(reference)

    ref_scores = score_by_reference(
        adata,
        cluster_key,
        reference=ref_adata,
        ref_label_key=ref_label_key,
        model_path=model_path,
        backend=backend,
    )
    console.print(f"  Reference: {len(ref_scores['cluster'].unique())} clusters scored")

    # Step 4: Ensemble fusion
    if no_ensemble:
        console.print("[bold blue]Step 4/5:[/bold blue] Using reference scores only (no ensemble)")
        final_df = ref_scores
    else:
        console.print("[bold blue]Step 4/5:[/bold blue] Ensemble fusion...")
        final_df = ensemble_scores(
            marker_scores,
            ref_scores,
            marker_weight=marker_weight,
            adaptive=True,
        )
        ensemble_summary = generate_ensemble_summary(final_df)

        # Show agreement stats
        n_agree = sum(1 for _, r in ensemble_summary.iterrows() if r.get("agreement", True))
        n_total = len(ensemble_summary)
        console.print(f"  Agreement: {n_agree}/{n_total} clusters agree between methods")

    # Step 5: Detect transitional states
    console.print("[bold blue]Step 5/5:[/bold blue] Detecting transitional states...")
    transitions = detect_transitional_states(ref_scores, marker_scores)
    disagreements = analyze_disagreements(final_df)

    n_transitional = transitions["is_transitional"].sum() if not transitions.empty else 0
    console.print(f"  Transitional candidates: {n_transitional}")

    if not disagreements.empty:
        console.print("\n  [bold]Disagreements (potential novel/transitional states):[/bold]")
        for _, row in disagreements.head(5).iterrows():
            console.print(
                f"    Cluster {row['cluster']}: "
                f"marker→{row['marker_type']} ({row['marker_score']:.2f}) "
                f"vs ref→{row['ref_type']} ({row['ref_score']:.2f})"
            )
            console.print(f"      → {row['interpretation'][:80]}...")

    # Save outputs
    ensemble_path = output_path / "ensemble_scores.csv"
    final_df.to_csv(ensemble_path, index=False)
    console.print(f"\n  Ensemble scores: {ensemble_path}")

    if not transitions.empty:
        trans_path = output_path / "transitional_states.csv"
        transitions.to_csv(trans_path, index=False)
        console.print(f"  Transitional states: {trans_path}")

    if not disagreements.empty:
        disagree_path = output_path / "disagreements.csv"
        disagreements.to_csv(disagree_path, index=False)
        console.print(f"  Disagreements: {disagree_path}")

    if json_output:
        output_json = {
            "ensemble": final_df.to_dict(orient="records"),
            "transitional": transitions.to_dict(orient="records") if not transitions.empty else [],
            "disagreements": disagreements.to_dict(orient="records")
            if not disagreements.empty
            else [],
        }
        console.print(json.dumps(output_json, indent=2, default=str))

    console.print("\n[bold green]Done![/bold green] Ensemble annotation complete.")


# ──────────────────────────────────────────────
# backends command (check reference scoring support)
# ──────────────────────────────────────────────
@app.command()
def backends():
    """Show available reference scoring backends."""
    from .native_backends import check_native_backend_runtimes
    from .reference_scorer import check_reference_backends

    status = check_reference_backends()

    console.print("[bold]Reference Scoring Backends:[/bold]\n")
    for name, available in status.items():
        icon = "[green]✓[/green]" if available else "[red]✗[/red]"
        console.print(f"  {icon} {name}")

    console.print("\n[bold]Install backends:[/bold]")
    console.print("  pip install celltypist       # Pre-trained models (recommended)")
    console.print("  pip install scvi-tools        # Custom reference atlas (scANVI)")
    console.print("  Or provide --reference .h5ad  # KNN/correlation (no extra deps)")
    console.print("\n[bold]Native annotate runners (static preflight):[/bold]")
    for name, details in check_native_backend_runtimes().items():
        icon = "[green]available[/green]" if details["available"] else "[yellow]missing[/yellow]"
        console.print(f"  {name}: {icon} ({details['check']})")
    console.print("  Static availability is not a successful biological or end-to-end run.")


# ──────────────────────────────────────────────
# license command
# ──────────────────────────────────────────────
@app.command("domain-validation-plan")
def domain_validation_plan(
    registry: Path = typer.Option(
        Path("benchmarks/public_v1/registry.json"),
        "--registry",
        help="Frozen public cohort registry JSON",
    ),
    output_dir: Path = typer.Option(
        Path("benchmarks/domain_depth_v1"), "--output", "-o", help="Plan/output root"
    ),
):
    """Inventory the three depth domains and lock an executable evidence plan."""
    from .domain_validation_pipeline import (
        DomainValidationPipelineError,
        build_domain_validation_plan,
    )

    try:
        result = build_domain_validation_plan(registry, output_dir)
    except DomainValidationPipelineError as exc:
        console.print(f"[red]Domain plan failed: {exc}[/red]")
        raise typer.Exit(1) from exc
    console.print(f"Domain validation plan: {result['plan_path']}")
    for domain_id, summary in result["plan"]["domains"].items():
        console.print(
            f"  {domain_id}: cohorts={summary['registered_cohorts']} "
            f"ready={summary['ready_cohorts']} claim_ready=false"
        )


@app.command("lineage-coverage-audit")
def lineage_coverage_audit(
    predictions: Path = typer.Option(..., "--predictions", help="OOF predictions CSV"),
    cluster_map: Path = typer.Option(..., "--cluster-map", help="Truth-free cluster map CSV"),
    domain: str = typer.Option(..., "--domain", help="lung, gut_ibd, or tumor_microenvironment"),
    output_dir: Path = typer.Option(
        Path("lineage_coverage_audit"), "--output", "-o", help="Audit artifact directory"
    ),
    methods: str = typer.Option(
        "celltypist,popv,singler", "--methods", help="Comma-separated candidate backends"
    ),
):
    """Audit multi-lineage addressability without reading benchmark truth."""
    from .lineage_coverage import LineageCoverageError, build_selector_lineage_audit

    selected_methods = tuple(value.strip() for value in methods.split(",") if value.strip())
    try:
        result = build_selector_lineage_audit(
            predictions,
            cluster_map,
            output_dir,
            domain_id=domain,
            methods=selected_methods,
        )
    except LineageCoverageError as exc:
        console.print(f"[red]Lineage coverage audit failed: {exc}[/red]")
        raise typer.Exit(1) from exc
    manifest = result["manifest"]
    console.print(f"Lineage coverage audit: {result['manifest_path']}")
    console.print(
        f"  status={manifest['status']} observed={','.join(manifest['observed_lineages'])} "
        f"accepted_clusters={manifest['n_accepted_clusters']}"
    )


@app.command("calibration-split")
def calibration_split(
    registry: Path = typer.Option(
        Path("benchmarks/public_v1/registry.json"),
        "--registry",
        help="Frozen cohort registry JSON",
    ),
    output_dir: Path = typer.Option(
        Path("benchmarks/calibration_v1"), "--output", "-o", help="Locked role directory"
    ),
    fraction: float = typer.Option(0.2, "--fraction", help="Calibration donor fraction"),
    seed: str = typer.Option(
        "celltypepilot-calibration-v1", "--seed", help="Predeclared deterministic seed"
    ),
):
    """Lock donor-disjoint calibration and evaluation roles without reading truth."""
    from .calibration_split import CalibrationSplitError, build_calibration_split

    try:
        result = build_calibration_split(
            registry, output_dir, calibration_fraction=fraction, seed=seed
        )
    except (CalibrationSplitError, OSError, json.JSONDecodeError) as exc:
        console.print(f"[red]Calibration split failed: {exc}[/red]")
        raise typer.Exit(1) from exc
    manifest = result["manifest"]
    console.print(f"Calibration split: {result['manifest_path']}")
    console.print(
        f"  calibration_donors={manifest['n_calibration_donors']} "
        f"evaluation_donors={manifest['n_evaluation_donors']}"
    )


@app.command("governance-freeze-verify")
def governance_freeze_verify(
    freeze: Path = typer.Option(
        Path(__file__).parent / "data" / "governance_freeze_v1.json",
        "--freeze",
        help="Frozen governance manifest",
    ),
):
    """Verify that release-governed code and data still match the freeze."""
    from .governance_freeze import GovernanceFreezeError, verify_governance_freeze

    try:
        result = verify_governance_freeze(freeze)
    except (GovernanceFreezeError, OSError, json.JSONDecodeError) as exc:
        console.print(f"[red]Governance freeze verification failed: {exc}[/red]")
        raise typer.Exit(1) from exc
    console.print(
        f"Governance freeze verified: {result['release_id']} "
        f"({result['n_files']} files, {result['freeze_sha256']})"
    )


@app.command("calibrate-locked-donors")
def calibrate_locked_donors(
    registry: Path = typer.Option(..., "--registry", help="Frozen cohort registry"),
    assignments: Path = typer.Option(..., "--assignments", help="Locked donor roles CSV"),
    cohort: str = typer.Option(..., "--cohort", help="Calibration cohort ID"),
    predictions: Path = typer.Option(..., "--predictions", help="Selector cell predictions"),
    label_map: Path = typer.Option(..., "--label-map", help="Predeclared label map"),
    output: Path = typer.Option(..., "--output", "-o", help="Abstention policy JSON"),
    max_error: float = typer.Option(0.25, "--max-error"),
    min_coverage: float = typer.Option(0.2, "--min-coverage"),
):
    """Fit a downgrade-only policy on donors locked to the calibration role."""
    from .calibration import CalibrationError
    from .calibration_split import CalibrationSplitError, fit_policy_from_locked_donors

    try:
        result = fit_policy_from_locked_donors(
            registry,
            assignments,
            cohort,
            predictions,
            label_map,
            output,
            max_selective_error=max_error,
            min_coverage=min_coverage,
        )
    except (CalibrationSplitError, CalibrationError, OSError, json.JSONDecodeError) as exc:
        console.print(f"[red]Locked-donor calibration failed: {exc}[/red]")
        raise typer.Exit(1) from exc
    policy = result["policy"]
    console.print(f"Locked-donor calibration policy: {result['policy_path']}")
    console.print(
        f"  threshold={policy['threshold']:.6g} donors={policy['n_calibration_donors']} "
        f"cells={policy['n_calibration_cells']}"
    )


@app.command("domain-validation-run")
def domain_validation_run(
    plan: Path = typer.Option(..., "--plan", help="Locked domain_validation_plan.json"),
    domain: list[str] = typer.Option([], "--domain", help="Optional domain filter; repeatable"),
    cohort: list[str] = typer.Option([], "--cohort", help="Optional cohort_id filter; repeatable"),
):
    """Execute ready depth-domain cohorts with resumable method/fold checkpoints."""
    from .domain_validation_pipeline import (
        DomainValidationPipelineError,
        execute_domain_validation_plan,
    )

    try:
        result = execute_domain_validation_plan(
            plan,
            domain_ids=set(domain) or None,
            cohort_ids=set(cohort) or None,
        )
    except DomainValidationPipelineError as exc:
        console.print(f"[red]Domain execution failed: {exc}[/red]")
        raise typer.Exit(1) from exc
    console.print(f"Domain run status: {result['status_path']}")
    console.print(f"Domain run manifest: {result['manifest_path']}")


@app.command()
def benchmark(
    input: str = typer.Option(..., "--input", "-i", help="Path to benchmark .h5ad"),
    truth_key: str = typer.Option(..., "--truth-key", help="Ground-truth label column in obs"),
    study_key: str = typer.Option(..., "--study-key", help="Study identifier column in obs"),
    donor_key: str = typer.Option(..., "--donor-key", help="Donor identifier column in obs"),
    output_dir: str = typer.Option("benchmark", "--output", "-o", help="Output directory"),
    strategy: str = typer.Option("study", "--strategy", help="study or donor holdout"),
    predictions: str | None = typer.Option(
        None,
        "--predictions",
        help="Out-of-fold long CSV: cell_id,fold_id,method,predicted_label",
    ),
    cluster_key: str | None = typer.Option(
        None, "--cluster-key", help="Required for cluster or both evaluation units"
    ),
    evaluation_unit: str = typer.Option("cell", "--evaluation-unit", help="cell, cluster, or both"),
):
    """Lock study/donor holdouts and optionally evaluate comparator predictions."""
    import pandas as pd

    from .benchmark import (
        build_calibration_artifacts,
        build_cluster_level_track,
        build_holdout_assignments,
        evaluate_holdout_predictions,
        save_benchmark_plan,
    )
    from .data_adapter import load_h5ad

    adata = load_h5ad(input)
    if truth_key not in adata.obs:
        raise typer.BadParameter(f"truth key '{truth_key}' not found in obs")
    assignments = build_holdout_assignments(adata.obs, study_key, donor_key, strategy)
    paths = save_benchmark_plan(assignments, output_dir, study_key, donor_key, strategy)
    console.print(f"Locked {assignments['fold_id'].nunique()} independent test folds")

    if predictions is not None:
        if evaluation_unit not in {"cell", "cluster", "both"}:
            raise typer.BadParameter("--evaluation-unit must be cell, cluster, or both")
        if evaluation_unit in {"cluster", "both"} and (
            cluster_key is None or cluster_key not in adata.obs
        ):
            raise typer.BadParameter("cluster/both evaluation requires a valid --cluster-key")
        prediction_table = pd.read_csv(predictions, dtype={"cell_id": str, "fold_id": str})
        if evaluation_unit in {"cell", "both"}:
            aggregate, per_fold = evaluate_holdout_predictions(
                adata.obs[truth_key], assignments, prediction_table
            )
            aggregate_path = Path(output_dir) / "benchmark_results.csv"
            per_fold_path = Path(output_dir) / "benchmark_results_by_fold.csv"
            aggregate.to_csv(aggregate_path, index=False)
            per_fold.to_csv(per_fold_path, index=False)
            paths.update({"results": aggregate_path, "per_fold": per_fold_path})
            if "confidence" in prediction_table.columns:
                bins, risk = build_calibration_artifacts(adata.obs[truth_key], prediction_table)
                bins_path = Path(output_dir) / "calibration_bins.csv"
                risk_path = Path(output_dir) / "risk_coverage.csv"
                bins.to_csv(bins_path, index=False)
                risk.to_csv(risk_path, index=False)
                paths.update({"calibration_bins": bins_path, "risk_coverage": risk_path})
            console.print(aggregate.to_string(index=False))
        if evaluation_unit in {"cluster", "both"}:
            cluster_truth, cluster_assignments, cluster_predictions, diagnostics = (
                build_cluster_level_track(
                    adata.obs[truth_key],
                    assignments,
                    prediction_table,
                    adata.obs[cluster_key],
                )
            )
            cluster_results, cluster_folds = evaluate_holdout_predictions(
                cluster_truth,
                cluster_assignments,
                cluster_predictions,
            )
            for name, frame in (
                ("cluster_track_results", cluster_results),
                ("cluster_track_results_by_fold", cluster_folds),
                ("cluster_track_predictions", cluster_predictions),
                ("cluster_track_diagnostics", diagnostics),
            ):
                path = Path(output_dir) / f"{name}.csv"
                frame.to_csv(path, index=False)
                paths[name] = path
    else:
        console.print(
            "No predictions supplied; comparator status remains not_provided until "
            "out-of-fold predictions are generated."
        )

    for label, path in paths.items():
        console.print(f"  {label}: {path}")


@app.command()
def calibrate(
    input: str = typer.Option(..., "--input", "-i", help="Calibration .h5ad"),
    truth_key: str = typer.Option(..., "--truth-key", help="Ground-truth label column"),
    predictions: str = typer.Option(
        ..., "--predictions", help="Calibration CSV with cell_id,method,predicted_label,confidence"
    ),
    method: str = typer.Option("celltypepilot", "--method", help="Method to calibrate"),
    output: str = typer.Option(
        "abstention_policy.json", "--output", "-o", help="Output policy JSON"
    ),
    max_selective_error: float = typer.Option(
        0.1, "--max-selective-error", help="Maximum empirical error among retained predictions"
    ),
    min_coverage: float = typer.Option(
        0.2, "--min-coverage", help="Minimum fraction of predictions retained"
    ),
):
    """Fit an abstention threshold on a dataset explicitly reserved for calibration."""
    import pandas as pd

    from .calibration import CalibrationError, fit_abstention_policy, save_abstention_policy
    from .data_adapter import load_h5ad

    adata = load_h5ad(input)
    if truth_key not in adata.obs:
        raise typer.BadParameter(f"truth key '{truth_key}' not found in obs")
    prediction_table = pd.read_csv(predictions, dtype={"cell_id": str})
    try:
        policy = fit_abstention_policy(
            adata.obs[truth_key],
            prediction_table,
            method=method,
            max_selective_error=max_selective_error,
            min_coverage=min_coverage,
            dataset_role="calibration",
        )
    except CalibrationError as exc:
        console.print(f"[red]Calibration failed: {exc}[/red]")
        raise typer.Exit(1) from exc
    path = save_abstention_policy(policy, output)
    console.print(f"Saved abstention policy: {path}")
    console.print(
        f"  threshold={policy['threshold']:.4f}  "
        f"empirical_error={policy['empirical_selective_error']:.4f}  "
        f"coverage={policy['empirical_coverage']:.4f}"
    )


@app.command("benchmark-run")
def benchmark_run(
    input: str = typer.Option(..., "--input", "-i", help="Benchmark .h5ad"),
    truth_key: str = typer.Option(..., "--truth-key", help="Locked canonical truth column"),
    study_key: str = typer.Option(..., "--study-key", help="Study identifier column"),
    constant_study_id: str | None = typer.Option(
        None,
        "--constant-study-id",
        help="Truth-blind constant used when an immutable single-study h5ad lacks study_key",
    ),
    donor_key: str = typer.Option(..., "--donor-key", help="Globally unique donor column"),
    cluster_key: str = typer.Option(..., "--cluster-key", "-k", help="Predeclared cluster column"),
    cluster_map: str | None = typer.Option(
        None,
        "--cluster-map",
        help="Optional locked CSV with cell_id,cluster when the source h5ad is immutable",
    ),
    species: str = typer.Option(..., "--species", "-s", help="Explicit species"),
    tissue: str = typer.Option(..., "--tissue", "-t", help="Explicit tissue"),
    output_dir: str = typer.Option("benchmark", "--output", "-o", help="Output directory"),
    strategy: str = typer.Option("study", "--strategy", help="study or donor holdout"),
    methods: str = typer.Option(
        "celltypepilot,celltypist",
        "--methods",
        help="Comma-separated methods; external methods require comparator configs",
    ),
    comparator_config: list[str] = typer.Option(
        [],
        "--comparator-config",
        help="Repeatable JSON argv adapter for SingleR, Azimuth, or popV",
    ),
    label_map: str | None = typer.Option(
        None,
        "--label-map",
        help="Predeclared CSV: method,raw_label,canonical_label",
    ),
    continue_on_unavailable: bool = typer.Option(
        False,
        "--continue-on-unavailable",
        help="Record unavailable comparators instead of failing the benchmark run",
    ),
    fold_id: list[str] = typer.Option(
        [],
        "--fold-id",
        help="Optional: only run these fold_id values (repeatable; for donor-fold workers)",
    ),
    no_aggregate_tables: bool = typer.Option(
        False,
        "--no-aggregate-tables",
        help=(
            "Write only atomic per-fold checkpoints; do not rewrite global "
            "out_of_fold_predictions.csv (required for multi-node workers)"
        ),
    ),
    worker_id: str | None = typer.Option(
        None, "--worker-id", help="Optional worker identity recorded in checkpoint status"
    ),
    batch_id: str | None = typer.Option(
        None, "--batch-id", help="Optional frozen batch id recorded in checkpoint status"
    ),
    evaluation_unit: str = typer.Option(
        "both", "--evaluation-unit", help="cell, cluster, or both; both keeps endpoints separate"
    ),
):
    """Execute locked-fold comparator runs; never expose test truth to adapters."""
    import pandas as pd

    from .benchmark import (
        BenchmarkValidationError,
        apply_truth_label_map,
        build_calibration_artifacts,
        build_cluster_level_track,
        build_holdout_assignments,
        evaluate_holdout_predictions,
        save_benchmark_plan,
    )
    from .benchmark_runner import (
        CommandComparator,
        configure_benchmark_runtime,
        run_benchmark_comparators,
    )
    from .data_adapter import compute_data_hash, load_h5ad

    runtime_config = configure_benchmark_runtime(output_dir)
    adata = load_h5ad(input)
    if study_key not in adata.obs:
        if not constant_study_id:
            raise typer.BadParameter(
                f"study key '{study_key}' not found; provide --constant-study-id"
            )
        adata.obs[study_key] = str(constant_study_id)
    if cluster_key not in adata.obs:
        if not cluster_map:
            raise typer.BadParameter(
                f"cluster key '{cluster_key}' not found; provide --cluster-map for immutable data"
            )
        cluster_frame = pd.read_csv(cluster_map, dtype=str)
        required_cluster_columns = {"cell_id", "cluster"}
        missing_cluster_columns = required_cluster_columns - set(cluster_frame)
        if missing_cluster_columns:
            raise typer.BadParameter(
                f"cluster map missing columns: {sorted(missing_cluster_columns)}"
            )
        if cluster_frame["cell_id"].duplicated().any():
            raise typer.BadParameter("cluster map contains duplicate cell_id values")
        cluster_values = cluster_frame.set_index("cell_id")["cluster"]
        cluster_values.index = cluster_values.index.astype(str)
        expected_cells = pd.Index(adata.obs_names.astype(str))
        missing_cells = expected_cells.difference(cluster_values.index)
        extra_cells = cluster_values.index.difference(expected_cells)
        if len(missing_cells) or len(extra_cells):
            raise typer.BadParameter(
                f"cluster map cell mismatch: missing={len(missing_cells)}, extra={len(extra_cells)}"
            )
        adata.obs[cluster_key] = cluster_values.reindex(expected_cells).to_numpy()
    requested = tuple(value.strip().lower() for value in methods.split(",") if value.strip())
    if not requested:
        raise typer.BadParameter("At least one method is required")
    if evaluation_unit not in {"cell", "cluster", "both"}:
        raise typer.BadParameter("--evaluation-unit must be cell, cluster, or both")
    assignments = build_holdout_assignments(adata.obs, study_key, donor_key, strategy)
    paths = save_benchmark_plan(assignments, output_dir, study_key, donor_key, strategy)
    specs = tuple(CommandComparator.from_json(path) for path in comparator_config)
    map_frame = pd.read_csv(label_map, dtype=str) if label_map else None
    evaluation_truth = apply_truth_label_map(adata.obs[truth_key], map_frame)
    benchmark_manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    benchmark_manifest["execution"] = {
        "input_sha256": compute_data_hash(input),
        "truth_exposure_policy": "test_h5ad_strips_truth_and_label_like_obs_columns",
        "cluster_key": cluster_key,
        "cluster_map_sha256": compute_data_hash(cluster_map) if cluster_map else None,
        "species": species,
        "tissue": tissue,
        "constant_study_id": constant_study_id,
        "methods": list(requested),
        "comparator_config_sha256": {
            str(path): compute_data_hash(path) for path in comparator_config
        },
        "label_map_sha256": compute_data_hash(label_map) if label_map else None,
        "continue_on_unavailable": continue_on_unavailable,
        "runtime_config": runtime_config,
        "checkpoint_contract": "atomic_per_method_fold_v1",
        "fold_ids": list(fold_id) if fold_id else None,
        "write_aggregate_tables": not no_aggregate_tables,
        "worker_id": worker_id,
        "batch_id": batch_id,
        "evaluation_unit": evaluation_unit,
    }
    paths["manifest"].write_text(json.dumps(benchmark_manifest, indent=2), encoding="utf-8")
    try:
        predictions, status = run_benchmark_comparators(
            adata,
            assignments,
            truth_key,
            cluster_key,
            output_dir,
            species,
            tissue,
            methods=requested,
            command_specs=specs,
            label_map=map_frame,
            continue_on_unavailable=continue_on_unavailable,
            fold_ids=tuple(fold_id) if fold_id else None,
            write_aggregate_tables=not no_aggregate_tables,
            worker_id=worker_id,
            batch_id=batch_id,
        )
    except BenchmarkValidationError as exc:
        console.print(f"[red]Benchmark execution failed: {exc}[/red]")
        raise typer.Exit(1) from exc

    status_path = Path(output_dir) / "comparator_status.csv"
    status.to_csv(status_path, index=False)
    paths["comparator_status"] = status_path
    if not predictions.empty:
        prediction_path = Path(output_dir) / "out_of_fold_predictions.csv"
        predictions.to_csv(prediction_path, index=False)
        paths["predictions"] = prediction_path
        if evaluation_unit in {"cell", "both"}:
            aggregate, per_fold = evaluate_holdout_predictions(
                evaluation_truth, assignments, predictions, expected_methods=requested
            )
            result_path = Path(output_dir) / "benchmark_results.csv"
            fold_path = Path(output_dir) / "benchmark_results_by_fold.csv"
            aggregate.to_csv(result_path, index=False)
            per_fold.to_csv(fold_path, index=False)
            bins, risk = build_calibration_artifacts(evaluation_truth, predictions)
            bins_path = Path(output_dir) / "calibration_bins.csv"
            risk_path = Path(output_dir) / "risk_coverage.csv"
            bins.to_csv(bins_path, index=False)
            risk.to_csv(risk_path, index=False)
            paths.update(
                {
                    "results": result_path,
                    "per_fold": fold_path,
                    "calibration_bins": bins_path,
                    "risk_coverage": risk_path,
                }
            )
            console.print(aggregate.to_string(index=False))
        if evaluation_unit in {"cluster", "both"}:
            cluster_truth, cluster_assignments, cluster_predictions, diagnostics = (
                build_cluster_level_track(
                    evaluation_truth,
                    assignments,
                    predictions,
                    adata.obs[cluster_key],
                )
            )
            cluster_results, cluster_folds = evaluate_holdout_predictions(
                cluster_truth,
                cluster_assignments,
                cluster_predictions,
                expected_methods=requested,
            )
            for name, frame in (
                ("cluster_track_results", cluster_results),
                ("cluster_track_results_by_fold", cluster_folds),
                ("cluster_track_predictions", cluster_predictions),
                ("cluster_track_diagnostics", diagnostics),
            ):
                path = Path(output_dir) / f"{name}.csv"
                frame.to_csv(path, index=False)
                paths[name] = path
    for label, path in paths.items():
        console.print(f"  {label}: {path}")


@app.command("benchmark-release")
def benchmark_release(
    registry: str = typer.Option(
        ...,
        "--registry",
        "-r",
        help="Locked public multi-cohort registry JSON",
    ),
    output_dir: str = typer.Option(
        "benchmark_release",
        "--output",
        "-o",
        help="Release artifact directory",
    ),
    allow_incomplete: bool = typer.Option(
        False,
        "--allow-incomplete",
        help="Write an explicitly non-claim-ready release when data or predictions are missing",
    ),
    n_boot: int = typer.Option(
        2000,
        "--n-boot",
        min=100,
        help="Donor/study hierarchical bootstrap replicates",
    ),
    seed: int = typer.Option(42, "--seed", help="Locked resampling seed"),
):
    """Build a donor-aware public benchmark release and retain negative results."""
    from .benchmark import BenchmarkValidationError
    from .benchmark_release import build_public_benchmark_release

    try:
        paths = build_public_benchmark_release(
            registry,
            output_dir,
            allow_incomplete=allow_incomplete,
            n_boot=n_boot,
            seed=seed,
        )
    except BenchmarkValidationError as exc:
        console.print(f"[red]Benchmark release failed: {exc}[/red]")
        raise typer.Exit(1) from exc
    manifest = json.loads(paths["release_manifest"].read_text(encoding="utf-8"))
    readiness = manifest["readiness"]
    colour = "green" if readiness["status"] == "claim_ready" else "yellow"
    console.print(f"[{colour}]Release status: {readiness['status']}[/{colour}]")
    console.print(
        f"  cohorts={readiness['n_evaluated_cohorts']}/{readiness['n_registered_cohorts']} "
        f"donors={readiness['n_independent_donors']} "
        f"blocking={readiness['blocking_findings']}"
    )
    for label, path in paths.items():
        console.print(f"  {label}: {path}")


@app.command()
def license(
    action: str = typer.Argument(..., help="Action: status, activate, deactivate"),
    key: str | None = typer.Option(None, "--key", "-k", help="License key (for activate)"),
    holder: str = typer.Option("", "--holder", help="License holder name"),
    email: str = typer.Option("", "--email", help="License holder email"),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """Manage CellTypePilot license."""
    from .license_manager import (
        ACADEMIC_FEATURES,
        COMMERCIAL_FEATURES,
        FREE_FEATURES,
        TRIAL_FEATURES,
        LicenseTier,
        activate_license,
        load_license,
    )

    if action == "status":
        lic = load_license()
        if json_output:
            from dataclasses import asdict

            data = asdict(lic)
            data["tier"] = lic.tier.value
            console.print(json.dumps(data, indent=2))
        else:
            available_features = {
                LicenseTier.FREE: FREE_FEATURES,
                LicenseTier.ACADEMIC: ACADEMIC_FEATURES,
                LicenseTier.COMMERCIAL: COMMERCIAL_FEATURES,
                LicenseTier.TRIAL: TRIAL_FEATURES,
            }[lic.tier]
            console.print("[bold]CellTypePilot License[/bold]")
            console.print(f"  Tier:      {lic.tier.value}")
            console.print(f"  Holder:    {lic.holder or 'N/A'}")
            console.print(f"  Email:     {lic.email or 'N/A'}")
            console.print(f"  Expires:   {lic.expires_at or 'Never'}")
            console.print(f"  Features:  {len(available_features)} available")
            if lic.is_expired():
                console.print("  [red]EXPIRED[/red]")
            console.print()
            console.print("[bold]Tier comparison:[/bold]")
            console.print(
                f"  Free:      {len(FREE_FEATURES)} features (all bundled MIT-licensed atlases)"
            )
            console.print(
                f"  Academic:  {len(ACADEMIC_FEATURES)} features "
                "(literature and review workflow services)"
            )
            console.print(
                f"  Commercial:{len(COMMERCIAL_FEATURES)} features "
                "(custom-panel services, team sharing, API)"
            )

    elif action == "activate":
        if not key:
            console.print("[red]--key is required for activation[/red]")
            raise typer.Exit(1)
        success, message = activate_license(key, holder=holder, email=email)
        if success:
            console.print(f"[green]{message}[/green]")
        else:
            console.print(f"[red]{message}[/red]")
            raise typer.Exit(1)

    elif action == "deactivate":
        from .license_manager import LicenseInfo, save_license

        save_license(LicenseInfo(tier=LicenseTier.FREE))
        console.print("[yellow]License deactivated. Reverted to free tier.[/yellow]")

    else:
        console.print(f"[red]Unknown action: {action}. Use: status, activate, deactivate[/red]")
        raise typer.Exit(1)


# ──────────────────────────────────────────────
# curate command (literature co-occurrence sweep)
# ──────────────────────────────────────────────
evidence_app = typer.Typer(
    name="evidence",
    help="Gene/cell identity contracts and human-gated marker-edge evidence",
    add_completion=False,
)
app.add_typer(evidence_app, name="evidence")


@evidence_app.command("propose-promotion")
def evidence_propose_promotion(
    atlas_path: str = typer.Option(..., "--atlas"),
    tissue: str = typer.Option(..., "--tissue"),
    cell_path: str = typer.Option(..., "--cell-path"),
    gene: str = typer.Option(..., "--gene"),
    polarity: str = typer.Option(..., "--polarity"),
    target_status: str = typer.Option(..., "--target-status"),
    evidence_json: str = typer.Option(..., "--evidence-json"),
    requested_by: str = typer.Option(..., "--requested-by"),
    output: str = typer.Option(..., "--output", "-o"),
    proposal_origin: str = typer.Option("human_curator", "--origin"),
):
    """Draft an immutable promotion proposal; this never changes the atlas."""
    from .evidence_promotion import build_promotion_proposal, write_promotion_proposal

    atlas = json.loads(Path(atlas_path).read_text(encoding="utf-8"))
    evidence = json.loads(Path(evidence_json).read_text(encoding="utf-8"))
    proposal = build_promotion_proposal(
        atlas,
        tissue=tissue,
        cell_path=cell_path,
        gene=gene,
        polarity=polarity,
        target_status=target_status,
        evidence=evidence,
        requested_by=requested_by,
        proposal_origin=proposal_origin,
    )
    path = write_promotion_proposal(proposal, output)
    console.print(f"Promotion proposal written: {path} (status={proposal['status']})")


@evidence_app.command("review-promotion")
def evidence_review_promotion(
    proposal_path: str = typer.Option(..., "--proposal"),
    reviewer: str = typer.Option(..., "--reviewer"),
    decision: str = typer.Option(..., "--decision"),
    notes: str = typer.Option(..., "--notes"),
    output: str | None = typer.Option(None, "--output", "-o"),
):
    """Append one independent human review to a promotion proposal."""
    from .evidence_promotion import add_promotion_review, write_promotion_proposal

    source = Path(proposal_path)
    proposal = json.loads(source.read_text(encoding="utf-8"))
    reviewed = add_promotion_review(proposal, reviewer=reviewer, decision=decision, notes=notes)
    path = write_promotion_proposal(reviewed, output or source)
    console.print(f"Promotion review written: {path} (status={reviewed['status']})")


@evidence_app.command("apply-promotion")
def evidence_apply_promotion(
    atlas_path: str = typer.Option(..., "--atlas"),
    proposal_path: str = typer.Option(..., "--proposal"),
    new_version: str = typer.Option(..., "--new-version"),
    output: str = typer.Option(..., "--output", "-o"),
):
    """Create a new atlas from a two-reviewer-approved proposal."""
    from .evidence_promotion import apply_approved_promotion

    atlas_source = Path(atlas_path).resolve()
    target = Path(output).resolve()
    if atlas_source == target:
        raise typer.BadParameter("--output must not overwrite the source atlas")
    atlas = json.loads(atlas_source.read_text(encoding="utf-8"))
    proposal = json.loads(Path(proposal_path).read_text(encoding="utf-8"))
    promoted = apply_approved_promotion(atlas, proposal, new_version=new_version)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(promoted, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    console.print(f"Promoted atlas written: {target}")


@app.command()
def curate(
    atlas_path: str | None = typer.Option(
        None,
        "--atlas",
        help="Atlas JSON to sweep (default: bundled marker atlas)",
    ),
    tissue: str | None = typer.Option(None, "--tissue", help="Restrict sweep to one tissue"),
    limit: int | None = typer.Option(None, "--limit", help="Max edges to sweep (pilot runs)"),
    delay: float = typer.Option(
        0.34, "--delay", help="Seconds between PubMed requests (rate limit)"
    ),
    output: str | None = typer.Option(
        None, "--output", "-o", help="Path for the sweep report JSON"
    ),
    apply: bool = typer.Option(
        False,
        "--apply",
        help="Write supported upgrades back into the atlas file (requires --new-version)",
    ),
    new_version: str | None = typer.Option(
        None, "--new-version", help="New atlas version when applying upgrades"
    ),
    json_output: bool = typer.Option(False, "--json", help="Output summary as JSON"),
):
    """Sweep aggregate-status marker edges for PubMed co-occurrence evidence.

    Supported edges are reported; with --apply they are upgraded to
    literature_cooccurrence_supported with a full evidence locator.
    """
    from .atlas_curation import (
        CurationError,
        apply_sweep_results,
        sweep_edges,
        write_sweep_report,
    )
    from .constants import ATLAS_PATH

    path = Path(atlas_path) if atlas_path else ATLAS_PATH
    try:
        atlas = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        console.print(f"[red]Cannot read atlas {path}: {exc}[/red]")
        raise typer.Exit(1) from exc

    def _progress(step: int, total: int, message: str):
        console.print(f"[dim]\rSweeping edge {step}/{total}: {message}[/dim]", end="")

    sweep = sweep_edges(atlas, tissue=tissue, delay_seconds=delay, limit=limit, progress=_progress)
    console.print()

    report_path = Path(output) if output else path.with_name(f"{path.stem}.sweep.json")
    write_sweep_report(sweep, report_path)

    applied = 0
    if apply:
        if not new_version:
            console.print("[red]--apply requires --new-version[/red]")
            raise typer.Exit(1)
        try:
            updated, applied = apply_sweep_results(atlas, sweep["results"], new_version)
        except CurationError as exc:
            console.print(f"[red]Apply failed (fail closed): {exc}[/red]")
            raise typer.Exit(1) from exc
        path.write_text(json.dumps(updated, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    summary = {
        "atlas": str(path),
        "swept": sweep["swept"],
        "supported": sweep["supported"],
        "errors": sweep["errors"],
        "report": str(report_path),
        "applied_upgrades": applied,
        "new_version": new_version if apply else None,
    }
    if json_output:
        print(json.dumps(summary, indent=2))
    else:
        console.print(
            f"Swept [cyan]{sweep['swept']}[/cyan] edges; "
            f"[green]{sweep['supported']}[/green] supported by literature co-occurrence; "
            f"[red]{sweep['errors']}[/red] search errors"
        )
        console.print(f"Report: {report_path}")
        if apply:
            console.print(
                f"Applied [green]{applied}[/green] upgrades → atlas version {new_version}"
            )


# ──────────────────────────────────────────────
# ontology commands (live Cell Ontology checks)
# ──────────────────────────────────────────────
ontology_app = typer.Typer(
    name="ontology",
    help="Cell Ontology cache and live identifier validation",
    add_completion=False,
)
app.add_typer(ontology_app, name="ontology")


@ontology_app.command("update")
def ontology_update(
    force: bool = typer.Option(False, "--force", help="Re-download even if cached"),
):
    """Download and cache the Cell Ontology (cl.obo) from the OBO Foundry."""
    from .ontology import OntologyError, download_ontology, ontology_cache_status

    try:
        metadata = download_ontology(force=force)
    except OntologyError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(1) from exc
    console.print(f"[green]Cell Ontology cache ready[/green] ({metadata.get('path')})")
    print(json.dumps(ontology_cache_status(), indent=2))


@ontology_app.command("check")
def ontology_check(
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """Validate bundled and installed atlas CL identifiers against the cache."""
    from .constants import ATLAS_PATH, FIRST_PARTY_PACKS_DIR
    from .ontology import (
        OntologyError,
        check_atlas_ontology,
        load_ontology,
        ontology_cache_status,
        summarize_findings,
    )

    status = ontology_cache_status()
    if not status.get("cached"):
        console.print(f"[yellow]{status['detail']}[/yellow]")
        raise typer.Exit(1)
    try:
        service = load_ontology()
    except OntologyError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(1) from exc

    targets = [("marker_atlas", ATLAS_PATH)]
    if FIRST_PARTY_PACKS_DIR.is_dir():
        for child in sorted(FIRST_PARTY_PACKS_DIR.iterdir()):
            atlas_file = child / "marker_atlas.json"
            if atlas_file.is_file():
                targets.append((f"pack:{child.name}", atlas_file))
    packs_root = None
    try:
        from .pack_manager import packs_dir

        packs_root = packs_dir()
    except Exception:
        pass
    if packs_root is not None and packs_root.is_dir():
        for child in sorted(packs_root.iterdir()):
            atlas_file = child / "marker_atlas.json"
            if atlas_file.is_file():
                targets.append((f"pack:{child.name}", atlas_file))

    report = {}
    exit_code = 0
    for label, atlas_path in targets:
        atlas = json.loads(atlas_path.read_text(encoding="utf-8"))
        findings = check_atlas_ontology(service, atlas)
        summary = summarize_findings(findings)
        report[label] = {"summary": summary, "findings": findings}
        if not summary["ok"]:
            exit_code = 1

    if json_output:
        print(json.dumps(report, indent=2, ensure_ascii=False))
    else:
        for label, entry in report.items():
            summary = entry["summary"]
            badge = "[green]OK[/green]" if summary["ok"] else "[red]ERRORS[/red]"
            console.print(
                f"{label}: {badge} | {summary['errors']} errors, "
                f"{summary['warnings']} label warnings"
            )
            for finding in entry["findings"]:
                if finding["severity"] == "error":
                    console.print(
                        f"  [red]- {finding['path']} [{finding['cl_id']}]: {finding['issue']}[/red]"
                    )
    if exit_code:
        raise typer.Exit(exit_code)


# ──────────────────────────────────────────────
# pack commands (domain extension packs)
# ──────────────────────────────────────────────
pack_app = typer.Typer(
    name="pack",
    help="Manage domain extension packs (data-only marker/state atlases)",
    add_completion=False,
)
app.add_typer(pack_app, name="pack")


@pack_app.command("install")
def pack_install(
    source: str = typer.Argument(
        ..., help="Local pack directory or git URL (e.g. https://.../celltypepilot-tme-pack.git)"
    ),
    trust: str = typer.Option(
        "atlas",
        "--trust",
        help="atlas: require full provenance validation; hypothesis: accept as draft evidence",
    ),
    force: bool = typer.Option(False, "--force", help="Reinstall over an existing pack"),
):
    """Install a domain extension pack from a local path or git URL."""
    from .pack_manager import PackError, install_pack

    try:
        summary = install_pack(source, trust=trust, force=force)
    except PackError as exc:
        console.print(f"[red]Pack install failed: {exc}[/red]")
        raise typer.Exit(1) from exc
    console.print(f"[green]Installed pack {summary['name']!r} v{summary['version']}[/green]")
    console.print(f"  Trust:   {summary['trust']}")
    console.print(f"  Path:    {summary['path']}")
    console.print(f"  Tissues: {', '.join(summary['tissues']) or '(none declared)'}")
    if summary["validation_issues"]:
        console.print("[yellow]  Installed at hypothesis trust with issues:[/yellow]")
        for issue in summary["validation_issues"][:5]:
            console.print(f"[yellow]  - {issue}[/yellow]")


@pack_app.command("list")
def pack_list(
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """List first-party and installed extension packs."""
    from .pack_manager import list_installed_packs

    entries = list_installed_packs()
    if json_output:
        console.print(json.dumps(entries, indent=2))
        return
    if not entries:
        console.print("No extension packs installed. Use: celltypepilot pack install <source>")
        return
    for entry in entries:
        console.print(f"[cyan]{entry['name']}[/cyan] v{entry['version']} ({entry['origin']})")
        console.print(f"  Trust: {entry['trust']} | License: {entry['license_tier']}")
        console.print(f"  Species: {', '.join(entry['species'])}")
        console.print(f"  Tissues: {', '.join(entry['tissues']) or '(none declared)'}")
        if entry["description"]:
            console.print(f"  [dim]{entry['description']}[/dim]")
        console.print()


@pack_app.command("validate")
def pack_validate(
    name: str = typer.Argument(..., help="Installed pack name"),
):
    """Re-validate an installed pack against schema and provenance gates."""
    from pathlib import Path

    from .pack_manager import list_installed_packs, validate_pack

    index = {entry["name"]: entry for entry in list_installed_packs()}
    if name not in index:
        console.print(f"[red]Pack {name!r} is not installed[/red]")
        raise typer.Exit(1)
    issues = validate_pack(Path(index[name]["path"]))
    if issues:
        console.print(f"[red]Pack {name!r} FAILED validation:[/red]")
        for issue in issues:
            console.print(f"  - {issue}")
        raise typer.Exit(1)
    console.print(f"[green]Pack {name!r} passed schema and provenance validation[/green]")


@pack_app.command("remove")
def pack_remove(
    name: str = typer.Argument(..., help="Installed pack name"),
):
    """Remove a user-installed extension pack."""
    from .pack_manager import PackError, remove_pack

    try:
        remove_pack(name)
    except PackError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(1) from exc
    console.print(f"[green]Removed pack {name!r}[/green]")


@pack_app.command("scaffold")
def pack_scaffold(
    name: str = typer.Argument(..., help="Pack name [a-z0-9_-]"),
    output: Path = typer.Option(..., "--output", "-o", help="Directory to create"),
    tissue: list[str] = typer.Option(["general"], "--tissue", help="Tissue keys (repeatable)"),
    disease: list[str] = typer.Option([], "--disease", help="Disease/context labels (repeatable)"),
    pack_kind: str = typer.Option("evidence", "--kind", help="evidence | reference | mixed"),
    version: str = typer.Option("0.1.0", "--version"),
    license_spdx: str = typer.Option("CC-BY-4.0", "--license"),
):
    """Scaffold a signed-ready, data-only evidence/reference pack (no code)."""
    from .pack_manager import PackError
    from .pack_signing import scaffold_pack

    try:
        root = scaffold_pack(
            output,
            name=name,
            version=version,
            tissues=list(tissue),
            diseases=list(disease),
            pack_kind=pack_kind,
            license_spdx=license_spdx,
        )
    except PackError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(1) from exc
    console.print(f"[green]Scaffolded data-only pack at {root}[/green]")
    console.print(
        "  Fill marker_atlas.json edges with provenance, then: celltypepilot pack sign <dir>"
    )


@pack_app.command("sign")
def pack_sign(
    pack_dir: Path = typer.Argument(..., help="Pack directory containing pack.json"),
    signer: str = typer.Option("local-curator", "--signer"),
    hmac_secret: str | None = typer.Option(
        None, "--hmac-secret", help="Optional HMAC secret (dev/community); else default dev secret"
    ),
    private_key: Path | None = typer.Option(
        None, "--private-key", help="Optional RSA PEM private key for pack.sig.json"
    ),
):
    """Sign a data-only pack (binds content hashes + license/provenance/ontology)."""
    from .pack_manager import PackError
    from .pack_signing import sign_pack

    pem = private_key.read_text(encoding="utf-8") if private_key else None
    try:
        payload = sign_pack(
            pack_dir,
            private_key_pem=pem,
            hmac_secret=hmac_secret,
            signer=signer,
        )
    except PackError as exc:
        console.print(f"[red]Sign failed: {exc}[/red]")
        raise typer.Exit(1) from exc
    console.print(f"[green]Signed pack[/green] algorithm={payload['algorithm']}")
    console.print(f"  fingerprint: {payload['fingerprint_sha256'][:16]}…")
    console.print(f"  wrote: {pack_dir / 'pack.sig.json'}")


@pack_app.command("verify")
def pack_verify(
    pack_dir: Path = typer.Argument(..., help="Pack directory"),
    require_signature: bool = typer.Option(
        False, "--require-signature", help="Fail if pack.sig.json missing/invalid"
    ),
    json_output: bool = typer.Option(False, "--json"),
):
    """Validate data-only policy, provenance, and signature for a pack."""
    from .pack_signing import validate_pack_ecosystem, verify_pack_signature

    issues = validate_pack_ecosystem(pack_dir, require_signature=require_signature)
    sig = verify_pack_signature(pack_dir)
    report = {
        "issues": issues,
        "signature": sig,
        "ok": not issues and (sig.get("valid") or not require_signature),
    }
    if json_output:
        print(json.dumps(report, indent=2))
    else:
        if issues:
            console.print("[red]Pack validation FAILED[/red]")
            for issue in issues:
                console.print(f"  - {issue}")
        else:
            console.print("[green]Pack data-only + provenance OK[/green]")
        console.print(f"  signature: {sig.get('status')} ({sig.get('reason')})")
    if issues or (require_signature and not sig.get("valid")):
        raise typer.Exit(1)


@app.command("review-resign")
def review_resign_cmd(
    output_dir: Path = typer.Option(..., "--output", "-o", help="Annotation output directory"),
    signer: str = typer.Option("cli-reviewer", "--signer"),
    no_regenerate: bool = typer.Option(
        False, "--no-regenerate", help="Only re-hash existing artifacts without HTML refresh"
    ),
    json_output: bool = typer.Option(False, "--json"),
):
    """Regenerate derived artifacts after human edits and re-sign content hashes.

    Append-only audit is preserved. Clears stale flags only after successful resign.
    """
    from .review_resign import ReviewResignError, resign_review_outputs

    try:
        result = resign_review_outputs(output_dir, signer=signer, regenerate=not no_regenerate)
    except ReviewResignError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(1) from exc
    if json_output:
        print(json.dumps(result, indent=2, default=str))
    else:
        sig = result["signature"]
        console.print("[green]Review outputs re-signed[/green]")
        console.print(f"  signature: {sig['signature_sha256'][:16]}…")
        console.print(f"  regenerated: {sig.get('regenerated')}")
        console.print(f"  artifact_status: {result['artifact_status']['review_state']}")


# ──────────────────────────────────────────────
# assets commands (immutable object-store catalog)
# ──────────────────────────────────────────────
assets_app = typer.Typer(
    name="assets",
    help="Immutable asset catalog (CELLxGENE, Azimuth refs, label maps, Docker images)",
    add_completion=False,
)
app.add_typer(assets_app, name="assets")

_DEFAULT_ASSET_CATALOG = (
    Path(__file__).resolve().parents[2] / "benchmarks" / "assets" / "catalog.json"
)
_DEFAULT_ASSET_POLICY = (
    Path(__file__).resolve().parents[2] / "benchmarks" / "assets" / "storage_policy.json"
)


def _resolve_asset_paths(catalog: Path | None, policy: Path | None) -> tuple[Path, Path]:
    catalog_path = Path(catalog) if catalog is not None else _DEFAULT_ASSET_CATALOG
    policy_path = Path(policy) if policy is not None else _DEFAULT_ASSET_POLICY
    # Prefer package-adjacent repo layout; fall back to CWD-relative benchmarks/assets.
    if not catalog_path.is_file():
        cwd_catalog = Path("benchmarks/assets/catalog.json")
        if cwd_catalog.is_file():
            catalog_path = cwd_catalog
    if not policy_path.is_file():
        cwd_policy = Path("benchmarks/assets/storage_policy.json")
        if cwd_policy.is_file():
            policy_path = cwd_policy
    return catalog_path, policy_path


@assets_app.command("list")
def assets_list(
    kind: str | None = typer.Option(None, "--kind", help="Filter by asset kind"),
    availability: str | None = typer.Option(
        None, "--availability", help="Filter by availability status"
    ),
    catalog: Path | None = typer.Option(None, "--catalog", help="Catalog JSON path"),
    policy: Path | None = typer.Option(None, "--policy", help="Storage policy JSON path"),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """List immutable assets with URL, version, SHA-256, license, and status."""
    from .asset_catalog import (
        AssetCatalogError,
        filter_assets,
        load_asset_catalog,
        load_storage_policy,
        summarize_catalog,
    )

    catalog_path, policy_path = _resolve_asset_paths(catalog, policy)
    try:
        payload, _ = load_asset_catalog(catalog_path)
        storage, _ = load_storage_policy(policy_path)
    except (OSError, json.JSONDecodeError, AssetCatalogError) as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(1) from exc

    kinds = [kind] if kind else None
    avails = [availability] if availability else None
    rows = filter_assets(payload, kinds=kinds, availability=avails)
    view = {**payload, "assets": rows}
    summary = summarize_catalog(view, policy=storage)
    if json_output:
        print(json.dumps(summary, indent=2))
        return
    console.print(
        f"[bold]{payload.get('catalog_id')}[/bold]  ({len(rows)} assets; never writes fold runs/)"
    )
    for asset in summary["assets"]:
        console.print(
            f"  [{asset['kind']}] {asset['asset_id']}@{asset['version']}  "
            f"{asset['availability']}  sha256={asset['sha256'][:12]}…  "
            f"{asset['species']}/{asset['tissue']}  license={asset['license']}"
        )
        console.print(f"      url={asset['url']}")


@assets_app.command("verify")
def assets_verify(
    catalog: Path | None = typer.Option(None, "--catalog", help="Catalog JSON path"),
    policy: Path | None = typer.Option(None, "--policy", help="Storage policy JSON path"),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """Validate catalog schema and verify local file:/object-cache digests."""
    from .asset_catalog import (
        AssetCatalogError,
        load_asset_catalog,
        load_storage_policy,
        summarize_catalog,
    )

    catalog_path, policy_path = _resolve_asset_paths(catalog, policy)
    try:
        payload, resolved_catalog = load_asset_catalog(catalog_path)
        storage, _ = load_storage_policy(policy_path)
        summary = summarize_catalog(
            payload,
            policy=storage,
            catalog_root=resolved_catalog.parent,
            verify_local=True,
        )
    except (OSError, json.JSONDecodeError, AssetCatalogError) as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(1) from exc

    if json_output:
        print(json.dumps(summary, indent=2))
    else:
        console.print(f"[green]Catalog OK[/green]: {summary['catalog_id']}")
        console.print(f"  assets: {summary['n_assets']}")
        console.print(f"  availability: {summary['by_availability']}")
        console.print(
            f"  local verified: {summary.get('local_verified_count', 0)}; "
            f"hard failures: {summary.get('local_failure_count', 0)}"
        )
        for row in summary.get("local_verification", []):
            colour = (
                "green"
                if row["status"] == "verified" or row["status"].startswith("ok_")
                else "yellow"
            )
            if row["status"] in {"sha256_mismatch", "byte_size_mismatch"}:
                colour = "red"
            console.print(f"  [{colour}]{row['asset_id']}[/{colour}]: {row['status']}")

    hard = summary.get("local_failure_count", 0)
    available_missing = sum(
        1
        for row in summary.get("local_verification", [])
        if row["status"] == "missing_local" and row["availability"] == "available"
    )
    if hard or available_missing:
        raise typer.Exit(1)


@assets_app.command("materialize")
def assets_materialize(
    catalog: Path | None = typer.Option(None, "--catalog", help="Catalog JSON path"),
    policy: Path | None = typer.Option(None, "--policy", help="Storage policy JSON path"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Report without writing"),
    kind: str | None = typer.Option(None, "--kind", help="Limit to one asset kind"),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """Copy file: sources into content-addressed objects/ (never under runs/)."""
    from .asset_catalog import (
        AssetCatalogError,
        filter_assets,
        load_asset_catalog,
        load_storage_policy,
        materialize_source_to_object_cache,
    )

    catalog_path, policy_path = _resolve_asset_paths(catalog, policy)
    try:
        payload, resolved_catalog = load_asset_catalog(catalog_path)
        storage, _ = load_storage_policy(policy_path)
        assets = filter_assets(payload, kinds=[kind] if kind else None)
        rows = [
            materialize_source_to_object_cache(
                asset,
                resolved_catalog.parent,
                storage,
                dry_run=dry_run,
            )
            for asset in assets
        ]
    except (OSError, json.JSONDecodeError, AssetCatalogError) as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(1) from exc

    if json_output:
        print(json.dumps({"materialize": rows}, indent=2))
        return
    for row in rows:
        console.print(f"  {row['asset_id']}: {row['status']}")
    failures = {
        "source_missing",
        "source_sha256_mismatch",
        "source_byte_size_mismatch",
        "target_exists_with_mismatch",
        "write_verify_failed",
    }
    if any(row["status"] in failures for row in rows):
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
