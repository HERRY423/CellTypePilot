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
    help="CellTypePilot — Single-cell annotation intelligence layer",
    add_completion=False,
)


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
    """CellTypePilot — Single-cell annotation intelligence layer."""
    pass


# ──────────────────────────────────────────────
# doctor command
# ──────────────────────────────────────────────
@app.command()
def doctor():
    """Check environment: Python version, dependencies, capability level."""
    from .doctor import print_doctor

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
        console.print(json.dumps(report, indent=2))
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
):
    """Run the full annotation pipeline: marker scoring → critic → report."""
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
        console.print(json.dumps(output_json, indent=2, default=str))

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

    atlas = load_marker_atlas(species)
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
        console.print(f"  Score:       {row.get('combined_score', 0):.3f}")
        console.print(f"  Confidence:  {row.get('critic_confidence', 'N/A')}")
        console.print(f"  Flags:       {row.get('critic_flags', 'PASS')}")
        console.print(f"  Evidence:    {row.get('critic_evidence', '')}")
        console.print(f"  Notes:       {row.get('critic_notes', '')}")

    # Show top-5 candidates
    console.print(f"\n[bold]Top 5 Candidates for Cluster {focus}:[/bold]")
    cluster_scores = scores[scores["cluster"] == focus].head(5)
    for _, row in cluster_scores.iterrows():
        console.print(
            f"  #{int(row['rank'])} {row['cell_type']} (score={row['combined_score']:.3f}, "
            f"overlap={row['pct_overlap']:.0%}, neg_conflict={row['neg_conflict']:.0%})"
        )


# ──────────────────────────────────────────────
# markers command (list available markers)
# ──────────────────────────────────────────────
@app.command()
def markers(
    tissue: str | None = typer.Option(None, "--tissue", "-t", help="Tissue to list markers for"),
    species: str = typer.Option("human", "--species", "-s", help="Species: human/mouse"),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """List available cell types and markers in the knowledge graph."""
    from .data_adapter import get_all_markers_for_tissue, load_marker_atlas

    atlas = load_marker_atlas(species)
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
        console.print(json.dumps(markers_dict, indent=2))
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
        help="Path to CellTypePilot output directory",
    ),
    host: str = typer.Option("127.0.0.1", "--host", help="Host to bind to"),
    port: int = typer.Option(8765, "--port", "-p", help="Port to listen on"),
):
    """Launch the Web Inspector — interactive annotation review panel."""
    from .web_inspector import run_inspector

    if not output_dir.exists():
        console.print(f"[red]Output directory not found: {output_dir}[/red]")
        console.print("Run 'celltypepilot annotate' first to generate output.")
        raise typer.Exit(1)

    console.print("[bold]Launching Web Inspector...[/bold]")
    console.print(f"  Output dir: {output_dir}")
    console.print(f"  URL: http://{host}:{port}")
    console.print()

    try:
        run_inspector(output_dir, host=host, port=port)
    except KeyboardInterrupt:
        console.print("\n[yellow]Web Inspector stopped.[/yellow]")
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1) from e


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
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """Annotate using reference embedding + marker ensemble.

    Combines marker-based scoring with deep learning reference mapping
    (CellTypist / scANVI / correlation) to resolve continuous
    differentiation trajectories and rare transitional states.

    Examples:
        # Auto-detect backend with reference atlas
        celltypepilot annotate-embedding -i data.h5ad -k leiden -r reference.h5ad

        # Use CellTypist pre-trained model
        celltypepilot annotate-embedding -i data.h5ad -k leiden -m models/Immune_All_Low.pkl

        # Force correlation backend (no extra deps)
        celltypepilot annotate-embedding -i data.h5ad -k leiden -r ref.h5ad -b correlation
    """
    from .data_adapter import (
        detect_species,
        detect_tissue,
        get_all_markers_for_tissue,
        load_h5ad,
        load_marker_atlas,
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


# ──────────────────────────────────────────────
# license command
# ──────────────────────────────────────────────
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
            console.print("[bold]CellTypePilot License[/bold]")
            console.print(f"  Tier:      {lic.tier.value}")
            console.print(f"  Holder:    {lic.holder or 'N/A'}")
            console.print(f"  Email:     {lic.email or 'N/A'}")
            console.print(f"  Expires:   {lic.expires_at or 'Never'}")
            console.print(f"  Features:  {len(lic.features)} enabled")
            if lic.is_expired():
                console.print("  [red]EXPIRED[/red]")
            console.print()
            console.print("[bold]Tier comparison:[/bold]")
            console.print(f"  Free:      {len(FREE_FEATURES)} features (basic atlas, 11 tissues)")
            console.print(
                f"  Academic:  {len(ACADEMIC_FEATURES)} features (extended atlas, disease states)"
            )
            console.print(
                f"  Commercial:{len(COMMERCIAL_FEATURES)} features (full atlas, custom panels, API)"
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


if __name__ == "__main__":
    app()
