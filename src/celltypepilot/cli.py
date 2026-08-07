"""CellTypePilot CLI — command-line interface."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

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
        False, "--version", "-v", callback=version_callback,
        is_eager=True, help="Show version and exit.",
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
    cluster_key: Optional[str] = typer.Option(None, "--cluster-key", "-k", help="Cluster key in obs"),
    embedding_key: Optional[str] = typer.Option(None, "--embedding-key", "-e", help="Embedding key in obsm"),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """Inspect an h5ad file: detect species, tissue, clusters, embeddings, layers."""
    from .data_adapter import inspect_adata, format_inspect_report

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
    species: Optional[str] = typer.Option(None, "--species", "-s", help="Species: human/mouse (auto-detect if omitted)"),
    tissue: Optional[str] = typer.Option(None, "--tissue", "-t", help="Tissue context (e.g., blood, lung, brain)"),
    embedding_key: Optional[str] = typer.Option(None, "--embedding-key", "-e", help="Embedding key in obsm"),
    layer: Optional[str] = typer.Option(None, "--layer", help="Layer to use for expression (default: X)"),
    json_output: bool = typer.Option(False, "--json", help="Output results as JSON"),
    no_figures: bool = typer.Option(False, "--no-figures", help="Skip figure generation"),
):
    """Run the full annotation pipeline: marker scoring → critic → report."""
    from .data_adapter import (
        load_h5ad, compute_data_hash, detect_species, detect_tissue,
        load_marker_atlas, get_all_markers_for_tissue,
    )
    from .marker_scorer import compute_marker_scores, generate_annotation_summary
    from .critic import run_critic, generate_critic_summary
    from .visualizer import generate_all_figures
    from .reporter import save_evidence_table, generate_html_report, generate_methodology_text
    from .provenance import create_manifest, update_manifest_outputs, save_manifest, format_manifest_summary

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Step 1: Load data
    console.print("[bold blue]Step 1/6:[/bold blue] Loading data...")
    adata = load_h5ad(input)
    data_hash = compute_data_hash(input)

    # Step 2: Detect/auto-set parameters
    console.print("[bold blue]Step 2/6:[/bold blue] Detecting parameters...")
    if species is None:
        species = detect_species(adata)
        console.print(f"  Detected species: [cyan]{species}[/cyan]")
    if tissue is None:
        tissue = detect_tissue(adata)
        if tissue:
            console.print(f"  Detected tissue: [cyan]{tissue}[/cyan]")
        else:
            tissue = "general"
            console.print(f"  Tissue not detected, using [cyan]general[/cyan] marker set")
    if embedding_key is None:
        from .data_adapter import find_embedding_keys
        candidates = find_embedding_keys(adata)
        if candidates:
            embedding_key = candidates[0]
            console.print(f"  Using embedding: [cyan]{embedding_key}[/cyan]")
        else:
            console.print("[yellow]  Warning: No embedding found. Figures will be skipped.[/yellow]")

    # Validate cluster key
    if cluster_key not in adata.obs.columns:
        console.print(f"[red]Error: cluster key '{cluster_key}' not found in obs.[/red]")
        console.print(f"Available columns: {list(adata.obs.columns)}")
        raise typer.Exit(1)

    # Step 3: Load marker atlas and score
    console.print("[bold blue]Step 3/6:[/bold blue] Computing marker scores...")
    atlas = load_marker_atlas(species)
    markers = get_all_markers_for_tissue(atlas, tissue)
    console.print(f"  Using {len(markers)} cell types from '{tissue}' tissue atlas")

    scores = compute_marker_scores(adata, cluster_key, markers, layer=layer)
    summary = generate_annotation_summary(scores, cluster_key)

    if summary.empty:
        console.print("[red]Error: No annotations generated. Check marker gene overlap with your data.[/red]")
        raise typer.Exit(1)

    console.print(f"  Annotated {len(summary)} clusters")

    # Step 4: Run critic
    console.print("[bold blue]Step 4/6:[/bold blue] Running Annotation Critic...")
    critic_results = run_critic(adata, cluster_key, summary, atlas, tissue)
    critic_summary = generate_critic_summary(critic_results)

    console.print(f"  Passed: [green]{critic_summary['pass']}[/green] | "
                  f"Flagged: [red]{critic_summary['flagged']}[/red]")

    # Step 5: Generate figures
    figure_paths = []
    if not no_figures and embedding_key:
        console.print("[bold blue]Step 5/6:[/bold blue] Generating figures...")
        figure_paths = generate_all_figures(
            adata, cluster_key, embedding_key, critic_results, output_path, tissue
        )
        console.print(f"  Generated {len(figure_paths)} figures")
    else:
        console.print("[bold blue]Step 5/6:[/bold blue] Skipping figures")

    # Step 6: Save outputs
    console.print("[bold blue]Step 6/6:[/bold blue] Saving outputs...")

    # Evidence table
    evidence_path = save_evidence_table(critic_results, output_path)
    console.print(f"  Evidence table: {evidence_path}")

    # Manifest
    manifest = create_manifest(
        input_path=input,
        data_hash=data_hash,
        cluster_key=cluster_key,
        species=species,
        tissue=tissue,
        parameters={
            "embedding_key": embedding_key,
            "layer": layer,
        },
        output_dir=output_path,
    )

    # HTML report
    report_path = generate_html_report(
        critic_results, critic_results, critic_summary, manifest, figure_paths, output_path
    )
    console.print(f"  HTML report: {report_path}")

    # Methodology text
    method_text = generate_methodology_text(manifest, critic_summary, critic_results)
    method_path = output_path / "methodology_draft.txt"
    with open(method_path, "w") as f:
        f.write(method_text)
    console.print(f"  Methodology draft: {method_path}")

    # Update and save manifest
    manifest = update_manifest_outputs(manifest, output_path)
    manifest_path = save_manifest(manifest, output_path)
    console.print(f"  Manifest: {manifest_path}")

    # Write annotations back to adata
    _write_annotations_to_adata(adata, critic_results, cluster_key, output_path)

    # JSON output
    if json_output:
        output_json = {
            "annotations": critic_results.to_dict(orient="records"),
            "critic_summary": critic_summary,
            "manifest": manifest,
        }
        console.print(json.dumps(output_json, indent=2, default=str))

    console.print("\n[bold green]Done![/bold green] CellTypePilot annotation complete.")
    console.print(f"Output directory: {output_path.resolve()}")


def _write_annotations_to_adata(
    adata, critic_results: "pd.DataFrame", cluster_key: str, output_dir: Path
):
    """Write annotation results back into adata obs and save."""
    import anndata as ad

    # Map cluster → annotation
    cluster_to_ct = dict(zip(critic_results["cluster"], critic_results["cell_type"]))
    cluster_to_cl = dict(zip(critic_results["cluster"], critic_results.get("cl_id", [""] * len(critic_results))))
    cluster_to_conf = dict(zip(critic_results["cluster"], critic_results.get("critic_confidence", [""] * len(critic_results))))

    adata.obs["ctp_cell_type"] = adata.obs[cluster_key].map(cluster_to_ct).fillna("Unknown")
    adata.obs["ctp_cl_id"] = adata.obs[cluster_key].map(cluster_to_cl).fillna("")
    adata.obs["ctp_confidence"] = adata.obs[cluster_key].map(cluster_to_conf).fillna("unknown")

    output_path = output_dir / "data.annotated.h5ad"
    adata.write(output_path)


# ──────────────────────────────────────────────
# critic command (re-review a specific cluster)
# ──────────────────────────────────────────────
@app.command()
def critic(
    input: str = typer.Option(..., "--input", "-i", help="Path to .h5ad file"),
    cluster_key: str = typer.Option(..., "--cluster-key", "-k", help="Cluster key in obs"),
    focus: str = typer.Option(..., "--focus", "-f", help="Cluster ID to deep-review"),
    species: Optional[str] = typer.Option(None, "--species", "-s", help="Species"),
    tissue: Optional[str] = typer.Option(None, "--tissue", "-t", help="Tissue context"),
):
    """Deep-review a specific cluster flagged by the critic."""
    from .data_adapter import load_h5ad, detect_species, detect_tissue, load_marker_atlas
    from .marker_scorer import compute_marker_scores, generate_annotation_summary
    from .critic import run_critic

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
        console.print(f"  #{int(row['rank'])} {row['cell_type']} (score={row['combined_score']:.3f}, "
                      f"overlap={row['pct_overlap']:.0%}, neg_conflict={row['neg_conflict']:.0%})")


# ──────────────────────────────────────────────
# markers command (list available markers)
# ──────────────────────────────────────────────
@app.command()
def markers(
    tissue: Optional[str] = typer.Option(None, "--tissue", "-t", help="Tissue to list markers for"),
    species: str = typer.Option("human", "--species", "-s", help="Species: human/mouse"),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """List available cell types and markers in the knowledge graph."""
    from .data_adapter import load_marker_atlas, get_all_markers_for_tissue

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
    markers: Optional[str] = typer.Option(None, "--markers", "-m", help="Comma-separated marker genes"),
    max_refs: int = typer.Option(5, "--max-refs", help="Max references per query"),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """Search literature for marker validation (PubMed/bioRxiv)."""
    from .literature import (
        validate_marker_in_literature,
        validate_annotation_with_literature,
        check_mcp_availability,
        generate_mcp_search_queries,
    )

    # Check MCP availability
    mcp_status = check_mcp_availability()
    if not mcp_status.get("pubmed_direct"):
        console.print("[yellow]Warning: PubMed direct access not available. Check network.[/yellow]")

    marker_list = [m.strip() for m in markers.split(",")] if markers else []

    if marker_list:
        # Validate specific markers for a cell type
        results = validate_annotation_with_literature(cell_type, marker_list, max_refs_per_marker=max_refs)

        if json_output:
            console.print(json.dumps(results, indent=2))
        else:
            console.print(f"\n[bold]Literature Validation for '{cell_type}':[/bold]\n")
            console.print(f"  Positive markers checked: {results['positive_markers_checked']}")
            console.print(f"  Markers supported: {results['positive_markers_supported']}")
            console.print(f"  Total refs found: {results['total_literature_refs']}")
            console.print(f"  Assessment: {results['overall_assessment']}\n")

            for ev in results.get("positive_evidence", []):
                status = "[green]OK[/green]" if ev["consensus"] == "supported" else "[yellow]?[/yellow]"
                console.print(f"  {status} {ev['gene']}: {ev['total_refs']} refs")
                for hit in ev.get("top_hits", []):
                    console.print(f"      - {hit['authors']} ({hit['year']}). {hit['title'][:60]}...")
    else:
        # Generate search queries for manual MCP use
        queries = generate_mcp_search_queries(cell_type, [])
        if json_output:
            console.print(json.dumps({"queries": queries, "mcp_status": mcp_status}, indent=2))
        else:
            console.print(f"\n[bold]Suggested search queries for '{cell_type}':[/bold]\n")
            for i, q in enumerate(queries, 1):
                console.print(f"  {i}. {q}")
            console.print(f"\n[bold]MCP Status:[/bold]")
            for tool, available in mcp_status.items():
                status = "[green]available[/green]" if available else "[red]not available[/red]"
                console.print(f"  {tool}: {status}")


# ──────────────────────────────────────────────
# inspect-web command (Web Inspector)
# ──────────────────────────────────────────────
@app.command()
def inspect_web(
    output_dir: Path = typer.Option(
        "./ctp_output", "--output", "-o",
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

    console.print(f"[bold]Launching Web Inspector...[/bold]")
    console.print(f"  Output dir: {output_dir}")
    console.print(f"  URL: http://{host}:{port}")
    console.print()

    try:
        run_inspector(output_dir, host=host, port=port)
    except KeyboardInterrupt:
        console.print("\n[yellow]Web Inspector stopped.[/yellow]")
    except Exception as e:
        console.print(f"[red]Error: {e}[/red]")
        raise typer.Exit(1)


# ──────────────────────────────────────────────
# convert-rds command (Seurat support)
# ──────────────────────────────────────────────
@app.command()
def convert_rds(
    input_rds: Path = typer.Option(..., "--input", "-i", help="Path to Seurat .rds file"),
    output_h5ad: Optional[Path] = typer.Option(None, "--output", "-o", help="Output .h5ad path"),
):
    """Convert Seurat .rds to .h5ad for CellTypePilot annotation."""
    from .seurat_adapter import seurat_to_h5ad, check_seurat_support

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

    console.print(f"[bold]Converting Seurat .rds to .h5ad...[/bold]")
    console.print(f"  Input:  {input_rds}")
    console.print(f"  Output: {output_h5ad}")

    try:
        result_path = seurat_to_h5ad(input_rds, output_h5ad)
        console.print(f"[green]Conversion complete: {result_path}[/green]")
        console.print(f"\nNow run: celltypepilot annotate --input {result_path}")
    except Exception as e:
        console.print(f"[red]Conversion failed: {e}[/red]")
        raise typer.Exit(1)


# ──────────────────────────────────────────────
# apply-overrides command (write corrections to .h5ad)
# ──────────────────────────────────────────────
@app.command()
def apply_overrides(
    output_dir: Path = typer.Option(
        "./ctp_output", "--output", "-o",
        help="CellTypePilot output directory containing data.annotated.h5ad",
    ),
    overrides_file: Path = typer.Option(
        ..., "--overrides", "-f",
        help="Path to annotation_overrides.json",
    ),
    regenerate: bool = typer.Option(
        False, "--regenerate", "-r",
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
    import shutil
    from datetime import datetime

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
        raise typer.Exit(1)

    if not overrides:
        console.print("[yellow]No overrides to apply.[/yellow]")
        return

    console.print(f"[bold]Applying {len(overrides)} override(s)...[/bold]")

    # Backup original
    backup_name = f"data.annotated.backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}.h5ad"
    backup_path = output_dir / backup_name
    shutil.copy2(h5ad_path, backup_path)
    console.print(f"  Backup: {backup_path}")

    # Load and modify
    import scanpy as sc
    adata = sc.read_h5ad(h5ad_path)
    obs = adata.obs

    applied = 0
    skipped = 0

    for cluster_id, override in overrides.items():
        new_type = override.get("new_type", "")
        reason = override.get("reason", "")
        if not new_type:
            skipped += 1
            continue

        # Find cluster column
        cluster_col = None
        for col in ["ctp_cl_id", "cluster", "cl_id"]:
            if col in obs.columns:
                cluster_col = col
                break

        if cluster_col is None:
            console.print(f"  [yellow]Warning: No cluster column found, skipping[/yellow]")
            skipped += 1
            continue

        mask = obs[cluster_col].astype(str) == str(cluster_id)
        n_cells = mask.sum()

        if n_cells == 0:
            console.print(f"  [yellow]Cluster {cluster_id}: no cells found, skipped[/yellow]")
            skipped += 1
            continue

        old_type = obs.loc[mask, "ctp_cell_type"].iloc[0] if "ctp_cell_type" in obs.columns else "Unknown"

        # Apply
        if "ctp_cell_type" in obs.columns:
            adata.obs.loc[mask, "ctp_cell_type"] = new_type
        if "ctp_override_reason" not in obs.columns:
            adata.obs["ctp_override_reason"] = ""
        adata.obs.loc[mask, "ctp_override_reason"] = reason
        if "ctp_overridden" not in obs.columns:
            adata.obs["ctp_overridden"] = False
        adata.obs.loc[mask, "ctp_overridden"] = True

        applied += 1
        console.print(f"  [green]OK[/green] Cluster {cluster_id}: {old_type} → {new_type} ({n_cells} cells)")

    # Save
    adata.write(h5ad_path)
    console.print(f"\n[bold green]Applied {applied} override(s)[/bold green], {skipped} skipped")

    # Regenerate figures if requested
    if regenerate:
        console.print("\n[bold blue]Regenerating figures...[/bold blue]")
        try:
            from .visualizer import generate_all_figures
            embedding_key = None
            from .data_adapter import find_embedding_keys
            candidates = find_embedding_keys(adata)
            if candidates:
                embedding_key = candidates[0]
            if embedding_key:
                tissue = "general"
                figure_paths = generate_all_figures(
                    adata, cluster_col, embedding_key, None, output_dir, tissue
                )
                console.print(f"  Regenerated {len(figure_paths)} figures")
            else:
                console.print("  [yellow]No embedding found, skipping figure regeneration[/yellow]")
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
    reference: Optional[str] = typer.Option(None, "--reference", "-r", help="Reference .h5ad with cell type labels"),
    ref_label_key: str = typer.Option("cell_type", "--ref-label", help="Cell type column in reference.obs"),
    model_path: Optional[str] = typer.Option(None, "--model", "-m", help="CellTypist model path (.pkl)"),
    backend: str = typer.Option("auto", "--backend", "-b", help="Backend: auto/celltypist/scanvi/correlation"),
    output_dir: str = typer.Option(".", "--output", "-o", help="Output directory"),
    species: Optional[str] = typer.Option(None, "--species", "-s", help="Species: human/mouse"),
    tissue: Optional[str] = typer.Option(None, "--tissue", "-t", help="Tissue context"),
    marker_weight: float = typer.Option(0.5, "--marker-weight", help="Base marker weight (0-1)"),
    no_ensemble: bool = typer.Option(False, "--no-ensemble", help="Skip ensemble, use reference only"),
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
    from .data_adapter import load_h5ad, detect_species, detect_tissue, load_marker_atlas, get_all_markers_for_tissue
    from .marker_scorer import compute_marker_scores, generate_annotation_summary
    from .reference_scorer import score_by_reference, detect_transitional_states, check_reference_backends
    from .ensemble_scorer import ensemble_scores, generate_ensemble_summary, analyze_disagreements

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
        adata, cluster_key,
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
            marker_scores, ref_scores,
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
        console.print(f"\n  [bold]Disagreements (potential novel/transitional states):[/bold]")
        for _, row in disagreements.head(5).iterrows():
            console.print(f"    Cluster {row['cluster']}: "
                          f"marker→{row['marker_type']} ({row['marker_score']:.2f}) "
                          f"vs ref→{row['ref_type']} ({row['ref_score']:.2f})")
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
            "disagreements": disagreements.to_dict(orient="records") if not disagreements.empty else [],
        }
        console.print(json.dumps(output_json, indent=2, default=str))

    console.print(f"\n[bold green]Done![/bold green] Ensemble annotation complete.")


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
    console.print("  Or provide --reference .h5ad  # Correlation (no extra deps)")


# ──────────────────────────────────────────────
# license command
# ──────────────────────────────────────────────
@app.command()
def license(
    action: str = typer.Argument(..., help="Action: status, activate, deactivate"),
    key: Optional[str] = typer.Option(None, "--key", "-k", help="License key (for activate)"),
    holder: str = typer.Option("", "--holder", help="License holder name"),
    email: str = typer.Option("", "--email", help="License holder email"),
    json_output: bool = typer.Option(False, "--json", help="Output as JSON"),
):
    """Manage CellTypePilot license."""
    from .license_manager import (
        load_license, activate_license, LicenseTier,
        FREE_FEATURES, ACADEMIC_FEATURES, COMMERCIAL_FEATURES,
    )

    if action == "status":
        lic = load_license()
        if json_output:
            from dataclasses import asdict
            data = asdict(lic)
            data["tier"] = lic.tier.value
            console.print(json.dumps(data, indent=2))
        else:
            console.print(f"[bold]CellTypePilot License[/bold]")
            console.print(f"  Tier:      {lic.tier.value}")
            console.print(f"  Holder:    {lic.holder or 'N/A'}")
            console.print(f"  Email:     {lic.email or 'N/A'}")
            console.print(f"  Expires:   {lic.expires_at or 'Never'}")
            console.print(f"  Features:  {len(lic.features)} enabled")
            if lic.is_expired():
                console.print(f"  [red]EXPIRED[/red]")
            console.print()
            console.print(f"[bold]Tier comparison:[/bold]")
            console.print(f"  Free:      {len(FREE_FEATURES)} features (basic atlas, 11 tissues)")
            console.print(f"  Academic:  {len(ACADEMIC_FEATURES)} features (extended atlas, disease states)")
            console.print(f"  Commercial:{len(COMMERCIAL_FEATURES)} features (full atlas, custom panels, API)")

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
        from .license_manager import save_license, LicenseInfo
        save_license(LicenseInfo(tier=LicenseTier.FREE))
        console.print("[yellow]License deactivated. Reverted to free tier.[/yellow]")

    else:
        console.print(f"[red]Unknown action: {action}. Use: status, activate, deactivate[/red]")
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
