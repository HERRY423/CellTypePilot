"""Visualization layer — UMAP, dotplot, confidence plots."""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt

from .constants import CB_PALETTE, OUTPUT_FIGURES_DIR


def generate_all_figures(
    adata,
    cluster_key: str,
    embedding_key: str,
    annotations: pd.DataFrame,
    output_dir: str | Path,
    tissue: str = "",
) -> list[str]:
    """Generate all standard figures.

    Returns list of generated file paths.
    """
    output_dir = Path(output_dir)
    figures_dir = output_dir / OUTPUT_FIGURES_DIR
    figures_dir.mkdir(parents=True, exist_ok=True)

    generated = []

    # 1. UMAP colored by cluster
    path = plot_umap_clusters(adata, cluster_key, embedding_key, figures_dir)
    if path:
        generated.append(path)

    # 2. UMAP colored by cell type annotation
    if "cell_type" in annotations.columns:
        path = plot_umap_celltype(adata, annotations, cluster_key, embedding_key, figures_dir)
        if path:
            generated.append(path)

    # 3. Confidence visualization
    if "critic_confidence" in annotations.columns:
        path = plot_umap_confidence(adata, annotations, cluster_key, embedding_key, figures_dir)
        if path:
            generated.append(path)

    # 4. Marker dotplot
    path = plot_marker_dotplot(adata, annotations, cluster_key, figures_dir)
    if path:
        generated.append(path)

    # 5. Confidence bar chart
    if "critic_confidence" in annotations.columns:
        path = plot_confidence_bar(annotations, figures_dir)
        if path:
            generated.append(path)

    # 6. UMAP colored by Identity x State composite label
    if "display_label" in annotations.columns:
        path = plot_umap_identity_state(adata, annotations, cluster_key, embedding_key, figures_dir)
        if path:
            generated.append(path)

    # 7. Identity x State distribution heatmap
    if "cell_state_candidate" in annotations.columns:
        path = plot_identity_state_distribution(annotations, figures_dir)
        if path:
            generated.append(path)

    return generated


def plot_umap_clusters(adata, cluster_key: str, embedding_key: str, output_dir: Path) -> str | None:
    """UMAP colored by cluster ID."""
    if embedding_key not in adata.obsm:
        return None

    coords = adata.obsm[embedding_key]
    clusters = adata.obs[cluster_key]
    unique_clusters = sorted(clusters.unique(), key=str)
    n_clusters = len(unique_clusters)

    colors = _generate_colors(n_clusters)
    cluster_color_map = {str(c): colors[i % len(colors)] for i, c in enumerate(unique_clusters)}

    fig, ax = plt.subplots(figsize=(8, 6))
    for cluster in unique_clusters:
        mask = clusters == cluster
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=[cluster_color_map[str(cluster)]],
            label=str(cluster),
            s=1,
            alpha=0.6,
            rasterized=True,
        )

    ax.set_title(f"Clusters ({cluster_key})", fontsize=12)
    ax.set_xlabel(f"{embedding_key}1")
    ax.set_ylabel(f"{embedding_key}2")
    ax.legend(
        markerscale=4,
        fontsize=7,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        frameon=False,
    )
    ax.set_aspect("equal")

    path = output_dir / "umap_clusters.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(path)


def plot_umap_celltype(
    adata, annotations: pd.DataFrame, cluster_key: str, embedding_key: str, output_dir: Path
) -> str | None:
    """UMAP colored by cell type annotation."""
    if embedding_key not in adata.obsm:
        return None

    coords = adata.obsm[embedding_key]

    # Map cluster → cell_type
    cluster_to_ct = dict(zip(annotations["cluster"], annotations["cell_type"], strict=True))
    cell_types = [cluster_to_ct.get(str(c), "Unknown") for c in adata.obs[cluster_key]]

    unique_ct = sorted(set(cell_types))
    colors = _generate_colors(len(unique_ct))
    ct_color_map = {ct: colors[i % len(colors)] for i, ct in enumerate(unique_ct)}

    fig, ax = plt.subplots(figsize=(8, 6))
    for ct in unique_ct:
        mask = np.array([c == ct for c in cell_types])
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=[ct_color_map[ct]],
            label=ct,
            s=1,
            alpha=0.6,
            rasterized=True,
        )

    ax.set_title("Cell Type Annotation", fontsize=12)
    ax.set_xlabel(f"{embedding_key}1")
    ax.set_ylabel(f"{embedding_key}2")
    ax.legend(
        markerscale=4,
        fontsize=7,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        frameon=False,
    )
    ax.set_aspect("equal")

    path = output_dir / "umap_celltype.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(path)


def plot_umap_confidence(
    adata, annotations: pd.DataFrame, cluster_key: str, embedding_key: str, output_dir: Path
) -> str | None:
    """UMAP colored by critic confidence level."""
    if embedding_key not in adata.obsm:
        return None

    coords = adata.obsm[embedding_key]

    cluster_to_conf = dict(
        zip(
            annotations["cluster"],
            annotations.get("critic_confidence", ["unknown"] * len(annotations)),
            strict=True,
        )
    )
    conf_levels = [cluster_to_conf.get(str(c), "unknown") for c in adata.obs[cluster_key]]

    conf_colors = {
        "high": "#009E73",
        "medium": "#E69F00",
        "low": "#D55E00",
        "needs_review": "#CC0000",
        "unknown": "#999999",
    }

    fig, ax = plt.subplots(figsize=(8, 6))
    for level, color in conf_colors.items():
        mask = np.array([c == level for c in conf_levels])
        if mask.any():
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                c=[color],
                label=level,
                s=2,
                alpha=0.6,
                rasterized=True,
            )

    ax.set_title("Annotation Confidence", fontsize=12)
    ax.set_xlabel(f"{embedding_key}1")
    ax.set_ylabel(f"{embedding_key}2")
    ax.legend(fontsize=8, frameon=False)
    ax.set_aspect("equal")

    path = output_dir / "umap_confidence.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(path)


def plot_marker_dotplot(
    adata, annotations: pd.DataFrame, cluster_key: str, output_dir: Path
) -> str | None:
    """Marker gene dotplot: cell types × key markers."""
    # Collect top markers per cell type
    ct_markers = {}
    for _, row in annotations.iterrows():
        ct = row.get("cell_type", "")
        # We'll use a fixed set of common markers for the dotplot
        ct_markers[ct] = []

    # Get all unique markers from the atlas data embedded in annotations
    if "n_pos_markers" not in annotations.columns:
        return None

    # Use scanpy's dotplot if available
    try:
        import scanpy as sc

        # Collect markers from the atlas
        from .data_adapter import get_all_markers_for_tissue, load_marker_atlas

        atlas = load_marker_atlas()

        # Get tissue from annotations or use general
        tissue = "general"
        markers_dict = get_all_markers_for_tissue(atlas, tissue)

        # Select top markers for annotated cell types
        selected_markers = {}
        for _, row in annotations.iterrows():
            ct = row.get("cell_type", "")
            if ct in markers_dict:
                pos = markers_dict[ct].get("positive_markers", [])[:5]  # top 5
                selected_markers[ct] = pos

        if not selected_markers:
            return None

        # Flatten and deduplicate markers
        all_marker_genes = []
        for genes in selected_markers.values():
            for g in genes:
                if g in adata.var_names and g not in all_marker_genes:
                    all_marker_genes.append(g)

        if not all_marker_genes:
            return None

        # Create dotplot
        dp = sc.pl.dotplot(
            adata,
            var_names=all_marker_genes,
            groupby=cluster_key,
            show=False,
            ax=None,
            return_fig=True,
        )
        fig = dp.figure
        fig.set_size_inches(
            max(10, len(all_marker_genes) * 0.5), max(6, adata.obs[cluster_key].nunique() * 0.4)
        )

        path = output_dir / "marker_dotplot.png"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return str(path)

    except Exception:
        return None


def plot_confidence_bar(annotations: pd.DataFrame, output_dir: Path) -> str | None:
    """Bar chart showing confidence distribution across clusters."""
    if "critic_confidence" not in annotations.columns:
        return None

    conf_counts = annotations["critic_confidence"].value_counts()

    conf_colors = {
        "high": "#009E73",
        "medium": "#E69F00",
        "low": "#D55E00",
        "needs_review": "#CC0000",
    }

    fig, ax = plt.subplots(figsize=(6, 4))
    colors = [conf_colors.get(c, "#999999") for c in conf_counts.index]
    bars = ax.bar(conf_counts.index, conf_counts.values, color=colors, edgecolor="white")

    ax.set_xlabel("Confidence Level")
    ax.set_ylabel("Number of Clusters")
    ax.set_title("Annotation Confidence Distribution")

    for bar, val in zip(bars, conf_counts.values, strict=True):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.1,
            str(val),
            ha="center",
            va="bottom",
            fontsize=10,
        )

    path = output_dir / "confidence_distribution.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(path)


def _generate_colors(n: int) -> list[str]:
    """Generate n distinct colors using the colorblind-friendly palette."""
    if n <= len(CB_PALETTE):
        return CB_PALETTE[:n]

    # Extend with additional distinct colors
    extra_colors = [
        "#8B4513",
        "#20B2AA",
        "#9370DB",
        "#FF6347",
        "#4682B4",
        "#32CD32",
        "#FF8C00",
        "#8A2BE2",
        "#00CED1",
        "#DC143C",
        "#228B22",
        "#FF1493",
        "#1E90FF",
        "#B22222",
        "#3CB371",
    ]
    all_colors = CB_PALETTE + extra_colors
    return [all_colors[i % len(all_colors)] for i in range(n)]


def plot_umap_identity_state(
    adata, annotations: pd.DataFrame, cluster_key: str, embedding_key: str, output_dir: Path
) -> str | None:
    """UMAP colored by Identity x State composite display label."""
    if embedding_key not in adata.obsm or "display_label" not in annotations.columns:
        return None

    coords = adata.obsm[embedding_key]
    cluster_to_label = dict(
        zip(
            annotations["cluster"].astype(str),
            annotations["display_label"].astype(str),
            strict=True,
        )
    )
    labels = [cluster_to_label.get(str(c), "Unknown") for c in adata.obs[cluster_key]]

    unique_labels = sorted(set(labels))
    colors = _generate_colors(len(unique_labels))
    label_color_map = {lbl: colors[i % len(colors)] for i, lbl in enumerate(unique_labels)}

    fig, ax = plt.subplots(figsize=(8, 6))
    for lbl in unique_labels:
        mask = np.array([label == lbl for label in labels])
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=[label_color_map[lbl]],
            label=lbl,
            s=1,
            alpha=0.6,
            rasterized=True,
        )

    ax.set_title("Identity x Cell State Composite Annotation", fontsize=12)
    ax.set_xlabel(f"{embedding_key}1")
    ax.set_ylabel(f"{embedding_key}2")
    ax.legend(
        markerscale=4,
        fontsize=7,
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
        frameon=False,
    )
    ax.set_aspect("equal")

    path = output_dir / "umap_identity_state.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(path)


def plot_identity_state_distribution(annotations: pd.DataFrame, output_dir: Path) -> str | None:
    """Heatmap matrix showing Cell State distribution per Identity."""
    if "cell_type" not in annotations.columns or "cell_state_candidate" not in annotations.columns:
        return None

    valid = annotations[annotations["cell_state_candidate"] != "Unknown"]
    if valid.empty:
        return None

    ct_state_counts = pd.crosstab(valid["cell_type"], valid["cell_state_candidate"])
    if ct_state_counts.empty:
        return None

    n_rows, n_cols = ct_state_counts.shape
    fig, ax = plt.subplots(figsize=(max(6, n_cols * 1.2), max(4, n_rows * 0.6)))
    cax = ax.matshow(ct_state_counts.values, cmap="YlGnBu", alpha=0.85)

    ax.set_xticks(range(n_cols))
    ax.set_yticks(range(n_rows))
    ax.set_xticklabels(ct_state_counts.columns, rotation=45, ha="left", fontsize=9)
    ax.set_yticklabels(ct_state_counts.index, fontsize=9)

    fig.colorbar(cax, ax=ax, fraction=0.046, pad=0.04)

    for i in range(n_rows):
        for j in range(n_cols):
            val = ct_state_counts.values[i, j]
            if val > 0:
                ax.text(
                    j,
                    i,
                    str(val),
                    ha="center",
                    va="center",
                    color="black",
                    fontsize=9,
                    fontweight="bold",
                )

    ax.set_title("Cell State Distribution across Identities", fontsize=11, pad=40)
    fig.tight_layout()

    path = output_dir / "identity_state_distribution.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(path)
