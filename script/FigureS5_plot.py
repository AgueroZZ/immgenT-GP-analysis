"""Figure S5 (plot step): draw the EBMF vs RQVI level2-cluster heatmaps.

Reads the raw cluster-mean matrices written by script/FigureS5.R, orders EBMF
programs by hierarchical clustering (average linkage, correlation distance,
optimal leaf ordering), scales every program to [0, 1] across clusters, and draws
two Blues heatmaps (EBMF | corresponding RQVI) sharing row order, columns, and a
"Relative loading" colorbar.

Run from the repository root:
    python script/FigureS5_plot.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage

FIG_DIR = Path("figures/generated/Figure S5")
EBMF_MEANS = FIG_DIR / "S5_ebmf_raw_means_level2.csv"
# RQVI cluster means for the right panel: the matched programs from
# FigureS5_rematch.py. Pass --rqvi-means to plot a different matrix.
RQVI_MEANS = FIG_DIR / "S5_rqvi_rematched_raw_means_level2.csv"
CLUSTER_ORDER = FIG_DIR / "S5_cluster_order.csv"
LEVEL1_PALETTE = FIG_DIR / "S5_level1_palette.csv"
SUBFIG_DIR = FIG_DIR / "S5_subfigures"

# populated from S5_level1_palette.csv in main()
LEVEL1_COLORS: dict[str, str] = {}


def _zscore_columns(df: pd.DataFrame) -> pd.DataFrame:
    values = df.to_numpy(dtype=np.float64)
    means = values.mean(axis=0, keepdims=True)
    stds = values.std(axis=0, ddof=0, keepdims=True)
    informative = stds.ravel() > np.finfo(float).eps
    if not np.all(informative):
        raise ValueError(f"constant-profile factors: {df.columns[~informative].tolist()}")
    z = (values - means) / stds
    return pd.DataFrame(z, index=df.index, columns=df.columns)


def _scale_columns_to_unit_interval(df: pd.DataFrame) -> pd.DataFrame:
    values = df.to_numpy(dtype=np.float64)
    minima = values.min(axis=0, keepdims=True)
    ranges = values.max(axis=0, keepdims=True) - minima
    informative = ranges.ravel() > np.finfo(float).eps
    scaled = np.zeros_like(values)
    scaled[:, informative] = (values[:, informative] - minima[:, informative]) / ranges[:, informative]
    return pd.DataFrame(scaled, index=df.index, columns=df.columns)


def _group_spans(lineages: list[str]) -> list[tuple[str, int, int]]:
    spans: list[tuple[str, int, int]] = []
    start = 0
    for position in range(1, len(lineages) + 1):
        if position == len(lineages) or lineages[position] != lineages[start]:
            spans.append((lineages[start], start, position))
            start = position
    return spans


def _draw_lineage_strip(ax, cluster_lineages, spans) -> None:
    categories = list(dict.fromkeys(cluster_lineages))
    category_to_code = {label: index for index, label in enumerate(categories)}
    codes = np.asarray([[category_to_code[label] for label in cluster_lineages]])
    cmap = mcolors.ListedColormap([LEVEL1_COLORS.get(label, "#BDBDBD") for label in categories])
    ax.imshow(codes, aspect="auto", cmap=cmap, interpolation="none")
    ax.set_xlim(-0.5, len(cluster_lineages) - 0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for lineage, start, stop in spans:
        if stop - start >= 4:
            ax.text((start + stop - 1) / 2, -0.6, lineage, ha="center", va="bottom",
                    fontsize=7.5, clip_on=False)
        if start > 0:
            ax.axvline(start - 0.5, color="white", linewidth=1.0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def _style_heatmap(ax, spans) -> None:
    for _, start, _ in spans[1:]:
        ax.axvline(start - 0.5, color="#777777", linewidth=0.35)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color("#333333")
        spine.set_linewidth(0.55)


def _set_plot_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 300,
        "savefig.dpi": 300,
    })


def _plot_heatmap_subfigure(matrix, cluster_lineages, ylabel, ylabel_on_right, output_pdf) -> None:
    _set_plot_style()
    cmap = plt.get_cmap("Blues")
    norm = mcolors.Normalize(vmin=0.0, vmax=1.0)
    spans = _group_spans(cluster_lineages)
    fig = plt.figure(figsize=(5.1, 7.8), facecolor="white")
    grid = fig.add_gridspec(2, 1, height_ratios=[0.18, 7.6], hspace=0.02)
    ax_strip = fig.add_subplot(grid[0, 0])
    ax_heatmap = fig.add_subplot(grid[1, 0])
    _draw_lineage_strip(ax_strip, cluster_lineages, spans)
    ax_heatmap.imshow(matrix, aspect="auto", interpolation="none", cmap=cmap, norm=norm, rasterized=True)
    _style_heatmap(ax_heatmap, spans)
    ax_heatmap.set_ylabel(ylabel, fontsize=10, labelpad=8)
    if ylabel_on_right:
        ax_heatmap.yaxis.set_label_position("right")
        fig.subplots_adjust(left=0.04, right=0.86, top=0.95, bottom=0.04)
    else:
        fig.subplots_adjust(left=0.14, right=0.96, top=0.95, bottom=0.04)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_shared_colorbar(output_pdf) -> None:
    _set_plot_style()
    cmap = plt.get_cmap("Blues")
    norm = mcolors.Normalize(vmin=0.0, vmax=1.0)
    scalar_mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig = plt.figure(figsize=(2.1, 0.48), facecolor="white")
    colorbar_axis = fig.add_axes([0.08, 0.54, 0.84, 0.27])
    colorbar = fig.colorbar(scalar_mappable, cax=colorbar_axis, orientation="horizontal", ticks=[0.0, 0.5, 1.0])
    colorbar.set_label("Relative loading", fontsize=8, labelpad=2)
    colorbar.ax.tick_params(labelsize=7, length=2, pad=1)
    colorbar.outline.set_linewidth(0.45)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf, dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def _plot(ebmf_plot, rqvi_plot, cluster_lineages, output_pdf, output_png) -> None:
    _set_plot_style()
    cmap = plt.get_cmap("Blues")
    norm = mcolors.Normalize(vmin=0.0, vmax=1.0)
    fig = plt.figure(figsize=(10.8, 7.8), facecolor="white")
    grid = fig.add_gridspec(2, 3, height_ratios=[0.18, 7.6], width_ratios=[1.0, 0.14, 1.0],
                            hspace=0.02, wspace=0.0)
    ax_strip_ebmf = fig.add_subplot(grid[0, 0])
    ax_ebmf = fig.add_subplot(grid[1, 0])
    ax_strip_rqvi = fig.add_subplot(grid[0, 2])
    ax_rqvi = fig.add_subplot(grid[1, 2])
    spans = _group_spans(cluster_lineages)
    _draw_lineage_strip(ax_strip_ebmf, cluster_lineages, spans)
    _draw_lineage_strip(ax_strip_rqvi, cluster_lineages, spans)
    ax_ebmf.imshow(ebmf_plot, aspect="auto", interpolation="none", cmap=cmap, norm=norm, rasterized=True)
    ax_rqvi.imshow(rqvi_plot, aspect="auto", interpolation="none", cmap=cmap, norm=norm, rasterized=True)
    _style_heatmap(ax_ebmf, spans)
    _style_heatmap(ax_rqvi, spans)
    ax_ebmf.set_ylabel("EBMF factors", fontsize=10, labelpad=8)
    ax_rqvi.set_ylabel("Corresponding RQVI factors", fontsize=10, labelpad=8)
    ax_rqvi.yaxis.set_label_position("right")
    scalar_mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar_axis = fig.add_axes([0.80, 0.055, 0.12, 0.016])
    colorbar = fig.colorbar(scalar_mappable, cax=colorbar_axis, orientation="horizontal", ticks=[0.0, 0.5, 1.0])
    colorbar.set_label("Relative loading", fontsize=8, labelpad=2)
    colorbar.ax.tick_params(labelsize=7, length=2, pad=1)
    colorbar.outline.set_linewidth(0.45)
    fig.subplots_adjust(left=0.08, right=0.92, top=0.95, bottom=0.09)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rqvi-means", type=Path, default=RQVI_MEANS,
                        help="RQVI cluster-mean CSV to plot on the right panel.")
    parser.add_argument("--suffix", type=str, default="",
                        help="Suffix for output filenames (e.g. '_rematched').")
    args = parser.parse_args()

    global LEVEL1_COLORS
    palette = pd.read_csv(LEVEL1_PALETTE)
    LEVEL1_COLORS = dict(zip(palette["level1"].astype(str), palette["color"].astype(str)))

    cluster_info = pd.read_csv(CLUSTER_ORDER).sort_values("display_column")
    cluster_labels = cluster_info["level2_cluster"].astype(str).tolist()
    cluster_lineages = cluster_info["level1"].astype(str).tolist()

    ebmf_raw = pd.read_csv(EBMF_MEANS, index_col="level2_cluster")
    rqvi_raw = pd.read_csv(args.rqvi_means, index_col="level2_cluster")
    ebmf_raw.index = ebmf_raw.index.astype(str)
    rqvi_raw.index = rqvi_raw.index.astype(str)
    if set(cluster_labels) != set(ebmf_raw.index) or set(cluster_labels) != set(rqvi_raw.index):
        raise ValueError("cluster labels differ between order file and mean matrices")
    ebmf_raw = ebmf_raw.loc[cluster_labels]
    rqvi_raw = rqvi_raw.loc[cluster_labels]
    factors = [f"F{k}" for k in range(1, 201)]
    if ebmf_raw.columns.tolist() != factors or rqvi_raw.columns.tolist() != factors:
        raise ValueError("expected columns F1..F200 in both matrices")

    # row order: hierarchical clustering of z-scored EBMF profiles
    ebmf_z = _zscore_columns(ebmf_raw)
    tree = linkage(ebmf_z.to_numpy().T, method="average", metric="correlation", optimal_ordering=True)
    display_order = leaves_list(tree)

    ebmf_scaled = _scale_columns_to_unit_interval(ebmf_raw)
    rqvi_scaled = _scale_columns_to_unit_interval(rqvi_raw)
    ebmf_plot = ebmf_scaled.to_numpy().T[display_order]
    rqvi_plot = rqvi_scaled.to_numpy().T[display_order]

    sfx = args.suffix
    subfig_dir = SUBFIG_DIR if not sfx else SUBFIG_DIR.with_name(SUBFIG_DIR.name + sfx)
    _plot(ebmf_plot, rqvi_plot, cluster_lineages,
          FIG_DIR / f"S5_ebmf_rqvi_level2_comparison{sfx}.pdf",
          FIG_DIR / f"S5_ebmf_rqvi_level2_comparison{sfx}.png")
    _plot_heatmap_subfigure(ebmf_plot, cluster_lineages, "EBMF factors", False,
                            subfig_dir / "panel_A_ebmf_factors.pdf")
    _plot_heatmap_subfigure(rqvi_plot, cluster_lineages, "Corresponding RQVI factors", True,
                            subfig_dir / "panel_B_corresponding_rqvi_factors.pdf")
    _plot_shared_colorbar(subfig_dir / "shared_relative_loading_colorbar.pdf")

    ordered_factors = [factors[i] for i in display_order]
    pd.DataFrame(ebmf_plot, index=ordered_factors, columns=cluster_labels).to_csv(
        FIG_DIR / f"S5_ebmf_scaled_display{sfx}.csv")
    pd.DataFrame(rqvi_plot, index=ordered_factors, columns=cluster_labels).to_csv(
        FIG_DIR / f"S5_rqvi_scaled_display{sfx}.csv")

    print(f"rows: {len(display_order)} factors; cols: {len(cluster_labels)} level2 clusters")
    print(f"wrote {FIG_DIR / ('S5_ebmf_rqvi_level2_comparison'+sfx+'.pdf')} and .png")
    print(f"wrote subfigures to {subfig_dir}")


if __name__ == "__main__":
    main()
