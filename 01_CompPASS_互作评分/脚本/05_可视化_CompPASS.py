import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


COLORS = {
    "blue": "#0072B2",
    "orange": "#D55E00",
    "green": "#009E73",
    "sky": "#56B4E9",
    "gray": "#6C757D",
    "gold": "#E69F00",
}


def first_existing(df: pd.DataFrame, names):
    for name in names:
        if name in df.columns:
            return name
    return None


def numeric_series(df: pd.DataFrame, col: str) -> pd.Series:
    return pd.to_numeric(df[col], errors="coerce")


def safe_filename(text: str) -> str:
    return str(text).replace("/", "_").replace("\\", "_").replace(":", "_")


def save_wd_histogram(df, wd_col, out_dir: Path):
    wd = numeric_series(df, wd_col).dropna()
    if wd.empty:
        return

    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.hist(wd, bins=50, color=COLORS["blue"], edgecolor="white", alpha=0.9)

    median = wd.median()
    p95 = wd.quantile(0.95)
    ax.axvline(median, color=COLORS["orange"], linewidth=2.0, linestyle="--", label=f"Median = {median:.2f}")
    ax.axvline(p95, color=COLORS["green"], linewidth=2.0, linestyle=":", label=f"P95 = {p95:.2f}")

    wd_min, wd_max = wd.min(), wd.max()
    pad = max((wd_max - wd_min) * 0.05, 0.05)
    ax.set_xlim(wd_min - pad, wd_max + pad)
    ax.set_title("CompPASS WD Score Distribution", fontsize=16, pad=16)
    ax.set_xlabel(wd_col, fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.tick_params(axis="both", labelsize=10)
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    ax.grid(axis="y", alpha=0.18)

    fig.tight_layout()
    fig.savefig(out_dir / "01_wd_histogram.png", dpi=220, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)


def save_wd_z_scatter(df, wd_col, z_col, out_dir: Path, max_points=50000):
    plot_df = df[[wd_col, z_col]].copy()
    plot_df[wd_col] = numeric_series(plot_df, wd_col)
    plot_df[z_col] = numeric_series(plot_df, z_col)
    plot_df = plot_df.dropna()
    if plot_df.empty:
        return

    if len(plot_df) > max_points:
        plot_df = plot_df.sample(max_points, random_state=7)

    fig, ax = plt.subplots(figsize=(10, 7.5))
    ax.scatter(
        plot_df[wd_col],
        plot_df[z_col],
        s=12,
        alpha=0.18,
        color=COLORS["orange"],
        edgecolors="none",
        label=f"Sampled points = {len(plot_df):,}",
    )

    wd_min, wd_max = plot_df[wd_col].min(), plot_df[wd_col].max()
    z_min, z_max = plot_df[z_col].min(), plot_df[z_col].max()
    wd_pad = max((wd_max - wd_min) * 0.05, 0.05)
    z_pad = max((z_max - z_min) * 0.05, 0.05)
    ax.set_xlim(wd_min - wd_pad, wd_max + wd_pad)
    ax.set_ylim(z_min - z_pad, z_max + z_pad)

    ax.axhline(plot_df[z_col].quantile(0.95), color=COLORS["green"], linestyle=":", linewidth=1.8, alpha=0.9)
    ax.axvline(plot_df[wd_col].quantile(0.95), color=COLORS["sky"], linestyle=":", linewidth=1.8, alpha=0.9)
    ax.set_title("CompPASS WD vs Z", fontsize=16, pad=16)
    ax.set_xlabel(wd_col, fontsize=12)
    ax.set_ylabel(z_col, fontsize=12)
    ax.tick_params(axis="both", labelsize=10)
    ax.grid(alpha=0.16)
    ax.legend(frameon=False, fontsize=10, loc="upper left")

    fig.tight_layout()
    fig.savefig(out_dir / "02_wd_z_scatter.png", dpi=220, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)


def save_top_prey_bars(df, bait_col, prey_col, wd_col, out_dir: Path, top_baits=8, top_prey=10):
    bait_order = (
        df.groupby(bait_col)[wd_col]
        .max()
        .sort_values(ascending=False)
        .head(top_baits)
        .index.tolist()
    )
    for bait in bait_order:
        subset = df.loc[df[bait_col] == bait].copy()
        subset[wd_col] = numeric_series(subset, wd_col)
        subset = subset.dropna(subset=[wd_col]).sort_values(wd_col, ascending=False).head(top_prey)
        if subset.empty:
            continue

        fig, ax = plt.subplots(figsize=(11, 7))
        bars = ax.barh(
            subset[prey_col].astype(str),
            subset[wd_col],
            color=COLORS["green"],
            alpha=0.85,
            edgecolor="white",
        )
        ax.invert_yaxis()
        ax.set_title(f"Top prey for {bait}", fontsize=15, pad=14)
        ax.set_xlabel(wd_col, fontsize=12)
        ax.set_ylabel(prey_col, fontsize=12)
        ax.tick_params(axis="x", labelsize=10)
        ax.tick_params(axis="y", labelsize=9)
        ax.grid(axis="x", alpha=0.15)

        max_value = subset[wd_col].max()
        for bar in bars:
            value = bar.get_width()
            ax.text(
                value + max_value * 0.015,
                bar.get_y() + bar.get_height() / 2,
                f"{value:.2f}",
                va="center",
                fontsize=8,
                color="#333333",
            )

        fig.tight_layout()
        fig.savefig(out_dir / f"03_top_prey_{safe_filename(bait)}.png", dpi=220, bbox_inches="tight", pad_inches=0.35)
        plt.close(fig)


def save_heatmap(df, bait_col, prey_col, wd_col, out_dir: Path, max_baits=12, max_prey=25):
    bait_order = df.groupby(bait_col)[wd_col].max().sort_values(ascending=False).head(max_baits).index
    prey_order = df.groupby(prey_col)[wd_col].max().sort_values(ascending=False).head(max_prey).index
    plot_df = df.loc[df[bait_col].isin(bait_order) & df[prey_col].isin(prey_order)].copy()
    if plot_df.empty:
        return

    matrix = plot_df.pivot_table(index=prey_col, columns=bait_col, values=wd_col, aggfunc="max", fill_value=0)
    fig, ax = plt.subplots(figsize=(13, 10))
    im = ax.imshow(matrix.values, aspect="auto", cmap="cividis")
    ax.set_title("Top bait-prey WD heatmap", fontsize=16, pad=16)
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=75, ha="right", fontsize=9)
    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=8)
    colorbar = fig.colorbar(im, ax=ax, label=wd_col, fraction=0.035, pad=0.02)
    colorbar.ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "04_top_bait_prey_wd_heatmap.png", dpi=220, bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)


def make_bait_prey_table(df, bait_col, prey_col, wd_col, prey_label_col=None, top_baits=18, top_prey=10):
    cols = [bait_col, prey_col, wd_col]
    if prey_label_col is not None and prey_label_col not in cols:
        cols.append(prey_label_col)
    plot_df = df[cols].copy()
    plot_df[wd_col] = numeric_series(plot_df, wd_col)
    plot_df = plot_df.dropna(subset=[wd_col])
    if plot_df.empty:
        return plot_df

    bait_order = (
        plot_df.groupby(bait_col)[wd_col]
        .max()
        .sort_values(ascending=False)
        .head(top_baits)
        .index.tolist()
    )

    parts = []
    for bait in bait_order:
        subset = plot_df.loc[plot_df[bait_col] == bait].sort_values(wd_col, ascending=False).head(top_prey)
        parts.append(subset)

    if not parts:
        return plot_df.iloc[0:0].copy()
    return pd.concat(parts, ignore_index=True)


def build_bipartite_graph(plot_df, bait_col, prey_col, wd_col, prey_label_col=None):
    graph = nx.Graph()
    for _, row in plot_df.iterrows():
        bait = str(row[bait_col])
        prey = str(row[prey_col])
        weight = float(row[wd_col])
        prey_label = str(row[prey_label_col]) if prey_label_col is not None and pd.notna(row[prey_label_col]) and str(row[prey_label_col]).strip() else prey
        graph.add_node(bait, node_type="bait", label=bait)
        graph.add_node(prey, node_type="prey", label=prey_label)
        graph.add_edge(bait, prey, weight=weight)
    return graph


def node_strength(graph, node):
    return sum(graph[node][nbr].get("weight", 0.0) for nbr in graph.neighbors(node))


def draw_overview_network(graph, out_path: Path, title: str):
    if graph.number_of_nodes() == 0 or graph.number_of_edges() == 0:
        return

    bait_nodes = [node for node, data in graph.nodes(data=True) if data.get("node_type") == "bait"]
    prey_nodes = [node for node, data in graph.nodes(data=True) if data.get("node_type") == "prey"]
    bait_count = max(len(bait_nodes), 1)

    fixed_positions = {}
    radius = 7.0 + 0.18 * bait_count
    for idx, bait in enumerate(sorted(bait_nodes)):
        angle = 2 * np.pi * idx / bait_count
        fixed_positions[bait] = np.array([radius * np.cos(angle), radius * np.sin(angle)])

    init_positions = dict(fixed_positions)
    for prey in prey_nodes:
        neighbors = list(graph.neighbors(prey))
        if neighbors:
            anchor = np.mean([fixed_positions[n] for n in neighbors if n in fixed_positions], axis=0)
        else:
            anchor = np.array([0.0, 0.0])
        jitter_seed = sum(ord(ch) for ch in prey)
        jitter = np.array([
            ((jitter_seed % 11) - 5) * 0.12,
            (((jitter_seed // 11) % 11) - 5) * 0.12,
        ])
        init_positions[prey] = anchor + jitter

    pos = nx.spring_layout(
        graph,
        seed=7,
        pos=init_positions,
        fixed=bait_nodes,
        k=max(0.6, 2.8 / max(graph.number_of_nodes(), 4) ** 0.5),
        iterations=350,
        weight="weight",
    )

    edge_weights = [graph[u][v].get("weight", 0.0) for u, v in graph.edges()]
    max_weight = max(edge_weights) if edge_weights else 1.0
    min_weight = min(edge_weights) if edge_weights else 0.0

    edge_widths = [1.0 + 4.0 * (weight / max_weight) for weight in edge_weights]
    edge_alphas = [0.2 + 0.55 * (weight / max_weight) for weight in edge_weights]

    bait_sizes = [420 + 70 * graph.degree(node) + 18 * node_strength(graph, node) / max_weight for node in bait_nodes]
    prey_sizes = [110 + 45 * max(graph.degree(node) - 1, 0) + 12 * node_strength(graph, node) / max_weight for node in prey_nodes]

    fig, ax = plt.subplots(figsize=(20, 16))
    for (u, v), width, alpha in zip(graph.edges(), edge_widths, edge_alphas):
        nx.draw_networkx_edges(
            graph,
            pos,
            edgelist=[(u, v)],
            width=width,
            alpha=alpha,
            edge_color=COLORS["gray"],
            ax=ax,
        )

    nx.draw_networkx_nodes(
        graph,
        pos,
        nodelist=prey_nodes,
        node_color=COLORS["sky"],
        node_size=prey_sizes,
        alpha=0.88,
        linewidths=0.8,
        edgecolors="white",
        ax=ax,
    )
    nx.draw_networkx_nodes(
        graph,
        pos,
        nodelist=bait_nodes,
        node_color=COLORS["orange"],
        node_size=bait_sizes,
        alpha=0.95,
        linewidths=0.8,
        edgecolors="white",
        ax=ax,
    )

    shared_prey = [node for node in prey_nodes if graph.degree(node) > 1]
    ranked_shared = sorted(shared_prey, key=lambda node: (graph.degree(node), node_strength(graph, node)), reverse=True)[:28]
    label_nodes = bait_nodes + ranked_shared
    label_map = {node: graph.nodes[node].get("label", node) for node in label_nodes}
    nx.draw_networkx_labels(
        graph,
        pos,
        labels=label_map,
        font_size=10,
        font_family="sans-serif",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.78, pad=0.2),
        ax=ax,
    )

    legend_handles = [
        Line2D([0], [0], marker="o", color="w", label="Bait", markerfacecolor=COLORS["orange"], markersize=12, alpha=0.95),
        Line2D([0], [0], marker="o", color="w", label="Prey", markerfacecolor=COLORS["sky"], markersize=10, alpha=0.88),
        Line2D([0], [0], marker="o", color="w", label="Shared prey (labeled)", markerfacecolor=COLORS["sky"], markersize=10, alpha=0.88),
        Line2D([0], [0], color=COLORS["gray"], lw=3, alpha=0.6, label=f"Edge WD range: {min_weight:.2f} to {max_weight:.2f}"),
    ]

    ax.legend(handles=legend_handles, frameon=False, fontsize=11, loc="upper left", bbox_to_anchor=(0.01, 0.99))
    ax.set_title(
        f"{title}\nBaits={len(bait_nodes)}  Preys={len(prey_nodes)}  Edges={graph.number_of_edges()}  Shared preys={len(shared_prey)}",
        fontsize=18,
        pad=20,
    )
    ax.margins(0.18)
    ax.set_axis_off()
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=240, bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)


def draw_focus_grid(plot_df, bait_col, prey_col, wd_col, out_path: Path, prey_label_col=None, top_baits=8, top_prey=10):
    if plot_df.empty:
        return

    bait_order = (
        plot_df.groupby(bait_col)[wd_col]
        .max()
        .sort_values(ascending=False)
        .head(top_baits)
        .index.tolist()
    )
    if not bait_order:
        return

    ncols = 2
    nrows = int(np.ceil(len(bait_order) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(18, max(8, 5.7 * nrows)))
    axes = np.atleast_1d(axes).ravel()

    for ax, bait in zip(axes, bait_order):
        subset = plot_df.loc[plot_df[bait_col] == bait].sort_values(wd_col, ascending=False).head(top_prey).copy()
        subset[wd_col] = numeric_series(subset, wd_col)
        subset = subset.dropna(subset=[wd_col])
        if subset.empty:
            ax.axis("off")
            continue

        prey_nodes = subset[prey_col].astype(str).tolist()
        prey_labels = {}
        if prey_label_col is not None and prey_label_col in subset.columns:
            prey_labels = dict(zip(subset[prey_col].astype(str), subset[prey_label_col].fillna("").astype(str)))

        max_wd = subset[wd_col].max()
        angles = np.linspace(0, 2 * np.pi, len(prey_nodes), endpoint=False)
        center = np.array([0.0, 0.0])
        radius = 2.6 + 0.11 * len(prey_nodes)

        ax.scatter([center[0]], [center[1]], s=1200, color=COLORS["orange"], alpha=0.96, edgecolors="white", linewidths=1.2, zorder=4)
        ax.text(center[0], center[1], str(bait), ha="center", va="center", fontsize=10, color="white", weight="bold", zorder=5)

        for angle, prey, (_, row) in zip(angles, prey_nodes, subset.iterrows()):
            point = np.array([radius * np.cos(angle), radius * np.sin(angle)])
            value = float(row[wd_col])
            size = 180 + 700 * (value / max_wd)
            width = 1.0 + 4.0 * (value / max_wd)
            alpha = 0.25 + 0.55 * (value / max_wd)
            ax.plot([center[0], point[0]], [center[1], point[1]], color=COLORS["gray"], linewidth=width, alpha=alpha, zorder=1)
            ax.scatter([point[0]], [point[1]], s=size, color=COLORS["sky"], alpha=0.9, edgecolors="white", linewidths=0.8, zorder=3)

            label = prey_labels.get(prey, "") or prey
            label = label if label and label != "nan" else prey
            align = "left" if point[0] >= 0 else "right"
            offset = 0.22 if point[0] >= 0 else -0.22
            ax.text(
                point[0] + offset,
                point[1],
                f"{label}\nWD={value:.1f}",
                ha=align,
                va="center",
                fontsize=8,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=0.18),
                zorder=6,
            )

        lim = radius + 1.5
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_title(f"{bait}: bait-centered high-confidence prey", fontsize=12, pad=10)
        ax.axis("off")

    for ax in axes[len(bait_order):]:
        ax.axis("off")

    legend_handles = [
        Line2D([0], [0], marker="o", color="w", label="Bait", markerfacecolor=COLORS["orange"], markersize=12, alpha=0.96),
        Line2D([0], [0], marker="o", color="w", label="Prey", markerfacecolor=COLORS["sky"], markersize=10, alpha=0.9),
        Line2D([0], [0], color=COLORS["gray"], lw=3, alpha=0.5, label="Edge width ~ WD"),
    ]
    fig.legend(handles=legend_handles, frameon=False, fontsize=11, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.995))
    fig.suptitle("Bait-centered interaction maps", fontsize=18, y=1.02)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=240, bbox_inches="tight", pad_inches=0.35)
    plt.close(fig)


def save_ppi_global_network(df, bait_col, prey_col, wd_col, out_dir: Path, prey_label_col=None):
    plot_df = make_bait_prey_table(df, bait_col, prey_col, wd_col, prey_label_col=prey_label_col, top_baits=18, top_prey=10)
    if plot_df.empty:
        return
    graph = build_bipartite_graph(plot_df, bait_col, prey_col, wd_col, prey_label_col=prey_label_col)
    draw_overview_network(graph, out_dir / "05_ppi_global_network.png", "Bait-prey overview network")


def save_ppi_focus_grid(df, bait_col, prey_col, wd_col, out_dir: Path, prey_label_col=None):
    plot_df = make_bait_prey_table(df, bait_col, prey_col, wd_col, prey_label_col=prey_label_col, top_baits=8, top_prey=10)
    draw_focus_grid(plot_df, bait_col, prey_col, wd_col, out_dir / "06_ppi_core_module.png", prey_label_col=prey_label_col, top_baits=8, top_prey=10)


def main():
    parser = argparse.ArgumentParser(description="Create improved CompPASS visualizations.")
    parser.add_argument("--input", default="6.CompPASS评分结果/原始结果/comppass_results.tsv")
    parser.add_argument("--output-dir", default="7.可视化/图形输出")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    bait_col = first_existing(df, ["idBait", "bait", "Bait"])
    prey_col = first_existing(df, ["idPrey", "prey", "Prey"])
    wd_col = first_existing(df, ["scoreWD", "WD", "wd", "NWD", "nwd"])
    z_col = first_existing(df, ["scoreZ", "Z", "z", "ZScore", "zscore"])
    prey_label_col = first_existing(df, ["preyGene", "Gene", "gene"])

    if bait_col is None or prey_col is None or wd_col is None:
        raise ValueError("结果表至少需要 bait、prey 和 WD/NWD 分数字段。")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    save_wd_histogram(df, wd_col, out_dir)
    if z_col is not None:
        save_wd_z_scatter(df, wd_col, z_col, out_dir)
    save_top_prey_bars(df, bait_col, prey_col, wd_col, out_dir)
    save_heatmap(df, bait_col, prey_col, wd_col, out_dir)
    save_ppi_global_network(df, bait_col, prey_col, wd_col, out_dir, prey_label_col=prey_label_col)
    save_ppi_focus_grid(df, bait_col, prey_col, wd_col, out_dir, prey_label_col=prey_label_col)

    print(f"Figures written to: {out_dir}")


if __name__ == "__main__":
    main()
