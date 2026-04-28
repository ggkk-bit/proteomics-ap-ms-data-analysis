#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PPIrank 可视化脚本。

支持的通用图形：
- 排名分布图
- Top prey 条形图
- bait 局部网络图
- 方法间重叠矩阵图

文献风格：
- bait 节点突出显示
- 网络图优先服务于“补充排序”结果解释
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Set, Tuple


def _sniff_delimiter(path: Path) -> str:
    sample = path.read_text(encoding="utf-8", errors="ignore")[:4096]
    return "\t" if "\t" in sample else ","


def _read_rows(path: Path) -> List[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=_sniff_delimiter(path))
        return list(reader)


def plot_rank_distribution(rows: Sequence[dict], score_col: str, output: Path) -> None:
    import matplotlib.pyplot as plt

    scores = [float(r[score_col]) for r in rows if str(r.get(score_col, "")).strip()]
    plt.figure(figsize=(7, 4))
    plt.hist(scores, bins=30, color="#2f6fed", edgecolor="white")
    plt.xlabel("PPIRank score")
    plt.ylabel("频数")
    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output, dpi=300)
    plt.close()


def plot_top_pairs(rows: Sequence[dict], prey_col: str, score_col: str, top_n: int, output: Path) -> None:
    import matplotlib.pyplot as plt

    data = sorted(
        ((r.get(prey_col, ""), float(r.get(score_col, 0) or 0)) for r in rows),
        key=lambda x: x[1],
        reverse=True,
    )[:top_n]
    labels = [x[0] for x in data][::-1]
    values = [x[1] for x in data][::-1]
    plt.figure(figsize=(8, max(4, top_n * 0.35)))
    plt.barh(labels, values, color="#3c9d5c")
    plt.xlabel("PPIRank score")
    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output, dpi=300)
    plt.close()


def plot_network(rows: Sequence[dict], bait_col: str, prey_col: str, score_col: str, output: Path, top_n: int = 50) -> None:
    import matplotlib.pyplot as plt
    import networkx as nx

    selected = sorted(rows, key=lambda r: float(r.get(score_col, 0) or 0), reverse=True)[:top_n]
    graph = nx.Graph()
    for row in selected:
        bait = str(row.get(bait_col, "")).strip()
        prey = str(row.get(prey_col, "")).strip()
        score = float(row.get(score_col, 0) or 0)
        if not bait or not prey:
            continue
        graph.add_edge(bait, prey, weight=score)

    if not graph.nodes:
        raise ValueError("没有可用于绘图的网络边")

    pos = nx.spring_layout(graph, seed=42)
    bait_nodes = {str(row.get(bait_col, "")).strip() for row in selected}
    prey_nodes = set(graph.nodes) - bait_nodes
    plt.figure(figsize=(10, 8))
    nx.draw_networkx_nodes(graph, pos, nodelist=list(prey_nodes), node_color="#b9c3d1", node_size=350)
    nx.draw_networkx_nodes(graph, pos, nodelist=list(bait_nodes), node_color="#3c9d5c", node_size=500)
    nx.draw_networkx_edges(graph, pos, width=1.0, alpha=0.5)
    nx.draw_networkx_labels(graph, pos, font_size=8)
    plt.axis("off")
    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output, dpi=300)
    plt.close()


def plot_overlap_matrix(method_to_edges: Dict[str, Set[Tuple[str, str]]], output: Path) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    methods = list(method_to_edges)
    matrix = np.zeros((len(methods), len(methods)), dtype=int)
    for i, left in enumerate(methods):
        for j, right in enumerate(methods):
            matrix[i, j] = len(method_to_edges[left] & method_to_edges[right])

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(matrix, cmap="Blues")
    ax.set_xticks(range(len(methods)))
    ax.set_yticks(range(len(methods)))
    ax.set_xticklabels(methods, rotation=45, ha="right")
    ax.set_yticklabels(methods)
    for i in range(len(methods)):
        for j in range(len(methods)):
            ax.text(j, i, matrix[i, j], ha="center", va="center", color="black", fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.8)
    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output, dpi=300)
    plt.close()


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="PPIrank 常用可视化")
    parser.add_argument("--input", required=True, help="PPIrank 结果 TSV")
    parser.add_argument("--output-dir", required=True, help="图片输出目录")
    parser.add_argument("--bait-col", default="bait_name")
    parser.add_argument("--prey-col", default="prey_id")
    parser.add_argument("--score-col", default="ppirank_score")
    parser.add_argument("--top-n", type=int, default=50)
    args = parser.parse_args()

    rows = _read_rows(Path(args.input))
    output_dir = Path(args.output_dir)
    plot_rank_distribution(rows, args.score_col, output_dir / "01_分数分布图.png")
    plot_top_pairs(rows, args.prey_col, args.score_col, min(args.top_n, len(rows)), output_dir / "02_Top互作条形图.png")
    plot_network(rows, args.bait_col, args.prey_col, args.score_col, output_dir / "03_网络图.png", top_n=args.top_n)


if __name__ == "__main__":
    main()
