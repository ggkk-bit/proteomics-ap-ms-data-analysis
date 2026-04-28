#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""CS Score 可视化脚本。"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def save_figure(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="绘制 CS Score 可视化图")
    parser.add_argument("--workspace", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    workspace = args.workspace
    input_file = args.input or workspace / "6.CS_Score评分结果" / "CS_Score_筛选后候选.tsv"
    output_dir = args.output_dir or workspace / "7.可视化"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_file.exists():
        raise FileNotFoundError(f"找不到输入文件: {input_file}")

    df = pd.read_csv(input_file, sep="\t", dtype=str).fillna("")
    for col in ["observed_cooccurrence", "cs_score", "score", "z_score"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    score_col = "cs_score" if "cs_score" in df.columns else ("score" if "score" in df.columns else ("z_score" if "z_score" in df.columns else None))
    cooc_col = "observed_cooccurrence" if "observed_cooccurrence" in df.columns else None

    sns.set_theme(style="whitegrid")

    if score_col:
        plt.figure(figsize=(8, 5))
        sns.histplot(df[score_col].dropna(), bins=40, kde=True, color="#1f77b4")
        plt.xlabel("CS score")
        plt.ylabel("频数")
        plt.title("CS Score 分布")
        save_figure(output_dir / "01_CS_Score分布.png")

    if score_col and cooc_col:
        plt.figure(figsize=(8, 5))
        sns.scatterplot(data=df, x=cooc_col, y=score_col, s=30, color="#d62728")
        plt.xlabel("观察共现次数")
        plt.ylabel("CS score")
        plt.title("共现次数与 CS Score 关系")
        save_figure(output_dir / "02_共现次数_vs_CS_Score.png")

    if score_col and {"protein_a", "protein_b"}.issubset(df.columns):
        top = df.sort_values(score_col, ascending=False).head(20).copy()
        top["pair"] = top["protein_a"].astype(str) + " - " + top["protein_b"].astype(str)
        plt.figure(figsize=(10, 7))
        sns.barplot(data=top, y="pair", x=score_col, color="#2ca02c")
        plt.xlabel("CS score")
        plt.ylabel("蛋白对")
        plt.title("Top 20 候选蛋白对")
        save_figure(output_dir / "03_Top20候选.png")

    if cooc_col:
        plt.figure(figsize=(8, 5))
        sns.histplot(df[cooc_col].dropna(), bins=20, color="#9467bd")
        plt.xlabel("观察共现次数")
        plt.ylabel("频数")
        plt.title("共现次数分布")
        save_figure(output_dir / "04_共现次数分布.png")

    if score_col and {"protein_a", "protein_b"}.issubset(df.columns):
        top = df.sort_values(score_col, ascending=False).head(25)
        proteins = sorted(set(top["protein_a"].astype(str)).union(set(top["protein_b"].astype(str))))
        heat = pd.DataFrame(0.0, index=proteins, columns=proteins)
        for _, row in top.iterrows():
            a = str(row["protein_a"])
            b = str(row["protein_b"])
            value = float(row[score_col]) if pd.notna(row[score_col]) else 0.0
            heat.loc[a, b] = value
            heat.loc[b, a] = value
        plt.figure(figsize=(10, 8))
        sns.heatmap(heat, cmap="viridis")
        plt.title("Top 候选蛋白对热图")
        save_figure(output_dir / "05_Top候选热图.png")

    note = output_dir / "可视化说明.md"
    note.write_text(
        "# CS Score 可视化说明\n\n"
        "- 分数分布图：看整体是否偏正、是否长尾。\n"
        "- 共现次数 vs 分数散点图：看高分候选是否同时有足够的共现支撑。\n"
        "- Top 候选条形图：汇报最强的蛋白对。\n"
        "- 共现次数分布图：检查阈值是否过松或过紧。\n"
        "- Top 候选热图：看局部模块是否成团。\n",
        encoding="utf-8",
    )
    print(f"已输出: {note}")


if __name__ == "__main__":
    main()

