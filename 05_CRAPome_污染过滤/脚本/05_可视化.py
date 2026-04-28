"""CRAPome 通用可视化脚本。

可复用图形思路：
- 背景频率分布
- bait / control 谱数对比散点
- top 候选条形图
- 重复一致性热图

输入建议是 CRAPome 输出表或本地筛选汇总表，只要列名足够即可。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def _pick_column(df: pd.DataFrame, candidates: list[str]) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"找不到候选列: {', '.join(candidates)}")


def plot_frequency_hist(df: pd.DataFrame, outdir: Path) -> None:
    col = _pick_column(df, ["CRAPome Frequency", "Background Frequency", "频率", "frequency"])
    plt.figure(figsize=(7, 4))
    sns.histplot(pd.to_numeric(df[col], errors="coerce").dropna(), bins=20, kde=True)
    plt.xlabel(col)
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(outdir / "01_背景频率分布.png", dpi=300)
    plt.close()


def plot_top_bars(df: pd.DataFrame, outdir: Path) -> None:
    score_col = _pick_column(df, ["PPIRank Score", "FC-A", "FC-B", "Score", "score"])
    label_col = _pick_column(df, ["Prey Name", "Gene names", "Protein", "protein"])
    top = df[[label_col, score_col]].copy()
    top[score_col] = pd.to_numeric(top[score_col], errors="coerce")
    top = top.dropna().sort_values(score_col, ascending=False).head(20)
    plt.figure(figsize=(10, 6))
    sns.barplot(data=top, x=score_col, y=label_col, color="#3b82f6")
    plt.tight_layout()
    plt.savefig(outdir / "02_top候选条形图.png", dpi=300)
    plt.close()


def plot_bait_control_scatter(df: pd.DataFrame, outdir: Path) -> None:
    bait_col = _pick_column(df, ["Bait Count", "bait_count", "bait", "Bait Spectral Count"])
    ctrl_col = _pick_column(df, ["Control Count", "control_count", "control", "Control Spectral Count"])
    label_col = _pick_column(df, ["Prey Name", "Gene names", "Protein", "protein"])
    plot_df = df[[label_col, bait_col, ctrl_col]].copy()
    plot_df[bait_col] = pd.to_numeric(plot_df[bait_col], errors="coerce")
    plot_df[ctrl_col] = pd.to_numeric(plot_df[ctrl_col], errors="coerce")
    plot_df = plot_df.dropna()
    plt.figure(figsize=(6, 6))
    sns.scatterplot(data=plot_df, x=ctrl_col, y=bait_col, s=40)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel(ctrl_col)
    plt.ylabel(bait_col)
    plt.tight_layout()
    plt.savefig(outdir / "03_bait对照散点图.png", dpi=300)
    plt.close()


def plot_repeat_heatmap(df: pd.DataFrame, outdir: Path) -> None:
    row_col = _pick_column(df, ["Prey Name", "Gene names", "Protein", "protein"])
    rep_col = _pick_column(df, ["AP Name", "Replicate", "replicate", "Sample"])
    val_col = _pick_column(df, ["Spectral Count", "spectral_count", "Score", "score"])
    pivot = df.pivot_table(index=row_col, columns=rep_col, values=val_col, aggfunc="sum", fill_value=0)
    if pivot.empty:
        return
    plt.figure(figsize=(12, max(4, min(12, 0.28 * len(pivot)))))
    sns.heatmap(pivot, cmap="viridis")
    plt.tight_layout()
    plt.savefig(outdir / "04_重复一致性热图.png", dpi=300)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="CRAPome 可视化")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for fn in [plot_frequency_hist, plot_top_bars, plot_bait_control_scatter, plot_repeat_heatmap]:
        try:
            fn(df, args.output_dir)
        except Exception as exc:  # noqa: BLE001
            print(f"跳过 {fn.__name__}: {exc}")

    print(f"图形输出完成: {args.output_dir}")


if __name__ == "__main__":
    main()
