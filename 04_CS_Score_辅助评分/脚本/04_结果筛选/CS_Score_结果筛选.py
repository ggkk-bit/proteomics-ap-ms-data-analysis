#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""筛选 CS Score 结果。"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="筛选 CS Score 结果")
    parser.add_argument("--workspace", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--input", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--min-cooccurrence", type=int, default=2)
    parser.add_argument("--min-score", type=float, default=2.0)
    parser.add_argument("--keep-negative", action="store_true")
    args = parser.parse_args()

    workspace = args.workspace
    input_file = args.input or workspace / "6.CS_Score评分结果" / "CS_Score_结果.tsv"
    output_dir = args.output_dir or workspace / "6.CS_Score评分结果"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_file.exists():
        raise FileNotFoundError(f"找不到结果文件: {input_file}")

    df = pd.read_csv(input_file, sep="\t", dtype=str).fillna("")
    for col in ["observed_cooccurrence", "cs_score", "score", "z_score"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    cooc_col = "observed_cooccurrence" if "observed_cooccurrence" in df.columns else None
    score_col = "cs_score" if "cs_score" in df.columns else ("score" if "score" in df.columns else ("z_score" if "z_score" in df.columns else None))

    filtered = df.copy()
    if cooc_col:
        filtered = filtered[filtered[cooc_col] >= args.min_cooccurrence].copy()
    if score_col:
        filtered = filtered[filtered[score_col].notna()].copy()
        filtered = filtered[filtered[score_col] >= args.min_score].copy()

    if not args.keep_negative and score_col and score_col in filtered.columns:
        filtered = filtered[filtered[score_col] >= 0].copy()

    out_file = output_dir / "CS_Score_筛选后候选.tsv"
    filtered.to_csv(out_file, sep="\t", index=False)

    note = output_dir / "结果筛选说明.md"
    note.write_text(
        "# CS Score 结果筛选说明\n\n"
        "## 默认规则\n\n"
        f"- 观察共现次数 `>= {args.min_cooccurrence}`\n"
        f"- CS score `>= {args.min_score}`\n"
        "- 默认只保留正向高分候选\n\n"
        "## 这样筛的原因\n\n"
        "- 论文层面已经提示：共现次数太少的蛋白对不稳定，先用 `>=2` 做硬过滤更稳妥。\n"
        "- `CS score` 本质是 Z-score，`>=2` 是常见的高置信经验门槛，适合作为初筛。\n"
        "- 如果项目里已有对照库或已知负集，可以再收紧。\n",
        encoding="utf-8",
    )

    print(f"已输出: {out_file}")
    print(f"已输出: {note}")


if __name__ == "__main__":
    main()

