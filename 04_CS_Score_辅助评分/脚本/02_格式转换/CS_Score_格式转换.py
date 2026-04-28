#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""把成员表转换成 CS Score 可直接使用的输入文件。"""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="转换为 CS Score 输入文件")
    parser.add_argument("--workspace", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--member-table", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    workspace = args.workspace
    member_table = args.member_table or workspace / "3.数据清洗与预处理" / "CS_Score_成员表.tsv"
    output_dir = args.output_dir or workspace / "4.CS_Score输入文件"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not member_table.exists():
        raise FileNotFoundError(f"找不到成员表: {member_table}")

    df = pd.read_csv(member_table, sep="\t", dtype=str).fillna("")
    required = ["purification_id", "sample_id", "bait_name", "replicate_group", "prey_id", "presence"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"成员表缺少必要列: {', '.join(missing)}")

    df["presence"] = pd.to_numeric(df["presence"], errors="coerce").fillna(0).astype(int)
    df = df[df["presence"] > 0].copy()
    df["purification_id"] = df["purification_id"].replace("", df["sample_id"])

    member_out = output_dir / "CS_Score_成员表.tsv"
    matrix_out = output_dir / "CS_Score_二值矩阵.tsv"
    pair_out = output_dir / "CS_Score_配对候选.tsv"
    note_out = output_dir / "输入文件说明.md"

    df.to_csv(member_out, sep="\t", index=False)

    matrix = (
        pd.crosstab(df["purification_id"], df["prey_id"])
        .astype(int)
        .sort_index(axis=0)
        .sort_index(axis=1)
    )
    matrix.to_csv(matrix_out, sep="\t")

    pair_counts = []
    grouped = df.groupby("purification_id")["prey_id"].apply(lambda s: sorted(set(s.astype(str))))
    for proteins in grouped:
        if len(proteins) < 2:
            continue
        for a, b in combinations(proteins, 2):
            pair_counts.append((a, b))

    if pair_counts:
        pair_df = pd.DataFrame(pair_counts, columns=["protein_a", "protein_b"])
        pair_df = pair_df.value_counts().reset_index(name="observed_cooccurrence")
    else:
        pair_df = pd.DataFrame(columns=["protein_a", "protein_b", "observed_cooccurrence"])
    pair_df.to_csv(pair_out, sep="\t", index=False)

    note_out.write_text(
        "# CS Score 输入文件说明\n\n"
        "- `CS_Score_成员表.tsv`：推荐作为主输入的长表。\n"
        "- `CS_Score_二值矩阵.tsv`：每个 purification 一行、每个蛋白一列的 0/1 矩阵。\n"
        "- `CS_Score_配对候选.tsv`：仅做共现候选统计和质控，不替代正式评分结果。\n",
        encoding="utf-8",
    )

    print(f"已输出: {member_out}")
    print(f"已输出: {matrix_out}")
    print(f"已输出: {pair_out}")
    print(f"已输出: {note_out}")


if __name__ == "__main__":
    main()

