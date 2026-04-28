#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""检查 CS Score 所需输入是否齐全。"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="确认 CS Score 输入文件")
    parser.add_argument("--workspace", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--protein-groups", type=Path, default=None)
    parser.add_argument("--annotation", type=Path, default=None)
    parser.add_argument("--member-table", type=Path, default=None)
    args = parser.parse_args()

    workspace = args.workspace
    protein_groups = args.protein_groups or workspace / "1.原始输入" / "proteinGroups.txt"
    annotation = args.annotation or workspace / "2.样本信息与分组" / "sample_annotation_template.tsv"
    member_table = args.member_table or workspace / "3.数据清洗与预处理" / "CS_Score_成员表.tsv"

    report = []
    for label, path in [
        ("proteinGroups.txt", protein_groups),
        ("样本注释", annotation),
        ("成员表", member_table),
    ]:
        report.append((label, str(path), "存在" if path.exists() else "缺失"))

    print(pd.DataFrame(report, columns=["文件类型", "路径", "状态"]).to_string(index=False))

    if protein_groups.exists():
        pg = pd.read_csv(protein_groups, sep="\t", nrows=5, dtype=str)
        key_cols = ["Majority protein IDs", "Gene names", "Reverse", "Potential contaminant", "Only identified by site"]
        present = [c for c in key_cols if c in pg.columns]
        missing = [c for c in key_cols if c not in pg.columns]
        print("\nproteinGroups.txt 关键列：", ", ".join(present) if present else "无")
        if missing:
            print("缺少关键列：", ", ".join(missing))

    if annotation.exists():
        ann = pd.read_csv(annotation, sep="\t", dtype=str).fillna("")
        required = ["sample_id", "bait_name", "replicate_group", "keep"]
        missing = [c for c in required if c not in ann.columns]
        print("\n样本注释样本数：", ann["sample_id"].nunique() if "sample_id" in ann.columns else 0)
        if missing:
            print("样本注释缺少列：", ", ".join(missing))

    if member_table.exists():
        mem = pd.read_csv(member_table, sep="\t", dtype=str).fillna("")
        required = ["purification_id", "sample_id", "bait_name", "replicate_group", "prey_id", "presence"]
        missing = [c for c in required if c not in mem.columns]
        print("\n成员表行数：", len(mem))
        if missing:
            print("成员表缺少列：", ", ".join(missing))


if __name__ == "__main__":
    main()

