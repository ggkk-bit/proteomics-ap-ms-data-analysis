#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PPIrank 输入文件确认脚本。

职责：
- 检查目录里是否存在 PPIrank 所需文件
- 检查最小字段集
- 给出缺失项提示

这是一个确认脚本，不做分析。
"""

from __future__ import annotations

import csv
from pathlib import Path


REQUIRED_FILES = [
    "1.原始输入/proteinGroups.txt",
    "2.样本信息与分组/sample_annotation_template.tsv",
    "3.数据清洗与预处理/清洗后蛋白组.tsv",
    "4.PPIrank输入文件/PPIrank_输入_长表.tsv",
]

RECOMMENDED_FILES = [
    "4.PPIrank输入文件/PPIrank_输入_负对照.tsv",
    "4.PPIrank输入文件/PPIrank_输入_样本汇总.tsv",
    "4.PPIrank输入文件/PPIrank_输入_蛋白长度.tsv",
]

REQUIRED_COLUMNS = {"sample_id", "bait_name", "prey_id", "spectral_count", "is_control"}


def _sniff_delimiter(path: Path) -> str:
    sample = path.read_text(encoding="utf-8", errors="ignore")[:4096]
    return "\t" if "\t" in sample else ","


def check_files(root: Path) -> list[str]:
    missing = []
    for rel in REQUIRED_FILES:
        if not (root / rel).exists():
            missing.append(rel)
    return missing


def check_columns(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=_sniff_delimiter(path))
        columns = set(reader.fieldnames or [])
    return sorted(REQUIRED_COLUMNS - columns)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="确认 PPIrank 输入文件是否齐全")
    parser.add_argument("--root", required=True, help="PPIrank 项目根目录")
    args = parser.parse_args()

    root = Path(args.root)
    missing_files = check_files(root)
    if missing_files:
        print("缺失文件：")
        for item in missing_files:
            print(f"- {item}")
    else:
        print("文件齐全。")

    recommended_missing = [rel for rel in RECOMMENDED_FILES if not (root / rel).exists()]
    if recommended_missing:
        print("建议补齐文件：")
        for item in recommended_missing:
            print(f"- {item}")

    long_table = root / "4.PPIrank输入文件/PPIrank_输入_长表.tsv"
    if long_table.exists():
        missing_columns = check_columns(long_table)
        if missing_columns:
            print("PPIrank 长表缺失列：")
            for item in missing_columns:
                print(f"- {item}")
        else:
            print("PPIrank 长表字段齐全。")


if __name__ == "__main__":
    main()
