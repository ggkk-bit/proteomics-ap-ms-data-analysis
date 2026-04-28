#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PPIrank 数据清洗脚本。

职责：
- 删除 Reverse / Potential contaminant / Only identified by site
- 删除缺失 bait / prey / spectral count 的记录
- 生成清洗报告

说明：
- 这是针对 AP-MS / MaxQuant 导出的通用清洗器。
- 若原始表没有光谱计数列，请先从 evidence/msms 或上游汇总文件补齐。
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Tuple


TRUE_VALUES = {"+", "1", "true", "yes", "y"}


@dataclass
class CleanReport:
    input_rows: int = 0
    kept_rows: int = 0
    removed_reverse: int = 0
    removed_contaminant: int = 0
    removed_site_only: int = 0
    removed_missing_key: int = 0


def _sniff_delimiter(path: Path) -> str:
    with path.open("r", encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
    if "\t" in sample:
        return "\t"
    if "," in sample:
        return ","
    return "\t"


def _read_table(path: Path) -> Tuple[List[str], List[dict]]:
    delimiter = _sniff_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def _is_true(value: str | None) -> bool:
    return str(value or "").strip().lower() in TRUE_VALUES


def clean_rows(rows: Iterable[dict], bait_col: str, prey_col: str, count_col: str) -> Tuple[List[dict], CleanReport]:
    report = CleanReport()
    kept: List[dict] = []
    for row in rows:
        report.input_rows += 1
        reverse = _is_true(row.get("Reverse"))
        contaminant = _is_true(row.get("Potential contaminant"))
        site_only = _is_true(row.get("Only identified by site"))
        bait = str(row.get(bait_col, "")).strip()
        prey = str(row.get(prey_col, "")).strip()
        count = str(row.get(count_col, "")).strip()

        if reverse:
            report.removed_reverse += 1
            continue
        if contaminant:
            report.removed_contaminant += 1
            continue
        if site_only:
            report.removed_site_only += 1
            continue
        if not bait or not prey or not count:
            report.removed_missing_key += 1
            continue

        kept.append(row)
        report.kept_rows += 1
    return kept, report


def write_table(path: Path, fieldnames: List[str], rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_report(path: Path, report: CleanReport) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        ("输入行数", report.input_rows),
        ("保留行数", report.kept_rows),
        ("删除 Reverse", report.removed_reverse),
        ("删除 Potential contaminant", report.removed_contaminant),
        ("删除 Only identified by site", report.removed_site_only),
        ("删除缺失 bait/prey/计数", report.removed_missing_key),
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["指标", "数量"])
        writer.writerows(lines)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="清洗 PPIrank 输入前的蛋白组结果")
    parser.add_argument("--input", required=True, help="proteinGroups.txt 或其导出的 TSV/CSV")
    parser.add_argument("--output", required=True, help="清洗后输出 TSV")
    parser.add_argument("--report", required=True, help="清洗报告 TSV")
    parser.add_argument("--bait-col", default="bait_name", help="bait 列名")
    parser.add_argument("--prey-col", default="Gene names", help="prey 列名")
    parser.add_argument("--count-col", default="spectral_count", help="光谱计数列名")
    args = parser.parse_args()

    input_path = Path(args.input)
    fieldnames, rows = _read_table(input_path)
    cleaned, report = clean_rows(rows, args.bait_col, args.prey_col, args.count_col)
    write_table(Path(args.output), fieldnames, cleaned)
    write_report(Path(args.report), report)


if __name__ == "__main__":
    main()
