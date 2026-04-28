#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PPIrank 输入格式转换脚本。

职责：
- 把清洗后的 bait/prey/count 数据转成 PPIrank 长表
- 统一样本编号、重复编号、对照标识

说明：
- PPIrank 论文强调 spectral counts、负对照和重复性。
- 若你只有 proteinGroups.txt，而没有能还原 spectral counts 的列，请先补上计数表。
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, List


def _sniff_delimiter(path: Path) -> str:
    sample = path.read_text(encoding="utf-8", errors="ignore")[:4096]
    return "\t" if "\t" in sample else ","


def _read_rows(path: Path) -> List[dict]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=_sniff_delimiter(path))
        return list(reader)


def _write_tsv(path: Path, fieldnames: List[str], rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def convert_rows(rows: Iterable[dict], sample_col: str, bait_col: str, prey_col: str, count_col: str, control_col: str) -> List[dict]:
    output: List[dict] = []
    for row in rows:
        count = row.get(count_col, "")
        if str(count).strip() == "":
            raise ValueError(f"缺少光谱计数列 {count_col} 的值：{row}")
        output.append(
            {
                "sample_id": str(row.get(sample_col, "")).strip(),
                "bait_name": str(row.get(bait_col, "")).strip(),
                "prey_id": str(row.get(prey_col, "")).strip(),
                "spectral_count": count,
                "is_control": str(row.get(control_col, "")).strip(),
            }
        )
    return output


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="将清洗后的数据转成 PPIrank 输入长表")
    parser.add_argument("--input", required=True, help="清洗后的 TSV/CSV")
    parser.add_argument("--output", required=True, help="PPIrank 输入长表")
    parser.add_argument("--control-output", help="PPIrank 负对照子表")
    parser.add_argument("--summary-output", help="样本汇总表")
    parser.add_argument("--sample-col", default="sample_id")
    parser.add_argument("--bait-col", default="bait_name")
    parser.add_argument("--prey-col", default="Gene names")
    parser.add_argument("--count-col", default="spectral_count")
    parser.add_argument("--control-col", default="is_control")
    args = parser.parse_args()

    rows = _read_rows(Path(args.input))
    converted = convert_rows(rows, args.sample_col, args.bait_col, args.prey_col, args.count_col, args.control_col)
    _write_tsv(Path(args.output), ["sample_id", "bait_name", "prey_id", "spectral_count", "is_control"], converted)

    if args.control_output:
        control_rows = [row for row in converted if str(row.get("is_control", "")).strip().lower() in {"1", "true", "yes", "+"}]
        _write_tsv(Path(args.control_output), ["sample_id", "bait_name", "prey_id", "spectral_count", "is_control"], control_rows)

    if args.summary_output:
        summary: dict[tuple[str, str], int] = {}
        for row in converted:
            key = (row["sample_id"], row["bait_name"])
            summary[key] = summary.get(key, 0) + 1
        summary_rows = [
            {"sample_id": sample_id, "bait_name": bait_name, "条目数": count}
            for (sample_id, bait_name), count in sorted(summary.items())
        ]
        _write_tsv(Path(args.summary_output), ["sample_id", "bait_name", "条目数"], summary_rows)


if __name__ == "__main__":
    main()
