#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
通用样本注释表生成脚本

作用：
1. 从 proteinGroups.txt 或清洗后的 TSV 表头中识别样本定量列
2. 支持 Intensity / LFQ intensity / iBAQ 等常见 MaxQuant 列前缀
3. 生成可人工补充的样本注释表初稿
4. 不绑定某个数据集，可复用于 MiST / CompPASS / CRAPome / PPIrank 等流程

示例：
python 00_生成样本注释表.py \
  --input ..\\3.数据清洗与预处理\\proteinGroups_已清洗.tsv \
  --output ..\\2.样本信息与分组\\sample_annotation_auto.tsv
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

SUPPORTED_PREFIXES = [
    "Intensity ",
    "LFQ intensity ",
    "iBAQ ",
    "MS/MS count ",
]

BASE_COLUMNS = [
    "sample_id",
    "raw_file_name",
    "quant_type",
    "source_column",
    "bait_name",
    "condition",
    "is_control",
    "replicate_group",
    "keep",
    "notes",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="从 proteinGroups 表头自动生成样本注释表")
    parser.add_argument("--input", default=r"..\3.数据清洗与预处理\proteinGroups_已清洗.tsv", help="输入表路径")
    parser.add_argument("--output", default=r"..\2.样本信息与分组\sample_annotation_auto.tsv", help="输出注释表路径")
    parser.add_argument(
        "--prefixes",
        default=";".join(SUPPORTED_PREFIXES),
        help="要识别的列前缀，用分号分隔，例如 'Intensity ;LFQ intensity ;iBAQ '"
    )
    parser.add_argument(
        "--prefer",
        default="Intensity",
        help="推荐定量类型标签，仅写入说明字段，不影响全部样本列的提取"
    )
    parser.add_argument(
        "--raw-suffix",
        default=".raw",
        help="自动补到 raw_file_name 后缀，默认 .raw"
    )
    return parser.parse_args()


def read_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            return next(reader)
        except StopIteration as exc:
            raise ValueError(f"输入文件为空: {path}") from exc


def detect_sample_columns(header: list[str], prefixes: list[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()

    for column in header:
        for prefix in prefixes:
            if not column.startswith(prefix):
                continue
            sample_id = column[len(prefix):].strip()
            if not sample_id:
                continue
            quant_type = prefix.strip()
            key = (quant_type.lower(), sample_id.lower())
            if key in seen:
                continue
            seen.add(key)
            rows.append(
                {
                    "sample_id": sample_id,
                    "raw_file_name": "",
                    "quant_type": quant_type,
                    "source_column": column,
                    "bait_name": "",
                    "condition": "",
                    "is_control": "",
                    "replicate_group": "",
                    "keep": "Y",
                    "notes": "",
                }
            )
            break
    return rows


def fill_raw_name(rows: list[dict[str, str]], suffix: str) -> None:
    suffix = suffix or ""
    for row in rows:
        row["raw_file_name"] = f"{row['sample_id']}{suffix}"


def sort_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(rows, key=lambda x: (x["quant_type"].lower(), x["sample_id"].lower()))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=BASE_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    prefixes = [p for p in (x.strip() for x in args.prefixes.split(";")) if p]
    prefixes = [p if p.endswith(" ") else p + " " for p in prefixes]

    if not input_path.exists():
        raise FileNotFoundError(f"找不到输入文件: {input_path}")

    header = read_header(input_path)
    rows = detect_sample_columns(header, prefixes)
    if not rows:
        raise ValueError(
            "未识别到任何样本定量列。请检查输入文件是否包含 Intensity / LFQ intensity / iBAQ / MS/MS count 等列。"
        )

    fill_raw_name(rows, args.raw_suffix)
    rows = sort_rows(rows)
    write_tsv(output_path, rows)

    print(f"已生成样本注释表: {output_path}")
    print(f"识别到样本列数: {len(rows)}")
    print(f"识别前缀: {', '.join(prefix.strip() for prefix in prefixes)}")
    print("下一步: 请人工补充 bait_name、condition、is_control、replicate_group 等实验设计信息。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
