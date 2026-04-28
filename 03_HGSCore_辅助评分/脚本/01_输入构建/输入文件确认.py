from __future__ import annotations

import argparse
import csv
from pathlib import Path


FILTER_COLUMNS = ["Reverse", "Potential contaminant", "Only identified by site"]
COUNT_PREFIXES = ["Peptides ", "Razor + unique peptides ", "LFQ intensity ", "Intensity ", "iBAQ "]
REQUIRED_SAMPLE_COLUMNS = ["sample_id", "bait_name", "replicate_group", "keep"]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        return next(csv.reader(fh, delimiter="\t"))


def detect_measure_columns(header: list[str], sample_ids: list[str]) -> dict[str, str]:
    found: dict[str, str] = {}
    for sid in sample_ids:
        for prefix in COUNT_PREFIXES:
            candidate = f"{prefix}{sid}"
            if candidate in header:
                found[sid] = candidate
                break
    return found


def collect_sample_ids(rows: list[dict[str, str]]) -> list[str]:
    sample_ids: list[str] = []
    for row in rows:
        sid = row.get("sample_id", "").strip()
        raw_name = row.get("raw_file_name", "").strip()
        if sid:
            sample_ids.append(sid)
        if raw_name:
            sample_ids.append(Path(raw_name).stem)
    seen = set()
    unique: list[str] = []
    for sid in sample_ids:
        if sid not in seen:
            seen.add(sid)
            unique.append(sid)
    return unique


def write_report(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="检查 HGSCore 输入文件是否齐全")
    parser.add_argument("--protein-groups", required=True, type=Path)
    parser.add_argument("--sample-annotation", required=True, type=Path)
    parser.add_argument("--length-map", type=Path, default=None)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()

    header = read_header(args.protein_groups)
    samples = read_tsv(args.sample_annotation)
    sample_header = read_header(args.sample_annotation)

    missing_sample_cols = [c for c in REQUIRED_SAMPLE_COLUMNS if c not in sample_header]
    sample_ids = collect_sample_ids(samples)
    measure_columns = detect_measure_columns(header, sample_ids)
    missing_measure = [sid for sid in sample_ids if sid not in measure_columns]

    lines = [
        "# HGSCore 输入文件检查报告",
        "",
        f"- proteinGroups 文件: `{args.protein_groups}`",
        f"- 样本注释文件: `{args.sample_annotation}`",
        f"- 长度映射文件: `{args.length_map}`" if args.length_map else "- 长度映射文件: 未提供",
        "",
        "## 检查结果",
        f"- proteinGroups 列数: {len(header)}",
        f"- 样本数量: {len(sample_ids)}",
        f"- 可匹配的计数列数量: {len(measure_columns)}",
        f"- 缺失的样本计数列: {len(missing_measure)}",
    ]

    if missing_sample_cols:
        lines += ["", "## 样本注释缺失列", *[f"- {c}" for c in missing_sample_cols]]
    if missing_measure:
        lines += ["", "## 未找到的样本计数列", *[f"- {sid}" for sid in missing_measure[:50]]]
        if len(missing_measure) > 50:
            lines.append(f"- 其余 {len(missing_measure) - 50} 个样本未列出")

    if "Sequence length" not in header and "Length" not in header and not args.length_map:
        lines += [
            "",
            "## 不确定点",
            "- 当前 proteinGroups 中未发现明显的长度列，且未提供长度映射文件",
            "- HGSCore 原始定义依赖长度信息，建议补一个长度表",
        ]

    report_path = args.out_dir / "输入文件检查报告.md"
    write_report(report_path, lines)
    print(report_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
