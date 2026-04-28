from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


FILTER_COLUMNS = ["Reverse", "Potential contaminant", "Only identified by site"]


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def first_nonempty(row: dict[str, str], keys: list[str]) -> str:
    for key in keys:
        value = row.get(key, "").strip()
        if value:
            return value
    return ""


def split_first(value: str) -> str:
    return value.split(";")[0].strip() if value else ""


def load_length_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        mapping: dict[str, str] = {}
        for row in reader:
            protein = first_nonempty(row, ["protein_id", "Protein ID", "Protein IDs", "Majority protein IDs"]).split(";")[0].strip()
            length = first_nonempty(row, ["protein_length", "Sequence length", "Length"])
            if protein and length:
                mapping[protein] = length
        return mapping


def main() -> int:
    parser = argparse.ArgumentParser(description="清洗 MaxQuant proteinGroups.txt")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--length-map", type=Path, default=None)
    args = parser.parse_args()

    header, rows = read_rows(args.input)
    length_map = load_length_map(args.length_map)
    output_rows: list[dict[str, str]] = []
    reasons = Counter()

    for row in rows:
        reason = []
        for col in FILTER_COLUMNS:
            if row.get(col, "").strip() in {"+", "True", "true", "1", "yes", "Y"}:
                reason.append(col)
        if not first_nonempty(row, ["Protein IDs", "Majority protein IDs"]):
            reason.append("Missing protein IDs")
        if reason:
            reasons.update(reason)
            continue

        protein_id = split_first(first_nonempty(row, ["Majority protein IDs", "Protein IDs"]))
        gene_name = split_first(first_nonempty(row, ["Gene names", "Gene name"]))
        length = first_nonempty(row, ["Sequence length", "Length"])
        if not length and protein_id in length_map:
            length = length_map[protein_id]

        cleaned = dict(row)
        cleaned["clean_protein_id"] = protein_id
        cleaned["clean_gene_name"] = gene_name
        cleaned["protein_length"] = length
        output_rows.append(cleaned)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=(list(output_rows[0].keys()) if output_rows else header), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    log_path = args.output.with_suffix(".清洗日志.md")
    lines = [
        "# HGSCore 数据清洗日志",
        "",
        f"- 输入行数: {len(rows)}",
        f"- 输出行数: {len(output_rows)}",
        f"- 去除行数: {len(rows) - len(output_rows)}",
        "",
        "## 过滤理由统计",
    ]
    for key, value in reasons.most_common():
        lines.append(f"- {key}: {value}")
    if not reasons:
        lines.append("- 无")
    if not args.length_map:
        lines += ["", "## 不确定点", "- 未提供长度映射文件，protein_length 只能依赖 proteinGroups 里已有列"]
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(args.output)
    print(log_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
