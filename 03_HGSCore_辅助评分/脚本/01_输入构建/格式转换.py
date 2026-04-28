from __future__ import annotations

import argparse
import csv
from pathlib import Path


COUNT_PREFIXES = ["Peptides ", "Razor + unique peptides ", "LFQ intensity ", "Intensity ", "iBAQ "]


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def pick_value_column(header: list[str], sample_id: str, preference: str) -> str | None:
    preference_order = {
        "count": ["Peptides ", "Razor + unique peptides ", "LFQ intensity ", "Intensity ", "iBAQ "],
        "intensity": ["LFQ intensity ", "Intensity ", "iBAQ ", "Peptides ", "Razor + unique peptides "],
    }.get(preference, COUNT_PREFIXES)
    for prefix in preference_order:
        candidate = f"{prefix}{sample_id}"
        if candidate in header:
            return candidate
    return None


def first_nonempty(row: dict[str, str], keys: list[str]) -> str:
    for key in keys:
        value = row.get(key, "").strip()
        if value:
            return value
    return ""


def collect_sample_ids(rows: list[dict[str, str]]) -> tuple[list[str], dict[str, str]]:
    sample_ids: list[str] = []
    sample_bait: dict[str, str] = {}
    for row in rows:
        bait = row.get("bait_name", "").strip()
        sid = row.get("sample_id", "").strip()
        raw_name = row.get("raw_file_name", "").strip()
        for candidate in [sid, Path(raw_name).stem if raw_name else ""]:
            if candidate:
                sample_ids.append(candidate)
                if bait:
                    sample_bait[candidate] = bait
    unique: list[str] = []
    seen = set()
    for sid in sample_ids:
        if sid not in seen:
            seen.add(sid)
            unique.append(sid)
    return unique, sample_bait


def parse_length(value: str) -> float | None:
    try:
        return float(value)
    except Exception:
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description="把清洗后的 proteinGroups 转成 HGSCore 输入文件")
    parser.add_argument("--cleaned", required=True, type=Path)
    parser.add_argument("--sample-annotation", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--value-mode", default="count", choices=["count", "intensity"])
    args = parser.parse_args()

    header, rows = read_table(args.cleaned)
    _, samples = read_table(args.sample_annotation)
    sample_ids, sample_bait = collect_sample_ids(samples)

    value_columns = {sid: pick_value_column(header, sid, args.value_mode) for sid in sample_ids}
    missing = [sid for sid, col in value_columns.items() if not col]
    if missing:
        raise SystemExit(f"未找到这些样本的计数/定量列: {', '.join(missing[:20])}")

    long_rows: list[dict[str, str]] = []
    matrix: dict[str, dict[str, float]] = {}
    normalized: dict[str, dict[str, float]] = {}
    lengths: dict[str, str] = {}

    for row in rows:
        prey = first_nonempty(row, ["clean_protein_id", "Majority protein IDs", "Protein IDs"])
        length_text = first_nonempty(row, ["protein_length", "Sequence length", "Length"])
        if prey and length_text:
            lengths[prey] = length_text
        for sid, col in value_columns.items():
            raw = row.get(col, "").strip()
            value = parse_length(raw)
            if value is None:
                continue
            matrix.setdefault(prey, {})[sid] = value
            length = parse_length(length_text) or 0.0
            norm = value / length if length > 0 else value
            normalized.setdefault(prey, {})[sid] = norm
            long_rows.append(
                {
                    "sample_id": sid,
                    "bait_name": sample_bait.get(sid, ""),
                    "prey_id": prey,
                    "value": str(value),
                    "length": length_text,
                    "length_normalized_value": str(norm),
                    "source_column": col,
                }
            )

    args.out_dir.mkdir(parents=True, exist_ok=True)

    long_path = args.out_dir / "HGSCore_长表.tsv"
    with long_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(long_rows[0].keys()) if long_rows else ["sample_id", "bait_name", "prey_id", "value", "length", "length_normalized_value", "source_column"], delimiter="\t")
        writer.writeheader()
        writer.writerows(long_rows)

    matrix_path = args.out_dir / "HGSCore_run_prey矩阵.tsv"
    prey_ids = sorted(matrix)
    with matrix_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["prey_id", *sample_ids])
        for prey in prey_ids:
            writer.writerow([prey, *[matrix.get(prey, {}).get(sid, 0) for sid in sample_ids]])

    norm_path = args.out_dir / "HGSCore_run_prey长度归一化矩阵.tsv"
    with norm_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["prey_id", *sample_ids])
        for prey in prey_ids:
            writer.writerow([prey, *[normalized.get(prey, {}).get(sid, 0) for sid in sample_ids]])

    length_path = args.out_dir / "HGSCore_长度表.tsv"
    with length_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["prey_id", "protein_length"])
        for prey, length in sorted(lengths.items()):
            writer.writerow([prey, length])

    note_path = args.out_dir / "HGSCore输入文件说明.md"
    note_path.write_text(
        "\n".join(
            [
                "# HGSCore 输入文件说明",
                "",
                "- `HGSCore_长表.tsv`: 样本-蛋白长表",
                "- `HGSCore_run_prey矩阵.tsv`: 原始计数矩阵",
                "- `HGSCore_run_prey长度归一化矩阵.tsv`: 长度归一化矩阵",
                "- `HGSCore_长度表.tsv`: 蛋白长度表",
                "",
                "不确定点：原始 HGSCore 是否接受长度归一化后的矩阵，取决于具体实现版本；本脚本同时输出原始矩阵和适配矩阵。",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(long_path)
    print(matrix_path)
    print(norm_path)
    print(length_path)
    print(note_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
