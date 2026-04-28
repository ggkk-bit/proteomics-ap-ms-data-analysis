#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

SUPPORTED_PREFIXES = [
    "MS/MS Count ",
    "MS/MS count ",
    "Intensity.",
    "Intensity ",
    "LFQ intensity ",
]

BASE_COLUMNS = [
    "sample_id",
    "raw_file_name",
    "plate",
    "well",
    "replicate_letter",
    "ip_name",
    "bait_name",
    "replicate_group",
    "keep",
]


def default_path(*parts: str) -> str:
    return str(PROJECT_DIR.joinpath(*parts))


def resolve_user_path(text: str) -> Path:
    path = Path(text)
    if path.is_absolute():
        return path
    return Path.cwd() / path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="从 proteinGroups 和 archiveSummary 自动生成 CompPASS 样本注释表。")
    parser.add_argument("--input", default=default_path("1.原始输入", "proteinGroups.txt"), help="proteinGroups.txt 路径")
    parser.add_argument("--archive", default=default_path("1.原始输入", "archiveSummary.tsv"), help="archiveSummary.tsv 路径")
    parser.add_argument("--output", default=default_path("2.样本信息与分组", "sample_annotation.tsv"), help="输出样本注释表路径")
    parser.add_argument(
        "--prefixes",
        default=";".join(SUPPORTED_PREFIXES),
        help="要识别的定量列前缀，用分号分隔",
    )
    return parser.parse_args()


def read_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            return next(reader)
        except StopIteration as exc:
            raise ValueError(f"输入文件为空: {path}") from exc


def normalize_prefixes(text: str) -> list[str]:
    parts = [item.strip() for item in text.split(";")]
    parts = [item for item in parts if item]
    return [item if item.endswith((" ", ".")) else item + " " for item in parts]


def detect_samples(header: list[str], prefixes: list[str]) -> list[str]:
    sample_ids: list[str] = []
    seen: set[str] = set()
    for column in header:
        for prefix in prefixes:
            if not column.startswith(prefix):
                continue
            sample_id = column[len(prefix):].strip()
            if not sample_id:
                continue
            key = sample_id.lower()
            if key not in seen:
                sample_ids.append(sample_id)
                seen.add(key)
            break
    return sample_ids


def parse_info(info: str) -> tuple[str, str, str]:
    text = info.strip()
    bait = ""
    plate = ""
    well = ""

    bait_match = re.search(r"Bait:([^()]+)", text)
    if bait_match:
        bait = bait_match.group(1).strip()

    paren_match = re.search(r"\(([^()]+)\)", text)
    if paren_match:
        parts = [item.strip() for item in paren_match.group(1).split(",")]
        if len(parts) >= 1:
            plate = parts[0]
        if len(parts) >= 2:
            well = parts[1]

    return bait, plate, well


def parse_replicate_letter(old_name: str) -> str:
    match = re.search(r"_([A-Z])_", old_name.strip())
    return match.group(1) if match else ""


def load_archive_rows(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"oldName", "newName", "info"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            return {}

        mapping: dict[str, dict[str, str]] = {}
        for row in reader:
            new_name = (row.get("newName") or "").strip()
            if not new_name:
                continue
            sample_id = re.sub(r"\.raw$", "", new_name, flags=re.IGNORECASE)
            if not sample_id:
                continue
            bait_name, plate, well = parse_info(row.get("info") or "")
            replicate_letter = parse_replicate_letter(row.get("oldName") or "")
            mapping.setdefault(
                sample_id.lower(),
                {
                    "sample_id": sample_id,
                    "raw_file_name": new_name,
                    "plate": plate,
                    "well": well,
                    "replicate_letter": replicate_letter,
                    "ip_name": sample_id,
                    "bait_name": bait_name,
                    "replicate_group": bait_name,
                    "keep": "Y",
                },
            )
        return mapping


def build_rows(sample_ids: list[str], archive_map: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for sample_id in sample_ids:
        archive = archive_map.get(sample_id.lower(), {})
        rows.append(
            {
                "sample_id": sample_id,
                "raw_file_name": archive.get("raw_file_name", f"{sample_id}.raw"),
                "plate": archive.get("plate", ""),
                "well": archive.get("well", ""),
                "replicate_letter": archive.get("replicate_letter", ""),
                "ip_name": archive.get("ip_name", sample_id),
                "bait_name": archive.get("bait_name", ""),
                "replicate_group": archive.get("replicate_group", archive.get("bait_name", "")),
                "keep": archive.get("keep", "Y"),
            }
        )
    return sorted(rows, key=lambda row: row["sample_id"].lower())


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=BASE_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_path = resolve_user_path(args.input)
    archive_path = resolve_user_path(args.archive)
    output_path = resolve_user_path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"找不到 proteinGroups.txt: {input_path}")

    header = read_header(input_path)
    sample_ids = detect_samples(header, normalize_prefixes(args.prefixes))
    if not sample_ids:
        raise ValueError("未识别到任何样本定量列。")

    archive_map = load_archive_rows(archive_path)
    rows = build_rows(sample_ids, archive_map)
    write_tsv(output_path, rows)

    matched = sum(1 for row in rows if row["bait_name"])
    print(f"已生成样本注释表: {output_path}")
    print(f"识别到样本数: {len(rows)}")
    print(f"结合 archiveSummary 自动补全的样本数: {matched}")
    print("仍建议抽查 bait_name、replicate_group、replicate_letter 是否符合真实实验设计。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
