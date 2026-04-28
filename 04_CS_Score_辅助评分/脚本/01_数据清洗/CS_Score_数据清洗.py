#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""CS Score 数据清洗。

输入：
- `1.原始输入/proteinGroups.txt`
- `2.样本信息与分组/sample_annotation_template.tsv`

输出：
- `3.数据清洗与预处理/清洗后proteinGroups.tsv`
- `3.数据清洗与预处理/CS_Score_成员表.tsv`
- `3.数据清洗与预处理/清洗统计.tsv`
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def first_token(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    return text.split(";")[0].strip()


def pick_column(columns: list[str], candidates: list[str]) -> str | None:
    for candidate in candidates:
        if candidate in columns:
            return candidate
    return None


def optional_series(df: pd.DataFrame, column: str, default: str = "") -> pd.Series:
    if column in df.columns:
        return df[column].astype(str).str.strip()
    return pd.Series([default] * len(df), index=df.index, dtype=str)


def load_annotation(path: Path, sample_ids: list[str]) -> pd.DataFrame:
    if path.exists():
        annotation = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        return annotation

    template = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "raw_file_name": [f"{sample_id}.raw" for sample_id in sample_ids],
            "plate": "",
            "well": "",
            "replicate_letter": "",
            "ip_name": sample_ids,
            "bait_name": "",
            "replicate_group": "",
            "keep": "Y",
        }
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    template.to_csv(path, sep="\t", index=False)
    raise FileNotFoundError(f"样本注释模板已生成，请先补全后再运行：{path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="清洗 proteinGroups.txt 并展开 CS Score 成员表")
    parser.add_argument("--workspace", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--protein-groups", type=Path, default=None)
    parser.add_argument("--annotation", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--count-prefix", default="Intensity ")
    args = parser.parse_args()

    workspace = args.workspace
    protein_groups = args.protein_groups or workspace / "1.原始输入" / "proteinGroups.txt"
    annotation_file = args.annotation or workspace / "2.样本信息与分组" / "sample_annotation_template.tsv"
    output_dir = args.output_dir or workspace / "3.数据清洗与预处理"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not protein_groups.exists():
        raise FileNotFoundError(f"找不到 proteinGroups.txt: {protein_groups}")

    pg = pd.read_csv(protein_groups, sep="\t", dtype=str).fillna("")
    original_n = len(pg)

    for column in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if column in pg.columns:
            pg = pg[pg[column].ne("+")].copy()

    after_marker_n = len(pg)

    majority_col = pick_column(list(pg.columns), ["Majority protein IDs", "Majority protein IDs ".strip()])
    gene_col = pick_column(list(pg.columns), ["Gene names", "Gene names ".strip()])
    seq_col = pick_column(list(pg.columns), ["Sequence length", "Sequence length ".strip()])

    if not majority_col:
        raise KeyError("proteinGroups.txt 中缺少 `Majority protein IDs` 列")

    count_cols = [c for c in pg.columns if c.startswith(args.count_prefix)]
    sample_ids = [c[len(args.count_prefix):].strip() for c in count_cols]
    if not count_cols:
        raise KeyError(f"没有找到以 `{args.count_prefix}` 开头的定量列")

    annotation = load_annotation(annotation_file, sample_ids)
    annotation = annotation.copy()
    required = ["sample_id", "bait_name", "replicate_group", "keep"]
    missing = [c for c in required if c not in annotation.columns]
    if missing:
        raise KeyError(f"样本注释缺少必要列: {', '.join(missing)}")

    annotation["keep"] = annotation["keep"].astype(str).str.upper().str.strip()
    annotation = annotation[annotation["keep"].isin(["Y", "YES", "1", "TRUE"])].copy()
    annotation["sample_id"] = annotation["sample_id"].astype(str).str.strip()
    annotation["ip_name"] = annotation["sample_id"]
    if "ip_name" in annotation.columns:
        annotation["ip_name"] = optional_series(annotation, "ip_name", "")

    annotation["raw_file_name"] = optional_series(annotation, "raw_file_name", "")
    annotation.loc[annotation["raw_file_name"].eq(""), "raw_file_name"] = annotation.loc[
        annotation["raw_file_name"].eq(""), "sample_id"
    ] + ".raw"
    annotation["bait_name"] = optional_series(annotation, "bait_name", "")
    annotation["replicate_group"] = optional_series(annotation, "replicate_group", "")
    annotation["is_control"] = optional_series(annotation, "is_control", "").str.upper()

    present_samples = set(sample_ids)
    missing_samples = sorted(set(annotation["sample_id"]) - present_samples)
    if missing_samples:
        raise KeyError(f"这些 sample_id 不在 proteinGroups 定量列里: {', '.join(missing_samples)}")

    keep_rows = []
    member_rows = []
    for _, row in annotation.iterrows():
        sample_id = row["sample_id"]
        quant_col = f"{args.count_prefix}{sample_id}"
        if quant_col not in pg.columns:
            continue

        counts = pd.to_numeric(pg[quant_col], errors="coerce").fillna(0)
        prey_id = pg[majority_col].map(first_token)
        gene_name = pg[gene_col].map(first_token) if gene_col else prey_id
        gene_name = gene_name.where(gene_name.ne(""), prey_id)
        keep_idx = counts.gt(0) & prey_id.ne("")

        subset = pg.loc[keep_idx].copy()
        subset["sample_id"] = sample_id
        subset["bait_name"] = row["bait_name"]
        subset["replicate_group"] = row["replicate_group"]
        subset["purification_id"] = row.get("ip_name", sample_id) or sample_id
        subset["raw_file_name"] = row["raw_file_name"]
        subset["prey_id"] = prey_id[keep_idx].values
        subset["gene_name"] = gene_name[keep_idx].values
        subset["intensity"] = counts[keep_idx].values
        subset["presence"] = 1
        subset["is_control"] = row.get("is_control", "")
        subset["plate"] = row.get("plate", "")
        subset["well"] = row.get("well", "")
        subset["replicate_letter"] = row.get("replicate_letter", "")
        keep_rows.append(subset)

        member_rows.append(
            pd.DataFrame(
                {
                    "purification_id": subset["purification_id"],
                    "sample_id": subset["sample_id"],
                    "raw_file_name": subset["raw_file_name"],
                    "bait_name": subset["bait_name"],
                    "replicate_group": subset["replicate_group"],
                    "prey_id": subset["prey_id"],
                    "gene_name": subset["gene_name"],
                    "intensity": subset["intensity"],
                    "presence": subset["presence"],
                    "is_control": subset["is_control"],
                    "plate": subset["plate"],
                    "well": subset["well"],
                    "replicate_letter": subset["replicate_letter"],
                }
            )
        )

    cleaned_pg = pd.concat(keep_rows, ignore_index=True) if keep_rows else pd.DataFrame()
    member_table = pd.concat(member_rows, ignore_index=True) if member_rows else pd.DataFrame()

    cleaned_pg_path = output_dir / "清洗后proteinGroups.tsv"
    member_path = output_dir / "CS_Score_成员表.tsv"
    stats_path = output_dir / "清洗统计.tsv"

    cleaned_pg.to_csv(cleaned_pg_path, sep="\t", index=False)
    member_table.to_csv(member_path, sep="\t", index=False)

    stats = pd.DataFrame(
        [
            ["原始行数", original_n],
            ["去标记后行数", after_marker_n],
            ["成员表行数", len(member_table)],
            ["样本数", annotation["sample_id"].nunique()],
            ["bait 数", annotation["bait_name"].replace("", pd.NA).dropna().nunique()],
            ["重复组数", annotation["replicate_group"].replace("", pd.NA).dropna().nunique()],
        ],
        columns=["指标", "数值"],
    )
    stats.to_csv(stats_path, sep="\t", index=False)

    print(f"已输出: {cleaned_pg_path}")
    print(f"已输出: {member_path}")
    print(f"已输出: {stats_path}")


if __name__ == "__main__":
    main()
