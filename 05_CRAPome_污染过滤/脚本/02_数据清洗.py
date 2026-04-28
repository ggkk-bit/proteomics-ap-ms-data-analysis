"""CRAPome 数据清洗脚本。

默认只做 MaxQuant 常规三类硬过滤：
- Reverse/反向匹配
- Potential contaminant/潜在污染物
- Only identified by site/仅位点鉴定蛋白

可选附加规则：
- 去掉缺失主 ID 的蛋白组
- 按 Unique peptides 设最小阈值
这些规则需要结合项目具体质量和重复数决定，不建议无条件强加。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _flagged(series: pd.Series) -> pd.Series:
    values = series.fillna("").astype(str).str.strip()
    return values.isin({"+", "1", "True", "TRUE", "yes", "YES"}) | values.str.contains(
        "reverse|potential contaminant|only identified by site", case=False, na=False
    )


def clean_protein_groups(input_path: Path, output_path: Path, min_unique_peptides: int | None) -> pd.DataFrame:
    df = pd.read_csv(input_path, sep="\t", low_memory=False)

    for col in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if col in df.columns:
            df = df.loc[~_flagged(df[col])]

    if "Gene names" in df.columns:
        df = df.loc[df["Gene names"].fillna("").astype(str).str.strip() != ""]

    if min_unique_peptides is not None and "Unique peptides" in df.columns:
        df = df.loc[pd.to_numeric(df["Unique peptides"], errors="coerce").fillna(0) >= min_unique_peptides]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    return df


def main() -> None:
    parser = argparse.ArgumentParser(description="清洗 proteinGroups.txt")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-unique-peptides", type=int, default=None)
    args = parser.parse_args()

    df = clean_protein_groups(args.input, args.output, args.min_unique_peptides)
    print(f"已输出清洗结果: {args.output} ({len(df)} 行)")


if __name__ == "__main__":
    main()

