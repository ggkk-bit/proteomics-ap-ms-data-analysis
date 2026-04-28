from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


PREY_COL_CANDIDATES = ["Prey Name", "Prey", "Protein IDs", "Gene names", "Protein"]
ANNOTATION_COLUMN_CANDIDATES = [
    ("AP Name", "Bait Name"),
    ("ip_name", "bait_name"),
    ("sample_id", "bait_name"),
    ("样本ID", "bait_name"),
]


def _load_standard_counts_wide(path: Path) -> pd.DataFrame | None:
    counts = pd.read_csv(path, sep="\t", low_memory=False)
    prey_col = next((c for c in PREY_COL_CANDIDATES if c in counts.columns), None)
    if prey_col is None:
        return None

    counts = counts.rename(columns={prey_col: "Prey Name"})
    keep_cols = ["Prey Name"] + [c for c in counts.columns if c != "Prey Name"]
    counts = counts[keep_cols].copy()
    counts["Prey Name"] = counts["Prey Name"].fillna("").astype(str).str.strip()
    counts = counts.loc[counts["Prey Name"] != ""].copy()
    return counts


def _load_mist_matrix(path: Path) -> pd.DataFrame | None:
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        first_line = handle.readline().rstrip("\n").split("\t")
        second_line = handle.readline().rstrip("\n").split("\t")
        third_line = handle.readline().rstrip("\n").split("\t")

    if len(first_line) < 5 or len(second_line) < 5 or len(third_line) < 5:
        return None
    if first_line[:4] != ["a", "b", "c", "IP"]:
        return None
    if second_line[:4] != ["a", "b", "c", "Bait"]:
        return None
    if third_line[0] != "Preys":
        return None

    raw = pd.read_csv(path, sep="\t", header=None, low_memory=False)
    sample_ids = raw.iloc[0, 4:].fillna("").astype(str).str.strip()
    valid_positions = [idx for idx, sample_id in enumerate(sample_ids.tolist(), start=4) if sample_id]

    counts = raw.iloc[3:, [0, *valid_positions]].copy()
    counts.columns = ["Prey Name", *sample_ids[sample_ids != ""].tolist()]
    counts["Prey Name"] = counts["Prey Name"].fillna("").astype(str).str.strip()
    counts = counts.loc[counts["Prey Name"] != ""].copy()
    return counts


def _load_counts_wide(path: Path) -> tuple[pd.DataFrame, str]:
    if not path.exists():
        raise FileNotFoundError(
            f"找不到 counts-wide 文件: {path}\n"
            "如果你没有 sample_counts_wide.tsv，可以改用 MiST 生成的 preprocessed_MAT.txt。"
        )

    standard = _load_standard_counts_wide(path)
    if standard is not None:
        return standard, "standard"

    mist = _load_mist_matrix(path)
    if mist is not None:
        return mist, "mist_preprocessed_mat"

    raise ValueError(
        "无法识别 counts-wide 文件格式。\n"
        "支持两种输入：\n"
        "1. 标准宽表：第一列为 prey，后面每列为一个 AP 样本；\n"
        "2. MiST 的 preprocessed_MAT.txt。"
    )


def _load_annotation(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"找不到样本注释文件: {path}\n"
            "可使用 sample_id + bait_name 两列，或 AP Name + Bait Name 两列。"
        )

    annot = pd.read_csv(path, sep="\t", low_memory=False)
    for ap_col, bait_col in ANNOTATION_COLUMN_CANDIDATES:
        if {ap_col, bait_col}.issubset(annot.columns):
            annot_view = annot[[ap_col, bait_col]].rename(
                columns={ap_col: "AP Name", bait_col: "Bait Name"}
            )
            annot_view["AP Name"] = annot_view["AP Name"].fillna("").astype(str).str.strip()
            annot_view["Bait Name"] = annot_view["Bait Name"].fillna("").astype(str).str.strip()

            if "is_control" in annot.columns:
                is_control = (
                    annot["is_control"]
                    .fillna("")
                    .astype(str)
                    .str.strip()
                    .str.upper()
                    .isin({"Y", "YES", "TRUE", "1"})
                )
                annot_view.loc[is_control & (annot_view["Bait Name"] == ""), "Bait Name"] = "CONTROL"

            annot_view = annot_view.loc[annot_view["AP Name"] != ""].drop_duplicates()
            return annot_view

    raise ValueError(
        "样本注释缺少可识别列名。至少需要以下任一组合：\n"
        "- AP Name + Bait Name\n"
        "- ip_name + bait_name\n"
        "- sample_id + bait_name\n"
        "- 样本ID + bait_name"
    )


def wide_to_long(counts_wide: Path, sample_annotation: Path, output_path: Path) -> pd.DataFrame:
    counts, counts_format = _load_counts_wide(counts_wide)
    annot_view = _load_annotation(sample_annotation)

    value_vars = [c for c in counts.columns if c != "Prey Name"]
    overlap = sorted(set(value_vars) & set(annot_view["AP Name"]))
    if not overlap:
        raise ValueError(
            "counts-wide 的样本列名与样本注释中的 AP 名称完全对不上。\n"
            f"counts-wide 示例列: {', '.join(value_vars[:5])}\n"
            f"样本注释示例 AP: {', '.join(annot_view['AP Name'].head(5).tolist())}"
        )

    long_df = counts.melt(
        id_vars=["Prey Name"],
        value_vars=value_vars,
        var_name="AP Name",
        value_name="Spectral Count",
    )
    long_df = long_df.merge(annot_view, on="AP Name", how="left")
    long_df["Bait Name"] = long_df["Bait Name"].fillna("").astype(str).str.strip()
    long_df.loc[long_df["Bait Name"] == "", "Bait Name"] = "CONTROL"
    long_df["Spectral Count"] = pd.to_numeric(long_df["Spectral Count"], errors="coerce").fillna(0).astype(int)
    long_df = long_df.loc[long_df["Spectral Count"] > 0].copy()
    long_df = long_df[["Bait Name", "AP Name", "Prey Name", "Spectral Count"]]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    long_df.to_csv(output_path, sep="\t", index=False)
    print(f"已识别 counts-wide 格式: {counts_format}")
    return long_df


def main() -> None:
    parser = argparse.ArgumentParser(description="生成 CRAPome 上传表")
    parser.add_argument("--counts-wide", required=True, type=Path)
    parser.add_argument("--sample-annotation", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    df = wide_to_long(args.counts_wide, args.sample_annotation, args.output)
    print(f"已输出 CRAPome 上传表: {args.output} ({len(df)} 行)")


if __name__ == "__main__":
    main()
