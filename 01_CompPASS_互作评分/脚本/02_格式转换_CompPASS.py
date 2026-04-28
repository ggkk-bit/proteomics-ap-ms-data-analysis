import argparse
from pathlib import Path

import pandas as pd


def truthy(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "t", "+"}


def first_existing_column(columns, names):
    for name in names:
        if name in columns:
            return name
    return None


def first_token_series(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.split(";", n=1)
        .str[0]
        .fillna("")
        .str.strip()
    )


def locate_measurement_column(columns, sample_id):
    patterns = [
        (f"MS/MS Count {sample_id}", "MS/MS Count"),
        (f"MS/MS count {sample_id}", "MS/MS count"),
        (f"Intensity.{sample_id}", "Intensity"),
        (f"Intensity {sample_id}", "Intensity"),
        (f"LFQ intensity {sample_id}", "LFQ intensity"),
    ]
    for column_name, count_type in patterns:
        if column_name in columns:
            return column_name, count_type
    return None, None


def main():
    parser = argparse.ArgumentParser(description="把清洗后的 proteinGroups 转为 CompPASS 长表输入。")
    parser.add_argument("--cleaned", default="3.数据清洗与预处理/中间结果/proteinGroups.cleaned.tsv")
    parser.add_argument("--annotation", default="2.样本信息与分组/sample_annotation.tsv")
    parser.add_argument("--output-dir", default="4.CompPASS输入文件/最终输入")
    args = parser.parse_args()

    cleaned = pd.read_csv(args.cleaned, sep="\t", low_memory=False)
    annotation = pd.read_csv(args.annotation, sep="\t", dtype=str).fillna("")

    required = ["sample_id", "ip_name", "bait_name"]
    missing = [col for col in required if col not in annotation.columns]
    if missing:
        raise ValueError(f"样本注释缺少必要字段: {missing}")

    if "keep" in annotation.columns:
        annotation = annotation.loc[annotation["keep"].map(truthy)].copy()

    columns = set(cleaned.columns)
    majority_id_col = first_existing_column(columns, ["Majority protein IDs", "Protein IDs"])
    gene_col = first_existing_column(columns, ["Gene names", "Gene names (primary)"])

    if majority_id_col is None:
        raise ValueError("cleaned 表缺少 Majority protein IDs 或 Protein IDs，无法生成 idPrey。")

    prey_id_series = first_token_series(cleaned[majority_id_col])
    prey_gene_series = first_token_series(cleaned[gene_col]) if gene_col is not None else pd.Series("", index=cleaned.index)

    base_df = pd.DataFrame(
        {
            "idPrey": prey_id_series,
            "preyGene": prey_gene_series,
        },
        index=cleaned.index,
    )
    base_df = base_df.loc[base_df["idPrey"].ne("")].copy()

    interaction_parts = []
    bait_rows = []

    for _, sample in annotation.iterrows():
        sample_id = sample["sample_id"]
        id_run = sample["ip_name"] or sample_id
        id_bait = sample["bait_name"]
        column_name, count_type = locate_measurement_column(cleaned.columns, sample_id)
        if column_name is None:
            continue

        values = pd.to_numeric(cleaned[column_name], errors="coerce").fillna(0)
        positive_mask = values.gt(0) & base_df.index.isin(base_df.index)
        if positive_mask.any():
            subset = base_df.loc[positive_mask].copy()
            subset["idRun"] = id_run
            subset["idBait"] = id_bait
            subset["countPrey"] = values.loc[positive_mask].astype(float)
            subset["countType"] = count_type
            subset["sourceColumn"] = column_name
            interaction_parts.append(
                subset[["idRun", "idBait", "idPrey", "countPrey", "countType", "sourceColumn", "preyGene"]]
            )

        bait_rows.append(
            {
                "idRun": id_run,
                "idBait": id_bait,
                "sample_id": sample_id,
                "replicate_group": sample.get("replicate_group", ""),
                "replicate_letter": sample.get("replicate_letter", ""),
                "raw_file_name": sample.get("raw_file_name", ""),
            }
        )

    if not interaction_parts:
        raise ValueError("没有生成任何 CompPASS 输入记录。请检查 sample_annotation.tsv 和样本列名是否匹配。")

    interaction_df = pd.concat(interaction_parts, ignore_index=True)
    bait_df = pd.DataFrame(bait_rows).drop_duplicates()
    prey_df = (
        interaction_df[["idPrey", "preyGene"]]
        .drop_duplicates()
        .sort_values(["idPrey", "preyGene"])
        .reset_index(drop=True)
    )

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    interaction_df.to_csv(out_dir / "comppass_input.tsv", sep="\t", index=False)
    bait_df.to_csv(out_dir / "bait_run_map.tsv", sep="\t", index=False)
    prey_df.to_csv(out_dir / "prey_reference.tsv", sep="\t", index=False)

    print(f"已输出: {out_dir / 'comppass_input.tsv'}")
    print(f"已输出: {out_dir / 'bait_run_map.tsv'}")
    print(f"已输出: {out_dir / 'prey_reference.tsv'}")


if __name__ == "__main__":
    main()
