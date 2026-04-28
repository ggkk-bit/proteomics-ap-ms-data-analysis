import argparse
import json
from pathlib import Path

import pandas as pd


TRUTHY = {"+", "1", "true", "yes", "y", "t"}


def as_bool(series: pd.Series) -> pd.Series:
    if series is None:
        return pd.Series(False, index=[])
    values = series.fillna("").astype(str).str.strip().str.lower()
    return values.isin(TRUTHY)


def find_first_column(df: pd.DataFrame, aliases):
    lookup = {str(c).strip().lower(): c for c in df.columns}
    for alias in aliases:
        key = alias.strip().lower()
        if key in lookup:
            return lookup[key]
    return None


def load_kept_samples(annotation_path: Path):
    if not annotation_path.exists():
        return []
    ann = pd.read_csv(annotation_path, sep="\t", dtype=str).fillna("")
    if "sample_id" not in ann.columns:
        return []
    if "keep" in ann.columns:
        keep_mask = as_bool(ann["keep"])
        ann = ann.loc[keep_mask]
    return ann["sample_id"].astype(str).tolist()


def detect_quant_columns(columns, sample_ids):
    candidates = []
    all_columns = list(columns)
    prefixes = ["MS/MS Count ", "MS/MS count ", "Intensity.", "Intensity ", "LFQ intensity "]
    if sample_ids:
        for sample_id in sample_ids:
            for prefix in prefixes:
                name = f"{prefix}{sample_id}"
                if name in all_columns:
                    candidates.append(name)
                    break
    if candidates:
        return candidates
    fallback = []
    for col in all_columns:
        if col.startswith("MS/MS Count ") or col.startswith("MS/MS count ") or col.startswith("Intensity.") or col.startswith("Intensity ") or col.startswith("LFQ intensity "):
            fallback.append(col)
    return fallback


def main():
    parser = argparse.ArgumentParser(description="清洗 proteinGroups.txt，生成 CompPASS 预处理输入。")
    parser.add_argument("--input", default="1.原始输入/proteinGroups.txt")
    parser.add_argument("--annotation", default="2.样本信息与分组/sample_annotation.tsv")
    parser.add_argument("--output", default="3.数据清洗与预处理/中间结果/proteinGroups.cleaned.tsv")
    parser.add_argument("--log", default="3.数据清洗与预处理/过滤记录/过滤统计.json")
    args = parser.parse_args()

    input_path = Path(args.input)
    annotation_path = Path(args.annotation)
    output_path = Path(args.output)
    log_path = Path(args.log)

    df = pd.read_csv(input_path, sep="\t", low_memory=False)
    initial_rows = len(df)

    reverse_col = find_first_column(df, ["Reverse", "反向匹配"])
    contaminant_col = find_first_column(df, ["Potential contaminant", "潜在污染物"])
    site_only_col = find_first_column(df, ["Only identified by site", "仅位点鉴定蛋白"])
    majority_id_col = find_first_column(df, ["Majority protein IDs", "Protein IDs"])

    removed = {}

    if reverse_col:
        mask = as_bool(df[reverse_col])
        removed["reverse"] = int(mask.sum())
        df = df.loc[~mask].copy()
    else:
        removed["reverse"] = 0

    if contaminant_col:
        mask = as_bool(df[contaminant_col])
        removed["potential_contaminant"] = int(mask.sum())
        df = df.loc[~mask].copy()
    else:
        removed["potential_contaminant"] = 0

    if site_only_col:
        mask = as_bool(df[site_only_col])
        removed["only_identified_by_site"] = int(mask.sum())
        df = df.loc[~mask].copy()
    else:
        removed["only_identified_by_site"] = 0

    if majority_id_col:
        missing_id = df[majority_id_col].fillna("").astype(str).str.strip().eq("")
        removed["missing_prey_id"] = int(missing_id.sum())
        df = df.loc[~missing_id].copy()
    else:
        removed["missing_prey_id"] = 0

    kept_samples = load_kept_samples(annotation_path)
    quant_cols = detect_quant_columns(df.columns, kept_samples)

    if quant_cols:
        quant = df[quant_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        zero_mask = quant.sum(axis=1).eq(0)
        removed["all_zero_or_missing_in_selected_quant_columns"] = int(zero_mask.sum())
        df = df.loc[~zero_mask].copy()
    else:
        removed["all_zero_or_missing_in_selected_quant_columns"] = 0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)

    summary = {
        "input": str(input_path),
        "annotation": str(annotation_path),
        "initial_rows": int(initial_rows),
        "remaining_rows": int(len(df)),
        "removed": removed,
        "selected_quant_columns": quant_cols,
        "selected_sample_count": len(kept_samples),
    }
    log_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(f"已输出清洗结果: {output_path}")
    print(f"已输出过滤统计: {log_path}")


if __name__ == "__main__":
    main()
