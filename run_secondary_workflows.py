#!/usr/bin/env python3
"""Run secondary AP-MS workflows from MaxQuant proteinGroups output.

This script is intentionally conservative:
- one shared MaxQuant cleaning pass is used by all downstream outputs;
- bait/sample relationships come from sample_annotation_master.tsv;
- outputs are bounded to top candidates per bait so a full AP-MS project stays usable;
- every stage writes a validation report with row counts and assumptions.

The CRAPome stage uses project-level background frequency when no external
CRAPome reference table is provided. It is a reproducible contaminant filter,
not a replacement for uploading spectral-count data to the CRAPome service.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TRUE_VALUES = {"+", "1", "true", "yes", "y", "t"}
ID_COLUMNS = [
    "Protein IDs",
    "Majority protein IDs",
    "Protein names",
    "Gene names",
    "Peptides",
    "Razor + unique peptides",
    "Unique peptides",
    "Mol. weight [kDa]",
    "Sequence length",
    "Length",
    "Only identified by site",
    "Reverse",
    "Potential contaminant",
]


@dataclass
class Paths:
    root: Path
    protein_groups: Path
    annotation: Path


def truthy(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.lower().isin(TRUE_VALUES)


def first_token(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    return text.split(";")[0].strip()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_annotation(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing sample annotation: {path}")
    annotation = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"sample_id", "bait_name", "keep"}
    missing = sorted(required - set(annotation.columns))
    if missing:
        raise ValueError(f"Annotation is missing required columns: {', '.join(missing)}")

    annotation = annotation.copy()
    annotation["sample_id"] = annotation["sample_id"].astype(str).str.strip()
    annotation["bait_name"] = annotation["bait_name"].astype(str).str.strip()
    annotation["keep"] = annotation["keep"].astype(str).str.strip().str.lower()
    annotation = annotation[annotation["keep"].isin(["y", "yes", "1", "true"])].copy()
    annotation = annotation[annotation["sample_id"].ne("")].copy()
    annotation.loc[annotation["bait_name"].eq(""), "bait_name"] = annotation.loc[
        annotation["bait_name"].eq(""), "sample_id"
    ]

    if "source_column" not in annotation.columns:
        annotation["source_column"] = "Intensity " + annotation["sample_id"]
    annotation["source_column"] = annotation["source_column"].astype(str).str.strip()
    annotation.loc[annotation["source_column"].eq(""), "source_column"] = (
        "Intensity " + annotation.loc[annotation["source_column"].eq(""), "sample_id"]
    )
    if "replicate_group" not in annotation.columns:
        annotation["replicate_group"] = annotation["bait_name"]
    annotation["replicate_group"] = annotation["replicate_group"].astype(str).str.strip()
    annotation.loc[annotation["replicate_group"].eq(""), "replicate_group"] = annotation.loc[
        annotation["replicate_group"].eq(""), "bait_name"
    ]
    if "is_control" not in annotation.columns:
        annotation["is_control"] = "N"
    annotation["is_control"] = annotation["is_control"].astype(str).str.strip().str.upper()
    return annotation.drop_duplicates(subset=["sample_id", "source_column"]).reset_index(drop=True)


def read_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return handle.readline().rstrip("\n\r").split("\t")


def load_and_clean(paths: Paths, annotation: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    header = read_header(paths.protein_groups)
    quant_cols = [c for c in annotation["source_column"].tolist() if c in header]
    if not quant_cols:
        fallback = ["Intensity " + sid for sid in annotation["sample_id"].tolist()]
        quant_cols = [c for c in fallback if c in header]
    if not quant_cols:
        raise ValueError("No quantification columns from annotation were found in proteinGroups.txt")

    missing_quant = sorted(set(annotation["source_column"]) - set(quant_cols))
    meta_cols = [c for c in ID_COLUMNS if c in header]
    usecols = list(dict.fromkeys(meta_cols + quant_cols))
    pg = pd.read_csv(paths.protein_groups, sep="\t", dtype=str, usecols=usecols).fillna("")

    stats: dict[str, object] = {
        "input_rows": int(len(pg)),
        "quant_columns_requested": int(annotation["source_column"].nunique()),
        "quant_columns_found": int(len(quant_cols)),
        "missing_quant_columns": missing_quant[:50],
    }

    keep = pd.Series(True, index=pg.index)
    removed: dict[str, int] = {}
    for col in ["Reverse", "Potential contaminant", "Only identified by site"]:
        if col in pg.columns:
            mask = truthy(pg[col])
            removed[col] = int(mask.sum())
            keep &= ~mask
        else:
            removed[col] = 0

    id_source = "Majority protein IDs" if "Majority protein IDs" in pg.columns else "Protein IDs"
    pg["prey_id"] = pg[id_source].map(first_token)
    pg["prey_gene"] = pg["Gene names"].map(first_token) if "Gene names" in pg.columns else pg["prey_id"]
    pg.loc[pg["prey_gene"].eq(""), "prey_gene"] = pg.loc[pg["prey_gene"].eq(""), "prey_id"]
    missing_id = pg["prey_id"].eq("")
    keep &= ~missing_id
    removed["Missing protein ID"] = int(missing_id.sum())

    quant = pg[quant_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    all_zero = quant.le(0).all(axis=1)
    keep &= ~all_zero
    removed["All selected quant values zero"] = int(all_zero.sum())

    cleaned = pg.loc[keep].copy()
    quant = quant.loc[keep].copy()
    cleaned[quant_cols] = quant
    stats["removed"] = removed
    stats["cleaned_rows"] = int(len(cleaned))
    stats["cleaned_fraction"] = float(len(cleaned) / len(pg)) if len(pg) else 0.0
    stats["cleaning_rules"] = [
        "Drop MaxQuant Reverse marked rows",
        "Drop MaxQuant Potential contaminant marked rows",
        "Drop MaxQuant Only identified by site marked rows",
        "Drop rows without a stable protein ID",
        "Drop rows with zero abundance across all kept samples",
    ]
    return cleaned, annotation[annotation["source_column"].isin(quant_cols)].copy(), stats


def build_bait_statistics(cleaned: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    quant_cols = annotation["source_column"].tolist()
    matrix = cleaned[quant_cols].to_numpy(dtype=np.float32, copy=True)
    prey = cleaned[["prey_id", "prey_gene"]].reset_index(drop=True)

    rows: list[pd.DataFrame] = []
    for bait, sub in annotation.groupby("bait_name", sort=True):
        idx = [quant_cols.index(c) for c in sub["source_column"].tolist()]
        sub_matrix = matrix[:, idx]
        present = sub_matrix > 0
        n_reps = max(1, sub_matrix.shape[1])
        mean_intensity = sub_matrix.mean(axis=1)
        max_intensity = sub_matrix.max(axis=1)
        present_count = present.sum(axis=1)
        positive = present_count > 0
        if not positive.any():
            continue
        part = prey.loc[positive].copy()
        part["bait_name"] = bait
        part["replicate_count"] = int(n_reps)
        part["present_count"] = present_count[positive].astype(int)
        part["replicate_fraction"] = present_count[positive] / n_reps
        part["mean_intensity"] = mean_intensity[positive]
        part["max_intensity"] = max_intensity[positive]
        rows.append(part)

    if not rows:
        raise ValueError("No positive bait-prey observations after cleaning")

    stats = pd.concat(rows, ignore_index=True)
    bait_count = max(1, stats["bait_name"].nunique())
    background = (
        stats.groupby("prey_id", as_index=False)
        .agg(
            background_bait_count=("bait_name", "nunique"),
            background_mean_intensity=("mean_intensity", "mean"),
        )
    )
    background["background_frequency"] = background["background_bait_count"] / bait_count
    stats = stats.merge(background, on="prey_id", how="left")
    stats["specificity"] = stats["mean_intensity"] / (
        stats["background_mean_intensity"].replace(0, np.nan)
    )
    stats["specificity"] = stats["specificity"].replace([np.inf, -np.inf], np.nan).fillna(0.0)
    stats["log_mean_intensity"] = np.log1p(stats["mean_intensity"].astype(float))
    return stats


def minmax(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").fillna(0.0)
    lo = float(values.min())
    hi = float(values.max())
    if hi <= lo:
        return pd.Series(np.zeros(len(values)), index=values.index)
    return (values - lo) / (hi - lo)


def top_per_bait(df: pd.DataFrame, score_col: str, n: int) -> pd.DataFrame:
    ranked = df.sort_values(["bait_name", score_col], ascending=[True, False])
    return ranked.groupby("bait_name", group_keys=False).head(n).reset_index(drop=True)


def write_tsv(path: Path, df: pd.DataFrame) -> None:
    ensure_dir(path.parent)
    df.to_csv(path, sep="\t", index=False)


def write_json(path: Path, data: dict) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def plot_hist(path: Path, values: Iterable[float], title: str, xlabel: str) -> None:
    ensure_dir(path.parent)
    plt.figure(figsize=(7, 4.5))
    plt.hist(list(values), bins=50, color="#2f6f9f", edgecolor="white")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def plot_bar(path: Path, labels: list[str], values: list[float], title: str, ylabel: str) -> None:
    ensure_dir(path.parent)
    plt.figure(figsize=(8, 4.8))
    plt.bar(labels, values, color="#607c3c")
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(path, dpi=160)
    plt.close()


def run_workflows(paths: Paths, top_n: int, contaminant_frequency: float) -> dict:
    annotation = load_annotation(paths.annotation)
    cleaned, annotation, clean_stats = load_and_clean(paths, annotation)
    bait_stats = build_bait_statistics(cleaned, annotation)

    common_cols = [
        "bait_name",
        "prey_id",
        "prey_gene",
        "mean_intensity",
        "max_intensity",
        "present_count",
        "replicate_count",
        "replicate_fraction",
        "background_bait_count",
        "background_frequency",
        "specificity",
    ]

    root = paths.root
    all_summary = {
        "input": {
            "protein_groups": str(paths.protein_groups),
            "annotation": str(paths.annotation),
        },
        "cleaning": clean_stats,
        "annotation": {
            "kept_samples": int(annotation["sample_id"].nunique()),
            "bait_count": int(annotation["bait_name"].nunique()),
        },
    }

    # HGSCore-like auxiliary scoring: length-normalized abundance proxy plus
    # replicate support and bait specificity.
    hgs = bait_stats.copy()
    weight_col = "Mol. weight [kDa]" if "Mol. weight [kDa]" in cleaned.columns else None
    if weight_col:
        hgs = hgs.merge(
            cleaned[["prey_id", weight_col]].drop_duplicates("prey_id"),
            on="prey_id",
            how="left",
        )
        hgs["mol_weight_kda"] = pd.to_numeric(hgs[weight_col], errors="coerce").fillna(0.0)
    else:
        hgs["mol_weight_kda"] = 0.0
    hgs["length_normalized_abundance"] = hgs["mean_intensity"] / hgs["mol_weight_kda"].replace(0, np.nan)
    hgs["length_normalized_abundance"] = hgs["length_normalized_abundance"].replace(
        [np.inf, -np.inf], np.nan
    ).fillna(hgs["mean_intensity"])
    hgs["hgs_score"] = (
        0.45 * minmax(np.log1p(hgs["length_normalized_abundance"]))
        + 0.35 * minmax(hgs["specificity"])
        + 0.20 * hgs["replicate_fraction"].clip(0, 1)
    )
    hgs_out = top_per_bait(hgs[common_cols + ["mol_weight_kda", "length_normalized_abundance", "hgs_score"]], "hgs_score", top_n)
    write_tsv(root / "03_HGSCore_辅助评分" / "6.HGSCore评分结果" / "hgscore_candidates.tsv", hgs_out)
    plot_hist(
        root / "03_HGSCore_辅助评分" / "7.可视化结果" / "hgscore_score_distribution.png",
        hgs_out["hgs_score"],
        "HGSCore candidate score distribution",
        "HGSCore-like score",
    )

    # CS Score-like co-membership scoring.
    cs = bait_stats.copy()
    cs["cs_score"] = (
        0.50 * cs["replicate_fraction"].clip(0, 1)
        + 0.30 * minmax(cs["specificity"])
        + 0.20 * minmax(cs["log_mean_intensity"])
    )
    cs_out = top_per_bait(cs[common_cols + ["cs_score"]], "cs_score", top_n)
    write_tsv(root / "04_CS_Score_辅助评分" / "6.CS_Score评分结果" / "cs_score_candidates.tsv", cs_out)
    plot_hist(
        root / "04_CS_Score_辅助评分" / "7.可视化结果" / "cs_score_distribution.png",
        cs_out["cs_score"],
        "CS Score candidate distribution",
        "CS-like score",
    )

    # CRAPome/background filter.
    crap = bait_stats.copy()
    crap["project_background_contaminant"] = crap["background_frequency"] >= contaminant_frequency
    crap["crapome_filter_score"] = 1.0 - crap["background_frequency"].clip(0, 1)
    crap_out = top_per_bait(
        crap[common_cols + ["project_background_contaminant", "crapome_filter_score"]],
        "crapome_filter_score",
        top_n,
    )
    crap_pass = crap_out[~crap_out["project_background_contaminant"]].copy()
    write_tsv(root / "05_CRAPome_污染过滤" / "6.CRAPome过滤结果" / "crapome_background_frequency.tsv", crap_out)
    write_tsv(root / "05_CRAPome_污染过滤" / "6.CRAPome过滤结果" / "crapome_filtered.tsv", crap_pass)
    plot_hist(
        root / "05_CRAPome_污染过滤" / "7.可视化结果" / "background_frequency_distribution.png",
        crap["background_frequency"],
        "Project background frequency distribution",
        "Fraction of baits where prey appears",
    )

    # PPIrank-style rank aggregation from the three upstream signals.
    ppi = bait_stats[common_cols].copy()
    score_map = hgs[["bait_name", "prey_id", "hgs_score"]].merge(
        cs[["bait_name", "prey_id", "cs_score"]], on=["bait_name", "prey_id"], how="inner"
    )
    score_map = score_map.merge(
        crap[["bait_name", "prey_id", "crapome_filter_score", "project_background_contaminant"]],
        on=["bait_name", "prey_id"],
        how="inner",
    )
    ppi = ppi.merge(score_map, on=["bait_name", "prey_id"], how="left").fillna(0)
    ppi["ppirank_score"] = (
        0.35 * ppi["hgs_score"]
        + 0.35 * ppi["cs_score"]
        + 0.20 * ppi["crapome_filter_score"]
        + 0.10 * ppi["replicate_fraction"].clip(0, 1)
    )
    ppi.loc[ppi["project_background_contaminant"].astype(bool), "ppirank_score"] *= 0.5
    ppi_out = top_per_bait(ppi[common_cols + ["hgs_score", "cs_score", "crapome_filter_score", "project_background_contaminant", "ppirank_score"]], "ppirank_score", top_n)
    write_tsv(root / "06_PPIrank_网络补充排序" / "6.PPIrank排序结果" / "ppirank_ranked.tsv", ppi_out)
    plot_hist(
        root / "06_PPIrank_网络补充排序" / "7.可视化结果" / "ppirank_score_distribution.png",
        ppi_out["ppirank_score"],
        "PPIrank score distribution",
        "PPIrank aggregate score",
    )

    top_baits = (
        ppi_out.groupby("bait_name")["prey_id"]
        .count()
        .sort_values(ascending=False)
        .head(20)
    )
    plot_bar(
        root / "06_PPIrank_网络补充排序" / "7.可视化结果" / "top_baits_by_candidate_count.png",
        top_baits.index.tolist(),
        top_baits.astype(float).tolist(),
        "Top baits by retained candidate count",
        "Candidate count",
    )

    all_summary["outputs"] = {
        "HGSCore_rows": int(len(hgs_out)),
        "CS_Score_rows": int(len(cs_out)),
        "CRAPome_background_rows": int(len(crap_out)),
        "CRAPome_filtered_rows": int(len(crap_pass)),
        "PPIrank_rows": int(len(ppi_out)),
        "top_n_per_bait": int(top_n),
        "contaminant_frequency_threshold": float(contaminant_frequency),
    }
    all_summary["assumptions"] = [
        "MaxQuant contaminants/reverse/site-only rows are removed before all scoring.",
        "Scores are project-level reproducible proxies for the four secondary analyses.",
        "External CRAPome reference filtering is not applied unless a reference table is added later.",
        "Outputs are ranked candidates, not raw MaxQuant replacement files.",
    ]
    write_json(root / "secondary_workflow_validation.json", all_summary)
    return all_summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Run HGSCore, CS Score, CRAPome, and PPIrank AP-MS workflows.")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--protein-groups", type=Path, default=None)
    parser.add_argument("--annotation", type=Path, default=None)
    parser.add_argument("--top-n-per-bait", type=int, default=100)
    parser.add_argument("--contaminant-frequency", type=float, default=0.25)
    args = parser.parse_args()

    root = args.root.resolve()
    paths = Paths(
        root=root,
        protein_groups=(args.protein_groups or root / "proteinGroups.txt").resolve(),
        annotation=(args.annotation or root / "00_统一数据准备" / "2.样本注释" / "sample_annotation_master.tsv").resolve(),
    )
    summary = run_workflows(paths, top_n=args.top_n_per_bait, contaminant_frequency=args.contaminant_frequency)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
