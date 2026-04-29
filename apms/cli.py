from __future__ import annotations

import argparse
import html
import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from run_secondary_workflows import (
    Paths,
    build_bait_statistics,
    find_prefixed_dir,
    load_and_clean,
    load_annotation,
    run_workflows,
)


def read_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Config must be a mapping: {path}")
    return data


def project_root_from_config(config_path: Path) -> Path:
    if config_path.parent.name == "config":
        return config_path.parent.parent.resolve()
    return Path.cwd().resolve()


def resolve_project_path(root: Path, value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (root / path).resolve()


def git_commit(root: Path) -> str:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()
    except Exception:
        return "unknown"


def make_run_id(project_name: str) -> str:
    safe = "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in project_name).strip("_")
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{stamp}_{safe or 'apms'}"


def copy_if_exists(src: Path, dest: Path) -> None:
    if not src.exists():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)


def copy_algorithm_figures(root: Path, figures_dir: Path) -> dict[str, str]:
    figure_specs = [
        ("CompPASS", "01_CompPASS"),
        ("HGSCore", "03_HGSCore"),
        ("CS Score", "04_CS_Score"),
        ("CRAPome", "05_CRAPome"),
        ("PPIrank", "06_PPIrank"),
    ]
    copied: dict[str, str] = {}
    for display_name, dir_prefix in figure_specs:
        try:
            algorithm_dir = find_prefixed_dir(root, dir_prefix)
            source_dir = find_prefixed_dir(algorithm_dir, "7.")
        except FileNotFoundError:
            continue
        safe_prefix = display_name.lower().replace(" ", "_")
        for src in sorted(source_dir.rglob("*.png")):
            dest = figures_dir / f"{safe_prefix}_{src.name}"
            copy_if_exists(src, dest)
            if dest.exists():
                copied[f"{display_name}: {src.stem}"] = str(dest)
    return copied


def write_standard_table(path: Path, cleaned: pd.DataFrame, annotation: pd.DataFrame) -> dict[str, int]:
    path.parent.mkdir(parents=True, exist_ok=True)
    quant_cols = annotation["source_column"].tolist()
    row_count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(
            "\t".join(
                [
                    "sample_id",
                    "bait_name",
                    "prey_id",
                    "prey_gene",
                    "quant_type",
                    "quant_value",
                    "replicate_group",
                    "condition",
                    "is_control",
                    "source_column",
                ]
            )
            + "\n"
        )
        base = cleaned[["prey_id", "prey_gene"]].copy()
        for _, sample in annotation.iterrows():
            col = sample["source_column"]
            if col not in cleaned.columns:
                continue
            values = pd.to_numeric(cleaned[col], errors="coerce").fillna(0)
            positive = values.gt(0)
            if not positive.any():
                continue
            sample_rows = base.loc[positive].copy()
            sample_rows.insert(0, "sample_id", sample["sample_id"])
            sample_rows.insert(1, "bait_name", sample["bait_name"])
            sample_rows["quant_type"] = col.split(" ", 1)[0] if " " in col else "Intensity"
            sample_rows["quant_value"] = values.loc[positive].to_numpy()
            sample_rows["replicate_group"] = sample.get("replicate_group", sample["bait_name"])
            sample_rows["condition"] = sample.get("condition", "")
            sample_rows["is_control"] = sample.get("is_control", "")
            sample_rows["source_column"] = col
            sample_rows.to_csv(handle, sep="\t", index=False, header=False)
            row_count += len(sample_rows)
    return {"standard_table_rows": int(row_count), "quant_columns": int(len(quant_cols))}


def summarize_table(path: Path, score_col: str | None = None) -> dict[str, Any]:
    if not path.exists():
        return {"exists": False}
    df = pd.read_csv(path, sep="\t")
    summary: dict[str, Any] = {"exists": True, "rows": int(len(df)), "columns": list(df.columns)}
    if score_col and score_col in df.columns and len(df):
        summary["score_min"] = float(df[score_col].min())
        summary["score_max"] = float(df[score_col].max())
    return summary


def build_report(run_dir: Path, run_summary: dict[str, Any], copied: dict[str, str]) -> None:
    tables = run_summary.get("tables", {})
    figures = run_summary.get("figures", {})
    cleaning = run_summary.get("secondary_workflow", {}).get("cleaning", {})
    outputs = run_summary.get("secondary_workflow", {}).get("outputs", {})
    annotation = run_summary.get("secondary_workflow", {}).get("annotation", {})

    table_rows = "\n".join(
        f"<tr><td>{html.escape(name)}</td><td>{info.get('rows', 'missing')}</td><td><a href='tables/{html.escape(Path(path).name)}'>download</a></td></tr>"
        for name, (path, info) in tables.items()
    )
    figure_cards = "\n".join(
        f"<figure><img src='figures/{html.escape(Path(path).name)}' alt='{html.escape(name)}'><figcaption>{html.escape(name)}</figcaption></figure>"
        for name, path in figures.items()
    )
    removed = cleaning.get("removed", {})
    removed_items = "".join(f"<li>{html.escape(str(k))}: {v}</li>" for k, v in removed.items())
    output_items = "".join(f"<li>{html.escape(str(k))}: {v}</li>" for k, v in outputs.items())

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>AP-MS report - {html.escape(run_summary['run_id'])}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 32px; color: #1f2933; }}
    h1, h2 {{ color: #102a43; }}
    table {{ border-collapse: collapse; width: 100%; margin: 16px 0; }}
    th, td {{ border: 1px solid #d9e2ec; padding: 8px; text-align: left; }}
    th {{ background: #f0f4f8; }}
    img {{ max-width: 720px; width: 100%; border: 1px solid #d9e2ec; }}
    figure {{ margin: 20px 0; }}
    code {{ background: #f0f4f8; padding: 2px 4px; }}
  </style>
</head>
<body>
  <h1>AP-MS Analysis Report</h1>
  <p><b>Run ID:</b> <code>{html.escape(run_summary['run_id'])}</code></p>
  <p><b>Project:</b> {html.escape(run_summary['project_name'])}</p>
  <p><b>Git commit:</b> <code>{html.escape(run_summary['git_commit'])}</code></p>
  <p><b>Created at:</b> {html.escape(run_summary['created_at'])}</p>

  <h2>Inputs</h2>
  <ul>
    <li>proteinGroups: <code>{html.escape(run_summary['inputs']['protein_groups'])}</code></li>
    <li>sample annotation: <code>{html.escape(run_summary['inputs']['sample_annotation'])}</code></li>
  </ul>

  <h2>Cleaning</h2>
  <p>Input rows: {cleaning.get('input_rows')} | Cleaned rows: {cleaning.get('cleaned_rows')}</p>
  <p>Kept samples: {annotation.get('kept_samples')} | Baits: {annotation.get('bait_count')}</p>
  <ul>{removed_items}</ul>

  <h2>Outputs</h2>
  <ul>{output_items}</ul>
  <table>
    <thead><tr><th>Table</th><th>Rows</th><th>File</th></tr></thead>
    <tbody>{table_rows}</tbody>
  </table>

  <h2>Figures</h2>
  {figure_cards}

  <h2>Notes</h2>
  <p>CRAPome currently uses project-level background frequency unless an external CRAPome reference table is added later.</p>
</body>
</html>
"""
    (run_dir / "report.html").write_text(html_text, encoding="utf-8")


def run_command(args: argparse.Namespace) -> int:
    config_path = Path(args.config).resolve()
    config = read_config(config_path)
    root = project_root_from_config(config_path)

    project_name = config.get("project", {}).get("name", "apms")
    run_id = args.run_id or make_run_id(project_name)
    runs_dir = resolve_project_path(root, config.get("output", {}).get("runs_dir", "runs"))
    run_dir = runs_dir / run_id
    tables_dir = run_dir / "tables"
    figures_dir = run_dir / "figures"
    logs_dir = run_dir / "logs"
    for directory in [tables_dir, figures_dir, logs_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    input_cfg = config.get("input", {})
    paths = Paths(
        root=root,
        protein_groups=resolve_project_path(root, input_cfg.get("protein_groups", "proteinGroups.txt")),
        annotation=resolve_project_path(
            root,
            input_cfg.get(
                "sample_annotation",
                "00_统一数据准备/2.样本注释/sample_annotation_master.tsv",
            ),
        ),
    )

    params = config.get("parameters", {})
    top_n = int(args.top_n_per_bait or params.get("top_n_per_bait", 100))
    contaminant_frequency = float(
        args.contaminant_frequency
        if args.contaminant_frequency is not None
        else params.get("contaminant_frequency", 0.25)
    )

    secondary_summary = run_workflows(paths, top_n=top_n, contaminant_frequency=contaminant_frequency)

    annotation = load_annotation(paths.annotation)
    cleaned, annotation, clean_stats = load_and_clean(paths, annotation)
    standard_summary: dict[str, int] = {}
    if bool(params.get("write_standard_table", True)):
        standard_summary = write_standard_table(tables_dir / "standard_apms_table.tsv", cleaned, annotation)

    copied_tables = {
        "Standard interactions": (
            find_prefixed_dir(root, "07_") / "结果输出" / "standard_interactions.tsv",
            "score",
        ),
        "CompPASS candidates": (
            root / "01_CompPASS_互作评分" / "6.CompPASS评分结果" / "compass_candidates.tsv",
            "compass_scorewd",
        ),
        "HGSCore candidates": (
            root / "03_HGSCore_辅助评分" / "6.HGSCore评分结果" / "hgscore_candidates.tsv",
            "hgs_score",
        ),
        "CS Score candidates": (
            root / "04_CS_Score_辅助评分" / "6.CS_Score评分结果" / "cs_score_candidates.tsv",
            "cs_score",
        ),
        "CRAPome background": (
            root / "05_CRAPome_污染过滤" / "6.CRAPome过滤结果" / "crapome_background_frequency.tsv",
            "crapome_filter_score",
        ),
        "CRAPome filtered": (
            root / "05_CRAPome_污染过滤" / "6.CRAPome过滤结果" / "crapome_filtered.tsv",
            "crapome_filter_score",
        ),
        "PPIrank ranked": (
            root / "06_PPIrank_网络补充排序" / "6.PPIrank排序结果" / "ppirank_ranked.tsv",
            "ppirank_score",
        ),
    }
    table_summaries: dict[str, tuple[str, dict[str, Any]]] = {}
    for name, (src, score_col) in copied_tables.items():
        dest = tables_dir / src.name
        copy_if_exists(src, dest)
        table_summaries[name] = (str(dest), summarize_table(dest, score_col))

    copied_figures = {
        "HGSCore score distribution": root / "03_HGSCore_辅助评分" / "7.可视化结果" / "hgscore_score_distribution.png",
        "CS Score distribution": root / "04_CS_Score_辅助评分" / "7.可视化结果" / "cs_score_distribution.png",
        "Background frequency distribution": root / "05_CRAPome_污染过滤" / "7.可视化结果" / "background_frequency_distribution.png",
        "PPIrank score distribution": root / "06_PPIrank_网络补充排序" / "7.可视化结果" / "ppirank_score_distribution.png",
        "Top baits by candidate count": root / "06_PPIrank_网络补充排序" / "7.可视化结果" / "top_baits_by_candidate_count.png",
    }
    figure_summaries: dict[str, str] = {}
    for name, src in copied_figures.items():
        dest = figures_dir / src.name
        copy_if_exists(src, dest)
        if dest.exists():
            figure_summaries[name] = str(dest)
    figure_summaries.update(copy_algorithm_figures(root, figures_dir))

    shutil.copy2(config_path, run_dir / "config.yaml")
    provenance = {
        "run_id": run_id,
        "project_name": project_name,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "git_commit": git_commit(root),
        "inputs": {
            "protein_groups": str(paths.protein_groups),
            "sample_annotation": str(paths.annotation),
        },
        "parameters": {
            "top_n_per_bait": top_n,
            "contaminant_frequency": contaminant_frequency,
        },
        "standard_table": standard_summary,
        "cleaning": clean_stats,
        "secondary_workflow": secondary_summary,
    }
    (run_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    run_summary = {
        **provenance,
        "tables": table_summaries,
        "figures": figure_summaries,
    }
    build_report(run_dir, run_summary, {})
    print(json.dumps({"run_id": run_id, "run_dir": str(run_dir), "report": str(run_dir / "report.html")}, indent=2, ensure_ascii=False))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="AP-MS analyzer command line interface")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run the reproducible AP-MS workflow")
    run_parser.add_argument("--config", default="config/project.yaml", help="Project YAML config")
    run_parser.add_argument("--run-id", default=None, help="Optional run id")
    run_parser.add_argument("--top-n-per-bait", type=int, default=None, help="Override top candidates per bait")
    run_parser.add_argument("--contaminant-frequency", type=float, default=None, help="Override contaminant frequency threshold")
    run_parser.set_defaults(func=run_command)

    args = parser.parse_args(argv)
    return args.func(args)
