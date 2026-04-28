from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def to_float(value: str) -> float | None:
    try:
        return float(value)
    except Exception:
        return None


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def save_hist(scores: list[float], out: Path, title: str) -> None:
    plt.figure(figsize=(8, 5))
    plt.hist(scores, bins=min(40, max(10, int(math.sqrt(len(scores))))), color="#1f77b4", edgecolor="white")
    plt.title(title)
    plt.xlabel("score")
    plt.ylabel("count")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def save_scatter(xs: list[float], ys: list[float], out: Path, title: str, xlabel: str, ylabel: str) -> None:
    plt.figure(figsize=(7, 5))
    plt.scatter(xs, ys, s=20, alpha=0.7, color="#d62728")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def save_bar(labels: list[str], values: list[float], out: Path, title: str) -> None:
    plt.figure(figsize=(10, 6))
    plt.barh(labels[::-1], values[::-1], color="#2ca02c")
    plt.title(title)
    plt.xlabel("score")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def save_heatmap(matrix_path: Path, out: Path, title: str) -> None:
    header, rows = read_rows(matrix_path)
    if len(header) < 2:
        return
    sample_ids = header[1:]
    data = []
    labels = []
    for row in rows:
        labels.append(row.get(header[0], ""))
        data.append([to_float(row.get(sid, "")) or 0.0 for sid in sample_ids])
    if not data:
        return
    plt.figure(figsize=(max(8, len(sample_ids) * 0.4), max(6, len(labels) * 0.2)))
    plt.imshow(data, aspect="auto", interpolation="nearest", cmap="viridis")
    plt.colorbar(label="value")
    plt.title(title)
    plt.xticks(range(len(sample_ids)), sample_ids, rotation=90, fontsize=8)
    plt.yticks(range(len(labels)), labels, fontsize=8)
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()


def main() -> int:
    parser = argparse.ArgumentParser(description="生成 HGSCore 常用图表")
    parser.add_argument("--score-table", required=True, type=Path)
    parser.add_argument("--matrix", type=Path, default=None)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--score-column", default="score")
    parser.add_argument("--bait-column", default="bait_name")
    parser.add_argument("--support-column", default="support_count")
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()

    header, rows = read_rows(args.score_table)
    if args.score_column not in header:
        raise SystemExit(f"未找到得分列: {args.score_column}")

    ensure_dir(args.output_dir)
    scored = []
    for row in rows:
        score = to_float(row.get(args.score_column, ""))
        if score is None:
            continue
        support = to_float(row.get(args.support_column, "")) or 0.0
        bait = row.get(args.bait_column, "") or "NA"
        scored.append((score, support, bait, row))

    scores = [x[0] for x in scored]
    supports = [x[1] for x in scored]
    if scores:
        save_hist(scores, args.output_dir / "HGSCore_分数分布.png", "HGSCore 分数分布")
    if scores and supports:
        save_scatter(scores, supports, args.output_dir / "HGSCore_分数_vs_支持度.png", "HGSCore 分数与支持度", "score", "support")

    top = sorted(scored, key=lambda x: x[0], reverse=True)[: args.top_n]
    if top:
        labels = [f"{row.get(args.bait_column, 'NA')}::{row.get('prey_id', row.get('prey', ''))}" for _, _, _, row in top]
        values = [x[0] for x in top]
        save_bar(labels, values, args.output_dir / "HGSCore_TopHit.png", "HGSCore Top Hit")

    if args.matrix and args.matrix.exists():
        save_heatmap(args.matrix, args.output_dir / "HGSCore_矩阵热图.png", "HGSCore 矩阵热图")

    note = args.output_dir / "可视化说明.md"
    note.write_text(
        "\n".join(
            [
                "# HGSCore 可视化说明",
                "",
                "- `HGSCore_分数分布.png`: 看分数是否长尾",
                "- `HGSCore_分数_vs_支持度.png`: 看高分是否依赖单次或低支持证据",
                "- `HGSCore_TopHit.png`: 看前列候选",
                "- `HGSCore_矩阵热图.png`: 看样本间模式和重复一致性",
                "",
                "不确定点：若你的结果表没有 `score/support_count/bait_name` 这些列名，需要通过参数改列名。",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
