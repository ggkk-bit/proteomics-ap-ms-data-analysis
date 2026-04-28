from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def to_float(value: str) -> float | None:
    try:
        return float(value)
    except Exception:
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description="筛选 HGSCore 结果")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--score-column", default="score")
    parser.add_argument("--quantile", type=float, default=0.95)
    parser.add_argument("--support-column", default="support_count")
    parser.add_argument("--min-support", type=int, default=2)
    parser.add_argument("--crapome-column", default="crapome_frequency")
    parser.add_argument("--max-crapome", type=float, default=0.1)
    args = parser.parse_args()

    header, rows = read_rows(args.input)
    if args.score_column not in header:
        raise SystemExit(f"未找到得分列: {args.score_column}")
    has_support = args.support_column in header
    has_crapome = args.crapome_column in header

    scored = []
    for row in rows:
        score = to_float(row.get(args.score_column, ""))
        if score is not None:
            scored.append((score, row))

    if not scored:
        raise SystemExit("没有可用于筛选的分数")

    scores = sorted(score for score, _ in scored)
    idx = max(0, min(len(scores) - 1, int(round((len(scores) - 1) * args.quantile))))
    threshold = scores[idx]

    kept = []
    for score, row in scored:
        if score < threshold:
            continue
        if has_support:
            support = to_float(row.get(args.support_column, "")) or 0
            if support < args.min_support:
                continue
        if has_crapome:
            crapome = to_float(row.get(args.crapome_column, ""))
            if crapome is not None and crapome > args.max_crapome:
                continue
        kept.append(row)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept)

    log_path = args.output.with_suffix(".说明.md")
    log_path.write_text(
        "\n".join(
            [
                "# HGSCore 结果筛选说明",
                "",
                f"- 输入: `{args.input}`",
                f"- 输出: `{args.output}`",
                f"- 分位数阈值: {args.quantile}",
                f"- 实际分数阈值: {threshold}",
                f"- 最小支持数: {args.min_support}",
                f"- CRAPome 上限: {args.max_crapome}",
                f"- 是否存在支持度列: {has_support}",
                f"- 是否存在 CRAPome 列: {has_crapome}",
                "",
                "## 筛选逻辑",
                "- 先按分数分位数截断",
                "- 再要求至少有一定重复/支持度",
                "- 若存在 CRAPome 频率列，则继续排除高频污染物",
                "",
                "## 依据",
                "- HGSCore 没有统一通用阈值时，分位数筛选比硬编码绝对分数更稳妥",
                "- AP-MS 文献习惯通常会结合重复支持和污染物背景一起筛选",
                "",
                "## 不确定点",
                "- 如果你的结果表不是 `score/support_count/crapome_frequency` 这套列名，需要在参数里改列名",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    print(args.output)
    print(log_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
