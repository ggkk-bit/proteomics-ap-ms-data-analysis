"""CRAPome 输入文件确认与说明。

用途：
- 检查原始输入是否齐全
- 验证 proteinGroups.txt 的关键列是否存在
- 输出一份中文检查报告，提醒用户 CRAPome 真正需要的上传格式

说明：
- 这个脚本不做任何结果筛选，只做输入确认。
- 如果只有 proteinGroups.txt，它仍然无法替代 CRAPome 的四列表格式。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REQUIRED_PG_COLUMNS = [
    "Protein IDs",
    "Majority protein IDs",
    "Gene names",
    "Reverse",
    "Potential contaminant",
    "Only identified by site",
]

ANNOTATION_ALIASES = [
    ("AP Name", "Bait Name"),
    ("ip_name", "bait_name"),
    ("sample_id", "bait_name"),
    ("样本ID", "bait_name"),
]


def build_report(protein_groups: Path, sample_annotation: Path | None) -> str:
    df = pd.read_csv(protein_groups, sep="\t", low_memory=False)
    missing = [c for c in REQUIRED_PG_COLUMNS if c not in df.columns]

    lines = [
        "CRAPome 输入确认报告",
        f"文件: {protein_groups}",
        f"总行数: {len(df)}",
    ]

    if missing:
        lines.append("缺失列: " + ", ".join(missing))
    else:
        lines.append("proteinGroups 关键列检查: 通过")

    lines.extend(
        [
            "",
            "CRAPome 实际需要的上传表:",
            "Bait Name / AP Name / Prey Name / Spectral Count",
            "",
            "当前结论:",
            "proteinGroups.txt 只能作为上游过滤与 ID 统一来源，不能单独生成完整 CRAPome 上传表。",
            "如果要转换成上传表，还需要样本级谱数信息，例如 evidence.txt 或等价的样本级计数矩阵。",
        ]
    )

    if sample_annotation is not None:
        annot = pd.read_csv(sample_annotation, sep="\t", low_memory=False)
        matched = next((pair for pair in ANNOTATION_ALIASES if set(pair).issubset(annot.columns)), None)
        lines.extend(
            [
                "",
                f"样本注释: {sample_annotation}",
                f"可识别映射: {matched[0]} -> AP Name, {matched[1]} -> Bait Name" if matched else "样本注释列名暂时无法映射到 CRAPome 需要的字段",
            ]
        )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="检查 CRAPome 输入文件")
    parser.add_argument("--protein-groups", required=True, type=Path)
    parser.add_argument("--sample-annotation", type=Path)
    parser.add_argument("--output-report", type=Path)
    args = parser.parse_args()

    if not args.protein_groups.exists():
        raise FileNotFoundError(f"找不到 proteinGroups 文件: {args.protein_groups}")

    if args.sample_annotation is not None and not args.sample_annotation.exists():
        raise FileNotFoundError(
            "找不到样本注释文件: "
            f"{args.sample_annotation}\n"
            "如果还没准备好正式样本注释表，请先去掉 --sample-annotation 参数；"
            "如果需要模板，请复制 ..\\2.样本信息与分组\\sample_annotation_template.tsv 后再填写。"
        )

    report = build_report(args.protein_groups, args.sample_annotation)
    print(report)
    if args.output_report:
        args.output_report.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
