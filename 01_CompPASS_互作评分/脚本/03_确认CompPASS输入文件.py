import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="检查 CompPASS 输入文件是否满足基本格式要求。")
    parser.add_argument("--input", default="4.CompPASS输入文件/最终输入/comppass_input.tsv")
    parser.add_argument("--report", default="5.质控/输入检查/CompPASS输入检查报告.txt")
    args = parser.parse_args()

    input_path = Path(args.input)
    report_path = Path(args.report)

    df = pd.read_csv(input_path, sep="\t")
    required = ["idRun", "idBait", "idPrey", "countPrey"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"缺少核心字段: {missing}")

    lines = []
    lines.append(f"输入文件: {input_path}")
    lines.append(f"总记录数: {len(df)}")
    lines.append(f"唯一 run 数: {df['idRun'].nunique()}")
    lines.append(f"唯一 bait 数: {df['idBait'].nunique()}")
    lines.append(f"唯一 prey 数: {df['idPrey'].nunique()}")

    numeric = pd.to_numeric(df["countPrey"], errors="coerce")
    if numeric.isna().any():
        raise ValueError("countPrey 中存在无法转为数值的记录。")
    if (numeric <= 0).any():
        raise ValueError("countPrey 中存在小于等于 0 的记录。")

    duplicated = df.duplicated(subset=["idRun", "idBait", "idPrey"]).sum()
    lines.append(f"重复的 run-bait-prey 记录数: {int(duplicated)}")

    bait_per_run = df.groupby("idRun")["idBait"].nunique()
    bad_runs = bait_per_run[bait_per_run > 1]
    if not bad_runs.empty:
        raise ValueError("发现同一个 idRun 对应多个 idBait，这不符合 CompPASS 输入要求。")

    replicate_summary = df[["idRun", "idBait"]].drop_duplicates().groupby("idBait").size().sort_values(ascending=False)
    lines.append("")
    lines.append("各 bait 的 run 数:")
    for bait, run_count in replicate_summary.items():
        lines.append(f"- {bait}: {int(run_count)}")

    warnings = []
    if df["idBait"].nunique() < 5:
        warnings.append("bait 总数少于 5，CompPASS 的背景频率估计可能不稳定。")
    singletons = replicate_summary[replicate_summary == 1]
    if not singletons.empty:
        warnings.append(f"有 {len(singletons)} 个 bait 只有 1 个 run，重复支持较弱。")
    if duplicated > 0:
        warnings.append("存在重复的 run-bait-prey 行，建议在正式评分前去重或合并。")

    lines.append("")
    lines.append("警告:")
    if warnings:
        for warning in warnings:
            lines.append(f"- {warning}")
    else:
        lines.append("- 未发现明显结构性警告。")

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"已输出检查报告: {report_path}")


if __name__ == "__main__":
    main()
