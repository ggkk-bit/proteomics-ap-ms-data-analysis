"""CRAPome 结果筛选说明与阈值建议。

这个脚本只负责把筛选口径写清楚，便于后续统一执行。
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path


@dataclass(frozen=True)
class Thresholds:
    min_replicates: int = 2
    max_crapome_frequency: float = 0.05
    min_spectral_count_sum: int = 2
    require_bait_enrichment: bool = True


DEFAULT_THRESHOLDS = Thresholds()


def build_explain_text() -> str:
    t = asdict(DEFAULT_THRESHOLDS)
    return "\n".join(
        [
            "CRAPome 结果筛选说明",
            "",
            "推荐阈值:",
            f"- 最少重复命中数: {t['min_replicates']}",
            f"- CRAPome 背景频率上限: {t['max_crapome_frequency']}",
            f"- bait 端谱数总和下限: {t['min_spectral_count_sum']}",
            f"- 是否要求 bait 富集于对照: {t['require_bait_enrichment']}",
            "",
            "筛选逻辑:",
            "- 先去掉公开库里高频出现的常见污染物",
            "- 再看本项目 bait 是否在多个重复中稳定出现",
            "- 最后结合对照中的谱数是否明显偏高来决定是否保留",
            "",
            "阈值依据:",
            "- CRAPome 的核心价值是背景频率，不适合只看单一绝对分数",
            "- AP-MS 结果对 bait、重复数、对照数高度敏感，阈值需要按项目调整",
            "- 经验上，重复一致性和对照低频比单次高谱数更重要",
        ]
    )


def main() -> None:
    print(build_explain_text())


if __name__ == "__main__":
    main()

