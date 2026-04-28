# 11 导出多 sheet Excel 报告
source("07_复合分析与整合/脚本/common_integration.R", encoding = "UTF-8")
load_pkgs(c("openxlsx"))
cfg <- parse_args(list(out = "07_复合分析与整合/8.导出结果/Excel报告/phase4_integration_report.xlsx"))
ensure_dir(dirname(cfg$out))
files <- list(
  标准化结果 = "07_复合分析与整合/1.算法结果汇总/标准化结果/standardize_summary.tsv",
  算法重叠矩阵 = "07_复合分析与整合/4.算法交叉验证/重叠分析/algorithm_overlap_matrix.tsv",
  一致性评分 = "07_复合分析与整合/4.算法交叉验证/一致性评分/pair_consistency_scores.tsv",
  STRING算法验证率 = "07_复合分析与整合/5.STRING验证对照/验证率统计/validation_rate_by_algorithm.tsv",
  STRING一致性验证率 = "07_复合分析与整合/5.STRING验证对照/验证率统计/validation_rate_by_support_count.tsv",
  新发现候选 = "07_复合分析与整合/5.STRING验证对照/新发现互作/candidate_novel_interactions.tsv",
  综合评分 = "07_复合分析与整合/6.综合评分/评分结果/integrated_interaction_scores.tsv",
  高置信互作 = "07_复合分析与整合/6.综合评分/高置信互作/high_confidence_interactions.tsv"
)
wb <- openxlsx::createWorkbook()
for (nm in names(files)) {
  openxlsx::addWorksheet(wb, nm)
  openxlsx::writeData(wb, nm, read_table_auto(files[[nm]]))
}
openxlsx::saveWorkbook(wb, cfg$out, overwrite = TRUE)
cat("完成 Excel 报告：", cfg$out, "\n")
