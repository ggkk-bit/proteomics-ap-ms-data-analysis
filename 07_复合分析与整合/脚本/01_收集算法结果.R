# 01 收集 6 个算法输出并复制到汇总目录
source("07_复合分析与整合/脚本/common_integration.R", encoding = "UTF-8")
cfg <- parse_args(list(out_dir = "07_复合分析与整合/1.算法结果汇总/原始结果"))
ensure_dir(cfg$out_dir)

# 如路径不存在，脚本会跳过并在日志中记录；后续可替换为真实输出路径
inputs <- data.frame(
  algorithm = c("CompPASS", "MiST", "HGSCore", "CS_Score", "CRAPome", "PPIrank"),
  path = c(
    "01_CompPASS_互作评分/6.CompPASS评分结果/原始结果/comppass_results.tsv",
    "02_MiST_互作评分/6.MiST评分结果/preprocessed_NoC_MAT_MIST.txt",
    "03_HGSCore_辅助评分/6.HGSCore评分结果/hgscore_filtered.tsv",
    "04_CS_Score_辅助评分/6.CS_Score评分结果/cs_score_filtered_pairs.tsv",
    "05_CRAPome_污染过滤/6.CRAPome过滤结果/crapome_filtered.tsv",
    "06_PPIrank_网络补充排序/6.PPIrank排序结果/ppirank_filtered.tsv"
  ),
  stringsAsFactors = FALSE
)

log <- inputs
log$exists <- file.exists(inputs$path)
log$copied_to <- file.path(cfg$out_dir, paste0(inputs$algorithm, "_raw.tsv"))
for (i in seq_len(nrow(inputs))) {
  if (log$exists[i]) file.copy(inputs$path[i], log$copied_to[i], overwrite = TRUE)
}
write_utf8_tsv(log, file.path(cfg$out_dir, "collect_log.tsv"))
cat("完成算法结果收集：", cfg$out_dir, "\n")
