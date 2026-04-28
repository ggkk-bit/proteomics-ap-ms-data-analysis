#!/usr/bin/env Rscript
# HGSCore 结果筛选脚本

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

# 参数配置
input_file <- get_arg("input", "./6.HGSCore评分结果/hgscore_results.tsv")
output_file <- get_arg("output", "./6.HGSCore评分结果/hgscore_filtered.tsv")
min_score <- as.numeric(get_arg("min_score", "3"))

cat("HGSCore 结果筛选脚本\n")
cat("====================\n\n")

# 读取结果
if (!file.exists(input_file)) {
  stop("找不到输入文件: ", input_file)
}
dat <- read.delim(input_file, stringsAsFactors = FALSE)
cat(sprintf("✓ 读取结果: %d 个蛋白对\n", nrow(dat)))

# 筛选高置信结果
dat_filtered <- dat[dat$HG >= min_score, ]
cat(sprintf("✓ 筛选阈值: HG >= %.2f\n", min_score))
cat(sprintf("✓ 筛选后: %d 个蛋白对 (%.1f%%)\n",
            nrow(dat_filtered),
            100 * nrow(dat_filtered) / nrow(dat)))

# 按 HG 分数降序排列
dat_filtered <- dat_filtered[order(-dat_filtered$HG), ]

# 保存结果
write.table(dat_filtered, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("✓ 筛选结果已保存: %s\n", output_file))

# 输出 Top 10
if (nrow(dat_filtered) > 0) {
  cat("\nTop 10 蛋白对:\n")
  top10 <- head(dat_filtered, 10)
  for (i in 1:nrow(top10)) {
    cat(sprintf("  %d. %s - %s (HG=%.2f)\n",
                i, top10[i, 1], top10[i, 2], top10$HG[i]))
  }
}
