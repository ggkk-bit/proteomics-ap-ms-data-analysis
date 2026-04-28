#!/usr/bin/env Rscript
# HGSCore 可视化脚本

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

# 参数配置
input_file <- get_arg("input", "./6.HGSCore评分结果/hgscore_filtered.tsv")
output_dir <- get_arg("output_dir", "./7.可视化")

cat("HGSCore 可视化脚本\n")
cat("==================\n\n")

# 读取筛选后的结果
if (!file.exists(input_file)) {
  stop("找不到输入文件: ", input_file)
}
dat <- read.delim(input_file, stringsAsFactors = FALSE)
cat(sprintf("✓ 读取数据: %d 个蛋白对\n", nrow(dat)))

# 确保输出目录存在
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1. HG 分数分布图
p1 <- ggplot(dat, aes(x = HG)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  labs(title = "HGScore 分布", x = "HG Score", y = "频数") +
  theme_minimal()

ggsave(file.path(output_dir, "01_HGScore_distribution.png"), p1, width = 8, height = 6)
cat("✓ 生成: 01_HGScore_distribution.png\n")

# 2. Top 20 蛋白对条形图
if (nrow(dat) > 0) {
  top20 <- head(dat, 20)
  top20$pair <- paste(top20[, 1], top20[, 2], sep = " - ")
  top20$pair <- factor(top20$pair, levels = rev(top20$pair))

  p2 <- ggplot(top20, aes(x = HG, y = pair)) +
    geom_col(fill = "coral") +
    labs(title = "Top 20 蛋白对 (按 HG Score)", x = "HG Score", y = "") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))

  ggsave(file.path(output_dir, "02_Top20_pairs.png"), p2, width = 10, height = 8)
  cat("✓ 生成: 02_Top20_pairs.png\n")
}

# 3. 统计摘要
summary_stats <- data.frame(
  指标 = c("蛋白对数量", "HG 最小值", "HG 最大值", "HG 中位数", "HG 平均值"),
  值 = c(
    nrow(dat),
    round(min(dat$HG, na.rm = TRUE), 2),
    round(max(dat$HG, na.rm = TRUE), 2),
    round(median(dat$HG, na.rm = TRUE), 2),
    round(mean(dat$HG, na.rm = TRUE), 2)
  )
)

write.table(summary_stats, file.path(output_dir, "00_summary_stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("✓ 生成: 00_summary_stats.tsv\n")

cat("\n可视化完成！\n")
