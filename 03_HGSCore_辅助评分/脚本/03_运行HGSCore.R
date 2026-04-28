#!/usr/bin/env Rscript
# HGSCore 运行脚本
# 调用 SMAD::HG() 函数计算 HGScore

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

# 参数配置
input_file <- get_arg("input", "./4.HGSCore输入文件/hgscore_input.tsv")
output_file <- get_arg("output", "./6.HGSCore评分结果/hgscore_results.tsv")
smad_path <- get_arg("smad_path", "./8.文档与参考/03_源码或软件/SMAD-master")

cat("HGSCore 运行脚本\n")
cat("==================\n\n")

# 加载 SMAD 包
if (dir.exists(smad_path)) {
  cat("✓ 从本地加载 SMAD 包:", smad_path, "\n")
  suppressPackageStartupMessages({
    library(devtools)
    load_all(smad_path, quiet = TRUE)
  })
} else {
  cat("✗ 未找到本地 SMAD 包，尝试从已安装包加载\n")
  if (!require("SMAD", quietly = TRUE)) {
    stop("错误: 找不到 SMAD 包\n  请安装: BiocManager::install('SMAD')\n  或指定本地路径: --smad_path=...")
  }
}

# 读取输入数据
if (!file.exists(input_file)) {
  stop("找不到输入文件: ", input_file)
}
dat_input <- read.delim(input_file, stringsAsFactors = FALSE)
cat(sprintf("✓ 读取输入数据: %d 行\n", nrow(dat_input)))

# 检查必需列
required_cols <- c("idRun", "idPrey", "countPrey", "lenPrey")
missing_cols <- setdiff(required_cols, names(dat_input))
if (length(missing_cols) > 0) {
  stop("输入数据缺少必需列: ", paste(missing_cols, collapse = ", "))
}

# 运行 HG 算法
cat("\n运行 HG 算法...\n")
dat_score <- HG(dat_input)
cat(sprintf("✓ 计算完成: %d 个蛋白对\n", nrow(dat_score)))

# 确保输出目录存在
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 保存结果
write.table(dat_score, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("✓ 结果已保存: %s\n", output_file))

# 输出统计信息
cat("\n统计信息:\n")
cat(sprintf("  HG 分数范围: %.2f - %.2f\n", min(dat_score$HG, na.rm = TRUE), max(dat_score$HG, na.rm = TRUE)))
cat(sprintf("  HG 分数中位数: %.2f\n", median(dat_score$HG, na.rm = TRUE)))
