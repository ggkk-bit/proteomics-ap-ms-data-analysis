#!/usr/bin/env Rscript
# HGSCore 格式转换脚本
# 将清洗后的数据转换为 SMAD::HG() 函数需要的格式

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

# 参数配置
cleaned_file <- get_arg("cleaned", "./3.数据清洗与预处理/proteinGroups_cleaned.tsv")
annotation_file <- get_arg("annotation", "../../00_统一数据准备/2.样本注释/sample_annotation_master.tsv")
output_file <- get_arg("output", "./4.HGSCore输入文件/hgscore_input.tsv")
quant_type <- get_arg("quant_type", "auto")

cat("HGSCore 格式转换脚本\n")
cat("====================\n\n")

# 读取清洗后的数据
if (!file.exists(cleaned_file)) {
  stop("找不到清洗后的数据文件: ", cleaned_file)
}
dat <- read.delim(cleaned_file, check.names = FALSE, stringsAsFactors = FALSE)
cat(sprintf("✓ 读取清洗数据: %d 行\n", nrow(dat)))

# 读取样本注释（可选）
if (file.exists(annotation_file)) {
  annot <- read.delim(annotation_file, stringsAsFactors = FALSE)
  cat(sprintf("✓ 读取样本注释: %d 个样本\n", nrow(annot)))
} else {
  cat("✗ 未找到样本注释文件，将使用所有定量列\n")
  annot <- NULL
}

# 自动检测定量类型
quant_prefixes <- c("Intensity ", "LFQ intensity ", "iBAQ ", "MS/MS count ")
if (quant_type == "auto") {
  for (prefix in quant_prefixes) {
    if (any(startsWith(names(dat), prefix))) {
      quant_type <- prefix
      break
    }
  }
}
cat(sprintf("✓ 使用定量类型: %s\n", quant_type))

# 提取定量列
quant_cols <- names(dat)[startsWith(names(dat), quant_type)]
if (length(quant_cols) == 0) {
  stop("未找到定量列，前缀: ", quant_type)
}
cat(sprintf("✓ 找到 %d 个定量列\n", length(quant_cols)))

# 提取蛋白ID和长度
protein_id_col <- if ("Majority protein IDs" %in% names(dat)) "Majority protein IDs" else "Protein IDs"
length_col <- if ("Sequence length" %in% names(dat)) "Sequence length" else "Length"

if (!protein_id_col %in% names(dat)) {
  stop("找不到蛋白ID列")
}
if (!length_col %in% names(dat)) {
  stop("找不到蛋白长度列")
}

# 转换为长表格式
long_dat <- dat %>%
  select(all_of(c(protein_id_col, length_col)), all_of(quant_cols)) %>%
  pivot_longer(
    cols = all_of(quant_cols),
    names_to = "sample",
    values_to = "intensity"
  ) %>%
  mutate(
    sample = sub(quant_type, "", sample, fixed = TRUE),
    intensity = as.numeric(intensity)
  ) %>%
  filter(!is.na(intensity) & intensity > 0)

# 转换为 HG() 需要的格式
hg_input <- long_dat %>%
  transmute(
    idRun = sample,
    idPrey = !!sym(protein_id_col),
    countPrey = intensity,
    lenPrey = !!sym(length_col)
  )

cat(sprintf("✓ 转换完成: %d 行\n", nrow(hg_input)))

# 确保输出目录存在
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 保存结果
write.table(hg_input, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("✓ 输出已保存: %s\n", output_file))
