#!/usr/bin/env Rscript
# PPIrank 数据清洗脚本
# 优先复用 00_统一数据准备 的清洗结果；若不存在，则从本地 proteinGroups.txt 执行基础清洗。

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) { prefix <- paste0("--", name, "="); hit <- args[startsWith(args, prefix)]; if (length(hit) == 0) return(default); sub(prefix, "", hit[[1]], fixed = TRUE) }
ensure_dir <- function(path) { dir_path <- dirname(path); if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE) }
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
write_tsv <- function(x, path) { ensure_dir(path); write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "") }

unified_cleaned <- get_arg("unified_cleaned", "../00_统一数据准备/3.清洗后数据/proteinGroups_cleaned.tsv")
input_file <- get_arg("input", "./1.原始输入/proteinGroups.txt")
output_file <- get_arg("output", "./3.数据清洗与预处理/proteinGroups_cleaned.tsv")
log_file <- get_arg("log", "./5.质控/01_数据清洗记录.txt")

cat("PPIrank 数据清洗\n================\n\n")
if (file.exists(unified_cleaned)) {
  dat <- read_tsv(unified_cleaned); source_type <- "统一数据准备"; n_original <- nrow(dat); n_cleaned <- nrow(dat)
  cat("发现统一准备层清洗结果，直接复用: ", unified_cleaned, "\n", sep = "")
} else if (file.exists(input_file)) {
  dat <- read_tsv(input_file); n_original <- nrow(dat)
  # MaxQuant 常见污染、反库和仅位点鉴定过滤；列不存在时自动跳过。
  if ("Reverse" %in% names(dat)) dat <- dat[is.na(dat$Reverse) | dat$Reverse == "", , drop = FALSE]
  if ("Potential contaminant" %in% names(dat)) dat <- dat[is.na(dat$`Potential contaminant`) | dat$`Potential contaminant` == "", , drop = FALSE]
  if ("Only identified by site" %in% names(dat)) dat <- dat[is.na(dat$`Only identified by site`) | dat$`Only identified by site` == "", , drop = FALSE]
  n_cleaned <- nrow(dat); source_type <- "本地 proteinGroups 清洗"
} else stop("找不到输入文件。请提供 --input= 或 --unified_cleaned=。")

write_tsv(dat, output_file)
writeLines(c("PPIrank 数据清洗记录", paste("时间:", Sys.time()), paste("数据来源:", source_type), paste("原始行数:", n_original), paste("清洗后行数:", n_cleaned), paste("删除行数:", n_original - n_cleaned), paste("输出:", output_file)), log_file, useBytes = TRUE)
cat("输出: ", output_file, "\n", sep = "")
