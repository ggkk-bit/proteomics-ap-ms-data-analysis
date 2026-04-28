#!/usr/bin/env Rscript
# CS_Score 数据清洗：优先复用统一准备层，否则从 proteinGroups.txt 清洗。

source(file.path(script_dir <- {
  cmd <- commandArgs(FALSE); fa <- cmd[startsWith(cmd, "--file=")]
  if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
}, "common_cs_score.R"))

args <- commandArgs(trailingOnly = TRUE)
message_header("CS_Score 数据清洗")

unified_cleaned <- get_arg(args, "unified_cleaned", "../../00_统一数据准备/3.清洗后数据/proteinGroups_cleaned.tsv")
input_file <- get_arg(args, "input", "./1.原始输入/proteinGroups.txt")
output_file <- get_arg(args, "output", "./3.数据清洗与预处理/proteinGroups_cleaned.tsv")
log_file <- get_arg(args, "log", "./5.质控/01_数据清洗记录.txt")

if (file.exists(unified_cleaned)) {
  dat <- read_tsv(unified_cleaned)
  source_type <- paste0("统一准备层: ", unified_cleaned)
  n_original <- nrow(dat)
  n_cleaned <- nrow(dat)
} else if (file.exists(input_file)) {
  dat <- read_tsv(input_file)
  source_type <- paste0("本地输入: ", input_file)
  n_original <- nrow(dat)
  if ("Reverse" %in% names(dat)) dat <- dat[is.na(dat$Reverse) | dat$Reverse == "", , drop = FALSE]
  if ("Potential contaminant" %in% names(dat)) dat <- dat[is.na(dat$`Potential contaminant`) | dat$`Potential contaminant` == "", , drop = FALSE]
  if ("Only identified by site" %in% names(dat)) dat <- dat[is.na(dat$`Only identified by site`) | dat$`Only identified by site` == "", , drop = FALSE]
  n_cleaned <- nrow(dat)
} else {
  stop("找不到输入文件。请提供 --input= 或先运行统一准备层。")
}

write_tsv(dat, output_file)
writeLines(c(
  "CS_Score 数据清洗记录",
  paste("时间:", Sys.time()),
  paste("来源:", source_type),
  paste("原始行数:", n_original),
  paste("清洗后行数:", n_cleaned),
  paste("删除行数:", n_original - n_cleaned)
), log_file, useBytes = TRUE)

cat("输出: ", output_file, "\n", sep = "")
