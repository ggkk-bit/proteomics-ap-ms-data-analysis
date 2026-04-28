#!/usr/bin/env Rscript
# HGSCore 数据清洗脚本
# 优先复用统一准备层的清洗结果，找不到时从本地清洗

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

# 参数配置
unified_cleaned <- get_arg("unified_cleaned", "../../00_统一数据准备/3.清洗后数据/proteinGroups_cleaned.tsv")
local_input <- get_arg("input", "./1.原始输入/proteinGroups.txt")
output_file <- get_arg("output", "./3.数据清洗与预处理/proteinGroups_cleaned.tsv")
log_file <- get_arg("log", "./5.质控/数据清洗记录.txt")

cat("HGSCore 数据清洗脚本\n")
cat("====================\n\n")

# 优先使用统一准备层的清洗结果
if (file.exists(unified_cleaned)) {
  cat("✓ 发现统一准备层清洗结果，直接复用:", unified_cleaned, "\n")
  dat <- read.delim(unified_cleaned, check.names = FALSE, stringsAsFactors = FALSE)
  source_type <- "统一准备层"
} else if (file.exists(local_input)) {
  cat("✗ 未找到统一准备层结果，从本地清洗:", local_input, "\n")
  dat <- read.delim(local_input, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")

  # 执行清洗
  n_original <- nrow(dat)

  # 删除 Reverse
  if ("Reverse" %in% names(dat)) {
    dat <- dat[is.na(dat$Reverse) | dat$Reverse == "", ]
  }

  # 删除 Potential contaminant
  if ("Potential contaminant" %in% names(dat)) {
    dat <- dat[is.na(dat$`Potential contaminant`) | dat$`Potential contaminant` == "", ]
  }

  # 删除 Only identified by site
  if ("Only identified by site" %in% names(dat)) {
    dat <- dat[is.na(dat$`Only identified by site`) | dat$`Only identified by site` == "", ]
  }

  n_cleaned <- nrow(dat)
  cat(sprintf("  原始行数: %d\n", n_original))
  cat(sprintf("  清洗后行数: %d\n", n_cleaned))
  cat(sprintf("  删除行数: %d (%.1f%%)\n", n_original - n_cleaned, 100 * (n_original - n_cleaned) / n_original))

  source_type <- "本地清洗"
} else {
  stop("错误: 找不到输入文件\n  统一准备层: ", unified_cleaned, "\n  本地输入: ", local_input)
}

# 确保输出目录存在
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

log_dir <- dirname(log_file)
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}

# 保存清洗结果
write.table(dat, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n✓ 清洗结果已保存:", output_file, "\n")

# 保存日志
log_content <- sprintf(
  "HGSCore 数据清洗记录\n时间: %s\n数据来源: %s\n输出行数: %d\n",
  Sys.time(), source_type, nrow(dat)
)
writeLines(log_content, log_file)
cat("✓ 清洗日志已保存:", log_file, "\n")
