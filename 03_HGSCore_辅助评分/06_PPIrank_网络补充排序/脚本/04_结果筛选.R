#!/usr/bin/env Rscript
# PPIrank 结果筛选：根据分数、重复数和对照证据筛选高置信 bait-prey 互作。

suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) { prefix <- paste0("--", name, "="); hit <- args[startsWith(args, prefix)]; if (length(hit) == 0) return(default); sub(prefix, "", hit[[1]], fixed = TRUE) }
ensure_dir <- function(path) { dir_path <- dirname(path); if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE) }
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
write_tsv <- function(x, path) { ensure_dir(path); write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "") }

input_file <- get_arg("input", "./6.PPIrank排序结果/ppirank_scores.tsv")
output_file <- get_arg("output", "./6.PPIrank排序结果/ppirank_high_confidence.tsv")
summary_file <- get_arg("summary", "./6.PPIrank排序结果/ppirank_filter_summary.tsv")
qc_file <- get_arg("qc", "./5.质控/04_结果筛选记录.txt")
min_score <- as.numeric(get_arg("min_score", "1"))
min_replicates <- as.integer(get_arg("min_replicates", "2"))
max_control_total <- as.numeric(get_arg("max_control_total", "5"))
keep_class <- get_arg("keep_class", "high,medium")
classes <- trimws(strsplit(keep_class, ",", fixed = TRUE)[[1]])

scores <- read_tsv(input_file)
filtered <- scores %>% filter(ppirank_score >= min_score, bait_replicates_observed >= min_replicates, control_total_spectral <= max_control_total, confidence_class %in% classes) %>% arrange(desc(ppirank_score), desc(bait_replicates_observed), desc(bait_total_spectral))
summary <- data.frame(metric = c("total_interactions", "kept_interactions", "kept_fraction", "min_score", "min_replicates", "max_control_total", "keep_class"), value = c(nrow(scores), nrow(filtered), ifelse(nrow(scores) > 0, nrow(filtered) / nrow(scores), NA), min_score, min_replicates, max_control_total, keep_class), stringsAsFactors = FALSE)
write_tsv(filtered, output_file); write_tsv(summary, summary_file)
writeLines(c("PPIrank 结果筛选记录", paste("时间:", Sys.time()), paste("输入:", input_file), paste("总互作数:", nrow(scores)), paste("保留互作数:", nrow(filtered)), paste("分数阈值:", min_score), paste("重复数阈值:", min_replicates), paste("对照总谱数上限:", max_control_total), paste("保留置信类别:", keep_class), paste("输出:", output_file)), qc_file, useBytes = TRUE)
cat("高置信结果: ", output_file, "\n", sep = "")
