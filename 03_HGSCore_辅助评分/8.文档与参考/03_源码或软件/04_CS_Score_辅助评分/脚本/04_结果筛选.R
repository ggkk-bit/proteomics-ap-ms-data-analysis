#!/usr/bin/env Rscript
# 按论文共现次数约束和用户阈值筛选 CS_Score 结果。

cmd <- commandArgs(FALSE)
fa <- cmd[startsWith(cmd, "--file=")]
.script_dir <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "common_cs_score.R"))
args <- commandArgs(trailingOnly = TRUE)
message_header("CS_Score 结果筛选")

input_file <- get_arg(args, "input", "./6.CS_Score评分结果/cs_score_all_pairs.tsv")
output_file <- get_arg(args, "output", "./6.CS_Score评分结果/cs_score_filtered_pairs.tsv")
min_cooccurrence <- as.integer(get_arg(args, "min_cooccurrence", "2"))
min_z <- as.numeric(get_arg(args, "min_z", "2"))

if (!file.exists(input_file)) stop("找不到评分结果: ", input_file)
res <- read_tsv(input_file)

filtered <- res %>%
  filter(observed_cooccurrence >= min_cooccurrence, !is.na(cs_score_z), cs_score_z >= min_z) %>%
  arrange(desc(cs_score_z), desc(observed_cooccurrence), empirical_p)

write_tsv(filtered, output_file)
cat("筛选前: ", nrow(res), "\n筛选后: ", nrow(filtered), "\n输出: ", output_file, "\n", sep = "")
