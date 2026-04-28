#!/usr/bin/env Rscript
# 按 CRAPome 常用策略筛选：高 FC + 低对照频率 = 高置信互作。

script_path <- commandArgs(FALSE)
file_arg <- script_path[startsWith(script_path, "--file=")]
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "common_crapome.R"))

args <- commandArgs(trailingOnly = TRUE)
message_header("CRAPome 结果筛选")

input_file <- get_arg(args, "input", "./6.CRAPome过滤结果/crapome_fc_all.tsv")
output_file <- get_arg(args, "output", "./6.CRAPome过滤结果/crapome_filtered.tsv")
summary_file <- get_arg(args, "summary", "./6.CRAPome过滤结果/crapome_filter_summary.tsv")
min_fca <- as.numeric(get_arg(args, "min_fca", "10"))
min_fcb <- as.numeric(get_arg(args, "min_fcb", "5"))
max_freq <- as.numeric(get_arg(args, "max_control_frequency", "0.2"))

dat <- read_tsv(input_file) %>% mutate(pass_fca = `FC-A` >= min_fca, pass_fcb = `FC-B` >= min_fcb, pass_frequency = Control_Frequency < max_freq, CRAPome_Pass = pass_fca & pass_fcb & pass_frequency, Confidence = case_when(CRAPome_Pass ~ "High", (`FC-A` >= min_fca | `FC-B` >= min_fcb) & Control_Frequency < max_freq ~ "Medium", TRUE ~ "Low")) %>% arrange(desc(CRAPome_Pass), desc(`FC-A`), desc(`FC-B`), Control_Frequency)
filtered <- dat %>% filter(CRAPome_Pass)
summary <- dat %>% group_by(`Bait Name`, `AP Name`) %>% summarise(Total_Prey = n(), Passed_Prey = sum(CRAPome_Pass), High = sum(Confidence == "High"), Medium = sum(Confidence == "Medium"), Low = sum(Confidence == "Low"), .groups = "drop")
write_tsv(filtered, output_file)
write_tsv(summary, summary_file)
write_tsv(dat, sub("\\.tsv$", "_with_flags.tsv", output_file))
cat("通过筛选: ", nrow(filtered), " / ", nrow(dat), "\n", sep = "")
cat("输出: ", output_file, "\n", sep = "")
