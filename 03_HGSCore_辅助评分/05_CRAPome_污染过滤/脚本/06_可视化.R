#!/usr/bin/env Rscript
# 生成对照频率分布、FC 散点图、Top 蛋白条形图。

script_path <- commandArgs(FALSE)
file_arg <- script_path[startsWith(script_path, "--file=")]
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "common_crapome.R"))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
message_header("CRAPome 可视化")

input_file <- get_arg(args, "input", "./6.CRAPome过滤结果/crapome_filtered_with_flags.tsv")
output_dir <- get_arg(args, "output_dir", "./7.可视化")
top_n <- as.integer(get_arg(args, "top_n", "30"))
ensure_dir(output_dir)

dat <- read_tsv(input_file)
if (!"Confidence" %in% names(dat)) dat$Confidence <- "All"

p1 <- ggplot(dat, aes(x = Control_Frequency)) + geom_histogram(bins = 40, fill = "#4C78A8", color = "white") + geom_vline(xintercept = 0.2, linetype = 2, color = "#E45756") + theme_bw(base_size = 12) + labs(title = "对照频率分布", x = "Control Frequency", y = "蛋白数量")
p2 <- ggplot(dat, aes(x = `FC-A`, y = `FC-B`, color = Confidence, size = 1 - Control_Frequency)) + geom_point(alpha = 0.75) + scale_x_log10() + scale_y_log10() + geom_vline(xintercept = 10, linetype = 2, color = "#E45756") + geom_hline(yintercept = 5, linetype = 2, color = "#E45756") + theme_bw(base_size = 12) + labs(title = "FC-A / FC-B 富集散点图", x = "FC-A (log10)", y = "FC-B (log10)", size = "低污染权重")
top_dat <- dat %>% arrange(desc(CRAPome_Pass), desc(`FC-A`), Control_Frequency) %>% slice_head(n = top_n) %>% mutate(`Prey Name` = factor(`Prey Name`, levels = rev(unique(`Prey Name`))))
p3 <- ggplot(top_dat, aes(x = `Prey Name`, y = `FC-A`, fill = Confidence)) + geom_col() + coord_flip() + scale_y_log10() + theme_bw(base_size = 12) + labs(title = "Top 蛋白 FC-A", x = "Prey", y = "FC-A (log10)")

ggsave(file.path(output_dir, "01_对照频率分布.png"), p1, width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, "02_FC散点图.png"), p2, width = 7, height = 5, dpi = 300)
ggsave(file.path(output_dir, "03_Top蛋白条形图.png"), p3, width = 8, height = 6, dpi = 300)
cat("图表输出目录: ", output_dir, "\n", sep = "")
