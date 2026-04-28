#!/usr/bin/env Rscript
# 输出共现次数分布、Z-score 分布、Top 蛋白对和共现热图。

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
cmd <- commandArgs(FALSE)
fa <- cmd[startsWith(cmd, "--file=")]
.script_dir <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "common_cs_score.R"))
args <- commandArgs(trailingOnly = TRUE)
message_header("CS_Score 可视化")

score_file <- get_arg(args, "score", "./6.CS_Score评分结果/cs_score_all_pairs.tsv")
filtered_file <- get_arg(args, "filtered", "./6.CS_Score评分结果/cs_score_filtered_pairs.tsv")
out_dir <- get_arg(args, "out_dir", "./7.可视化")
top_n <- as.integer(get_arg(args, "top_n", "30"))

if (!file.exists(score_file)) stop("找不到评分结果: ", score_file)
ensure_dir(out_dir)
score <- read_tsv(score_file)
plot_dat <- if (file.exists(filtered_file)) read_tsv(filtered_file) else score

p1 <- ggplot(score, aes(x = observed_cooccurrence)) +
  geom_histogram(binwidth = 1, fill = "#4D908E", color = "white") +
  theme_bw() +
  labs(x = "Observed co-occurrence", y = "Protein pairs", title = "Co-occurrence distribution")
ggsave(file.path(out_dir, "01_共现次数分布.png"), p1, width = 7, height = 5, dpi = 300)

p2 <- ggplot(score, aes(x = cs_score_z)) +
  geom_histogram(bins = 50, fill = "#577590", color = "white", na.rm = TRUE) +
  theme_bw() +
  labs(x = "CS_Score Z-score", y = "Protein pairs", title = "Z-score distribution")
ggsave(file.path(out_dir, "02_Zscore分布.png"), p2, width = 7, height = 5, dpi = 300)

top_pairs <- plot_dat %>% arrange(desc(cs_score_z)) %>% head(top_n) %>%
  mutate(pair_label = paste(protein_a, protein_b, sep = " - "),
         pair_label = factor(pair_label, levels = rev(pair_label)))
p3 <- ggplot(top_pairs, aes(x = pair_label, y = cs_score_z, fill = observed_cooccurrence)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient(low = "#90BE6D", high = "#F94144") +
  theme_bw() +
  labs(x = "Protein pair", y = "CS_Score Z-score", fill = "Co-occurrence", title = "Top protein pairs")
ggsave(file.path(out_dir, "03_Top蛋白对.png"), p3, width = 9, height = 7, dpi = 300)

heat_dat <- top_pairs %>%
  select(protein_a, protein_b, observed_cooccurrence) %>%
  bind_rows(top_pairs %>% transmute(protein_a = protein_b, protein_b = protein_a, observed_cooccurrence)) %>%
  mutate(protein_a = factor(protein_a), protein_b = factor(protein_b))
p4 <- ggplot(heat_dat, aes(x = protein_a, y = protein_b, fill = observed_cooccurrence)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#F9C74F", high = "#277DA1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Protein", y = "Protein", fill = "Co-occurrence", title = "Top-pair co-occurrence heatmap")
ggsave(file.path(out_dir, "04_Top蛋白对热图.png"), p4, width = 8, height = 7, dpi = 300)

cat("图表输出目录: ", out_dir, "\n", sep = "")
