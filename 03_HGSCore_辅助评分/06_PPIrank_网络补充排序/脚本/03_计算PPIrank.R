#!/usr/bin/env Rscript
# PPIrank 核心算法实现
# 按论文工作流整合谱数、生物学重复和 negative controls；过滤低谱数噪声并输出网络排序。

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) { prefix <- paste0("--", name, "="); hit <- args[startsWith(args, prefix)]; if (length(hit) == 0) return(default); sub(prefix, "", hit[[1]], fixed = TRUE) }
ensure_dir <- function(path) { dir_path <- dirname(path); if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE) }
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
write_tsv <- function(x, path) { ensure_dir(path); write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "") }
safe_div <- function(num, den, pseudo = 1e-6) (num + pseudo) / (den + pseudo)

compute_reproducibility <- function(x) {
  # 统计每个 bait-prey 在多少个非对照生物学重复中出现。
  x %>% filter(!is_control, spectral_count > 0) %>% group_by(bait_name, prey_id) %>% summarise(replicate_count = n_distinct(replicate_id), total_spectral_count = sum(spectral_count), mean_spectral_count = mean(spectral_count), max_spectral_count = max(spectral_count), .groups = "drop")
}

filter_noise <- function(x, min_replicates = 3, noise_count = 1) {
  # 论文要求：三个重复中都只有 1 个谱图计数的互作视为噪声。
  x %>% group_by(bait_name, prey_id) %>% mutate(appeared_replicates = n_distinct(replicate_id[spectral_count > 0]), max_observed_count = max(spectral_count, na.rm = TRUE), is_low_count_noise = appeared_replicates >= min_replicates & max_observed_count <= noise_count) %>% ungroup() %>% filter(!is_low_count_noise) %>% select(-appeared_replicates, -max_observed_count, -is_low_count_noise)
}

compute_ppirank <- function(x, min_biological_replicates = 3, control_weight = 1) {
  # 评分假设：谱数强度 * 重复出现比例 * 对照惩罚。分数越高，bait-prey 互作越可信。
  bait_reps <- x %>% filter(!is_control) %>% group_by(bait_name) %>% summarise(total_replicates = n_distinct(replicate_id), .groups = "drop")
  low_rep_baits <- bait_reps %>% filter(total_replicates < min_biological_replicates)
  if (nrow(low_rep_baits) > 0) warning("以下 bait 的生物学重复少于 ", min_biological_replicates, ": ", paste(low_rep_baits$bait_name, collapse = ", "))

  bait_summary <- x %>% filter(!is_control) %>% group_by(bait_name, prey_id) %>% summarise(bait_replicates_observed = n_distinct(replicate_id[spectral_count > 0]), bait_total_spectral = sum(spectral_count), bait_mean_spectral = mean(spectral_count), bait_max_spectral = max(spectral_count), .groups = "drop") %>% left_join(bait_reps, by = "bait_name") %>% mutate(replicate_fraction = safe_div(bait_replicates_observed, total_replicates, 0), reproducibility_score = replicate_fraction)
  control_summary <- x %>% filter(is_control) %>% group_by(prey_id) %>% summarise(control_replicates_observed = n_distinct(replicate_id[spectral_count > 0]), control_total_spectral = sum(spectral_count), control_mean_spectral = mean(spectral_count), .groups = "drop")

  bait_summary %>% left_join(control_summary, by = "prey_id") %>% mutate(across(c(control_replicates_observed, control_total_spectral, control_mean_spectral), ~ ifelse(is.na(.x), 0, .x)), spectral_component = log1p(bait_total_spectral), control_penalty = 1 / (1 + control_weight * log1p(control_total_spectral + control_replicates_observed)), ppirank_score = spectral_component * reproducibility_score * control_penalty, confidence_class = case_when(total_replicates >= min_biological_replicates & bait_replicates_observed >= 3 & ppirank_score >= 2 ~ "high", bait_replicates_observed >= 2 & ppirank_score >= 1 ~ "medium", TRUE ~ "low")) %>% arrange(desc(ppirank_score), desc(bait_replicates_observed), desc(bait_total_spectral))
}

build_network <- function(score_table) {
  # 生成 Cytoscape/igraph 可用的边表和节点表；bait 节点标记为绿色。
  edges <- score_table %>% transmute(bait = bait_name, prey = prey_id, score = ppirank_score, replicates = bait_replicates_observed, spectral_count = bait_total_spectral, confidence_class = confidence_class)
  bait_nodes <- edges %>% group_by(protein_id = bait) %>% summarise(total_spectral_count = sum(spectral_count), max_score = max(score), .groups = "drop") %>% mutate(type = "bait", color = "green", size = 18)
  prey_nodes <- edges %>% group_by(protein_id = prey) %>% summarise(total_spectral_count = sum(spectral_count), max_score = max(score), .groups = "drop") %>% mutate(type = "prey", color = "gray", size = 10)
  nodes <- bind_rows(bait_nodes, prey_nodes) %>% group_by(protein_id) %>% summarise(type = ifelse(any(type == "bait"), "bait", "prey"), color = ifelse(any(type == "bait"), "green", "gray"), total_spectral_count = sum(total_spectral_count), max_score = max(max_score), size = ifelse(any(type == "bait"), 18, 10), .groups = "drop")
  list(edges = edges, nodes = nodes)
}

input_file <- get_arg("input", "./4.PPIrank输入文件/ppirank_input_long.tsv")
output_file <- get_arg("output", "./6.PPIrank排序结果/ppirank_scores.tsv")
edge_file <- get_arg("edge_output", "./6.PPIrank排序结果/ppirank_edges.tsv")
node_file <- get_arg("node_output", "./6.PPIrank排序结果/ppirank_nodes.tsv")
qc_file <- get_arg("qc", "./5.质控/03_PPIrank计算记录.txt")
min_replicates <- as.integer(get_arg("min_replicates", "3"))
noise_count <- as.numeric(get_arg("noise_count", "1"))
control_weight <- as.numeric(get_arg("control_weight", "1"))

dat <- read_tsv(input_file)
required_cols <- c("bait_name", "replicate_id", "prey_id", "spectral_count", "is_control")
missing_cols <- setdiff(required_cols, names(dat)); if (length(missing_cols) > 0) stop("PPIrank 输入缺少列: ", paste(missing_cols, collapse = ", "))
dat <- dat %>% mutate(bait_name = as.character(bait_name), replicate_id = as.character(replicate_id), prey_id = as.character(prey_id), spectral_count = suppressWarnings(as.numeric(spectral_count)), spectral_count = ifelse(is.na(spectral_count), 0, spectral_count), is_control = tolower(as.character(is_control)) %in% c("true", "t", "1", "yes", "y", "control", "ctrl"))

repro_before <- compute_reproducibility(dat)
filtered <- filter_noise(dat, min_replicates = min_replicates, noise_count = noise_count)
scores <- compute_ppirank(filtered, min_biological_replicates = min_replicates, control_weight = control_weight)
network <- build_network(scores)
write_tsv(scores, output_file); write_tsv(network$edges, edge_file); write_tsv(network$nodes, node_file)
writeLines(c("PPIrank 计算记录", paste("时间:", Sys.time()), paste("输入:", input_file), paste("输入行数:", nrow(dat)), paste("过滤后行数:", nrow(filtered)), paste("bait 数:", length(unique(dat$bait_name[!dat$is_control]))), paste("control 行数:", sum(dat$is_control)), paste("评分互作数:", nrow(scores)), paste("节点数:", nrow(network$nodes)), paste("边数:", nrow(network$edges)), paste("最小重复数参数:", min_replicates), paste("噪声谱数阈值:", noise_count), paste("对照惩罚权重:", control_weight), paste("重复性统计行数:", nrow(repro_before)), paste("输出:", output_file)), qc_file, useBytes = TRUE)
cat("评分结果: ", output_file, "\n边表: ", edge_file, "\n节点表: ", node_file, "\n", sep = "")
