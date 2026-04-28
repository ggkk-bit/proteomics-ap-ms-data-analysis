#!/usr/bin/env Rscript
# PPIrank 格式转换：将 MaxQuant 宽表转换为样本级长表 bait_name, replicate_id, prey_id, spectral_count, is_control。

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) { prefix <- paste0("--", name, "="); hit <- args[startsWith(args, prefix)]; if (length(hit) == 0) return(default); sub(prefix, "", hit[[1]], fixed = TRUE) }
ensure_dir <- function(path) { dir_path <- dirname(path); if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE) }
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
write_tsv <- function(x, path) { ensure_dir(path); write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "") }
detect_protein_id_col <- function(dat) { candidates <- c("Majority protein IDs", "Protein IDs", "Protein.IDs", "Gene names", "Protein names", "prey_id"); hit <- candidates[candidates %in% names(dat)]; if (length(hit) == 0) stop("未找到蛋白标识列。"); hit[[1]] }
normalize_protein_id <- function(x) { x <- as.character(x); x <- sub(";.*$", "", x); trimws(x) }
detect_quant_prefix <- function(dat, requested = "auto") { if (!identical(requested, "auto")) return(requested); prefixes <- c("MS/MS count ", "Spectral count ", "Spectral Count ", "Intensity ", "LFQ intensity ", "iBAQ "); for (prefix in prefixes) if (any(startsWith(names(dat), prefix))) return(prefix); stop("未找到可识别的定量列前缀。") }

input_file <- get_arg("input", "./3.数据清洗与预处理/proteinGroups_cleaned.tsv")
sample_info <- get_arg("sample_info", "./2.样本信息与分组/sample_info.tsv")
output_file <- get_arg("output", "./4.PPIrank输入文件/ppirank_input_long.tsv")
matrix_file <- get_arg("matrix_output", "./4.PPIrank输入文件/spectral_count_matrix.tsv")
template_file <- get_arg("template_output", "./2.样本信息与分组/sample_info_template.tsv")
quant_type <- get_arg("quant_type", "auto")
min_count <- as.numeric(get_arg("min_count", "0"))

dat <- read_tsv(input_file)
id_col <- detect_protein_id_col(dat)
quant_prefix <- detect_quant_prefix(dat, quant_type)
quant_cols <- names(dat)[startsWith(names(dat), quant_prefix)]
sample_ids <- sub(paste0("^", gsub("([\\W])", "\\\\\\1", quant_prefix)), "", quant_cols)

if (!file.exists(sample_info)) {
  info <- data.frame(sample_id = sample_ids, bait_name = sample_ids, replicate_id = ave(sample_ids, sample_ids, FUN = seq_along), is_control = FALSE, stringsAsFactors = FALSE)
  write_tsv(info, template_file)
  cat("未找到样本注释文件，已生成模板: ", template_file, "\n", sep = "")
} else info <- read_tsv(sample_info)

required_cols <- c("sample_id", "bait_name", "replicate_id", "is_control")
missing_cols <- setdiff(required_cols, names(info)); if (length(missing_cols) > 0) stop("样本注释缺少列: ", paste(missing_cols, collapse = ", "))
info <- info %>% mutate(sample_id = as.character(sample_id), bait_name = as.character(bait_name), replicate_id = as.character(replicate_id), is_control = tolower(as.character(is_control)) %in% c("true", "t", "1", "yes", "y", "control", "ctrl"))

long <- dat %>%
  transmute(prey_id = normalize_protein_id(.data[[id_col]]), across(all_of(quant_cols))) %>%
  pivot_longer(cols = all_of(quant_cols), names_to = "raw_sample", values_to = "spectral_count") %>%
  mutate(sample_id = sub(paste0("^", gsub("([\\W])", "\\\\\\1", quant_prefix)), "", raw_sample), spectral_count = suppressWarnings(as.numeric(spectral_count)), spectral_count = ifelse(is.na(spectral_count), 0, spectral_count)) %>%
  left_join(info, by = "sample_id") %>%
  mutate(bait_name = ifelse(is.na(bait_name) | bait_name == "", sample_id, bait_name), replicate_id = ifelse(is.na(replicate_id) | replicate_id == "", sample_id, replicate_id), is_control = ifelse(is.na(is_control), FALSE, is_control)) %>%
  filter(!is.na(prey_id), prey_id != "", spectral_count > min_count) %>%
  select(bait_name, replicate_id, prey_id, spectral_count, is_control, sample_id)

write_tsv(long %>% select(bait_name, replicate_id, prey_id, spectral_count, is_control), output_file)
write_tsv(long %>% group_by(prey_id, sample_id) %>% summarise(spectral_count = sum(spectral_count), .groups = "drop") %>% pivot_wider(names_from = sample_id, values_from = spectral_count, values_fill = 0), matrix_file)
writeLines(c("PPIrank 格式转换记录", paste("输入:", input_file), paste("蛋白列:", id_col), paste("定量前缀:", quant_prefix), paste("样本数:", length(unique(long$sample_id))), paste("bait 数:", length(unique(long$bait_name[!long$is_control]))), paste("control 样本数:", length(unique(long$sample_id[long$is_control]))), paste("长表行数:", nrow(long)), paste("输出:", output_file)), "./5.质控/02_格式转换记录.txt", useBytes = TRUE)
cat("输出: ", output_file, "\n", sep = "")
