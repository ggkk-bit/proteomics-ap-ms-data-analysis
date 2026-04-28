#!/usr/bin/env Rscript
# 将 MaxQuant 宽表转换为 CRAPome 四列长表：Bait Name, AP Name, Prey Name, Spectral Count。

script_path <- commandArgs(FALSE)
file_arg <- script_path[startsWith(script_path, "--file=")]
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "common_crapome.R"))

args <- commandArgs(trailingOnly = TRUE)
message_header("CRAPome 格式转换")

input_file <- get_arg(args, "input", "./3.数据清洗与预处理/proteinGroups_cleaned.tsv")
sample_info <- get_arg(args, "sample_info", "./2.样本信息与分组/sample_info.tsv")
output_file <- get_arg(args, "output", "./4.CRAPome输入文件/crapome_input_long.tsv")
matrix_file <- get_arg(args, "matrix_output", "./4.CRAPome输入文件/spectral_count_matrix.tsv")
quant_prefix <- get_arg(args, "quant_type", "auto")
min_count <- as.numeric(get_arg(args, "min_count", "0"))

dat <- read_tsv(input_file)
id_col <- detect_protein_id_col(dat)
quant_prefix <- detect_quant_prefix(dat, quant_prefix)
quant_cols <- names(dat)[startsWith(names(dat), quant_prefix)]
if (length(quant_cols) == 0) stop("没有找到定量列，请检查 --quant_type= 参数。")

long <- dat %>%
  transmute(`Prey Name` = normalize_protein_id(.data[[id_col]]), across(all_of(quant_cols))) %>%
  pivot_longer(cols = all_of(quant_cols), names_to = "Sample", values_to = "Spectral Count") %>%
  mutate(`Spectral Count` = suppressWarnings(as.numeric(`Spectral Count`)), `Spectral Count` = ifelse(is.na(`Spectral Count`), 0, `Spectral Count`), `AP Name` = sub(paste0("^", gsub("([\\W])", "\\\\\\1", quant_prefix)), "", Sample), `Bait Name` = `AP Name`) %>%
  filter(`Spectral Count` > min_count) %>%
  select(`Bait Name`, `AP Name`, `Prey Name`, `Spectral Count`)

if (file.exists(sample_info)) {
  info <- read_tsv(sample_info)
  if (all(c("AP Name", "Bait Name") %in% names(info))) {
    long <- long %>% select(-`Bait Name`) %>% left_join(info[, c("AP Name", "Bait Name")], by = "AP Name") %>% mutate(`Bait Name` = ifelse(is.na(`Bait Name`) | `Bait Name` == "", `AP Name`, `Bait Name`)) %>% select(`Bait Name`, `AP Name`, `Prey Name`, `Spectral Count`)
  }
}

matrix_out <- long %>% group_by(`Prey Name`, `AP Name`) %>% summarise(`Spectral Count` = sum(`Spectral Count`), .groups = "drop") %>% pivot_wider(names_from = `AP Name`, values_from = `Spectral Count`, values_fill = 0)
write_tsv(long, output_file)
write_tsv(matrix_out, matrix_file)
writeLines(c("CRAPome 格式转换质控", paste("输入:", input_file), paste("蛋白列:", id_col), paste("定量前缀:", quant_prefix), paste("样本数:", length(unique(long$`AP Name`))), paste("蛋白数:", length(unique(long$`Prey Name`))), paste("长表行数:", nrow(long))), "./5.质控/02_格式转换记录.txt", useBytes = TRUE)
cat("输出: ", output_file, "\n", sep = "")
