#!/usr/bin/env Rscript
# 生成 CS_Score 输入：成员长表 purification x prey，以及 purification x protein 二值矩阵。

cmd <- commandArgs(FALSE)
fa <- cmd[startsWith(cmd, "--file=")]
.script_dir <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "common_cs_score.R"))
args <- commandArgs(trailingOnly = TRUE)
message_header("CS_Score 格式转换")

cleaned_file <- get_arg(args, "cleaned", "./3.数据清洗与预处理/proteinGroups_cleaned.tsv")
annotation_file <- get_arg(args, "annotation", "../../00_统一数据准备/2.样本注释/sample_annotation_master.tsv")
member_out <- get_arg(args, "member_out", "./4.CS_Score输入文件/cs_score_member_long.tsv")
matrix_out <- get_arg(args, "matrix_out", "./4.CS_Score输入文件/cs_score_binary_matrix.tsv")
qc_out <- get_arg(args, "qc", "./5.质控/02_格式转换统计.tsv")
quant_type <- get_arg(args, "quant_type", "auto")
presence_threshold <- as.numeric(get_arg(args, "presence_threshold", "0"))

if (!file.exists(cleaned_file)) stop("找不到清洗后数据: ", cleaned_file)
dat <- read_tsv(cleaned_file)
protein_col <- detect_protein_id_col(dat)
quant_prefix <- detect_quant_prefix(dat, quant_type)
quant_cols <- names(dat)[startsWith(names(dat), quant_prefix)]
if (length(quant_cols) == 0) stop("未找到定量列: ", quant_prefix)

annot <- if (file.exists(annotation_file)) read_tsv(annotation_file) else NULL

long <- dat %>%
  select(all_of(protein_col), all_of(quant_cols)) %>%
  pivot_longer(cols = all_of(quant_cols), names_to = "sample_id", values_to = "abundance") %>%
  mutate(
    sample_id = sub(quant_prefix, "", sample_id, fixed = TRUE),
    abundance = suppressWarnings(as.numeric(abundance)),
    prey_id = normalize_protein_id(.data[[protein_col]])
  ) %>%
  filter(!is.na(abundance), abundance > presence_threshold, !is.na(prey_id), prey_id != "")

if (!is.null(annot) && "sample_id" %in% names(annot)) {
  long <- long %>% left_join(annot, by = "sample_id")
}

if (!"bait_name" %in% names(long)) long$bait_name <- long$sample_id
if (!"replicate_group" %in% names(long)) long$replicate_group <- long$sample_id

member <- long %>%
  transmute(
    purification_id = sample_id,
    sample_id = sample_id,
    bait_name = bait_name,
    replicate_group = replicate_group,
    prey_id = prey_id,
    gene_name = prey_id,
    presence = 1
  ) %>%
  distinct()

binary_matrix <- member %>%
  select(purification_id, prey_id, presence) %>%
  distinct() %>%
  pivot_wider(names_from = prey_id, values_from = presence, values_fill = 0)

qc <- data.frame(
  metric = c("proteins_in_cleaned", "quant_columns", "member_rows", "purifications", "unique_prey"),
  value = c(nrow(dat), length(quant_cols), nrow(member), n_distinct(member$purification_id), n_distinct(member$prey_id))
)

write_tsv(member, member_out)
write_tsv(binary_matrix, matrix_out)
write_tsv(qc, qc_out)
cat("成员长表: ", member_out, "\n二值矩阵: ", matrix_out, "\n", sep = "")
