#!/usr/bin/env Rscript

# MiST 格式转换脚本
# 输出 MiST 所需 data.txt、keys.txt 和 mist.yml。

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}
safe_trim <- function(x) trimws(as.character(x))
truthy <- function(x) toupper(safe_trim(x)) %in% c("Y", "YES", "TRUE", "1")
first_col <- function(cols, aliases) {
  hit <- aliases[aliases %in% cols]
  if (length(hit)) hit[1] else NA_character_
}
clean_field <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  gsub("[\t\r\n]+", " ", x)
}
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE)) else "."
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

cleaned_file <- get_arg("cleaned", file.path(project_dir, "3.清洗后数据", "proteinGroups_cleaned.tsv"))
annotation_file <- get_arg("annotation", file.path(project_dir, "2.样本注释", "sample_annotation_master.tsv"))
output_dir <- get_arg("output_dir", file.path(project_dir, "4.算法特定输入", "mist_input"))
quant_type_arg <- safe_trim(get_arg("quant_type", "auto"))

if (!file.exists(cleaned_file)) stop("找不到清洗后数据: ", cleaned_file)
if (!file.exists(annotation_file)) stop("找不到样本注释表: ", annotation_file)
pg <- read.delim(cleaned_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
anno <- read.delim(annotation_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")

if ("keep" %in% names(anno)) anno <- anno[truthy(anno$keep), , drop = FALSE]
if (!all(c("sample_id", "bait_name") %in% names(anno))) stop("样本注释表缺少 sample_id 或 bait_name。")

if ("quant_type" %in% names(anno)) {
  available <- unique(safe_trim(anno$quant_type[nzchar(safe_trim(anno$quant_type))]))
  preferred <- c("Intensity", "LFQ intensity", "iBAQ", "MS/MS count", "MS/MS Count")
  if (tolower(quant_type_arg) == "auto") {
    hit <- preferred[tolower(preferred) %in% tolower(available)]
    quant_type <- if (length(hit)) hit[1] else available[1]
  } else {
    quant_type <- quant_type_arg
  }
  anno <- anno[tolower(safe_trim(anno$quant_type)) == tolower(quant_type), , drop = FALSE]
} else {
  quant_type <- if (tolower(quant_type_arg) == "auto") "Intensity" else quant_type_arg
}
if (!nrow(anno)) stop("样本注释表中没有可用于 MiST 的样本。")
if (any(!nzchar(safe_trim(anno$bait_name)))) stop("样本注释表存在空 bait_name，请先补全。")

normalize_token <- function(x) {
  x <- safe_trim(x)
  x <- sub("^(Intensity|LFQ intensity|iBAQ|MS/MS count|MS/MS Count)\\s+", "", x, ignore.case = TRUE)
  x <- sub("^Intensity\\.", "", x, ignore.case = TRUE)
  x <- sub("\\.raw$", "", x, ignore.case = TRUE)
  tolower(x)
}
quant_prefix <- paste0(quant_type, ifelse(quant_type == "Intensity.", "", " "))
sample_cols <- names(pg)[startsWith(names(pg), quant_prefix)]
if ("source_column" %in% names(anno)) sample_cols <- unique(c(sample_cols, anno$source_column[anno$source_column %in% names(pg)]))
if (!length(sample_cols)) stop("未找到 quant_type=", quant_type, " 对应的定量列。")
lut <- setNames(sample_cols, normalize_token(sample_cols))

anno$matched_column <- NA_character_
for (i in seq_len(nrow(anno))) {
  candidates <- normalize_token(c(anno$sample_id[i], if ("raw_file_name" %in% names(anno)) anno$raw_file_name[i] else "", if ("source_column" %in% names(anno)) anno$source_column[i] else ""))
  hit <- candidates[candidates %in% names(lut)]
  if (length(hit)) anno$matched_column[i] <- lut[hit[1]]
}
if (any(is.na(anno$matched_column))) {
  bad <- paste(head(anno$sample_id[is.na(anno$matched_column)], 20), collapse = ", ")
  stop("以下样本无法匹配到 proteinGroups 定量列: ", bad)
}

protein_col <- first_col(names(pg), c("Majority protein IDs", "Protein IDs"))
gene_col <- first_col(names(pg), c("Gene names"))
protein_name_col <- first_col(names(pg), c("Protein names"))
unique_pep_col <- first_col(names(pg), c("Unique peptides", "Razor + unique peptides"))
mw_col <- first_col(names(pg), c("Mol. weight [kDa]", "Molecular weight"))
score_col <- first_col(names(pg), c("Score"))
coverage_col <- first_col(names(pg), c("Sequence coverage [%]"))
if (is.na(protein_col)) stop("清洗表缺少 Majority protein IDs 或 Protein IDs。")

mw <- rep(55000, nrow(pg))
if (!is.na(mw_col)) {
  raw_mw <- suppressWarnings(as.numeric(pg[[mw_col]]))
  positive <- raw_mw[!is.na(raw_mw) & raw_mw > 0]
  if (length(positive)) {
    multiplier <- ifelse(grepl("kDa", mw_col, ignore.case = TRUE) || stats::median(positive) < 1000, 1000, 1)
    mw <- raw_mw * multiplier
    mw[is.na(mw) | mw <= 0] <- stats::median(positive) * multiplier
  }
}

long_rows <- list()
for (i in seq_len(nrow(anno))) {
  values <- suppressWarnings(as.numeric(pg[[anno$matched_column[i]]]))
  keep_idx <- which(!is.na(values) & values > 0)
  if (!length(keep_idx)) next
  keep_idx <- keep_idx[order(-values[keep_idx])]
  tmp <- data.frame(
    id = anno$sample_id[i],
    ms_rank = paste0("[", seq_along(keep_idx), "]"),
    ms_uniq_pep = if (!is.na(unique_pep_col)) pg[[unique_pep_col]][keep_idx] else NA,
    ms_uniprot_ac = pg[[protein_col]][keep_idx],
    ms_num_unique_peptide = if (!is.na(unique_pep_col)) pg[[unique_pep_col]][keep_idx] else NA,
    ms_x_cov = if (!is.na(coverage_col)) pg[[coverage_col]][keep_idx] else NA,
    ms_best_disc_score = if (!is.na(score_col)) pg[[score_col]][keep_idx] else NA,
    ms_best_expected_val = NA,
    ms_protein_mw = mw[keep_idx],
    ms_species = NA,
    ms_protein_name = if (!is.na(protein_name_col)) pg[[protein_name_col]][keep_idx] else if (!is.na(gene_col)) pg[[gene_col]][keep_idx] else NA,
    stringsAsFactors = FALSE
  )
  long_rows[[length(long_rows) + 1]] <- tmp
}
if (!length(long_rows)) stop("MiST data.txt 没有生成任何记录，请检查定量列是否有正值。")
data_df <- do.call(rbind, long_rows)
data_df[] <- lapply(data_df, clean_field)
keys_df <- unique(data.frame(id = anno$sample_id, bait_name = anno$bait_name, stringsAsFactors = FALSE))
keys_df[] <- lapply(keys_df, clean_field)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(data_df, file.path(output_dir, "data.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
write.table(keys_df, file.path(output_dir, "keys.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")
writeLines(c(
  "data: data.txt",
  "keys: keys.txt",
  "remove: remove.txt",
  "collapse: collapse.txt",
  "specificity_exclusions: specificity_exclusions.txt"
), file.path(output_dir, "mist.yml"), useBytes = TRUE)
for (nm in c("remove.txt", "collapse.txt", "specificity_exclusions.txt")) {
  p <- file.path(output_dir, nm)
  if (!file.exists(p)) writeLines(character(), p)
}

message("MiST 输入已生成: ", output_dir)
message("data.txt 行数: ", nrow(data_df), "; keys.txt 行数: ", nrow(keys_df))
