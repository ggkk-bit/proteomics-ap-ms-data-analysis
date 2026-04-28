#!/usr/bin/env Rscript

# CompPASS 格式转换脚本
# 将统一清洗表和统一样本注释表转换为 CompPASS 长表输入、bait-run 映射和 prey 参考表。

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}
safe_trim <- function(x) trimws(as.character(x))
truthy <- function(x) tolower(safe_trim(x)) %in% c("1", "true", "yes", "y", "t", "+")
first_col <- function(cols, aliases) {
  hit <- aliases[aliases %in% cols]
  if (length(hit)) hit[1] else NA_character_
}
first_token <- function(x) {
  x <- safe_trim(x)
  sub(";.*$", "", x)
}
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE)) else "."
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

cleaned_file <- get_arg("cleaned", file.path(project_dir, "3.清洗后数据", "proteinGroups_cleaned.tsv"))
annotation_file <- get_arg("annotation", file.path(project_dir, "2.样本注释", "sample_annotation_master.tsv"))
output_dir <- get_arg("output_dir", file.path(project_dir, "4.算法特定输入", "compass_input"))

if (!file.exists(cleaned_file)) stop("找不到清洗后数据: ", cleaned_file)
if (!file.exists(annotation_file)) stop("找不到样本注释表: ", annotation_file)

pg <- read.delim(cleaned_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
anno <- read.delim(annotation_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
required <- c("sample_id", "ip_name", "bait_name")
missing <- setdiff(required, names(anno))
if (length(missing)) stop("样本注释表缺少字段: ", paste(missing, collapse = ", "))
if ("keep" %in% names(anno)) anno <- anno[truthy(anno$keep), , drop = FALSE]

protein_col <- first_col(names(pg), c("Majority protein IDs", "Protein IDs"))
gene_col <- first_col(names(pg), c("Gene names", "Gene names (primary)"))
if (is.na(protein_col)) stop("清洗表缺少 Majority protein IDs 或 Protein IDs。")

base <- data.frame(
  idPrey = first_token(pg[[protein_col]]),
  preyGene = if (!is.na(gene_col)) first_token(pg[[gene_col]]) else "",
  stringsAsFactors = FALSE
)

locate_measurement_column <- function(cols, sample_id, source_column = "") {
  if (nzchar(source_column) && source_column %in% cols) return(c(source_column, sub("\\s+[^ ]+$", "", source_column)))
  patterns <- c(
    paste0("MS/MS Count ", sample_id),
    paste0("MS/MS count ", sample_id),
    paste0("Intensity.", sample_id),
    paste0("Intensity ", sample_id),
    paste0("LFQ intensity ", sample_id),
    paste0("iBAQ ", sample_id)
  )
  hit <- patterns[patterns %in% cols]
  if (!length(hit)) return(c(NA_character_, NA_character_))
  count_type <- sub(paste0("([ .])", sample_id, "$"), "", hit[1])
  c(hit[1], count_type)
}

interaction_parts <- list()
bait_rows <- list()
for (i in seq_len(nrow(anno))) {
  sample_id <- safe_trim(anno$sample_id[i])
  source_column <- if ("source_column" %in% names(anno)) safe_trim(anno$source_column[i]) else ""
  located <- locate_measurement_column(names(pg), sample_id, source_column)
  column_name <- located[1]
  count_type <- located[2]
  if (is.na(column_name)) next
  values <- suppressWarnings(as.numeric(pg[[column_name]]))
  values[is.na(values)] <- 0
  positive <- values > 0 & nzchar(base$idPrey)
  if (any(positive)) {
    subset <- base[positive, , drop = FALSE]
    subset$idRun <- if (nzchar(safe_trim(anno$ip_name[i]))) safe_trim(anno$ip_name[i]) else sample_id
    subset$idBait <- safe_trim(anno$bait_name[i])
    subset$countPrey <- values[positive]
    subset$countType <- count_type
    subset$sourceColumn <- column_name
    interaction_parts[[length(interaction_parts) + 1]] <- subset[, c("idRun", "idBait", "idPrey", "countPrey", "countType", "sourceColumn", "preyGene")]
  }
  bait_rows[[length(bait_rows) + 1]] <- data.frame(
    idRun = if (nzchar(safe_trim(anno$ip_name[i]))) safe_trim(anno$ip_name[i]) else sample_id,
    idBait = safe_trim(anno$bait_name[i]),
    sample_id = sample_id,
    replicate_group = if ("replicate_group" %in% names(anno)) anno$replicate_group[i] else "",
    replicate_letter = if ("replicate_letter" %in% names(anno)) anno$replicate_letter[i] else "",
    raw_file_name = if ("raw_file_name" %in% names(anno)) anno$raw_file_name[i] else "",
    stringsAsFactors = FALSE
  )
}

if (!length(interaction_parts)) stop("没有生成任何 CompPASS 输入记录，请检查样本注释和定量列匹配。")
interaction <- do.call(rbind, interaction_parts)
bait_map <- unique(do.call(rbind, bait_rows))
prey_ref <- unique(interaction[, c("idPrey", "preyGene")])
prey_ref <- prey_ref[order(prey_ref$idPrey, prey_ref$preyGene), , drop = FALSE]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(interaction, file.path(output_dir, "comppass_input.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
write.table(bait_map, file.path(output_dir, "bait_run_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
write.table(prey_ref, file.path(output_dir, "prey_reference.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")

message("CompPASS 输入已生成: ", output_dir)
message("comppass_input.tsv 行数: ", nrow(interaction))
