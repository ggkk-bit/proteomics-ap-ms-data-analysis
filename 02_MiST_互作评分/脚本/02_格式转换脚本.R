#!/usr/bin/env Rscript

# MiST 格式转换脚本
# 作用：
# 1. 读取清洗后的 proteinGroups
# 2. 结合样本注释表生成 MiST 官方输入 data.txt 和 keys.txt
# 3. 可选生成 remove.txt / collapse.txt / specificity_exclusions.txt 的模板
# 4. 默认优先使用自动生成的样本注释表，提升通用性

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

safe_trim <- function(x) trimws(as.character(x))

clean_field <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  gsub("[\t\r\n]+", " ", x)
}

nonempty_count <- function(x) {
  y <- trimws(as.character(x))
  sum(!is.na(x) & (y != ""))
}

get_script_dir <- function() {
  hit <- grep("^--file=", all_args, value = TRUE)
  if (!length(hit)) return("")
  normalizePath(dirname(sub("^--file=", "", hit[[1]])), winslash = "\\", mustWork = FALSE)
}

find_existing_file <- function(candidates) {
  candidates <- unique(candidates[nzchar(candidates)])
  for (cand in candidates) {
    if (file.exists(cand)) return(cand)
  }
  ""
}

find_sample_file <- function(user_file = "") {
  script_dir <- get_script_dir()
  project_dir <- if (nzchar(script_dir)) normalizePath(file.path(script_dir, ".."), winslash = "\\", mustWork = FALSE) else ""

  candidates <- c(
    user_file,
    "2.样本信息与分组/sample_annotation_auto_R.tsv",
    "2.样本信息与分组/sample_annotation_auto.tsv",
    "2.样本信息与分组/sample_annotation_template.tsv",
    "..\\2.样本信息与分组\\sample_annotation_auto_R.tsv",
    "..\\2.样本信息与分组\\sample_annotation_auto.tsv",
    "..\\2.样本信息与分组\\sample_annotation_template.tsv"
  )

  if (nzchar(script_dir)) {
    candidates <- c(
      candidates,
      file.path(script_dir, "..", "2.样本信息与分组", "sample_annotation_auto_R.tsv"),
      file.path(script_dir, "..", "2.样本信息与分组", "sample_annotation_auto.tsv"),
      file.path(script_dir, "..", "2.样本信息与分组", "sample_annotation_template.tsv")
    )
  }

  if (nzchar(project_dir)) {
    candidates <- c(
      candidates,
      file.path(project_dir, "2.样本信息与分组", "sample_annotation_auto_R.tsv"),
      file.path(project_dir, "2.样本信息与分组", "sample_annotation_auto.tsv"),
      file.path(project_dir, "2.样本信息与分组", "sample_annotation_template.tsv")
    )
  }

  out <- find_existing_file(candidates)
  if (nzchar(out)) normalizePath(out, winslash = "\\", mustWork = TRUE) else ""
}

input_file <- get_arg("input", "3.数据清洗与预处理/proteinGroups_已清洗.tsv")
sample_file <- find_sample_file(get_arg("samples", ""))
quant_type_arg <- safe_trim(get_arg("quant_type", "auto"))
data_out <- get_arg("data_out", "4.MiST输入文件/data.txt")
keys_out <- get_arg("keys_out", "4.MiST输入文件/keys.txt")
remove_out <- get_arg("remove_out", "4.MiST输入文件/remove.txt")
collapse_out <- get_arg("collapse_out", "4.MiST输入文件/collapse.txt")
spec_out <- get_arg("spec_out", "4.MiST输入文件/specificity_exclusions.txt")

if (!file.exists(input_file)) stop("找不到清洗后的输入文件: ", input_file)
if (!nzchar(sample_file) || !file.exists(sample_file)) {
  stop("找不到样本注释表。已尝试 sample_annotation_auto_R.tsv / sample_annotation_auto.tsv / sample_annotation_template.tsv")
}

pg <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
anno <- read.delim(sample_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")

if (!nrow(pg)) stop("清洗后的输入文件没有数据行: ", input_file)
if (!nrow(anno)) stop("样本注释表没有数据行: ", sample_file)

if ("quant_type" %in% names(anno)) {
  available_quant_types <- unique(safe_trim(anno$quant_type[nzchar(safe_trim(anno$quant_type))]))
  preferred_quant_types <- c("Intensity", "LFQ intensity", "iBAQ", "MS/MS count")
  if (tolower(quant_type_arg) == "auto") {
    preferred_hit <- preferred_quant_types[tolower(preferred_quant_types) %in% tolower(available_quant_types)]
    quant_type <- if (length(preferred_hit)) preferred_hit[[1]] else available_quant_types[[1]]
  } else {
    quant_type <- quant_type_arg
  }
  keep_quant_type <- tolower(safe_trim(anno$quant_type)) == tolower(quant_type)
  if (!any(keep_quant_type)) {
    stop(
      "样本注释表中没有 quant_type=", quant_type, " 的样本。可用 quant_type: ",
      paste(available_quant_types, collapse = ", "),
      "。如需指定其他定量列，请加参数 --quant_type='LFQ intensity'。"
    )
  }
  anno <- anno[keep_quant_type, , drop = FALSE]
} else {
  quant_type <- if (tolower(quant_type_arg) == "auto") "Intensity" else quant_type_arg
}

required_anno <- c("sample_id", "bait_name")
missing_anno <- setdiff(required_anno, names(anno))
if (length(missing_anno)) stop("样本表缺少字段: ", paste(missing_anno, collapse = ", "))

if (!("keep" %in% names(anno))) anno$keep <- "Y"
anno <- anno[toupper(safe_trim(anno$keep)) %in% c("Y", "YES", "TRUE", "1"), , drop = FALSE]
if (!nrow(anno)) stop("样本注释表中没有 keep=Y 的样本")

anno$sample_id <- safe_trim(anno$sample_id)
anno$bait_name <- safe_trim(anno$bait_name)

missing_bait <- which(!nzchar(anno$bait_name))
if (length(missing_bait)) {
  bad <- paste(head(anno$sample_id[missing_bait], 20), collapse = ", ")
  stop("样本注释表中存在空的 bait_name。请先补全。示例 sample_id: ", bad)
}

quant_prefix <- paste0(quant_type, " ")
sample_cols <- names(pg)[startsWith(names(pg), quant_prefix)]
if ("source_column" %in% names(anno)) {
  source_cols <- safe_trim(anno$source_column)
  source_cols <- source_cols[nzchar(source_cols) & (source_cols %in% names(pg))]
  sample_cols <- unique(c(sample_cols, source_cols))
}
if (length(sample_cols) == 0) {
  stop("未找到定量列，当前 quant_type=", quant_type, "。需要类似 '", quant_prefix, "<sample>' 的列。")
}

sample_names <- ifelse(
  startsWith(sample_cols, quant_prefix),
  substr(sample_cols, nchar(quant_prefix) + 1L, nchar(sample_cols)),
  sample_cols
)
normalize_token <- function(x) {
  x <- safe_trim(x)
  x <- sub("^(Intensity|LFQ intensity|iBAQ|MS/MS count)\\s+", "", x, ignore.case = TRUE)
  x <- sub("\\.raw$", "", x, ignore.case = TRUE)
  tolower(x)
}

sample_names_norm <- normalize_token(sample_names)
sample_names_lut <- setNames(sample_cols, sample_names_norm)

anno$sample_id_norm <- normalize_token(anno$sample_id)
anno$matched_sample <- anno$sample_id_norm
anno$matched_by <- ifelse(anno$sample_id_norm %in% names(sample_names_lut), "sample_id", NA_character_)
anno$matched_column <- ifelse(anno$sample_id_norm %in% names(sample_names_lut), sample_names_lut[anno$sample_id_norm], NA_character_)

if ("raw_file_name" %in% names(anno)) {
  raw_norm <- normalize_token(anno$raw_file_name)
  use_raw <- is.na(anno$matched_by) & raw_norm %in% names(sample_names_lut)
  anno$matched_sample[use_raw] <- raw_norm[use_raw]
  anno$matched_by[use_raw] <- "raw_file_name"
  anno$matched_column[use_raw] <- sample_names_lut[raw_norm[use_raw]]
}

if ("source_column" %in% names(anno)) {
  source_col <- safe_trim(anno$source_column)
  use_source <- is.na(anno$matched_by) & source_col %in% sample_cols
  anno$matched_by[use_source] <- "source_column"
  anno$matched_column[use_source] <- source_col[use_source]
  anno$matched_sample[use_source] <- normalize_token(source_col[use_source])
}

unmatched <- is.na(anno$matched_by)
if (any(unmatched)) {
  preview <- paste(head(sample_names, 12), collapse = ", ")
  bad_rows <- apply(anno[unmatched, intersect(c("sample_id", "raw_file_name", "bait_name"), names(anno)), drop = FALSE], 1, function(z) paste(z, collapse = " | "))
  stop(
    "样本注释表中的样本无法匹配 proteinGroups 强度列。\n",
    "请填写 sample_id=强度列名去掉 'Intensity ' 后的值，例如: ", preview, "\n",
    "如果样本表包含 raw_file_name，也可以填写如 b10494.raw，脚本会自动去掉 .raw 后匹配。\n",
    "未匹配行:\n",
    paste(bad_rows, collapse = "\n")
  )
}

pick_col <- function(candidates) {
  hit <- intersect(candidates, names(pg))
  if (length(hit)) hit[[1]] else NA_character_
}

protein_id_col <- pick_col(c("Majority protein IDs", "Protein IDs"))
gene_col <- pick_col(c("Gene names"))
protein_name_col <- pick_col(c("Protein names"))
unique_pep_col <- pick_col(c("Unique peptides"))
mw_col <- pick_col(c("Mol. weight [kDa]"))
score_col <- pick_col(c("Score"))
coverage_col <- pick_col(c("Sequence coverage [%]"))

if (is.na(protein_id_col)) stop("找不到蛋白 ID 字段（Majority protein IDs / Protein IDs）")
if (nonempty_count(pg[[protein_id_col]]) == 0) {
  stop("蛋白 ID 字段全为空: ", protein_id_col, "。请检查清洗后的输入文件是否被写坏。")
}

mw_multiplier <- 1
mw_fallback <- 55000
mw_values_all <- rep(mw_fallback, nrow(pg))
if (!is.na(mw_col)) {
  mw_raw <- suppressWarnings(as.numeric(pg[[mw_col]]))
  positive_mw <- mw_raw[!is.na(mw_raw) & (mw_raw > 0)]
  if (length(positive_mw)) {
    # MaxQuant 的 Mol. weight [kDa] 是 kDa；MiST 官方代码按 Da/110 估算蛋白长度。
    if (grepl("kDa", mw_col, ignore.case = TRUE) || stats::median(positive_mw) < 1000) {
      mw_multiplier <- 1000
    }
    mw_fallback <- stats::median(positive_mw) * mw_multiplier
    mw_values_all <- mw_raw * mw_multiplier
    mw_values_all[is.na(mw_values_all) | mw_values_all <= 0] <- mw_fallback
  }
}

find_sample_metric_col <- function(prefixes, sample_label, sample_id, raw_file_name = "") {
  ids <- unique(clean_field(c(sample_label, sample_id, sub("\\.raw$", "", raw_file_name, ignore.case = TRUE))))
  ids <- ids[nzchar(ids)]
  candidates <- as.vector(outer(prefixes, ids, paste, sep = " "))
  hit <- candidates[candidates %in% names(pg)]
  if (length(hit)) hit[[1]] else NA_character_
}

long_rows <- vector("list", length = 0)
samples_with_signal <- 0L
total_positive_values <- 0L

for (i in seq_len(nrow(anno))) {
  sid <- anno$sample_id[i]
  col_name <- anno$matched_column[i]
  sample_label <- if (startsWith(col_name, quant_prefix)) substr(col_name, nchar(quant_prefix) + 1L, nchar(col_name)) else sid
  x <- suppressWarnings(as.numeric(pg[[col_name]]))
  keep_idx <- which((!is.na(x)) & (x > 0))
  if (!length(keep_idx)) next
  samples_with_signal <- samples_with_signal + 1L
  total_positive_values <- total_positive_values + length(keep_idx)

  ord <- order(-x[keep_idx], na.last = TRUE)
  keep_idx <- keep_idx[ord]

  raw_file_name <- if ("raw_file_name" %in% names(anno)) anno$raw_file_name[i] else ""
  sample_pep_col <- find_sample_metric_col(
    c("Unique peptides", "Razor + unique peptides", "Peptides"),
    sample_label,
    sid,
    raw_file_name
  )
  pep_values <- if (!is.na(sample_pep_col)) pg[[sample_pep_col]][keep_idx] else if (!is.na(unique_pep_col)) pg[[unique_pep_col]][keep_idx] else NA

  tmp <- data.frame(
    id = sid,
    ms_rank = paste0("[", seq_along(keep_idx), "]"),
    ms_uniq_pep = pep_values,
    ms_uniprot_ac = pg[[protein_id_col]][keep_idx],
    ms_num_unique_peptide = pep_values,
    ms_x_cov = if (!is.na(coverage_col)) pg[[coverage_col]][keep_idx] else NA,
    ms_best_disc_score = if (!is.na(score_col)) pg[[score_col]][keep_idx] else NA,
    ms_best_expected_val = NA,
    ms_protein_mw = mw_values_all[keep_idx],
    ms_species = NA,
    ms_protein_name = if (!is.na(protein_name_col)) pg[[protein_name_col]][keep_idx] else NA,
    stringsAsFactors = FALSE
  )

  if (!is.na(gene_col)) {
    tmp$ms_protein_name <- ifelse(is.na(tmp$ms_protein_name) | (tmp$ms_protein_name == ""), pg[[gene_col]][keep_idx], tmp$ms_protein_name)
  }

  long_rows[[length(long_rows) + 1]] <- tmp
}

data_out_df <- if (length(long_rows)) do.call(rbind, long_rows) else data.frame(
  id = character(),
  ms_rank = character(),
  ms_uniq_pep = character(),
  ms_uniprot_ac = character(),
  ms_num_unique_peptide = character(),
  ms_x_cov = character(),
  ms_best_disc_score = character(),
  ms_best_expected_val = character(),
  ms_protein_mw = character(),
  ms_species = character(),
  ms_protein_name = character(),
  stringsAsFactors = FALSE
)

if (!nrow(data_out_df)) {
  stop(
    "格式转换后 data.txt 为 0 行。\n",
    "已匹配样品数: ", nrow(anno), "\n",
    "有正强度的样品数: ", samples_with_signal, "\n",
    "正强度值总数: ", total_positive_values, "\n",
    "常见原因：清洗后的 proteinGroups 表被写坏、所有强度为 0/NA、或 quant_type 选错。\n",
    "当前 quant_type: ", quant_type
  )
}

data_out_df[] <- lapply(data_out_df, clean_field)

dir.create(dirname(data_out), recursive = TRUE, showWarnings = FALSE)
write.table(data_out_df, data_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")

keys_df <- unique(anno[, c("sample_id", "bait_name")])
names(keys_df) <- c("id", "bait_name")
keys_df[] <- lapply(keys_df, clean_field)
write.table(keys_df, keys_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "")

if (!file.exists(remove_out)) writeLines(character(), remove_out)
if (!file.exists(collapse_out)) writeLines(character(), collapse_out)
if (!file.exists(spec_out)) writeLines(character(), spec_out)

message("MiST 格式转换完成")
message("使用样本注释表: ", sample_file)
message("使用定量类型: ", quant_type)
if (!is.na(mw_col)) {
  message("使用分子量字段: ", mw_col, ifelse(mw_multiplier == 1000, "（已从 kDa 转为 Da）", ""))
} else {
  message("未找到分子量字段，使用默认分子量: ", mw_fallback)
}
message("有正强度的样品数: ", samples_with_signal, "/", nrow(anno))
message("输出 data.txt 行数: ", nrow(data_out_df))
message("输出 keys.txt 样品/run 数: ", nrow(keys_df))
message("输出 keys.txt 唯一 bait 数: ", length(unique(keys_df$bait_name)))
