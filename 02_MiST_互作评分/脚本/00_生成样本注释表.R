#!/usr/bin/env Rscript

# 通用样本注释表生成脚本
# 作用：
# 1. 从 proteinGroups.txt 或清洗后的 TSV 表头中识别样本定量列
# 2. 支持 Intensity / LFQ intensity / iBAQ / MS/MS count 等常见 MaxQuant 列前缀
# 3. 可选结合 archiveSummary.tsv 自动补充 bait、批次、孔位、重复和常见对照标记
# 4. 生成可直接用于后续算法流程的样本注释表初稿

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

safe_trim <- function(x) trimws(as.character(x))

read_header <- function(path) {
  con <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  line <- readLines(con, n = 1, warn = FALSE)
  if (!length(line)) stop("输入文件为空: ", path)
  strsplit(line, "\t", fixed = TRUE)[[1]]
}

normalize_prefixes <- function(x) {
  pieces <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  pieces <- pieces[nzchar(pieces)]
  pieces <- ifelse(grepl(" $", pieces), pieces, paste0(pieces, " "))
  unique(pieces)
}

extract_replicate_from_name <- function(x) {
  x <- safe_trim(x)
  out <- rep("", length(x))
  hit1 <- regexpr("_[A-Z]_", x)
  ok1 <- hit1 > 0
  out[ok1] <- substring(x[ok1], hit1[ok1] + 1L, hit1[ok1] + 1L)
  out
}

extract_bait_from_info <- function(info) {
  info <- safe_trim(info)
  sub("^.*Bait:([^()]+).*$", "\\1", info)
}

extract_condition_from_info <- function(info) {
  info <- safe_trim(info)
  has_paren <- grepl("\\(", info)
  out <- rep("", length(info))
  out[has_paren] <- sub("^.*\\(([^,()]+).*$", "\\1", info[has_paren])
  out
}

extract_platewell_from_info <- function(info) {
  info <- safe_trim(info)
  has_two <- grepl(",", info)
  out <- rep("", length(info))
  out[has_two] <- sub("^.*,[[:space:]]*([^,()]+)\\).*$", "\\1", info[has_two])
  out
}

is_control_bait <- function(bait) {
  bait <- tolower(safe_trim(bait))
  bait %in% c("ctrl", "control", "mock", "empty", "vector", "beads", "igg", "gfp", "flag", "ha")
}

find_archive_file <- function(user_archive = "", input_path = "") {
  if (nzchar(user_archive) && file.exists(user_archive)) {
    return(normalizePath(user_archive, winslash = "\\", mustWork = TRUE))
  }

  candidates <- c(
    "archiveSummary.tsv",
    ".\\archiveSummary.tsv",
    "..\\archiveSummary.tsv",
    "..\\1.原始输入\\archiveSummary.tsv",
    "..\\..\\MIST\\1.原始输入\\archiveSummary.tsv",
    "..\\..\\..\\MIST\\1.原始输入\\archiveSummary.tsv"
  )

  if (nzchar(input_path)) {
    input_dir <- dirname(input_path)
    candidates <- c(
      candidates,
      file.path(input_dir, "archiveSummary.tsv"),
      file.path(input_dir, "..", "archiveSummary.tsv"),
      file.path(input_dir, "..", "1.原始输入", "archiveSummary.tsv"),
      file.path(input_dir, "..", "..", "MIST", "1.原始输入", "archiveSummary.tsv"),
      file.path(input_dir, "..", "..", "..", "MIST", "1.原始输入", "archiveSummary.tsv")
    )
  }

  candidates <- unique(candidates)
  for (cand in candidates) {
    if (file.exists(cand)) return(normalizePath(cand, winslash = "\\", mustWork = TRUE))
  }
  ""
}

input_file <- get_arg("input", ".\\1.原始输入\\proteinGroups.txt")
output_file <- get_arg("output", ".\\2.样本信息与分组\\sample_annotation_auto_R.tsv")
prefix_arg <- get_arg("prefixes", "Intensity ;LFQ intensity ;iBAQ ;MS/MS count ")
raw_suffix <- get_arg("raw_suffix", ".raw")
archive_file <- get_arg("archive", "")
archive_file <- find_archive_file(archive_file, input_file)

if (!file.exists(input_file)) stop("找不到输入文件: ", input_file)

header <- read_header(input_file)
prefixes <- normalize_prefixes(prefix_arg)

rows <- list()
seen <- character()

for (column in header) {
  for (prefix in prefixes) {
    if (!startsWith(column, prefix)) next
    sample_id <- safe_trim(substr(column, nchar(prefix) + 1L, nchar(column)))
    if (!nzchar(sample_id)) next
    key <- paste0(tolower(safe_trim(prefix)), "::", tolower(sample_id))
    if (key %in% seen) next
    seen <- c(seen, key)
    rows[[length(rows) + 1]] <- data.frame(
      sample_id = sample_id,
      raw_file_name = paste0(sample_id, raw_suffix),
      quant_type = safe_trim(prefix),
      source_column = column,
      bait_name = "",
      condition = "",
      is_control = "",
      replicate_group = "",
      keep = "Y",
      notes = "",
      stringsAsFactors = FALSE
    )
    break
  }
}

if (!length(rows)) {
  stop("未识别到任何样本定量列。请检查输入文件是否包含 Intensity / LFQ intensity / iBAQ / MS/MS count 等列。")
}

anno <- do.call(rbind, rows)

if (nzchar(archive_file) && file.exists(archive_file)) {
  archive <- read.delim(archive_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  required_archive <- c("newName", "oldName", "info")
  missing_archive <- setdiff(required_archive, names(archive))
  if (!length(missing_archive)) {
    archive$newName <- safe_trim(archive$newName)
    archive$oldName <- safe_trim(archive$oldName)
    archive$info <- safe_trim(archive$info)
    archive$sample_id <- sub("\\.raw$", "", archive$newName, ignore.case = TRUE)
    archive$bait_name_auto <- safe_trim(extract_bait_from_info(archive$info))
    archive$condition_auto <- safe_trim(extract_condition_from_info(archive$info))
    archive$platewell_auto <- safe_trim(extract_platewell_from_info(archive$info))
    archive$replicate_letter_auto <- safe_trim(extract_replicate_from_name(archive$oldName))
    archive$is_control_auto <- ifelse(is_control_bait(archive$bait_name_auto), "Y", "N")
    archive$replicate_group_auto <- ifelse(
      nzchar(archive$bait_name_auto) & nzchar(archive$condition_auto),
      paste0(archive$bait_name_auto, "__", archive$condition_auto),
      archive$bait_name_auto
    )

    archive_keep <- archive[
      !duplicated(tolower(archive$sample_id)) & nzchar(archive$sample_id),
      c("sample_id", "newName", "bait_name_auto", "condition_auto", "platewell_auto", "replicate_letter_auto", "is_control_auto", "replicate_group_auto"),
      drop = FALSE
    ]

    idx <- match(tolower(anno$sample_id), tolower(archive_keep$sample_id))
    matched <- !is.na(idx)

    anno$raw_file_name[matched] <- ifelse(
      nzchar(archive_keep$newName[idx[matched]]),
      archive_keep$newName[idx[matched]],
      anno$raw_file_name[matched]
    )
    anno$bait_name[matched] <- archive_keep$bait_name_auto[idx[matched]]
    anno$condition[matched] <- archive_keep$condition_auto[idx[matched]]
    anno$is_control[matched] <- archive_keep$is_control_auto[idx[matched]]
    anno$replicate_group[matched] <- archive_keep$replicate_group_auto[idx[matched]]

    anno$notes[matched] <- ifelse(
      nzchar(archive_keep$platewell_auto[idx[matched]]) | nzchar(archive_keep$replicate_letter_auto[idx[matched]]),
      paste0(
        ifelse(nzchar(archive_keep$platewell_auto[idx[matched]]), paste0("platewell=", archive_keep$platewell_auto[idx[matched]]), ""),
        ifelse(nzchar(archive_keep$platewell_auto[idx[matched]]) & nzchar(archive_keep$replicate_letter_auto[idx[matched]]), "; ", ""),
        ifelse(nzchar(archive_keep$replicate_letter_auto[idx[matched]]), paste0("replicate=", archive_keep$replicate_letter_auto[idx[matched]]), "")
      ),
      anno$notes[matched]
    )
  }
}

anno <- anno[order(tolower(anno$quant_type), tolower(anno$sample_id)), , drop = FALSE]

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(anno, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")

message("已生成样本注释表: ", output_file)
message("识别到样本列数: ", nrow(anno))
message("识别前缀: ", paste(safe_trim(prefixes), collapse = ", "))
if (nzchar(archive_file) && file.exists(archive_file)) {
  message("已尝试结合 archiveSummary 自动补充 bait 和分组信息: ", archive_file)
} else {
  message("未提供可用 archiveSummary.tsv，将只生成样本列骨架。")
}
message("仍建议抽查 bait_name、is_control、replicate_group 是否符合真实实验设计。")
