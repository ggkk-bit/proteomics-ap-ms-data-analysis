#!/usr/bin/env Rscript

# 统一样本注释表生成脚本
# 从 proteinGroups 表头识别定量列，并可从 archiveSummary.tsv 自动回填 bait、孔位、重复等信息。

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(FALSE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

safe_trim <- function(x) trimws(as.character(x))

script_file <- sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE)) else "."
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

read_header <- function(path) {
  con <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  line <- readLines(con, n = 1, warn = FALSE)
  if (!length(line)) stop("输入文件为空: ", path)
  strsplit(line, "\t", fixed = TRUE)[[1]]
}

normalize_prefixes <- function(text) {
  parts <- safe_trim(strsplit(text, ";", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  ifelse(grepl("[ .]$", parts), parts, paste0(parts, " "))
}

parse_info_bait <- function(info) {
  info <- safe_trim(info)
  out <- rep("", length(info))
  hit <- grepl("Bait:", info)
  out[hit] <- sub("^.*Bait:([^()]+).*$", "\\1", info[hit])
  safe_trim(out)
}

parse_info_plate <- function(info) {
  info <- safe_trim(info)
  out <- rep("", length(info))
  hit <- grepl("\\(", info)
  out[hit] <- sub("^.*\\(([^,()]+).*$", "\\1", info[hit])
  safe_trim(out)
}

parse_info_well <- function(info) {
  info <- safe_trim(info)
  out <- rep("", length(info))
  hit <- grepl(",", info)
  out[hit] <- sub("^.*,[[:space:]]*([^,()]+)\\).*$", "\\1", info[hit])
  safe_trim(out)
}

parse_replicate <- function(old_name) {
  old_name <- safe_trim(old_name)
  out <- rep("", length(old_name))
  m <- regexpr("_([A-Z])_", old_name)
  ok <- m > 0
  out[ok] <- substring(old_name[ok], m[ok] + 1L, m[ok] + 1L)
  out
}

is_control_bait <- function(bait) {
  tolower(safe_trim(bait)) %in% c("ctrl", "control", "mock", "empty", "vector", "beads", "igg", "gfp", "flag", "ha")
}

input_file <- get_arg("input", file.path(project_dir, "1.原始输入", "proteinGroups.txt"))
archive_file <- get_arg("archive", file.path(project_dir, "1.原始输入", "archiveSummary.tsv"))
output_file <- get_arg("output", file.path(project_dir, "2.样本注释", "sample_annotation_master.tsv"))
prefix_text <- get_arg("prefixes", "Intensity ;LFQ intensity ;iBAQ ;MS/MS count ;MS/MS Count ;Intensity.")
raw_suffix <- get_arg("raw_suffix", ".raw")

if (!file.exists(input_file)) stop("找不到输入文件: ", input_file)
header <- read_header(input_file)
prefixes <- normalize_prefixes(prefix_text)

rows <- list()
seen <- character()
for (column in header) {
  for (prefix in prefixes) {
    if (!startsWith(column, prefix)) next
    sample_id <- safe_trim(substr(column, nchar(prefix) + 1L, nchar(column)))
    if (!nzchar(sample_id)) next
    key <- paste(tolower(prefix), tolower(sample_id), sep = "::")
    if (key %in% seen) next
    seen <- c(seen, key)
    rows[[length(rows) + 1]] <- data.frame(
      sample_id = sample_id,
      raw_file_name = paste0(sample_id, raw_suffix),
      quant_type = safe_trim(prefix),
      source_column = column,
      plate = "",
      well = "",
      replicate_letter = "",
      ip_name = sample_id,
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

if (!length(rows)) stop("未识别到任何样本定量列，请检查 prefixes 或 proteinGroups 表头。")
anno <- do.call(rbind, rows)

if (file.exists(archive_file)) {
  archive <- read.delim(archive_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (all(c("oldName", "newName", "info") %in% names(archive))) {
    archive$sample_id <- sub("\\.raw$", "", safe_trim(archive$newName), ignore.case = TRUE)
    archive$raw_file_name <- safe_trim(archive$newName)
    archive$bait_name <- parse_info_bait(archive$info)
    archive$plate <- parse_info_plate(archive$info)
    archive$well <- parse_info_well(archive$info)
    archive$replicate_letter <- parse_replicate(archive$oldName)
    archive$condition <- archive$plate
    archive$is_control <- ifelse(is_control_bait(archive$bait_name), "Y", "N")
    archive$replicate_group <- ifelse(nzchar(archive$bait_name) & nzchar(archive$condition), paste0(archive$bait_name, "__", archive$condition), archive$bait_name)
    archive <- archive[nzchar(archive$sample_id) & !duplicated(tolower(archive$sample_id)), , drop = FALSE]
    idx <- match(tolower(anno$sample_id), tolower(archive$sample_id))
    matched <- !is.na(idx)
    for (col in c("raw_file_name", "plate", "well", "replicate_letter", "bait_name", "condition", "is_control", "replicate_group")) {
      anno[[col]][matched] <- archive[[col]][idx[matched]]
    }
    anno$ip_name[matched] <- anno$sample_id[matched]
    anno$notes[matched] <- "auto_filled_from_archiveSummary"
  } else {
    warning("archiveSummary.tsv 缺少 oldName/newName/info 字段，跳过自动回填。")
  }
}

anno <- anno[order(tolower(anno$quant_type), tolower(anno$sample_id)), , drop = FALSE]
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(anno, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")

message("统一样本注释表已生成: ", output_file)
message("识别到样本定量列数: ", nrow(anno))
message("请人工抽查 bait_name、is_control、replicate_group 是否符合实验设计。")
