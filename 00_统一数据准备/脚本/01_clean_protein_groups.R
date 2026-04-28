#!/usr/bin/env Rscript

# 统一 proteinGroups 清洗脚本
# 作用：
# 1. 删除 MaxQuant 标记的 Reverse / Potential contaminant / Only identified by site
# 2. 删除缺少蛋白 ID 的记录
# 3. 可选删除所有定量列均为 0 或 NA 的记录
# 4. 输出清洗后的 proteinGroups_cleaned.tsv 和 JSON 清洗报告

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(FALSE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

`%||%` <- function(x, y) {
  if (length(x) == 0 || is.na(x) || is.null(x)) y else x
}

script_file <- sub("^--file=", "", grep("^--file=", all_args, value = TRUE)[1])
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE)) else "."
project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)

json_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub('"', '\\"', x)
  x <- gsub("\r", "\\\\r", x)
  x <- gsub("\n", "\\\\n", x)
  x
}

to_json <- function(x, indent = 0) {
  pad <- paste(rep(" ", indent), collapse = "")
  pad2 <- paste(rep(" ", indent + 2), collapse = "")
  if (is.list(x) && is.null(names(x))) {
    return(paste0("[", paste(vapply(x, to_json, character(1), indent = indent + 2), collapse = ", "), "]"))
  }
  if (is.list(x)) {
    parts <- vapply(names(x), function(nm) {
      paste0(pad2, '"', json_escape(nm), '": ', to_json(x[[nm]], indent + 2))
    }, character(1))
    return(paste0("{\n", paste(parts, collapse = ",\n"), "\n", pad, "}"))
  }
  if (length(x) > 1) return(paste0("[", paste(vapply(x, to_json, character(1), indent = indent), collapse = ", "), "]"))
  if (is.logical(x)) return(ifelse(is.na(x), "null", ifelse(x, "true", "false")))
  if (is.numeric(x) || is.integer(x)) return(ifelse(is.na(x), "null", as.character(x)))
  if (is.na(x)) return("null")
  paste0('"', json_escape(x), '"')
}

truthy <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y %in% c("+", "1", "true", "yes", "y", "t")
}

first_col <- function(df, aliases) {
  lookup <- setNames(names(df), tolower(trimws(names(df))))
  hit <- lookup[tolower(trimws(aliases))]
  hit <- hit[!is.na(hit)]
  if (length(hit)) unname(hit[1]) else NA_character_
}

input_file <- get_arg("input", file.path(project_dir, "1.原始输入", "proteinGroups.txt"))
output_file <- get_arg("output", file.path(project_dir, "3.清洗后数据", "proteinGroups_cleaned.tsv"))
report_file <- get_arg("report", file.path(project_dir, "3.清洗后数据", "cleaning_report.json"))
drop_all_zero <- tolower(get_arg("drop_all_zero_quant", "true")) %in% c("true", "1", "yes", "y")
quant_prefixes <- strsplit(get_arg("quant_prefixes", "Intensity ;LFQ intensity ;iBAQ ;MS/MS count ;MS/MS Count ;Intensity."), ";", fixed = TRUE)[[1]]
quant_prefixes <- trimws(quant_prefixes)
quant_prefixes <- quant_prefixes[nzchar(quant_prefixes)]

if (!file.exists(input_file)) stop("找不到输入文件: ", input_file)

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
if (!nrow(dat)) stop("输入文件没有数据行: ", input_file)

initial_rows <- nrow(dat)
keep <- rep(TRUE, nrow(dat))
removed <- list()

filter_specs <- list(
  reverse = c("Reverse", "反向匹配"),
  potential_contaminant = c("Potential contaminant", "潜在污染物"),
  only_identified_by_site = c("Only identified by site", "仅位点鉴定蛋白")
)

for (filter_name in names(filter_specs)) {
  col <- first_col(dat, filter_specs[[filter_name]])
  if (is.na(col)) {
    removed[[filter_name]] <- 0L
    next
  }
  mask <- keep & truthy(dat[[col]])
  removed[[filter_name]] <- as.integer(sum(mask, na.rm = TRUE))
  keep <- keep & !truthy(dat[[col]])
}

protein_col <- first_col(dat, c("Majority protein IDs", "Protein IDs"))
if (is.na(protein_col)) stop("缺少蛋白 ID 字段：Majority protein IDs 或 Protein IDs")
missing_id <- keep & !nzchar(trimws(as.character(dat[[protein_col]])))
removed$missing_protein_id <- as.integer(sum(missing_id, na.rm = TRUE))
keep <- keep & !missing_id

quant_cols <- character()
for (prefix in quant_prefixes) {
  quant_cols <- unique(c(quant_cols, names(dat)[startsWith(names(dat), prefix)]))
}

if (drop_all_zero && length(quant_cols) > 0) {
  quant <- as.data.frame(lapply(dat[quant_cols], function(x) suppressWarnings(as.numeric(x))))
  quant[is.na(quant)] <- 0
  all_zero <- keep & rowSums(quant) == 0
  removed$all_zero_or_missing_quant <- as.integer(sum(all_zero, na.rm = TRUE))
  keep <- keep & !all_zero
} else {
  removed$all_zero_or_missing_quant <- 0L
}

keep[is.na(keep)] <- FALSE
cleaned <- dat[keep, , drop = FALSE]
if (!nrow(cleaned)) stop("清洗后没有剩余蛋白，请检查过滤条件和定量列。")

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)
write.table(cleaned, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")

report <- list(
  input = normalizePath(input_file, winslash = "/", mustWork = FALSE),
  output = normalizePath(output_file, winslash = "/", mustWork = FALSE),
  protein_id_column = protein_col,
  quant_column_count = length(quant_cols),
  drop_all_zero_quant = drop_all_zero,
  initial_rows = initial_rows,
  remaining_rows = nrow(cleaned),
  removed_rows = initial_rows - nrow(cleaned),
  removed = removed,
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
writeLines(to_json(report), report_file, useBytes = TRUE)

message("统一清洗完成: ", output_file)
message("清洗报告: ", report_file)
message("输入行数: ", initial_rows, "; 输出行数: ", nrow(cleaned), "; 删除行数: ", initial_rows - nrow(cleaned))
