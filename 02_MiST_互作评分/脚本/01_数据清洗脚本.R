#!/usr/bin/env Rscript

# MiST 数据清洗脚本
# 作用：
# 1. 读取 MaxQuant proteinGroups.txt 或同结构表格
# 2. 删除 Reverse / Potential contaminant / Only identified by site
# 3. 按可配置阈值过滤唯一肽段数
# 4. 可选删除所有样本强度都为 0/NA 的蛋白
# 5. 输出清洗后的通用表，供后续格式转换使用

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

pick_col <- function(candidates, names_vec) {
  hit <- intersect(candidates, names_vec)
  if (length(hit)) hit[[1]] else NA_character_
}

nonempty_count <- function(x) {
  y <- trimws(as.character(x))
  sum(!is.na(x) & (y != ""))
}

flagged_value <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y %in% c("true", "yes", "1", "+", "y")
}

input_file <- get_arg("input", ".\\1.原始输入\\proteinGroups.txt")
output_file <- get_arg("output", ".\\3.数据清洗与预处理\\proteinGroups_已清洗.tsv")
log_file <- get_arg("log", ".\\5.质控\\数据清洗记录.tsv")
min_unique_peptides <- as.integer(get_arg("min_unique_peptides", "1"))
drop_all_zero_intensity <- tolower(get_arg("drop_all_zero_intensity", "true")) %in% c("true", "1", "yes", "y")
require_protein_id <- tolower(get_arg("require_protein_id", "true")) %in% c("true", "1", "yes", "y")

if (!file.exists(input_file)) stop("找不到输入文件: ", input_file)

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
if (!nrow(dat)) stop("输入文件没有数据行: ", input_file)

protein_id_col <- pick_col(c("Majority protein IDs", "Protein IDs"), names(dat))
unique_pep_col <- pick_col(c("Unique peptides", "Razor + unique peptides"), names(dat))
filter_cols <- c("Reverse", "Potential contaminant", "Only identified by site")

if (require_protein_id && is.na(protein_id_col)) {
  stop("缺少蛋白 ID 字段。至少需要 'Majority protein IDs' 或 'Protein IDs' 之一")
}

keep <- rep(TRUE, nrow(dat))

for (col in filter_cols) {
  if (col %in% names(dat)) {
    is_flagged <- flagged_value(dat[[col]])
    is_flagged[is.na(is_flagged)] <- FALSE
    keep <- keep & !is_flagged
  }
}

if (!is.na(unique_pep_col) && !is.na(min_unique_peptides)) {
  uniq_val <- suppressWarnings(as.numeric(dat[[unique_pep_col]]))
  keep <- keep & (!is.na(uniq_val)) & (uniq_val >= min_unique_peptides)
}

if (!is.na(protein_id_col)) {
  pid <- trimws(as.character(dat[[protein_id_col]]))
  keep <- keep & (!is.na(pid)) & (pid != "")
}

intensity_cols <- grep("^Intensity ", names(dat), value = TRUE)
if (drop_all_zero_intensity && length(intensity_cols) > 0) {
  intensity_mat <- as.matrix(as.data.frame(lapply(dat[intensity_cols], function(x) suppressWarnings(as.numeric(x)))))
  zero_or_missing <- is.na(intensity_mat) | (intensity_mat == 0)
  zero_or_missing[is.na(zero_or_missing)] <- TRUE
  all_zero <- rowSums(zero_or_missing) == length(intensity_cols)
  keep <- keep & !all_zero
}

keep[is.na(keep)] <- FALSE
cleaned <- dat[keep, , drop = FALSE]
if (!nrow(cleaned)) {
  stop("清洗后没有剩余蛋白。请检查过滤条件、唯一肽段阈值或强度列是否全为 0/NA。")
}

if (!is.na(protein_id_col) && nonempty_count(cleaned[[protein_id_col]]) == 0) {
  stop("清洗后蛋白 ID 字段全为空: ", protein_id_col)
}

positive_intensity_count <- NA_integer_
if (length(intensity_cols) > 0) {
  positive_intensity_count <- sum(vapply(
    cleaned[intensity_cols],
    function(x) sum(suppressWarnings(as.numeric(x)) > 0, na.rm = TRUE),
    numeric(1)
  ))
  if (drop_all_zero_intensity && positive_intensity_count == 0) {
    stop("清洗后所有强度列都没有正数，无法生成 MiST 输入。")
  }
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

input_lines <- readLines(input_file, warn = FALSE, encoding = "UTF-8")
data_lines <- input_lines[-1]

if (length(data_lines) == nrow(dat)) {
  # 超宽 MaxQuant 表用 write.table() 在部分 Windows/R 组合下可能写成空字段。
  # 这里原样筛选原始行，最大限度保留 proteinGroups.txt 的列和值。
  writeLines(c(input_lines[1], data_lines[keep]), con = output_file, useBytes = TRUE)
} else {
  warning("输入物理行数与 R 解析行数不一致，改用带引号写出。请检查原始文件中是否有嵌入换行。")
  write.table(cleaned, output_file, sep = "\t", quote = TRUE, row.names = FALSE, col.names = TRUE, na = "")
}

log_df <- data.frame(
  input_file = normalizePath(input_file, winslash = "\\", mustWork = FALSE),
  output_file = normalizePath(output_file, winslash = "\\", mustWork = FALSE),
  protein_id_col = ifelse(is.na(protein_id_col), "", protein_id_col),
  unique_pep_col = ifelse(is.na(unique_pep_col), "", unique_pep_col),
  intensity_col_count = length(intensity_cols),
  positive_intensity_count = positive_intensity_count,
  input_rows = nrow(dat),
  output_rows = nrow(cleaned),
  removed_rows = nrow(dat) - nrow(cleaned),
  min_unique_peptides = min_unique_peptides,
  drop_all_zero_intensity = drop_all_zero_intensity,
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)
write.table(log_df, log_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")

message("MiST 数据清洗完成")
message("输入文件: ", normalizePath(input_file, winslash = "\\", mustWork = FALSE))
message("输出文件: ", normalizePath(output_file, winslash = "\\", mustWork = FALSE))
message("日志文件: ", normalizePath(log_file, winslash = "\\", mustWork = FALSE))
message("输入行数: ", nrow(dat))
message("输出行数: ", nrow(cleaned))
message("移除行数: ", nrow(dat) - nrow(cleaned))
if (!is.na(positive_intensity_count)) message("保留表中正强度值数量: ", positive_intensity_count)
