suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE)
})

args <- commandArgs(trailingOnly = TRUE)

default_input <- "C:/Users/dakey/Desktop/01_进行中/古老蛋白算法调研/下载的文献进行测试/proteinGroups.txt"
default_annotation <- "C:/Users/dakey/Desktop/01_进行中/古老蛋白算法调研/下载的文献进行测试/多算法分析工作区/00_公共准备_从proteinGroups开始/输出/样本注释模板.tsv"
default_output_dir <- "C:/Users/dakey/Desktop/01_进行中/古老蛋白算法调研/下载的文献进行测试/多算法分析工作区/00_公共准备_从proteinGroups开始/输出"

input_path <- if (length(args) >= 1) args[[1]] else default_input
annotation_path <- if (length(args) >= 2) args[[2]] else default_annotation
output_dir <- if (length(args) >= 3) args[[3]] else default_output_dir

if (!file.exists(input_path)) {
  stop("找不到 proteinGroups.txt: ", input_path)
}

if (!file.exists(annotation_path)) {
  stop("找不到样本注释文件: ", annotation_path)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pg <- read.delim(input_path, sep = "\t", check.names = FALSE)
ann <- read.delim(annotation_path, sep = "\t", check.names = FALSE)

required_ann <- c(
  "sample_id", "sample_type", "bait_name", "replicate_group",
  "replicate_id", "is_control", "keep"
)

missing_ann <- setdiff(required_ann, names(ann))
if (length(missing_ann) > 0) {
  stop("样本注释缺少字段: ", paste(missing_ann, collapse = ", "))
}

flag_column <- function(df, nm) {
  if (!nm %in% names(df)) {
    return(rep(FALSE, nrow(df)))
  }
  x <- trimws(as.character(df[[nm]]))
  x == "+"
}

pg <- pg[
  !flag_column(pg, "Reverse") &
    !flag_column(pg, "Potential contaminant") &
    !flag_column(pg, "Only identified by site"),
]

if (!"Majority protein IDs" %in% names(pg)) {
  stop("proteinGroups.txt 缺少 `Majority protein IDs` 列。")
}

extract_first_token <- function(x) {
  x <- as.character(x)
  x <- ifelse(is.na(x), "", x)
  sub(";.*$", "", x)
}

prey_id <- extract_first_token(pg[["Majority protein IDs"]])
gene_name <- if ("Gene names" %in% names(pg)) extract_first_token(pg[["Gene names"]]) else ""
protein_name <- if ("Protein names" %in% names(pg)) as.character(pg[["Protein names"]]) else ""
sequence_length <- if ("Sequence length" %in% names(pg)) pg[["Sequence length"]] else NA

pick_quant_column <- function(sample_id, header) {
  candidates <- c(
    paste0("LFQ intensity ", sample_id),
    paste0("Intensity ", sample_id),
    paste0("iBAQ ", sample_id)
  )
  hit <- candidates[candidates %in% header]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[[1]]
}

ann$quant_column <- vapply(ann$sample_id, pick_quant_column, character(1), header = names(pg))
ann$quant_type <- ifelse(
  startsWith(ann$quant_column, "LFQ intensity "), "LFQ intensity",
  ifelse(startsWith(ann$quant_column, "Intensity "), "Intensity",
         ifelse(startsWith(ann$quant_column, "iBAQ "), "iBAQ", ""))
)

usable_ann <- ann[ann$keep == "yes" & !is.na(ann$quant_column) & nzchar(ann$quant_column), ]

if (nrow(usable_ann) == 0) {
  stop("没有可用样本。请先补全样本注释模板中的 keep / bait_name / replicate 信息。")
}

long_tables <- lapply(seq_len(nrow(usable_ann)), function(i) {
  meta <- usable_ann[i, ]
  values <- suppressWarnings(as.numeric(pg[[meta$quant_column]]))
  keep_idx <- !is.na(values) & values > 0

  data.frame(
    sample_id = meta$sample_id,
    sample_type = meta$sample_type,
    bait_name = meta$bait_name,
    replicate_group = meta$replicate_group,
    replicate_id = meta$replicate_id,
    is_control = meta$is_control,
    quant_type = meta$quant_type,
    prey_id = prey_id[keep_idx],
    gene_name = gene_name[keep_idx],
    protein_name = protein_name[keep_idx],
    sequence_length = sequence_length[keep_idx],
    abundance = values[keep_idx],
    stringsAsFactors = FALSE
  )
})

long_table <- do.call(rbind, long_tables)

write.table(
  long_table,
  file = file.path(output_dir, "样本级长表.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

sample_summary <- aggregate(
  abundance ~ sample_id + bait_name + is_control + quant_type,
  data = long_table,
  FUN = length
)
names(sample_summary)[names(sample_summary) == "abundance"] <- "prey_count"

write.table(
  sample_summary,
  file = file.path(output_dir, "样本级长表_每样本统计.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("已生成样本级长表：", file.path(output_dir, "样本级长表.tsv"))
