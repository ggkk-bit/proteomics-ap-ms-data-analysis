suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE)
})

args <- commandArgs(trailingOnly = TRUE)

default_input <- "C:/Users/dakey/Desktop/01_进行中/古老蛋白算法调研/下载的文献进行测试/proteinGroups.txt"
default_output_dir <- "C:/Users/dakey/Desktop/01_进行中/古老蛋白算法调研/下载的文献进行测试/多算法分析工作区/00_公共准备_从proteinGroups开始/输出"

input_path <- if (length(args) >= 1) args[[1]] else default_input
output_dir <- if (length(args) >= 2) args[[2]] else default_output_dir

if (!file.exists(input_path)) {
  stop("找不到 proteinGroups.txt: ", input_path)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

header <- names(read.delim(input_path, nrows = 0, sep = "\t", check.names = FALSE))

extract_samples <- function(prefix_label, prefix_text, all_names) {
  hits <- all_names[startsWith(all_names, prefix_text)]
  if (length(hits) == 0) {
    return(data.frame())
  }

  data.frame(
    quant_type = prefix_label,
    column_name = hits,
    sample_id = substring(hits, nchar(prefix_text) + 1L),
    stringsAsFactors = FALSE
  )
}

sample_columns <- do.call(
  rbind,
  list(
    extract_samples("Intensity", "Intensity ", header),
    extract_samples("LFQ intensity", "LFQ intensity ", header),
    extract_samples("iBAQ", "iBAQ ", header)
  )
)

if (nrow(sample_columns) == 0) {
  stop("没有识别到可用于生成模板的样本列。")
}

sample_columns <- sample_columns[order(sample_columns$sample_id, sample_columns$quant_type), ]
sample_columns <- sample_columns[!duplicated(sample_columns$sample_id), ]

annotation_template <- data.frame(
  sample_id = sample_columns$sample_id,
  raw_file_name = "",
  sample_type = "experiment",
  bait_name = "",
  replicate_group = "",
  replicate_id = "",
  is_control = "no",
  keep = "yes",
  batch = "",
  note = "",
  stringsAsFactors = FALSE
)

write.table(
  annotation_template,
  file = file.path(output_dir, "样本注释模板.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

guide_lines <- c(
  "字段说明：",
  "sample_id\t必须与 proteinGroups.txt 中的样本列后缀一致",
  "raw_file_name\t原始 raw 文件名，可回填",
  "sample_type\texperiment / control / blank 等",
  "bait_name\t诱饵名称，同一 bait 的重复应一致",
  "replicate_group\t同一 bait 的重复组，例如 BaitA",
  "replicate_id\t重复编号，例如 R1/R2/R3",
  "is_control\tyes/no",
  "keep\tyes/no",
  "batch\t批次信息，可选",
  "note\t补充说明"
)

writeLines(guide_lines, con = file.path(output_dir, "样本注释模板说明.txt"))

message("已生成样本注释模板：", file.path(output_dir, "样本注释模板.tsv"))
