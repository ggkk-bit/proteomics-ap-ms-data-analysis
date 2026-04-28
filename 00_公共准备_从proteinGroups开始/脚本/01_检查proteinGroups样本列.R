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

prefix_map <- c(
  "Intensity " = "Intensity",
  "LFQ intensity " = "LFQ intensity",
  "iBAQ " = "iBAQ"
)

collect_by_prefix <- function(prefix_label, prefix_text, all_names) {
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

sample_tables <- Map(
  f = collect_by_prefix,
  prefix_label = names(prefix_map),
  prefix_text = unname(prefix_map),
  MoreArgs = list(all_names = header)
)

sample_columns <- do.call(rbind, sample_tables)

if (nrow(sample_columns) == 0) {
  warning("没有识别到 Intensity / LFQ intensity / iBAQ 样本列。")
  sample_columns <- data.frame(
    quant_type = character(),
    column_name = character(),
    sample_id = character(),
    stringsAsFactors = FALSE
  )
}

write.table(
  sample_columns,
  file = file.path(output_dir, "样本列清单.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

summary_lines <- c(
  paste0("输入文件\t", normalizePath(input_path, winslash = "/", mustWork = FALSE)),
  paste0("样本列总数\t", nrow(sample_columns)),
  paste0("Intensity列数\t", sum(sample_columns$quant_type == "Intensity")),
  paste0("LFQ intensity列数\t", sum(sample_columns$quant_type == "LFQ intensity")),
  paste0("iBAQ列数\t", sum(sample_columns$quant_type == "iBAQ"))
)

writeLines(summary_lines, con = file.path(output_dir, "运行说明.txt"))

message("已输出样本列清单：", file.path(output_dir, "样本列清单.tsv"))
