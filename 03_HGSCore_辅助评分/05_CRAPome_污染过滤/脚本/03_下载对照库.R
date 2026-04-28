#!/usr/bin/env Rscript
# CRAPome/REPRINT 对照库准备。在线库无法稳定自动下载时，生成可替换的示例库。

script_path <- commandArgs(FALSE)
file_arg <- script_path[startsWith(script_path, "--file=")]
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "common_crapome.R"))

args <- commandArgs(trailingOnly = TRUE)
message_header("CRAPome 对照库准备")

input_long <- get_arg(args, "input", "./4.CRAPome输入文件/crapome_input_long.tsv")
output_file <- get_arg(args, "output", "./4.CRAPome输入文件/control_library_example.tsv")
mode <- get_arg(args, "mode", "example")
url <- get_arg(args, "url", "")

if (mode == "download" && nzchar(url)) {
  suppressPackageStartupMessages(library(httr))
  resp <- httr::GET(url)
  if (httr::status_code(resp) >= 300) stop("下载失败，HTTP 状态码: ", httr::status_code(resp))
  ensure_dir(output_file)
  writeBin(httr::content(resp, "raw"), output_file)
  cat("已下载对照库: ", output_file, "\n", sep = "")
} else {
  if (file.exists(input_long)) prey <- sort(unique(read_tsv(input_long)$`Prey Name`)) else prey <- c("P02768", "P60709", "P62937", "P02787", "P68104", "Q9Y490", "P63261", "P00533", "P04637", "P38398")
  set.seed(1)
  lib <- data.frame(Prey_ID = prey, Control_Frequency = pmin(0.95, pmax(0.01, rbeta(length(prey), 1.5, 8))), Total_Controls = 411, Control_Spectral_Count = sample(1:5000, length(prey), replace = TRUE), Control_Total_Spectral_Count = 1000000, Control_Type = "all", stringsAsFactors = FALSE)
  write_tsv(lib, output_file)
  cat("已生成示例对照库: ", output_file, "\n", sep = "")
}

writeLines(c("CRAPome 对照库说明", "真实分析建议从 https://www.reprint-apms.org/ 手动导出对应物种/实验类型的对照统计表。", "本脚本示例库仅用于打通流程，不应用于正式结论。", "最低必需列: Prey_ID, Control_Frequency, Total_Controls。", "推荐附加列: Control_Spectral_Count, Control_Total_Spectral_Count, Control_Type。"), "./5.质控/03_对照库说明.txt", useBytes = TRUE)
