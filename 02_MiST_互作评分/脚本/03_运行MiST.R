#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

get_script_dir <- function() {
  hit <- grep("^--file=", all_args, value = TRUE)
  if (!length(hit)) return("")
  normalizePath(dirname(sub("^--file=", "", hit[[1]])), winslash = "/", mustWork = FALSE)
}

find_existing_file <- function(candidates) {
  candidates <- unique(candidates[nzchar(candidates)])
  for (cand in candidates) {
    if (file.exists(cand)) return(normalizePath(cand, winslash = "/", mustWork = TRUE))
  }
  ""
}

script_dir <- get_script_dir()
project_dir <- if (nzchar(script_dir)) normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)

config_file <- find_existing_file(c(
  get_arg("config", ""),
  file.path(project_dir, "4.MiST输入文件", "mist.yml"),
  file.path(project_dir, "8.文档与参考", "mist.yml")
))

main_file <- find_existing_file(c(
  get_arg("main", ""),
  file.path(project_dir, "8.文档与参考", "03_源码或软件", "mist_master", "mist-master", "main.R")
))

if (!nzchar(config_file)) stop("找不到 MiST 配置文件，请确认 4.MiST输入文件/mist.yml 是否存在")
if (!nzchar(main_file)) stop("找不到 MiST 主程序 main.R，请确认源码目录是否完整")

required_packages <- c("yaml", "getopt", "optparse", "reshape2", "pheatmap", "RColorBrewer", "ggplot2", "gridExtra", "MESS")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    paste0(
      "当前 R 环境缺少 MiST 依赖包: ", paste(missing_packages, collapse = ", "), "\n",
      "请先运行 脚本/03_安装MiST依赖.R"
    )
  )
}

expected_inputs <- c(
  file.path(project_dir, "4.MiST输入文件", "data.txt"),
  file.path(project_dir, "4.MiST输入文件", "keys.txt")
)
missing_inputs <- expected_inputs[!file.exists(expected_inputs)]
if (length(missing_inputs)) {
  stop(
    paste0(
      "找不到 MiST 输入文件:\n",
      paste(missing_inputs, collapse = "\n"), "\n",
      "请先运行 脚本/02_格式转换脚本.R"
    )
  )
}

rscript_bin <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript_bin)) rscript_bin <- file.path(R.home("bin"), "Rscript")

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(project_dir)

message("开始运行 MiST")
message("项目目录: ", project_dir)
message("配置文件: ", config_file)
message("主程序: ", main_file)

status <- system2(rscript_bin, c(main_file, "-c", config_file))
if (!identical(status, 0L)) stop("MiST 运行失败，退出状态码: ", status)

message("MiST 运行完成")
message("请查看目录: ", file.path(project_dir, "6.MiST评分结果"))