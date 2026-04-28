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
project_dir <- if (nzchar(script_dir)) {
  normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

config_file <- find_existing_file(c(
  get_arg("config", ""),
  file.path(project_dir, "4.MiST输入文件", "mist.yml")
))

main_file <- find_existing_file(c(
  get_arg("main", ""),
  file.path(project_dir, "8.文档与参考", "03_源码或软件", "mist_master", "mist-master", "main.R")
))

stage <- tolower(get_arg("stage", "mist"))
if (!stage %in% c("preprocess", "qc", "mist")) {
  stop("`--stage` 只能取 preprocess / qc / mist")
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("缺少 R 包 yaml，请先运行 脚本/03_安装MiST依赖.R")
}
if (!nzchar(config_file)) stop("找不到 MiST 配置文件 mist.yml")
if (!nzchar(main_file)) stop("找不到 MiST 主程序 main.R")

cfg <- yaml::read_yaml(config_file)
cfg$preprocess$enabled <- 0
cfg$qc$enabled <- 0
cfg$mist$enabled <- 0

if (stage == "preprocess") {
  cfg$preprocess$enabled <- 1
}

if (stage == "qc") {
  cfg$qc$enabled <- 1
  if (is.null(cfg$qc$matrix_file) || !nzchar(cfg$qc$matrix_file)) {
    cfg$qc$matrix_file <- file.path(project_dir, "6.MiST评分结果", "preprocessed_MAT.txt")
  }
}

if (stage == "mist") {
  cfg$mist$enabled <- 1
  if (is.null(cfg$mist$matrix_file) || !nzchar(cfg$mist$matrix_file)) {
    cfg$mist$matrix_file <- file.path(project_dir, "6.MiST评分结果", "preprocessed_MAT.txt")
  }
}

diag_dir <- file.path(project_dir, "6.MiST评分结果", "diagnostic")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

stage_config <- file.path(diag_dir, paste0("mist_", stage, ".yml"))
stage_log <- file.path(diag_dir, paste0("mist_", stage, ".log"))
rscript_bin <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript_bin)) rscript_bin <- file.path(R.home("bin"), "Rscript")

yaml::write_yaml(cfg, stage_config)
writeLines(character(), stage_log, useBytes = TRUE)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(project_dir)

message("MiST 诊断运行开始")
message("阶段: ", stage)
message("项目目录: ", project_dir)
message("配置文件: ", stage_config)
message("日志文件: ", stage_log)

status <- system2(
  rscript_bin,
  c(main_file, "-c", stage_config),
  stdout = stage_log,
  stderr = stage_log
)

if (!identical(status, 0L)) {
  stop("诊断运行失败，退出码: ", status, "\n请查看日志文件: ", stage_log)
}

message("MiST 诊断运行完成")
message("阶段: ", stage)
message("日志文件: ", stage_log)
