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

log_file <- get_arg("log", file.path(project_dir, "6.MiST评分结果", "mist_run.log"))
status_file <- get_arg("status", file.path(project_dir, "6.MiST评分结果", "mist_status.tsv"))

if (!nzchar(config_file)) stop("找不到 MiST 配置文件，请确认 4.MiST输入文件/mist.yml 是否存在")
if (!nzchar(main_file)) stop("找不到 MiST 主程序 main.R，请确认源码目录是否完整")

required_packages <- c("yaml", "getopt", "optparse", "reshape2", "pheatmap", "RColorBrewer", "ggplot2", "gridExtra", "MESS")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    paste0(
      "当前 R 环境缺少 MiST 依赖包: ",
      paste(missing_packages, collapse = ", "),
      "\n请先运行 脚本/03_安装MiST依赖.R"
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
      paste(missing_inputs, collapse = "\n"),
      "\n请先运行 脚本/02_格式转换脚本.R"
    )
  )
}

rscript_bin <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript_bin)) rscript_bin <- file.path(R.home("bin"), "Rscript")

dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

can_write_log <- function(path) {
  ok <- tryCatch({
    con <- suppressWarnings(file(path, open = "w"))
    close(con)
    TRUE
  }, error = function(e) FALSE)
  ok
}

if (!can_write_log(log_file)) {
  fallback_log <- file.path(
    dirname(log_file),
    paste0("mist_run_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
  )
  message("默认日志文件被占用，改用新日志文件: ", normalizePath(fallback_log, winslash = "/", mustWork = FALSE))
  log_file <- fallback_log
  writeLines(character(), log_file, useBytes = TRUE)
}

write_status <- function(status, exit_code = "", message_text = "") {
  clean_status_field <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    gsub("[\t\r\n]+", " ", x)
  }
  status_df <- data.frame(
    status = clean_status_field(status),
    exit_code = clean_status_field(exit_code),
    message = clean_status_field(message_text),
    start_time = clean_status_field(start_time),
    update_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    project_dir = clean_status_field(project_dir),
    config_file = clean_status_field(config_file),
    main_file = clean_status_field(main_file),
    log_file = clean_status_field(normalizePath(log_file, winslash = "/", mustWork = FALSE)),
    stringsAsFactors = FALSE
  )
  write.table(status_df, status_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
}

append_log <- function(line) {
  cat(line, file = log_file, sep = "\n", append = TRUE)
}

run_with_live_log <- function(rscript_bin, main_file, config_file) {
  quote_arg <- function(x) {
    if (.Platform$OS.type == "windows") {
      if (file.exists(x)) {
        utils::shortPathName(normalizePath(x, winslash = "\\", mustWork = TRUE))
      } else {
        x
      }
    } else {
      shQuote(x)
    }
  }
  cmd <- paste(quote_arg(rscript_bin), quote_arg(main_file), "-c", quote_arg(config_file), "2>&1")
  con <- pipe(cmd, open = "r")
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  last_status_update <- Sys.time()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE, encoding = "UTF-8")
    if (!length(line)) break
    append_log(line)

    now <- Sys.time()
    if (as.numeric(difftime(now, last_status_update, units = "secs")) >= 30) {
      write_status("RUNNING", "", paste0("最近日志: ", substr(line, 1, 120)))
      last_status_update <- now
    }
  }

  status <- close(con)
  on.exit(NULL, add = FALSE)
  if (is.null(status)) 0L else as.integer(status)
}

start_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
write_status("RUNNING", "", "MiST 已启动")

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(project_dir)

message("开始运行 MiST")
message("项目目录: ", project_dir)
message("配置文件: ", config_file)
message("主程序: ", main_file)
message("日志文件: ", normalizePath(log_file, winslash = "/", mustWork = FALSE))
message("状态文件: ", normalizePath(status_file, winslash = "/", mustWork = FALSE))
message("另开 PowerShell 可查看状态:")
message("& '", rscript_bin, "' '.\\脚本\\03_查看MiST状态.R'")

status <- run_with_live_log(rscript_bin, main_file, config_file)

if (!identical(status, 0L)) {
  write_status("FAILED", status, "MiST 运行失败")
  stop(
    "MiST 运行失败，退出状态码: ", status,
    "\n请查看日志文件: ", normalizePath(log_file, winslash = "/", mustWork = FALSE)
  )
}

write_status("FINISHED", 0, "MiST 运行完成")
message("MiST 运行完成")
message("结果目录: ", file.path(project_dir, "6.MiST评分结果"))
message("日志文件: ", normalizePath(log_file, winslash = "/", mustWork = FALSE))
