#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

as_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "y")
}

get_script_dir <- function() {
  hit <- grep("^--file=", all_args, value = TRUE)
  if (!length(hit)) return("")
  normalizePath(dirname(sub("^--file=", "", hit[[1]])), winslash = "/", mustWork = FALSE)
}

tail_lines <- function(path, n = 30) {
  if (!file.exists(path) || file.info(path)$size == 0) return(character())
  x <- tryCatch(
    suppressWarnings(readLines(path, warn = FALSE, encoding = "UTF-8")),
    error = function(e) {
      out <- character()
      attr(out, "read_error") <- conditionMessage(e)
      attr(out, "read_error_path") <- path
      out
    }
  )
  if (!is.null(attr(x, "read_error"))) return(x)
  utils::tail(x, n)
}

latest_existing_log <- function(log_file, result_dir) {
  candidates <- character()
  if (file.exists(log_file)) candidates <- c(candidates, log_file)
  extra <- list.files(result_dir, pattern = "^mist_run.*\\.log$", full.names = TRUE)
  candidates <- unique(c(candidates, extra))
  if (!length(candidates)) return(log_file)
  info <- file.info(candidates)
  candidates[which.max(info$mtime)]
}

log_has_scoring_finished <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(FALSE)
  lines <- tryCatch(
    suppressWarnings(readLines(path, warn = FALSE, encoding = "UTF-8")),
    error = function(e) character()
  )
  any(grepl("SCORING FINISHED", lines, ignore.case = TRUE))
}

find_completed_run <- function(result_dir) {
  result_file <- file.path(result_dir, "preprocessed_NoC_MAT_MIST.txt")
  logs <- list.files(result_dir, pattern = "^mist_run.*\\.log$", full.names = TRUE)
  if (!length(logs)) {
    return(list(found = FALSE, log_file = "", result_file = result_file, result_rows = NA_integer_))
  }

  info <- file.info(logs)
  logs <- logs[order(info$mtime, decreasing = TRUE)]
  finished <- vapply(logs, log_has_scoring_finished, logical(1))
  if (!any(finished)) {
    return(list(found = FALSE, log_file = "", result_file = result_file, result_rows = NA_integer_))
  }

  completed_log <- logs[which(finished)[1]]
  result_rows <- NA_integer_
  if (file.exists(result_file)) {
    result_rows <- tryCatch(
      max(0L, as.integer(length(readLines(result_file, warn = FALSE)) - 1L)),
      error = function(e) NA_integer_
    )
  }

  list(
    found = TRUE,
    log_file = completed_log,
    result_file = result_file,
    result_rows = result_rows
  )
}

detect_stage <- function(lines) {
  txt <- paste(lines, collapse = "\n")
  if (grepl("SCORING FINISHED", txt, ignore.case = TRUE)) return("FINISHED")
  if (grepl("Execution halted|Error", txt, ignore.case = TRUE)) return("ERROR")
  if (grepl(">> MIST", txt, fixed = TRUE)) return("MIST 打分")
  if (grepl(">> QUALITY CONTROL", txt, fixed = TRUE)) return("QUALITY CONTROL")
  if (grepl("CONVERTING TO MATRIX", txt, fixed = TRUE)) return("构建矩阵")
  if (grepl("MERGING KEYS WITH DATA", txt, fixed = TRUE)) return("合并 keys 和 data")
  if (grepl("READING FILES", txt, fixed = TRUE)) return("读取输入")
  if (grepl("LOADING DEPENDENCIES", txt, fixed = TRUE)) return("加载依赖")
  "未知/尚未写日志"
}

find_mist_processes <- function(project_dir) {
  if (.Platform$OS.type != "windows") return(character())
  cmd <- "Get-Process Rscript,Rterm -ErrorAction SilentlyContinue | Select-Object Id,ProcessName,CPU,WS | ConvertTo-Csv -NoTypeInformation"
  out <- suppressWarnings(system2("powershell", c("-NoProfile", "-Command", cmd), stdout = TRUE, stderr = FALSE))
  if (!length(out)) return(character())
  parsed <- tryCatch(read.csv(text = paste(out, collapse = "\n"), stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(parsed) || !nrow(parsed)) return(character())
  parsed <- parsed[parsed$Id != Sys.getpid(), , drop = FALSE]
  if (!nrow(parsed)) return(character())
  paste0(
    "PID=", parsed$Id,
    "\t", parsed$ProcessName,
    "\tCPU=", round(as.numeric(parsed$CPU), 1),
    "\t内存MB=", round(as.numeric(parsed$WS) / 1024 / 1024, 1)
  )
}

get_r_process_table <- function() {
  if (.Platform$OS.type != "windows") return(data.frame())
  cmd <- "Get-Process Rscript,Rterm -ErrorAction SilentlyContinue | Select-Object Id,ProcessName,CPU,WS | ConvertTo-Csv -NoTypeInformation"
  out <- suppressWarnings(system2("powershell", c("-NoProfile", "-Command", cmd), stdout = TRUE, stderr = FALSE))
  if (!length(out)) return(data.frame())
  parsed <- tryCatch(read.csv(text = paste(out, collapse = "\n"), stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(parsed) || !nrow(parsed)) return(data.frame())
  parsed <- parsed[parsed$Id != Sys.getpid(), , drop = FALSE]
  if (!nrow(parsed)) return(data.frame())
  parsed$CPU <- suppressWarnings(as.numeric(parsed$CPU))
  parsed$WS <- suppressWarnings(as.numeric(parsed$WS))
  parsed
}

sample_r_activity <- function(cpu_interval, active_cpu_delta) {
  p1 <- get_r_process_table()
  if (!nrow(p1)) {
    return(list(running = FALSE, active = FALSE, table = p1, max_delta = 0, note = "未发现 Rscript/Rterm 进程"))
  }
  Sys.sleep(cpu_interval)
  p2 <- get_r_process_table()
  if (!nrow(p2)) {
    return(list(running = FALSE, active = FALSE, table = p1, max_delta = 0, note = "第一次采样有 R 进程，第二次采样已消失"))
  }
  merged <- merge(
    p1[, c("Id", "ProcessName", "CPU", "WS"), drop = FALSE],
    p2[, c("Id", "CPU", "WS"), drop = FALSE],
    by = "Id",
    suffixes = c("_start", "_end"),
    all.y = TRUE
  )
  merged$CPU_delta <- merged$CPU_end - merged$CPU_start
  merged$CPU_delta[is.na(merged$CPU_delta)] <- 0
  merged$MemoryMB <- round(merged$WS_end / 1024 / 1024, 1)
  max_delta <- if (nrow(merged)) max(merged$CPU_delta, na.rm = TRUE) else 0
  list(
    running = nrow(merged) > 0,
    active = is.finite(max_delta) && max_delta >= active_cpu_delta,
    table = merged,
    max_delta = ifelse(is.finite(max_delta), max_delta, 0),
    note = ""
  )
}

make_conclusion <- function(status_state, stage, running, active, log_age_min, status_age_min, read_error, stale_minutes, completed_run) {
  if (isTRUE(completed_run$found) && !identical(stage, "FINISHED") && identical(status_state, "RUNNING") && !running) {
    return("结论：已有完成结果；当前 RUNNING 是后一次未完成/中断运行留下的状态残留。")
  }
  if (isTRUE(completed_run$found) && !identical(stage, "FINISHED") && !active && !is.na(log_age_min) && log_age_min > stale_minutes) {
    return("结论：已有完成结果；最近一次运行疑似中断或卡住，但可使用已完成的 MiST 结果。")
  }
  if (identical(status_state, "FINISHED") || identical(stage, "FINISHED")) {
    return("结论：已完成，不是卡住。")
  }
  if (identical(status_state, "FAILED") || identical(stage, "ERROR")) {
    return("结论：已失败，不是卡住；请看日志尾部的 Error。")
  }
  if (running && active) {
    return("结论：正在正常运行，没有卡住。CPU 仍在增长。")
  }
  if (running && !active && !is.na(log_age_min) && log_age_min <= stale_minutes) {
    return("结论：正在运行，暂时不像卡住。日志刚更新过，但本次 CPU 采样变化不明显。")
  }
  if (running && !active && !is.na(log_age_min) && log_age_min > stale_minutes) {
    return("结论：疑似卡住。Rscript 进程还在，但 CPU 没有明显增长，日志也超过阈值未更新。")
  }
  if (!running && identical(status_state, "RUNNING")) {
    return("结论：不是正常运行。状态停在 RUNNING，但没有检测到 Rscript/Rterm，疑似已异常退出或状态残留。")
  }
  if (!running) {
    return("结论：当前没有运行。未检测到 Rscript/Rterm。")
  }
  if (!is.null(read_error)) {
    return("结论：无法可靠判断。日志暂时不可读，请稍后重查。")
  }
  "结论：无法可靠判断。"
}

print_status_once <- function(project_dir, log_file, status_file, stale_minutes, tail_n, cpu_interval, active_cpu_delta) {
  now <- Sys.time()
  result_dir <- file.path(project_dir, "6.MiST评分结果")
  completed_run <- find_completed_run(result_dir)

  status_df <- NULL
  status_has_log <- FALSE
  if (file.exists(status_file)) {
    status_df <- tryCatch(
      read.delim(status_file, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (!is.null(status_df) && nrow(status_df) && "log_file" %in% names(status_df) && nzchar(status_df$log_file[1])) {
      log_file <- status_df$log_file[1]
      status_has_log <- TRUE
    }
  }
  if (!status_has_log) log_file <- latest_existing_log(log_file, result_dir)

  cat("\n================ MiST 状态 ================\n")
  cat("检查时间: ", format(now, "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  cat("项目目录: ", project_dir, "\n", sep = "")
  cat("日志文件: ", log_file, "\n", sep = "")
  if (isTRUE(completed_run$found)) {
    cat("已完成日志: ", completed_run$log_file, "\n", sep = "")
    cat("已完成结果: ", completed_run$result_file, "\n", sep = "")
    if (!is.na(completed_run$result_rows)) {
      cat("已完成结果行数: ", completed_run$result_rows, "\n", sep = "")
    }
  } else {
    cat("已完成结果: 未发现 SCORING FINISHED 日志\n")
  }

  status_state <- ""
  status_age_min <- NA_real_
  if (!is.null(status_df) && nrow(status_df)) {
    status_state <- status_df$status[1]
    status_time <- suppressWarnings(as.POSIXct(status_df$update_time[1], format = "%Y-%m-%d %H:%M:%S"))
    if (!is.na(status_time)) status_age_min <- as.numeric(difftime(now, status_time, units = "mins"))
    cat("状态文件: ", status_state, "\n", sep = "")
    cat("状态更新时间: ", status_df$update_time[1], "\n", sep = "")
    if (!is.na(status_age_min)) cat("状态距今未更新: ", sprintf("%.1f", status_age_min), " 分钟\n", sep = "")
    if (nzchar(status_df$message[1])) cat("状态说明: ", status_df$message[1], "\n", sep = "")
  } else {
    if (file.exists(status_file)) {
      cat("状态文件: 存在但读取失败\n")
    } else {
      cat("状态文件: 未生成\n")
    }
  }

  activity <- sample_r_activity(cpu_interval, active_cpu_delta)
  running <- activity$running
  active <- activity$active

  stage <- "未知/尚未写日志"
  read_error <- NULL
  age_min <- NA_real_
  if (file.exists(log_file)) {
    info <- file.info(log_file)
    age_min <- as.numeric(difftime(now, info$mtime, units = "mins"))
    lines <- tail_lines(log_file, max(tail_n, 200))
    read_error <- attr(lines, "read_error")
    stage <- detect_stage(lines)
    cat("日志大小: ", info$size, " bytes\n", sep = "")
    cat("日志最后更新: ", format(info$mtime, "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
    cat("距今未更新: ", sprintf("%.1f", age_min), " 分钟\n", sep = "")
    if (!is.null(read_error)) {
      cat("日志读取: 暂时失败，原因: ", read_error, "\n", sep = "")
      cat("推断阶段: 日志被占用，暂不可读\n")
    } else {
      cat("推断阶段: ", stage, "\n", sep = "")
    }

    if (is.null(read_error)) {
      cat("\n最近日志:\n")
      cat(paste(utils::tail(lines, tail_n), collapse = "\n"), "\n", sep = "")
    }
  } else {
    cat("日志文件: 不存在\n")
  }

  cat("\nCPU 采样: 间隔 ", cpu_interval, " 秒；活跃阈值 CPU 增量 >= ", active_cpu_delta, " 秒\n", sep = "")
  cat("Rscript/Rterm 进程: ", if (running) "检测到" else "未检测到", "\n", sep = "")
  if (running && nrow(activity$table)) {
    activity$table <- activity$table[order(activity$table$CPU_delta, decreasing = TRUE), , drop = FALSE]
    cat("最大 CPU 增量: ", sprintf("%.2f", activity$max_delta), " 秒\n", sep = "")
    cat("进程信息:\n")
    for (i in seq_len(nrow(activity$table))) {
      cat(
        "PID=", activity$table$Id[i],
        "\t", activity$table$ProcessName[i],
        "\tCPU增量=", sprintf("%.2f", activity$table$CPU_delta[i]),
        "\t内存MB=", activity$table$MemoryMB[i],
        "\n",
        sep = ""
      )
    }
  }

  cat("\n", make_conclusion(status_state, stage, running, active, age_min, status_age_min, read_error, stale_minutes, completed_run), "\n", sep = "")
}

script_dir <- get_script_dir()
project_dir <- if (nzchar(script_dir)) {
  normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

log_file <- get_arg("log", file.path(project_dir, "6.MiST评分结果", "mist_run.log"))
status_file <- get_arg("status", file.path(project_dir, "6.MiST评分结果", "mist_status.tsv"))
stale_minutes <- as.numeric(get_arg("stale_minutes", "10"))
tail_n <- as.integer(get_arg("tail", "30"))
cpu_interval <- as.numeric(get_arg("cpu_interval", "5"))
active_cpu_delta <- as.numeric(get_arg("active_cpu_delta", "1"))
watch <- as_bool(get_arg("watch", "false"))
interval <- as.numeric(get_arg("interval", "30"))

repeat {
  print_status_once(project_dir, log_file, status_file, stale_minutes, tail_n, cpu_interval, active_cpu_delta)
  if (!watch) break
  Sys.sleep(interval)
}
