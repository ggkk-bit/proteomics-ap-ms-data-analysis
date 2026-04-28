#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

get_arg <- function(args, name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

ensure_dir <- function(path) {
  dir_path <- if (grepl("\\.[A-Za-z0-9]+$", basename(path))) dirname(path) else path
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
}

write_tsv <- function(x, path) {
  ensure_dir(path)
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

script_dir <- function() {
  cmd <- commandArgs(FALSE)
  file_arg <- cmd[startsWith(cmd, "--file=")]
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)))
  getwd()
}

detect_quant_prefix <- function(dat, requested = "auto") {
  if (!identical(requested, "auto")) return(requested)
  prefixes <- c("LFQ intensity ", "Intensity ", "iBAQ ", "MS/MS count ")
  for (prefix in prefixes) if (any(startsWith(names(dat), prefix))) return(prefix)
  stop("未找到可识别的定量列前缀: LFQ intensity / Intensity / iBAQ / MS/MS count")
}

detect_protein_id_col <- function(dat) {
  candidates <- c("Majority protein IDs", "Protein IDs", "Protein.IDs", "protein_id", "prey_id")
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) == 0) stop("未找到蛋白 ID 列")
  hit[[1]]
}

normalize_protein_id <- function(x) {
  x <- as.character(x)
  x <- sub(";.*$", "", x)
  trimws(x)
}

message_header <- function(title) {
  cat(title, "\n", paste(rep("=", nchar(title)), collapse = ""), "\n\n", sep = "")
}
