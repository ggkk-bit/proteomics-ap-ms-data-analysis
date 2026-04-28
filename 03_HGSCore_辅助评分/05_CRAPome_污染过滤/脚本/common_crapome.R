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

message_header <- function(title) {
  cat(title, "\n", paste(rep("=", nchar(title)), collapse = ""), "\n\n", sep = "")
}

detect_protein_id_col <- function(dat) {
  candidates <- c("Majority protein IDs", "Protein IDs", "Protein.IDs", "Gene names", "Protein names", "Prey_ID", "prey_id")
  hit <- candidates[candidates %in% names(dat)]
  if (length(hit) == 0) stop("未找到蛋白标识列。可用列包括 Majority protein IDs / Protein IDs / Gene names。")
  hit[[1]]
}

normalize_protein_id <- function(x) {
  x <- as.character(x)
  x <- sub(";.*$", "", x)
  trimws(x)
}

detect_quant_prefix <- function(dat, requested = "auto") {
  if (!identical(requested, "auto")) return(requested)
  prefixes <- c("MS/MS count ", "Spectral count ", "Spectral Count ", "Intensity ", "LFQ intensity ", "iBAQ ")
  for (prefix in prefixes) if (any(startsWith(names(dat), prefix))) return(prefix)
  stop("未找到可识别的定量列前缀：MS/MS count / Spectral count / Intensity / LFQ intensity / iBAQ。")
}

safe_div <- function(num, den, pseudo = 1e-6) {
  (num + pseudo) / (den + pseudo)
}
