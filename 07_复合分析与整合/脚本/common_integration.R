# Phase 4 公共函数库：复合分析、STRING 验证和导出工具

parse_args <- function(defaults = list()) {
  args <- commandArgs(trailingOnly = TRUE)
  cfg <- defaults
  if (length(args) == 0) return(cfg)
  for (arg in args) {
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(kv) == 2) cfg[[kv[1]]] <- kv[2]
  }
  cfg
}

load_pkgs <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("缺少 R 包：%s；请先 install.packages('%s')", p, p), call. = FALSE)
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_utf8_tsv <- function(x, path) {
  ensure_dir(dirname(path))
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
}

read_table_auto <- function(path) {
  if (!file.exists(path)) {
    warning(sprintf("输入文件不存在，跳过：%s", path))
    return(data.frame())
  }
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8")
}

safe_col <- function(df, name, default = NA) {
  if (name %in% names(df)) df[[name]] else rep(default, nrow(df))
}

normalize_score <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(rep(0, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || diff(rng) == 0) return(ifelse(is.na(x), 0, 1))
  ifelse(is.na(x), 0, (x - rng[1]) / diff(rng))
}

make_pair_key <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  paste(pmin(a, b), pmax(a, b), sep = "||")
}

standardize_pairs <- function(df, algorithm, bait_col, prey_col, score_col, extra_cols = character()) {
  if (nrow(df) == 0) return(data.frame())
  out <- data.frame(
    bait = as.character(df[[bait_col]]),
    prey = as.character(df[[prey_col]]),
    raw_score = suppressWarnings(as.numeric(df[[score_col]])),
    algorithm = algorithm,
    stringsAsFactors = FALSE
  )
  out$norm_score <- normalize_score(out$raw_score)
  out$pair_key <- make_pair_key(out$bait, out$prey)
  for (cc in extra_cols) out[[cc]] <- safe_col(df, cc)
  out[!is.na(out$bait) & !is.na(out$prey) & out$bait != "" & out$prey != "", , drop = FALSE]
}

read_standardized_all <- function(path) {
  df <- read_table_auto(path)
  need <- c("bait", "prey", "raw_score", "norm_score", "algorithm", "pair_key")
  miss <- setdiff(need, names(df))
  if (length(miss) > 0) stop("标准化结果缺少列：", paste(miss, collapse = ", "))
  df
}

read_string_links <- function(path, min_score = 400) {
  df <- read_table_auto(path)
  if (nrow(df) == 0) return(data.frame())
  if ("combined_score" %in% names(df)) {
    out <- data.frame(string_a = df$protein1, string_b = df$protein2,
                      string_score = as.numeric(df$combined_score), stringsAsFactors = FALSE)
  } else if ("score" %in% names(df)) {
    out <- data.frame(string_a = df$stringId_A, string_b = df$stringId_B,
                      string_score = as.numeric(df$score), stringsAsFactors = FALSE)
  } else {
    stop("STRING 文件需包含 protein1/protein2/combined_score 或 stringId_A/stringId_B/score")
  }
  out <- out[out$string_score >= as.numeric(min_score), , drop = FALSE]
  out$string_pair_key <- make_pair_key(out$string_a, out$string_b)
  out
}

download_string_api <- function(identifiers, species = 9606, min_score = 400) {
  load_pkgs(c("httr", "jsonlite"))
  ids <- paste(unique(identifiers), collapse = "%0d")
  url <- "https://string-db.org/api/json/network"
  res <- httr::POST(url, body = list(identifiers = ids, species = species, required_score = min_score))
  httr::stop_for_status(res)
  jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))
}

summarise_algorithm_overlap <- function(std_df) {
  load_pkgs(c("dplyr", "tidyr"))
  std_df |>
    dplyr::distinct(pair_key, algorithm, bait, prey, norm_score) |>
    dplyr::group_by(pair_key) |>
    dplyr::summarise(
      bait = dplyr::first(bait),
      prey = dplyr::first(prey),
      support_n = dplyr::n_distinct(algorithm),
      algorithms = paste(sort(unique(algorithm)), collapse = ";"),
      mean_norm_score = mean(norm_score, na.rm = TRUE),
      .groups = "drop"
    )
}

calculate_validation <- function(std_df, mapping_df, string_df) {
  load_pkgs(c("dplyr"))
  map <- mapping_df |>
    dplyr::filter(!is.na(string_id), string_id != "") |>
    dplyr::distinct(protein_id, string_id)
  mapped <- std_df |>
    dplyr::left_join(map, by = c("bait" = "protein_id")) |>
    dplyr::rename(string_bait = string_id) |>
    dplyr::left_join(map, by = c("prey" = "protein_id")) |>
    dplyr::rename(string_prey = string_id)
  mapped$string_pair_key <- make_pair_key(mapped$string_bait, mapped$string_prey)
  validated <- mapped |>
    dplyr::left_join(string_df[, c("string_pair_key", "string_score")], by = "string_pair_key") |>
    dplyr::mutate(
      mapped_to_string = !is.na(string_bait) & !is.na(string_prey),
      string_validated = mapped_to_string & !is.na(string_score)
    )
  validated
}

validation_rates <- function(validated_df) {
  load_pkgs(c("dplyr"))
  by_algorithm <- validated_df |>
    dplyr::group_by(algorithm) |>
    dplyr::summarise(total_pairs = dplyr::n(), mapped_pairs = sum(mapped_to_string),
                     validated_pairs = sum(string_validated),
                     validation_rate = validated_pairs / pmax(mapped_pairs, 1), .groups = "drop")
  support <- summarise_algorithm_overlap(validated_df)
  pair_valid <- validated_df |>
    dplyr::group_by(pair_key) |>
    dplyr::summarise(string_validated = any(string_validated), mapped_to_string = any(mapped_to_string), .groups = "drop")
  by_support <- support |>
    dplyr::left_join(pair_valid, by = "pair_key") |>
    dplyr::group_by(support_n) |>
    dplyr::summarise(total_pairs = dplyr::n(), mapped_pairs = sum(mapped_to_string),
                     validated_pairs = sum(string_validated),
                     validation_rate = validated_pairs / pmax(mapped_pairs, 1), .groups = "drop")
  list(by_algorithm = by_algorithm, by_support = by_support)
}

integrated_score <- function(overlap_df, validated_pair_df, w1 = 0.4, w2 = 0.3, w3 = 0.3) {
  load_pkgs(c("dplyr"))
  max_support <- max(overlap_df$support_n, na.rm = TRUE)
  overlap_df |>
    dplyr::left_join(validated_pair_df, by = "pair_key") |>
    dplyr::mutate(
      string_bonus = ifelse(!is.na(string_validated) & string_validated == TRUE, 1, 0),
      integrated_score = as.numeric(w1) * support_n / max_support +
        as.numeric(w2) * mean_norm_score +
        as.numeric(w3) * string_bonus
    ) |>
    dplyr::arrange(dplyr::desc(integrated_score), dplyr::desc(support_n))
}
