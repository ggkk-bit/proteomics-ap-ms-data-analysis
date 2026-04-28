#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

script_path <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(), value = TRUE)[1]), winslash = '/', mustWork = FALSE)
script_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(script_dir, '..'), winslash = '/', mustWork = FALSE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0('--', name, '=')
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, '', hit[[1]], fixed = TRUE)
}

resolve_path <- function(path_text) {
  if (is.null(path_text) || !nzchar(path_text)) return(path_text)
  if (grepl('^[A-Za-z]:[/\\]|^/', path_text)) {
    return(normalizePath(path_text, winslash = '/', mustWork = FALSE))
  }
  normalizePath(file.path(project_dir, path_text), winslash = '/', mustWork = FALSE)
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop('缺少 R 包: ', pkg, '。请先安装后再运行。')
  }
}

input_file <- resolve_path(get_arg('input', '4.CompPASS输入文件/最终输入/comppass_input.tsv'))
bait_map_file <- resolve_path(get_arg('bait_map', '4.CompPASS输入文件/最终输入/bait_run_map.tsv'))
prey_ref_file <- resolve_path(get_arg('prey_ref', '4.CompPASS输入文件/最终输入/prey_reference.tsv'))
zip_file <- resolve_path(get_arg('zip', '8.文档与参考/03_源码或软件/cRomppass_master.zip'))
output_dir <- resolve_path(get_arg('output_dir', '6.CompPASS评分结果/原始结果'))
norm_factor <- as.numeric(get_arg('norm_factor', '0.98'))

if (!file.exists(input_file)) stop('找不到输入文件: ', input_file)
if (!file.exists(bait_map_file)) stop('找不到 bait_run_map.tsv: ', bait_map_file)
if (!file.exists(zip_file)) stop('找不到 cRomppass_master.zip: ', zip_file)

require_pkg('dplyr')
require_pkg('magrittr')

suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
})

input_df <- read.delim(input_file, sep = '\t', check.names = FALSE, stringsAsFactors = FALSE)
bait_map <- read.delim(bait_map_file, sep = '\t', check.names = FALSE, stringsAsFactors = FALSE)
prey_ref <- if (file.exists(prey_ref_file)) read.delim(prey_ref_file, sep = '\t', check.names = FALSE, stringsAsFactors = FALSE) else data.frame()

required_input <- c('idRun', 'idBait', 'idPrey', 'countPrey')
missing_input <- setdiff(required_input, names(input_df))
if (length(missing_input) > 0) stop('comppass_input.tsv 缺少字段: ', paste(missing_input, collapse = ', '))

replicate_df <- bait_map[, intersect(c('idRun', 'replicate_letter', 'sample_id'), names(bait_map)), drop = FALSE]
if (!'idRun' %in% names(replicate_df)) replicate_df$idRun <- character()
if (!'replicate_letter' %in% names(replicate_df)) replicate_df$replicate_letter <- ''
if (!'sample_id' %in% names(replicate_df)) replicate_df$sample_id <- ''
replicate_df <- replicate_df[!duplicated(replicate_df$idRun), , drop = FALSE]

score_input <- merge(input_df, replicate_df, by = 'idRun', all.x = TRUE, sort = FALSE)
score_input$replicate_letter <- trimws(as.character(score_input$replicate_letter))
score_input$sample_id <- trimws(as.character(score_input$sample_id))
score_input$Replicate <- ifelse(nzchar(score_input$replicate_letter), score_input$replicate_letter, score_input$sample_id)
score_input$Replicate[!nzchar(score_input$Replicate)] <- 'R1'
score_input$Experiment.ID <- score_input$idRun
score_input$Experiment.Type <- 'APMS'
score_input$Bait <- score_input$idBait
score_input$Prey <- score_input$idPrey
score_input$Spectral.Count <- suppressWarnings(as.numeric(score_input$countPrey))

if (any(is.na(score_input$Spectral.Count))) stop('countPrey 中存在无法转换为数值的记录。')
if (any(score_input$Spectral.Count <= 0)) stop('countPrey 中存在小于等于 0 的记录。')

official_input <- score_input[, c('Experiment.ID', 'Replicate', 'Experiment.Type', 'Bait', 'Prey', 'Spectral.Count')]

work_dir <- file.path(tempdir(), 'cRomppass_official_run')
if (dir.exists(work_dir)) unlink(work_dir, recursive = TRUE, force = TRUE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
unzip(zip_file, exdir = work_dir)
package_root <- file.path(work_dir, 'cRomppass-master')
source(file.path(package_root, 'R', 'cRomppass.R'), local = .GlobalEnv, encoding = 'UTF-8')
source(file.path(package_root, 'R', 'comppass.R'), local = .GlobalEnv, encoding = 'UTF-8')

results <- comppass(official_input, norm.factor = norm_factor)
results$idRun <- results$Experiment.ID
results$idBait <- results$Bait
results$idPrey <- results$Prey
results$scoreZ <- results$Z
results$scoreWD <- results$WD

if (nrow(prey_ref) > 0 && 'idPrey' %in% names(prey_ref)) {
  prey_keep <- prey_ref[, intersect(c('idPrey', 'preyGene'), names(prey_ref)), drop = FALSE]
  prey_keep <- prey_keep[!duplicated(prey_keep$idPrey), , drop = FALSE]
  results <- merge(results, prey_keep, by = 'idPrey', all.x = TRUE, sort = FALSE)
}

out_dir <- output_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(official_input, file.path(out_dir, 'comppass_official_input.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(results, file.path(out_dir, 'comppass_results.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)

cat('已输出官方格式输入: ', file.path(out_dir, 'comppass_official_input.tsv'), '\n', sep = '')
cat('已输出 CompPASS 结果: ', file.path(out_dir, 'comppass_results.tsv'), '\n', sep = '')
cat('提示: 05_可视化_CompPASS.py 现在可以直接读取 comppass_results.tsv。\n')