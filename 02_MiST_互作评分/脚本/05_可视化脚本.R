#!/usr/bin/env Rscript

# MiST 可视化脚本
# 作用：
# 1. 读取 MiST 评分结果或任意包含 bait/prey/score 的结果表
# 2. 自动识别评分列、bait 列、prey 列及可选 abundance / reproducibility 列
# 3. 输出通用图形和 Cytoscape 边表

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0) return(default)
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

pick_col <- function(candidates, names_vec) {
  hit <- intersect(candidates, names_vec)
  if (length(hit)) hit[[1]] else NA_character_
}

input_file <- get_arg("input", "..\\6.MiST评分结果\\MiST原始评分.tsv")
out_dir <- get_arg("out_dir", "..\\7.可视化")
threshold <- as.numeric(get_arg("threshold", "0.75"))
top_n <- as.integer(get_arg("top_n", "20"))

if (!file.exists(input_file)) stop("找不到 MiST 评分文件: ", input_file)

dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")

score_col <- pick_col(c("MiST", "mist", "score", "Score", "final_score"), names(dat))
abundance_col <- pick_col(c("Abundance", "abundance", "avg_abundance", "mean_intensity"), names(dat))
repro_col <- pick_col(c("Reproducibility", "reproducibility", "repro", "replicability"), names(dat))
bait_col <- pick_col(c("Bait", "bait", "bait_name", "BaitName"), names(dat))
prey_col <- pick_col(c("Prey", "prey", "prey_id", "prey_name", "PreyName"), names(dat))

if (is.na(score_col)) stop("评分表里找不到 MiST 分数字段，请至少包含 MiST/mist/score 之一")
if (is.na(bait_col) || is.na(prey_col)) stop("评分表里找不到 bait/prey 字段，请至少包含 bait_name 或 bait，以及 prey_id 或 prey")

dat[[score_col]] <- suppressWarnings(as.numeric(dat[[score_col]]))
dat <- dat[!is.na(dat[[score_col]]), , drop = FALSE]
if (!nrow(dat)) stop("评分表中没有可用的数值型评分")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png(file.path(out_dir, "01_MiST分数分布.png"), width = 1600, height = 1200, res = 200)
hist(dat[[score_col]], breaks = 30, col = "#4C78A8", border = "white",
     main = "MiST 分数分布", xlab = score_col, ylab = "Count")
abline(v = threshold, col = "#E45756", lwd = 3, lty = 2)
dev.off()

if (!is.na(abundance_col) && !is.na(repro_col)) {
  dat[[abundance_col]] <- suppressWarnings(as.numeric(dat[[abundance_col]]))
  dat[[repro_col]] <- suppressWarnings(as.numeric(dat[[repro_col]]))
  use <- !is.na(dat[[abundance_col]]) & !is.na(dat[[repro_col]])
  if (any(use)) {
    png(file.path(out_dir, "02_丰度_重复性散点图.png"), width = 1600, height = 1200, res = 200)
    plot(dat[[abundance_col]][use], dat[[repro_col]][use], pch = 16,
         col = ifelse(dat[[score_col]][use] >= threshold, "#54A24B", "#9D9D9D"),
         xlab = abundance_col, ylab = repro_col, main = "Abundance vs Reproducibility")
    abline(v = median(dat[[abundance_col]][use], na.rm = TRUE), lty = 3, col = "gray60")
    abline(h = median(dat[[repro_col]][use], na.rm = TRUE), lty = 3, col = "gray60")
    dev.off()
  }
}

top <- dat[order(-dat[[score_col]]), , drop = FALSE]
top <- head(top, top_n)
png(file.path(out_dir, "03_Top互作条形图.png"), width = 1800, height = 1200, res = 200)
par(mar = c(10, 5, 4, 2))
labs <- paste(top[[bait_col]], top[[prey_col]], sep = " -> ")
barplot(rev(top[[score_col]]), names.arg = rev(labs), las = 2, horiz = FALSE,
        col = "#72B7B2", main = paste0("Top ", nrow(top), " MiST 互作"), ylab = score_col)
abline(h = threshold, col = "#E45756", lty = 2)
dev.off()

edge <- data.frame(
  bait = dat[[bait_col]],
  prey = dat[[prey_col]],
  score = dat[[score_col]],
  pass_threshold = dat[[score_col]] >= threshold,
  stringsAsFactors = FALSE
)
write.table(edge, file.path(out_dir, "MiST_Cytoscape边表.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")

summary_df <- data.frame(
  score_col = score_col,
  bait_col = bait_col,
  prey_col = prey_col,
  abundance_col = ifelse(is.na(abundance_col), "", abundance_col),
  reproducibility_col = ifelse(is.na(repro_col), "", repro_col),
  threshold = threshold,
  total_rows = nrow(dat),
  passing_rows = sum(dat[[score_col]] >= threshold, na.rm = TRUE),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)
write.table(summary_df, file.path(out_dir, "可视化运行摘要.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")
