#!/usr/bin/env Rscript
# CS_Score 核心算法。
# 论文逻辑：统计蛋白对在 purification 中的共现次数；通过保持每个 purification 成员数不变的随机打乱
# 生成背景分布；以 Z = (observed - mean_random) / sd_random 表示共现显著性。

cmd <- commandArgs(FALSE)
fa <- cmd[startsWith(cmd, "--file=")]
.script_dir <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(.script_dir, "common_cs_score.R"))
args <- commandArgs(trailingOnly = TRUE)
message_header("CS_Score 计算")

member_file <- get_arg(args, "member", "./4.CS_Score输入文件/cs_score_member_long.tsv")
matrix_file <- get_arg(args, "matrix", "./4.CS_Score输入文件/cs_score_binary_matrix.tsv")
output_file <- get_arg(args, "output", "./6.CS_Score评分结果/cs_score_all_pairs.tsv")
random_summary_file <- get_arg(args, "random_summary", "./6.CS_Score评分结果/cs_score_random_summary.tsv")
qc_file <- get_arg(args, "qc", "./5.质控/03_CS_Score计算记录.txt")
n_perm <- as.integer(get_arg(args, "n_perm", "1000"))
seed <- as.integer(get_arg(args, "seed", "1"))
min_observed <- as.integer(get_arg(args, "min_observed", "2"))

if (file.exists(matrix_file)) {
  mat_df <- read_tsv(matrix_file)
  row_ids <- mat_df[[1]]
  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat) <- row_ids
  storage.mode(mat) <- "numeric"
} else if (file.exists(member_file)) {
  member <- read_tsv(member_file)
  mat_df <- member %>%
    select(purification_id, prey_id, presence) %>%
    distinct() %>%
    pivot_wider(names_from = prey_id, values_from = presence, values_fill = 0)
  row_ids <- mat_df[[1]]
  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat) <- row_ids
  storage.mode(mat) <- "numeric"
} else {
  stop("找不到 CS_Score 输入文件")
}

mat[mat > 0] <- 1
prey_ids <- colnames(mat)
if (ncol(mat) < 2) stop("蛋白数量少于 2，无法计算蛋白对共现")

observed_mat <- crossprod(mat)
diag(observed_mat) <- 0
pair_index <- which(upper.tri(observed_mat) & observed_mat >= min_observed, arr.ind = TRUE)
if (nrow(pair_index) == 0) stop("没有观察共现次数 >= ", min_observed, " 的蛋白对")

p1 <- rownames(observed_mat)[pair_index[, 1]]
p2 <- colnames(observed_mat)[pair_index[, 2]]
observed <- observed_mat[pair_index]
pair_key <- paste(p1, p2, sep = "||")

set.seed(seed)
row_sizes <- rowSums(mat)
n_prey <- ncol(mat)
random_values <- matrix(0, nrow = length(observed), ncol = n_perm)

for (i in seq_len(n_perm)) {
  rand <- matrix(0, nrow = nrow(mat), ncol = n_prey)
  for (r in seq_len(nrow(mat))) {
    k <- row_sizes[[r]]
    if (k > 0) rand[r, sample.int(n_prey, k)] <- 1
  }
  rand_co <- crossprod(rand)
  random_values[, i] <- rand_co[pair_index]
  if (i %% 100 == 0) cat("随机打乱完成: ", i, "/", n_perm, "\n", sep = "")
}

random_mean <- rowMeans(random_values)
random_sd <- apply(random_values, 1, sd)
z_score <- ifelse(random_sd > 0, (observed - random_mean) / random_sd, NA_real_)
p_empirical <- (rowSums(t(t(random_values) >= observed)) + 1) / (n_perm + 1)

result <- data.frame(
  protein_a = p1,
  protein_b = p2,
  pair_id = pair_key,
  observed_cooccurrence = observed,
  random_mean = random_mean,
  random_sd = random_sd,
  cs_score_z = z_score,
  empirical_p = p_empirical,
  n_random = n_perm,
  stringsAsFactors = FALSE
) %>% arrange(desc(cs_score_z), desc(observed_cooccurrence))

random_summary <- result %>%
  summarise(
    pairs_scored = n(),
    n_random = first(n_random),
    mean_observed = mean(observed_cooccurrence),
    mean_random_mean = mean(random_mean),
    sd_z = sd(cs_score_z, na.rm = TRUE)
  )

write_tsv(result, output_file)
write_tsv(random_summary, random_summary_file)
writeLines(c(
  "CS_Score 计算记录",
  paste("时间:", Sys.time()),
  paste("purification 数:", nrow(mat)),
  paste("prey 数:", ncol(mat)),
  paste("保留蛋白对阈值 observed >", min_observed - 1),
  paste("评分蛋白对数:", nrow(result)),
  paste("随机次数:", n_perm),
  paste("随机种子:", seed)
), qc_file, useBytes = TRUE)

cat("结果: ", output_file, "\n", sep = "")
