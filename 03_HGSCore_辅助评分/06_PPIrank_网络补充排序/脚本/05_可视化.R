#!/usr/bin/env Rscript
# PPIrank 可视化：输出网络图、统计图、Cytoscape 边表/节点表和 SIF 文件；bait 节点标绿。

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(igraph) })
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) { prefix <- paste0("--", name, "="); hit <- args[startsWith(args, prefix)]; if (length(hit) == 0) return(default); sub(prefix, "", hit[[1]], fixed = TRUE) }
ensure_dir <- function(path) { dir_path <- dirname(path); if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE) }
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
write_tsv <- function(x, path) { ensure_dir(path); write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "") }

score_file <- get_arg("scores", "./6.PPIrank排序结果/ppirank_scores.tsv")
edge_file <- get_arg("edges", "./6.PPIrank排序结果/ppirank_edges.tsv")
node_file <- get_arg("nodes", "./6.PPIrank排序结果/ppirank_nodes.tsv")
out_dir <- get_arg("out_dir", "./7.可视化")
top_n <- as.integer(get_arg("top_n", "100"))
min_score <- as.numeric(get_arg("min_score", "1"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scores <- read_tsv(score_file)
edges <- read_tsv(edge_file) %>% filter(score >= min_score) %>% arrange(desc(score)) %>% head(top_n)
nodes <- read_tsv(node_file) %>% filter(protein_id %in% unique(c(edges$bait, edges$prey)))
if (nrow(edges) == 0) stop("没有满足可视化阈值的边，请降低 --min_score=。")

# Cytoscape 标准导入表与 SIF。
cy_edges <- edges %>% transmute(source = bait, target = prey, interaction = "ppirank", score = score, replicates = replicates, confidence_class = confidence_class)
cy_nodes <- nodes %>% transmute(id = protein_id, type = type, color = color, size = size, total_spectral_count = total_spectral_count, max_score = max_score)
write_tsv(cy_edges, file.path(out_dir, "cytoscape_edges.tsv")); write_tsv(cy_nodes, file.path(out_dir, "cytoscape_nodes.tsv"))
write.table(cy_edges[, c("source", "interaction", "target")], file.path(out_dir, "ppirank_network.sif"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# igraph 网络图：bait 绿色且较大，prey 灰色；边宽反映 PPIrank 分数。
g <- graph_from_data_frame(d = edges %>% transmute(from = bait, to = prey, score = score, replicates = replicates, confidence_class = confidence_class), vertices = nodes %>% transmute(name = protein_id, type = type, color = color, size = size, total_spectral_count = total_spectral_count), directed = FALSE)
V(g)$color <- ifelse(V(g)$type == "bait", "#2ca25f", "#bdbdbd")
V(g)$size <- ifelse(V(g)$type == "bait", 9, 5) + log1p(V(g)$total_spectral_count)
V(g)$label.cex <- 0.55; V(g)$label.color <- "#222222"
E(g)$width <- pmax(1, 1 + 4 * E(g)$score / max(E(g)$score, na.rm = TRUE))
E(g)$color <- ifelse(E(g)$confidence_class == "high", "#3182bd", "#9ecae1")
png(file.path(out_dir, "ppirank_network.png"), width = 1800, height = 1400, res = 180); set.seed(1); plot(g, layout = layout_with_fr(g), main = "PPIrank bait-prey network"); dev.off()

p_score <- ggplot(scores, aes(x = ppirank_score, fill = confidence_class)) + geom_histogram(bins = 40, color = "white") + theme_bw(base_size = 12) + labs(x = "PPIrank score", y = "Interaction count", fill = "Confidence", title = "PPIrank score distribution")
ggsave(file.path(out_dir, "ppirank_score_distribution.png"), p_score, width = 7, height = 5, dpi = 300)
p_rep <- ggplot(scores, aes(x = factor(bait_replicates_observed), fill = confidence_class)) + geom_bar() + theme_bw(base_size = 12) + labs(x = "Observed biological replicates", y = "Interaction count", fill = "Confidence", title = "Reproducibility summary")
ggsave(file.path(out_dir, "ppirank_replicate_summary.png"), p_rep, width = 7, height = 5, dpi = 300)
write_tsv(scores %>% group_by(bait_name) %>% summarise(interactions = n(), high_confidence = sum(confidence_class == "high"), .groups = "drop") %>% arrange(desc(high_confidence), desc(interactions)), file.path(out_dir, "bait_interaction_summary.tsv"))
writeLines(c("PPIrank 可视化记录", paste("时间:", Sys.time()), paste("网络边数:", nrow(edges)), paste("网络节点数:", nrow(nodes)), paste("top_n:", top_n), paste("min_score:", min_score), paste("输出目录:", out_dir), "输出文件: ppirank_network.png, ppirank_score_distribution.png, ppirank_replicate_summary.png, cytoscape_edges.tsv, cytoscape_nodes.tsv, ppirank_network.sif"), "./5.质控/05_可视化记录.txt", useBytes = TRUE)
cat("可视化输出目录: ", out_dir, "\n", sep = "")
