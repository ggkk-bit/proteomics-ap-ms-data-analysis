# 08 可视化算法重叠：热图、Venn/UpSet 基础数据图
source("07_复合分析与整合/脚本/common_integration.R", encoding = "UTF-8")
load_pkgs(c("ggplot2", "pheatmap"))
cfg <- parse_args(list(matrix_file = "07_复合分析与整合/4.算法交叉验证/重叠分析/algorithm_overlap_matrix.tsv",
                       support_file = "07_复合分析与整合/4.算法交叉验证/一致性评分/pair_consistency_scores.tsv",
                       out_dir = "07_复合分析与整合/7.可视化/算法重叠图"))
ensure_dir(cfg$out_dir)
mat <- read_table_auto(cfg$matrix_file)
png(file.path(cfg$out_dir, "algorithm_overlap_heatmap.png"), width = 1600, height = 1400, res = 180)
pheatmap::pheatmap(as.matrix(mat), cluster_rows = FALSE, cluster_cols = FALSE, main = "算法互作重叠矩阵")
dev.off()
support <- read_table_auto(cfg$support_file)
p <- ggplot(support, aes(x = factor(support_n))) + geom_bar(fill = "#4E79A7") +
  labs(x = "支持算法数量", y = "蛋白对数量", title = "多算法一致性分布") + theme_bw()
ggsave(file.path(cfg$out_dir, "support_count_distribution.png"), p, width = 7, height = 5, dpi = 300)
cat("完成重叠分析可视化\n")
