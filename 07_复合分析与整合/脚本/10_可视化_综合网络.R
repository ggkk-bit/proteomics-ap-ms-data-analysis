# 10 使用 igraph 绘制高置信互作网络
source("07_复合分析与整合/脚本/common_integration.R", encoding = "UTF-8")
load_pkgs(c("igraph"))
cfg <- parse_args(list(input = "07_复合分析与整合/6.综合评分/高置信互作/high_confidence_interactions.tsv",
                       top_n = "200", out_dir = "07_复合分析与整合/7.可视化/综合网络图"))
ensure_dir(cfg$out_dir)
df <- read_table_auto(cfg$input)
df <- head(df[order(-df$integrated_score), ], as.numeric(cfg$top_n))
edges <- data.frame(from = df$bait, to = df$prey, weight = df$integrated_score, stringsAsFactors = FALSE)
g <- igraph::graph_from_data_frame(edges, directed = FALSE)
png(file.path(cfg$out_dir, "integrated_network.png"), width = 1800, height = 1600, res = 180)
plot(g, vertex.size = 5, vertex.label.cex = 0.6, edge.width = pmax(1, igraph::E(g)$weight * 4),
     main = "高置信综合互作网络")
dev.off()
cat("完成综合网络图\n")
