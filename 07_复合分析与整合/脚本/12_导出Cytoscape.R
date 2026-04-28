# 12 导出 Cytoscape .sif、边表和节点表
source("07_复合分析与整合/脚本/common_integration.R", encoding = "UTF-8")
load_pkgs(c("dplyr"))
cfg <- parse_args(list(input = "07_复合分析与整合/6.综合评分/高置信互作/high_confidence_interactions.tsv",
                       out_dir = "07_复合分析与整合/8.导出结果/Cytoscape文件"))
ensure_dir(cfg$out_dir)
df <- read_table_auto(cfg$input)
edges <- df |> dplyr::transmute(source = bait, interaction = "pp", target = prey,
                                integrated_score, support_n, algorithms, string_validated)
nodes <- data.frame(id = sort(unique(c(edges$source, edges$target))), stringsAsFactors = FALSE)
nodes$degree_in_export <- vapply(nodes$id, function(x) sum(edges$source == x | edges$target == x), numeric(1))
utils::write.table(edges[, c("source", "interaction", "target")], file.path(cfg$out_dir, "high_confidence_interactions.sif"),
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, fileEncoding = "UTF-8")
write_utf8_tsv(edges, file.path(cfg$out_dir, "cytoscape_edges.tsv"))
write_utf8_tsv(nodes, file.path(cfg$out_dir, "cytoscape_nodes.tsv"))
cat("完成 Cytoscape 导出\n")
