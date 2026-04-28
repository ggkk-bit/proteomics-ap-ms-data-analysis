# 09 STRING 验证率可视化：按算法和按一致性分层
source("07_复合分析与整合/脚本/common_integration.R", encoding = "UTF-8")
load_pkgs(c("ggplot2"))
cfg <- parse_args(list(rate_algorithm = "07_复合分析与整合/5.STRING验证对照/验证率统计/validation_rate_by_algorithm.tsv",
                       rate_support = "07_复合分析与整合/5.STRING验证对照/验证率统计/validation_rate_by_support_count.tsv",
                       out_dir = "07_复合分析与整合/7.可视化/STRING验证图"))
ensure_dir(cfg$out_dir)
ra <- read_table_auto(cfg$rate_algorithm)
p1 <- ggplot(ra, aes(x = algorithm, y = validation_rate)) + geom_col(fill = "#59A14F") +
  geom_text(aes(label = sprintf("%.1f%%", validation_rate * 100)), vjust = -0.3, size = 3) +
  labs(x = "算法", y = "STRING 验证率", title = "各算法 STRING 验证率") + theme_bw()
ggsave(file.path(cfg$out_dir, "validation_rate_by_algorithm.png"), p1, width = 8, height = 5, dpi = 300)
rs <- read_table_auto(cfg$rate_support)
p2 <- ggplot(rs, aes(x = support_n, y = validation_rate)) + geom_point(size = 3, color = "#E15759") +
  geom_line(color = "#E15759") + labs(x = "支持算法数量", y = "STRING 验证率", title = "一致性分层 STRING 验证率") + theme_bw()
ggsave(file.path(cfg$out_dir, "validation_rate_by_support_count.png"), p2, width = 7, height = 5, dpi = 300)
cat("完成 STRING 验证图\n")
