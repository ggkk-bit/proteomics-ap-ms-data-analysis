# 05 计算算法重叠矩阵和一致性评分
source("07_复合分析与整合/脚本/common_integration.R", encoding = "UTF-8")
load_pkgs(c("dplyr", "tidyr"))
cfg <- parse_args(list(std_file = "07_复合分析与整合/1.算法结果汇总/标准化结果/all_algorithms_standardized.tsv",
                       out_dir = "07_复合分析与整合/4.算法交叉验证"))
std <- read_standardized_all(cfg$std_file)
ensure_dir(file.path(cfg$out_dir, "重叠分析")); ensure_dir(file.path(cfg$out_dir, "一致性评分")); ensure_dir(file.path(cfg$out_dir, "Venn图数据"))

sets <- split(std$pair_key, std$algorithm)
algs <- sort(names(sets))
mat <- outer(algs, algs, Vectorize(function(a, b) length(intersect(unique(sets[[a]]), unique(sets[[b]])))))
dimnames(mat) <- list(algs, algs)
write_utf8_tsv(as.data.frame(mat), file.path(cfg$out_dir, "重叠分析", "algorithm_overlap_matrix.tsv"))

support <- summarise_algorithm_overlap(std)
write_utf8_tsv(support, file.path(cfg$out_dir, "一致性评分", "pair_consistency_scores.tsv"))
venn_data <- std |> dplyr::distinct(pair_key, algorithm) |> tidyr::pivot_wider(names_from = algorithm, values_from = algorithm, values_fn = length, values_fill = 0)
write_utf8_tsv(venn_data, file.path(cfg$out_dir, "Venn图数据", "venn_upset_binary.tsv"))
cat("完成算法交叉验证\n")
