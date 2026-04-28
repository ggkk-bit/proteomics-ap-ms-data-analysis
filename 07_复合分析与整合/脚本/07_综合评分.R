# 07 综合评分：w1*支持算法数量 + w2*归一化平均分数 + w3*STRING验证加分
source("07_复合分析与整合/脚本/common_integration.R", encoding = "UTF-8")
load_pkgs(c("dplyr"))
cfg <- parse_args(list(w1 = "0.4", w2 = "0.3", w3 = "0.3", high_conf = "0.7",
                       consistency_file = "07_复合分析与整合/4.算法交叉验证/一致性评分/pair_consistency_scores.tsv",
                       validation_file = "07_复合分析与整合/5.STRING验证对照/验证结果/pair_level_string_validation.tsv",
                       out_dir = "07_复合分析与整合/6.综合评分"))
ensure_dir(file.path(cfg$out_dir, "评分结果")); ensure_dir(file.path(cfg$out_dir, "高置信互作"))
overlap <- read_table_auto(cfg$consistency_file)
valid <- read_table_auto(cfg$validation_file) |> dplyr::select(pair_key, mapped_to_string, string_validated, max_string_score)
scored <- integrated_score(overlap, valid, cfg$w1, cfg$w2, cfg$w3)
write_utf8_tsv(scored, file.path(cfg$out_dir, "评分结果", "integrated_interaction_scores.tsv"))
write_utf8_tsv(dplyr::filter(scored, integrated_score >= as.numeric(cfg$high_conf)),
               file.path(cfg$out_dir, "高置信互作", "high_confidence_interactions.tsv"))
cat("完成综合评分\n")
