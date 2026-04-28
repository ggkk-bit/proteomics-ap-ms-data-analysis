# 02 将 6 个算法统一为 bait/prey/score/algorithm 格式
source("07_复合分析与整合/脚本/common_integration.R", encoding = "UTF-8")
load_pkgs(c("dplyr"))
cfg <- parse_args(list(in_dir = "07_复合分析与整合/1.算法结果汇总/原始结果",
                       out_dir = "07_复合分析与整合/1.算法结果汇总/标准化结果"))
ensure_dir(cfg$out_dir)

read_raw <- function(name) read_table_auto(file.path(cfg$in_dir, paste0(name, "_raw.tsv")))
all_std <- dplyr::bind_rows(
  standardize_pairs(read_raw("CompPASS"), "CompPASS", "Bait", "Prey", "scoreWD", c("scoreZ", "AvePSM")),
  standardize_pairs(read_raw("MiST"), "MiST", "Bait", "Prey", "MIST", c("Abundance", "Reproducibility", "Specificity")),
  standardize_pairs(read_raw("HGSCore"), "HGSCore", "Protein_A", "Protein_B", "HG"),
  standardize_pairs(read_raw("CS_Score"), "CS_Score", "protein_a", "protein_b", "cs_score_z", "observed_cooccurrence"),
  standardize_pairs(read_raw("CRAPome"), "CRAPome", "Bait Name", "Prey Name", "FC_A", c("Control_Frequency", "FC_B")),
  standardize_pairs(read_raw("PPIrank"), "PPIrank", "bait_name", "prey_id", "ppirank_score", "num_replicates")
)
write_utf8_tsv(all_std, file.path(cfg$out_dir, "all_algorithms_standardized.tsv"))
write_utf8_tsv(dplyr::count(all_std, algorithm, name = "n_pairs"), file.path(cfg$out_dir, "standardize_summary.tsv"))
cat("完成标准化：", nrow(all_std), "行\n")
