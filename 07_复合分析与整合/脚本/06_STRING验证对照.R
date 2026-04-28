# 06 核心脚本：与 STRING 数据库对照，统计验证率并识别新发现互作
source("07_复合分析与整合/脚本/common_integration.R", encoding = "UTF-8")
load_pkgs(c("dplyr"))
cfg <- parse_args(list(
  min_score = "400",
  std_file = "07_复合分析与整合/1.算法结果汇总/标准化结果/all_algorithms_standardized.tsv",
  mapping_file = "07_复合分析与整合/3.蛋白ID映射/映射表/protein_to_string_mapping.tsv",
  string_file = "07_复合分析与整合/2.STRING数据库/处理后数据/string_links_filtered.tsv",
  out_dir = "07_复合分析与整合/5.STRING验证对照"
))
ensure_dir(file.path(cfg$out_dir, "验证结果")); ensure_dir(file.path(cfg$out_dir, "验证率统计")); ensure_dir(file.path(cfg$out_dir, "新发现互作"))

std <- read_standardized_all(cfg$std_file)
mapping <- read_table_auto(cfg$mapping_file)
string_links <- read_string_links(cfg$string_file, min_score = as.numeric(cfg$min_score))

# 将每条算法互作映射到 STRING ID，并用无方向 pair_key 匹配 STRING 互作
validated <- calculate_validation(std, mapping, string_links)
write_utf8_tsv(validated, file.path(cfg$out_dir, "验证结果", "algorithm_pairs_string_validation.tsv"))

# 统计 1：每个算法各自有多少互作能被 STRING 证据验证
rates <- validation_rates(validated)
write_utf8_tsv(rates$by_algorithm, file.path(cfg$out_dir, "验证率统计", "validation_rate_by_algorithm.tsv"))

# 统计 2：按 1-6 个算法支持数分层，观察一致性越高是否越容易被 STRING 验证
write_utf8_tsv(rates$by_support, file.path(cfg$out_dir, "验证率统计", "validation_rate_by_support_count.tsv"))

# 汇总到蛋白对层面：STRING 未验证但多算法支持的结果作为候选新发现互作
overlap <- summarise_algorithm_overlap(validated)
pair_valid <- validated |>
  dplyr::group_by(pair_key) |>
  dplyr::summarise(mapped_to_string = any(mapped_to_string),
                   string_validated = any(string_validated),
                   max_string_score = suppressWarnings(max(string_score, na.rm = TRUE)),
                   .groups = "drop") |>
  dplyr::mutate(max_string_score = ifelse(is.infinite(max_string_score), NA, max_string_score))
pair_summary <- overlap |> dplyr::left_join(pair_valid, by = "pair_key")
write_utf8_tsv(pair_summary, file.path(cfg$out_dir, "验证结果", "pair_level_string_validation.tsv"))

novel <- pair_summary |>
  dplyr::filter(mapped_to_string, !string_validated) |>
  dplyr::arrange(dplyr::desc(support_n), dplyr::desc(mean_norm_score))
write_utf8_tsv(novel, file.path(cfg$out_dir, "新发现互作", "candidate_novel_interactions.tsv"))
cat("STRING 验证完成：算法验证率、支持数分层验证率、新发现候选均已输出\n")
