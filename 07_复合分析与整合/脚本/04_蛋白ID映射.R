# 04 蛋白 ID 映射到 STRING Ensembl Protein ID
source("07_复合分析与整合/脚本/common_integration.R", encoding = "UTF-8")
load_pkgs(c("dplyr"))
cfg <- parse_args(list(species = "9606",
                       std_file = "07_复合分析与整合/1.算法结果汇总/标准化结果/all_algorithms_standardized.tsv",
                       mapping_file = "",
                       out_dir = "07_复合分析与整合/3.蛋白ID映射/映射表",
                       log_dir = "07_复合分析与整合/3.蛋白ID映射/映射日志"))
ensure_dir(cfg$out_dir); ensure_dir(cfg$log_dir)
std <- read_standardized_all(cfg$std_file)
proteins <- data.frame(protein_id = sort(unique(c(std$bait, std$prey))), stringsAsFactors = FALSE)

if (nzchar(cfg$mapping_file) && file.exists(cfg$mapping_file)) {
  map <- read_table_auto(cfg$mapping_file)
  names(map)[1:2] <- c("protein_id", "string_id")
} else {
  # 占位映射：若输入本身已是 STRING ID（9606.ENSP...），直接保留；否则留空并记录
  map <- proteins |>
    dplyr::mutate(string_id = ifelse(grepl(paste0("^", cfg$species, "\\."), protein_id), protein_id, NA_character_))
}
map <- proteins |> dplyr::left_join(map, by = "protein_id") |> dplyr::distinct(protein_id, .keep_all = TRUE)
write_utf8_tsv(map, file.path(cfg$out_dir, "protein_to_string_mapping.tsv"))
write_utf8_tsv(dplyr::summarise(map, total = dplyr::n(), mapped = sum(!is.na(string_id)), unmapped = sum(is.na(string_id))),
               file.path(cfg$log_dir, "mapping_summary.tsv"))
cat("完成 ID 映射；如 unmapped 较多，请提供 --mapping_file=蛋白到STRING映射表.tsv\n")
