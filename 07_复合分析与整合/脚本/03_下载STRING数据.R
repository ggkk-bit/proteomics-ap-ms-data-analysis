# 03 下载 STRING v12.0 数据：支持 API 和批量文件 URL
source("07_复合分析与整合/脚本/common_integration.R", encoding = "UTF-8")
cfg <- parse_args(list(species = "9606", min_score = "400", mode = "url",
                       out_dir = "07_复合分析与整合/2.STRING数据库/原始下载",
                       std_file = "07_复合分析与整合/1.算法结果汇总/标准化结果/all_algorithms_standardized.tsv"))
ensure_dir(cfg$out_dir)

if (cfg$mode == "api") {
  std <- read_standardized_all(cfg$std_file)
  ids <- unique(c(std$bait, std$prey))
  api_df <- download_string_api(ids, species = cfg$species, min_score = cfg$min_score)
  write_utf8_tsv(api_df, file.path(cfg$out_dir, paste0("string_api_", cfg$species, ".tsv")))
} else {
  # STRING v12.0 批量下载地址；网络受限时请手动下载后放入此目录
  url <- sprintf("https://stringdb-downloads.org/download/protein.links.v12.0/%s.protein.links.v12.0.txt.gz", cfg$species)
  dest <- file.path(cfg$out_dir, sprintf("%s.protein.links.v12.0.txt.gz", cfg$species))
  ok <- tryCatch({ utils::download.file(url, dest, mode = "wb"); TRUE }, error = function(e) FALSE)
  log <- data.frame(species = cfg$species, url = url, dest = dest, downloaded = ok, stringsAsFactors = FALSE)
  write_utf8_tsv(log, file.path(cfg$out_dir, "download_log.tsv"))
}
cat("STRING 下载步骤完成\n")
