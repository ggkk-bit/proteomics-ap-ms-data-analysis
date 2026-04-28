#!/usr/bin/env Rscript

packages <- c("yaml", "getopt", "optparse", "reshape2", "pheatmap", "RColorBrewer", "ggplot2", "gridExtra", "MESS")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

if (!length(missing)) {
  message("MiST 依赖包已全部安装")
  quit(save = "no", status = 0)
}

message("准备安装缺失依赖: ", paste(missing, collapse = ", "))
install.packages(missing, repos = "https://cloud.r-project.org")

still_missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) {
  stop("仍有依赖未安装成功: ", paste(still_missing, collapse = ", "))
}

message("MiST 依赖包安装完成")