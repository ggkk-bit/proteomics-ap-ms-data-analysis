#!/usr/bin/env Rscript
# 核心算法：计算 CRAPome FC-A / FC-B。

script_path <- commandArgs(FALSE)
file_arg <- script_path[startsWith(script_path, "--file=")]
script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "common_crapome.R"))

args <- commandArgs(trailingOnly = TRUE)
message_header("CRAPome FC 计算")

input_long <- get_arg(args, "input", "./4.CRAPome输入文件/crapome_input_long.tsv")
control_file <- get_arg(args, "control", "./4.CRAPome输入文件/control_library_example.tsv")
output_file <- get_arg(args, "output", "./6.CRAPome过滤结果/crapome_fc_all.tsv")
experiment_type <- get_arg(args, "experiment_type", "all")
pseudo_count <- as.numeric(get_arg(args, "pseudo_count", "1"))

exp_long <- read_tsv(input_long) %>% mutate(`Spectral Count` = suppressWarnings(as.numeric(`Spectral Count`))) %>% filter(!is.na(`Spectral Count`), `Spectral Count` > 0)
ctrl <- read_tsv(control_file)
if (!all(c("Prey_ID", "Control_Frequency", "Total_Controls") %in% names(ctrl))) stop("对照库至少需要列: Prey_ID, Control_Frequency, Total_Controls。")
if (!"Control_Spectral_Count" %in% names(ctrl)) ctrl$Control_Spectral_Count <- ctrl$Control_Frequency * ctrl$Total_Controls
if (!"Control_Total_Spectral_Count" %in% names(ctrl)) ctrl$Control_Total_Spectral_Count <- max(sum(ctrl$Control_Spectral_Count), 1)
if (!"Control_Type" %in% names(ctrl)) ctrl$Control_Type <- "all"

exp_total <- exp_long %>% group_by(`AP Name`) %>% summarise(Experiment_Total_Spectral_Count = sum(`Spectral Count`), .groups = "drop")
exp_summary <- exp_long %>% group_by(`Bait Name`, `AP Name`, `Prey Name`) %>% summarise(Experiment_Spectral_Count = sum(`Spectral Count`), .groups = "drop") %>% left_join(exp_total, by = "AP Name") %>% mutate(Experiment_Ratio = safe_div(Experiment_Spectral_Count, Experiment_Total_Spectral_Count, pseudo_count))

ctrl_all <- ctrl %>% group_by(Prey_ID) %>% summarise(Control_Frequency = max(as.numeric(Control_Frequency), na.rm = TRUE), Total_Controls = max(as.numeric(Total_Controls), na.rm = TRUE), Control_Spectral_Count_A = sum(as.numeric(Control_Spectral_Count), na.rm = TRUE), Control_Total_Spectral_Count_A = max(sum(as.numeric(Control_Total_Spectral_Count), na.rm = TRUE), 1), .groups = "drop") %>% mutate(Control_Ratio_A = safe_div(Control_Spectral_Count_A, Control_Total_Spectral_Count_A, pseudo_count))
ctrl_b <- ctrl %>% filter(Control_Type == experiment_type | Control_Type == "all") %>% group_by(Prey_ID) %>% summarise(Control_Spectral_Count_B = sum(as.numeric(Control_Spectral_Count), na.rm = TRUE), Control_Total_Spectral_Count_B = max(sum(as.numeric(Control_Total_Spectral_Count), na.rm = TRUE), 1), .groups = "drop") %>% mutate(Control_Ratio_B = safe_div(Control_Spectral_Count_B, Control_Total_Spectral_Count_B, pseudo_count))

result <- exp_summary %>% left_join(ctrl_all, by = c("Prey Name" = "Prey_ID")) %>% left_join(ctrl_b, by = c("Prey Name" = "Prey_ID")) %>% mutate(Control_Frequency = ifelse(is.na(Control_Frequency), 0, Control_Frequency), Total_Controls = ifelse(is.na(Total_Controls), 0, Total_Controls), Control_Ratio_A = ifelse(is.na(Control_Ratio_A), safe_div(0, 1, pseudo_count), Control_Ratio_A), Control_Ratio_B = ifelse(is.na(Control_Ratio_B), Control_Ratio_A, Control_Ratio_B), `FC-A` = Experiment_Ratio / Control_Ratio_A, `FC-B` = Experiment_Ratio / Control_Ratio_B, log10_FC_A = log10(`FC-A`), log10_FC_B = log10(`FC-B`)) %>% arrange(desc(`FC-A`), Control_Frequency)

write_tsv(result, output_file)
writeLines(c("CRAPome FC 计算记录", paste("输入长表:", input_long), paste("对照库:", control_file), paste("实验类型:", experiment_type), paste("输出行数:", nrow(result))), "./5.质控/04_FC计算记录.txt", useBytes = TRUE)
cat("输出: ", output_file, "\n", sep = "")
