# 统一数据准备层

本目录是 6 个互作分析算法前面的统一入口。它把 MaxQuant `proteinGroups.txt` 先清洗成一份通用表，再生成统一样本注释，并按需转换为 CompPASS、MiST 等算法自己的输入格式。

现有 `01_CompPASS_互作评分/`、`02_MiST_互作评分/` 等目录不会被修改，原来的单独分析流程仍可继续使用。这里新增的是可复用的 Phase 2 数据准备层。

## 目录作用

- `1.原始输入/`：放置 `proteinGroups.txt` 和可选的 `archiveSummary.tsv`。也可以不复制文件，运行脚本时用参数指定路径。
- `2.样本注释/`：输出 `sample_annotation_master.tsv`，同时兼容 MiST 和 CompPASS 的字段。
- `3.清洗后数据/`：输出 `proteinGroups_cleaned.tsv` 和 `cleaning_report.json`。
- `4.算法特定输入/compass_input/`：输出 CompPASS 所需输入。
- `4.算法特定输入/mist_input/`：输出 MiST 所需输入。
- `脚本/`：统一 R 脚本。

## 使用方式 1：保持单独分析

不运行本目录脚本，继续使用各算法目录中原有脚本即可。此方式用于复现实验或避免改变既有结果。

## 使用方式 2：统一准备后再单独分析

在本目录生成统一输入，再把结果复制或作为参数传给各算法后续脚本。

```powershell
cd 00_统一数据准备
Rscript 脚本/02_generate_sample_annotation.R --input="../proteinGroups.txt" --archive="../01_CompPASS_互作评分/1.原始输入/archiveSummary.tsv"
Rscript 脚本/01_clean_protein_groups.R --input="../proteinGroups.txt"
Rscript 脚本/03_transform_for_compass.R
Rscript 脚本/04_transform_for_mist.R
```

## 使用方式 3：复合分析

复合分析时只维护一份 `sample_annotation_master.tsv` 和一份 `proteinGroups_cleaned.tsv`，然后为多个算法生成输入。

```powershell
cd 00_统一数据准备
Rscript 脚本/02_generate_sample_annotation.R --input="../proteinGroups.txt" --archive="../01_CompPASS_互作评分/1.原始输入/archiveSummary.tsv"
Rscript 脚本/01_clean_protein_groups.R --input="../proteinGroups.txt"
Rscript 脚本/03_transform_for_compass.R
Rscript 脚本/04_transform_for_mist.R --quant_type="Intensity"
```

## 配置参数

所有脚本都提供默认路径，也支持命令行覆盖，避免在脚本中硬编码项目路径。

- `01_clean_protein_groups.R`
  - `--input=...`
  - `--output=...`
  - `--report=...`
  - `--drop_all_zero_quant=true`
  - `--quant_prefixes="Intensity ;LFQ intensity ;iBAQ ;MS/MS count ;MS/MS Count ;Intensity."`
- `02_generate_sample_annotation.R`
  - `--input=...`
  - `--archive=...`
  - `--output=...`
  - `--prefixes=...`
- `03_transform_for_compass.R`
  - `--cleaned=...`
  - `--annotation=...`
  - `--output_dir=...`
- `04_transform_for_mist.R`
  - `--cleaned=...`
  - `--annotation=...`
  - `--output_dir=...`
  - `--quant_type=auto`

## 与现有流程的关系

统一层只新增文件，不修改现有 6 个算法目录。后续如果要逐步减少重复脚本，可以让各算法目录的清洗和格式转换步骤读取本目录的统一输出；在完全验证之前，建议保留原脚本作为可回退流程。
