# CRAPome 污染过滤 R 工作流

本目录提供从 MaxQuant `proteinGroups.txt` 到 CRAPome 污染过滤结果的完整 R 工作流。核心目标是把 AP-MS 定量结果转换为 CRAPome 四列长表，结合本地对照库计算 `FC-A` / `FC-B`，并按“高富集倍数 + 低对照频率”筛选高置信互作蛋白。

## 目录结构

```text
05_CRAPome_污染过滤/
├── 1.原始输入/
├── 2.样本信息与分组/
├── 3.数据清洗与预处理/
├── 4.CRAPome输入文件/
├── 5.质控/
├── 6.CRAPome过滤结果/
├── 7.可视化/
└── 脚本/
```

## 算法原理

CRAPome/REPRINT 的思想是使用大量 AP-MS 阴性或背景对照实验统计每个蛋白在对照中的出现频率。常见污染物会在对照中高频出现，真实互作蛋白则倾向于在目标 bait 实验中富集、在对照库中低频出现。

本地实现使用以下近似计算：

```text
FC-A = (实验谱数 / 实验总谱数) / (所有对照谱数 / 所有对照总谱数)
FC-B = (实验谱数 / 实验总谱数) / (相同类型对照谱数 / 相同类型对照总谱数)
```

默认筛选阈值：

```text
FC-A >= 10
FC-B >= 5
Control_Frequency < 0.2
```

## 输入格式

默认输入：

```text
1.原始输入/proteinGroups.txt
```

脚本 `02_格式转换_CRAPome.R` 会生成 CRAPome 四列长表：

```text
Bait Name    AP Name    Prey Name    Spectral Count
```

可选样本信息文件：

```text
2.样本信息与分组/sample_info.tsv
```

至少包含：

```text
AP Name    Bait Name
```

## 对照库格式

最低必需列：

```text
Prey_ID    Control_Frequency    Total_Controls
```

推荐附加列：

```text
Control_Spectral_Count
Control_Total_Spectral_Count
Control_Type
```

如果无法直接访问 CRAPome 在线数据库，可到 `https://www.reprint-apms.org/` 手动导出对应物种或实验类型的对照统计表，整理为上述 TSV 文件后传给 `04_计算FC.R --control=`。脚本 `03_下载对照库.R` 会生成示例对照库，只用于测试流程，不用于正式生物学结论。

本地参考资料：

```text
8.文档与参考/02_官方页面/CRAPome_REPRINT主页.html
8.文档与参考/01_论文原文/CRAPome_原始论文.html
```

## 执行命令

在 `05_CRAPome_污染过滤` 目录下运行：

```powershell
Rscript ./脚本/01_数据清洗.R --input=./1.原始输入/proteinGroups.txt
Rscript ./脚本/02_格式转换_CRAPome.R --quant_type=auto --min_count=0
Rscript ./脚本/03_下载对照库.R --mode=example
Rscript ./脚本/04_计算FC.R --control=./4.CRAPome输入文件/control_library_example.tsv --experiment_type=all
Rscript ./脚本/05_结果筛选.R --min_fca=10 --min_fcb=5 --max_control_frequency=0.2
Rscript ./脚本/06_可视化.R --top_n=30
```

## 参数说明

`01_数据清洗.R`

- `--input=`：MaxQuant `proteinGroups.txt`
- `--unified_cleaned=`：统一准备层清洗结果，存在时优先复用
- `--output=`：清洗后 TSV

`02_格式转换_CRAPome.R`

- `--quant_type=`：定量列前缀，默认 `auto`，自动识别 `MS/MS count `、`Spectral count `、`Intensity `、`LFQ intensity `、`iBAQ `
- `--sample_info=`：样本到 bait 的映射表
- `--min_count=`：进入长表的最小谱数，默认 `0`

`03_下载对照库.R`

- `--mode=example`：生成示例对照库
- `--mode=download --url=...`：尝试从指定 URL 下载
- `--output=`：对照库输出路径

`04_计算FC.R`

- `--input=`：四列长表
- `--control=`：本地 CRAPome 对照库
- `--experiment_type=`：FC-B 使用的对照类型，默认 `all`
- `--pseudo_count=`：伪计数，默认 `1`

`05_结果筛选.R`

- `--min_fca=`：FC-A 阈值，默认 `10`
- `--min_fcb=`：FC-B 阈值，默认 `5`
- `--max_control_frequency=`：对照频率上限，默认 `0.2`

`06_可视化.R`

- `--input=`：带筛选标记的结果表
- `--top_n=`：Top 蛋白数量，默认 `30`

## 输出文件

主要结果：

```text
3.数据清洗与预处理/proteinGroups_cleaned.tsv
4.CRAPome输入文件/crapome_input_long.tsv
4.CRAPome输入文件/spectral_count_matrix.tsv
4.CRAPome输入文件/control_library_example.tsv
6.CRAPome过滤结果/crapome_fc_all.tsv
6.CRAPome过滤结果/crapome_filtered.tsv
6.CRAPome过滤结果/crapome_filtered_with_flags.tsv
6.CRAPome过滤结果/crapome_filter_summary.tsv
7.可视化/01_对照频率分布.png
7.可视化/02_FC散点图.png
7.可视化/03_Top蛋白条形图.png
```

## 依赖

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "httr"))
```

`httr` 只在 `03_下载对照库.R --mode=download` 时需要。
