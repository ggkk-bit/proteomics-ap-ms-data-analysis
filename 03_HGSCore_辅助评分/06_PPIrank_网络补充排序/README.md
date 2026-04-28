# PPIrank 网络补充排序工作流

本目录提供 PPIrank 的完整 R 语言工作流，用于从 MaxQuant `proteinGroups.txt` 或统一准备层清洗结果生成 bait-prey 网络排序结果、筛选高置信互作，并输出网络图与 Cytoscape 文件。

## 目录结构

```text
06_PPIrank_网络补充排序/
├── 1.原始输入/
├── 2.样本信息与分组/
├── 3.数据清洗与预处理/
├── 4.PPIrank输入文件/
├── 5.质控/
├── 6.PPIrank排序结果/
├── 7.可视化/
├── 8.文档与参考/
└── 脚本/
    ├── 01_数据清洗.R
    ├── 02_格式转换_PPIrank.R
    ├── 03_计算PPIrank.R
    ├── 04_结果筛选.R
    └── 05_可视化.R
```

## 输入要求

核心输入为样本级长表，字段必须为：

```text
bait_name    replicate_id    prey_id    spectral_count    is_control
```

`02_格式转换_PPIrank.R` 可从 MaxQuant 宽表自动生成该长表。样本注释文件建议保存为 `2.样本信息与分组/sample_info.tsv`：

```text
sample_id    bait_name    replicate_id    is_control
BaitA_rep1   BaitA        rep1            FALSE
BaitA_rep2   BaitA        rep2            FALSE
BaitA_rep3   BaitA        rep3            FALSE
Ctrl_1       Control      ctrl1           TRUE
Ctrl_2       Control      ctrl2           TRUE
```

论文要求每个 bait 至少 3 个生物学重复。脚本会对不足 3 个重复的 bait 给出警告，但不会强制中止，便于调试旧数据。

## 算法原理

PPIrank 没有官方 R 包。本实现根据论文工作流和本地任务要求提取网络排序逻辑：

1. 数据整合：整合多个生物学重复的谱数，区分 bait 样本和 negative controls。
2. 噪声过滤：若 bait-prey 在至少 3 个重复中出现，但每个重复都只有 1 个谱图计数，则按低谱数噪声过滤。
3. 重复性评分：`compute_reproducibility()` 统计每个 prey 在同一 bait 的多少个重复中出现。
4. 对照惩罚：prey 在 controls 中出现越多、谱数越高，惩罚越强。
5. 网络排序：`compute_ppirank()` 使用谱数、重复性和对照惩罚生成 `ppirank_score`。
6. 网络构建：`build_network()` 输出 bait-prey 边表和节点表，bait 节点标记为绿色。

由于当前可访问目录中未找到 `8.文档与参考/01_论文原文/PPIrank_原始论文.html`，核心分数采用可追溯的组合假设：

```text
ppirank_score = log1p(bait_total_spectral) * replicate_fraction * control_penalty
control_penalty = 1 / (1 + control_weight * log1p(control_total_spectral + control_replicates_observed))
```

该分数保留论文强调的三个信息源：谱数强度、生物学重复稳定性、negative controls 背景证据。

## 执行命令

在 `06_PPIrank_网络补充排序` 目录下运行：

```bash
Rscript 脚本/01_数据清洗.R
Rscript 脚本/02_格式转换_PPIrank.R --sample_info=./2.样本信息与分组/sample_info.tsv
Rscript 脚本/03_计算PPIrank.R --min_replicates=3 --noise_count=1 --control_weight=1
Rscript 脚本/04_结果筛选.R --min_score=1 --min_replicates=2 --max_control_total=5
Rscript 脚本/05_可视化.R --min_score=1 --top_n=100
```

如果没有 `sample_info.tsv`，格式转换脚本会生成 `sample_info_template.tsv`。补全 bait、重复和对照信息后重新运行转换步骤。

## 参数说明

`01_数据清洗.R`

- `--unified_cleaned=`：统一准备层清洗结果，默认 `../00_统一数据准备/3.清洗后数据/proteinGroups_cleaned.tsv`
- `--input=`：本地 MaxQuant `proteinGroups.txt`
- `--output=`：清洗后输出文件

`02_格式转换_PPIrank.R`

- `--input=`：清洗后的 proteinGroups 表
- `--sample_info=`：样本注释文件
- `--quant_type=`：定量列前缀，默认 `auto`，支持 `MS/MS count `、`Spectral count `、`Intensity ` 等
- `--min_count=`：长表保留的最小谱数，默认 0

`03_计算PPIrank.R`

- `--input=`：PPIrank 样本级长表
- `--min_replicates=`：论文推荐最小生物学重复数，默认 3
- `--noise_count=`：低谱数噪声阈值，默认 1
- `--control_weight=`：对照惩罚权重，默认 1

`04_结果筛选.R`

- `--min_score=`：PPIrank 分数阈值，默认 1
- `--min_replicates=`：至少出现的 bait 重复数，默认 2
- `--max_control_total=`：control 总谱数上限，默认 5
- `--keep_class=`：保留置信类别，默认 `high,medium`

`05_可视化.R`

- `--min_score=`：进入网络图的最低分数
- `--top_n=`：绘图边数上限
- `--out_dir=`：可视化输出目录

## 输出文件

主要结果：

- `4.PPIrank输入文件/ppirank_input_long.tsv`
- `6.PPIrank排序结果/ppirank_scores.tsv`
- `6.PPIrank排序结果/ppirank_high_confidence.tsv`
- `6.PPIrank排序结果/ppirank_edges.tsv`
- `6.PPIrank排序结果/ppirank_nodes.tsv`

可视化与 Cytoscape：

- `7.可视化/ppirank_network.png`
- `7.可视化/ppirank_score_distribution.png`
- `7.可视化/ppirank_replicate_summary.png`
- `7.可视化/cytoscape_edges.tsv`
- `7.可视化/cytoscape_nodes.tsv`
- `7.可视化/ppirank_network.sif`

## 依赖

必需 R 包：

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "igraph"))
```

可选：如果后续要直接连接 Cytoscape，可安装 `RCy3`。当前工作流不强制依赖 Cytoscape，默认输出可导入 Cytoscape 的 `.sif` 和属性表。
