# Phase 4：复合分析与 STRING 数据库验证工作流

本目录用于整合 CompPASS、MiST、HGSCore、CS_Score、CRAPome、PPIrank 六类算法结果，并与 STRING v12.0 数据库进行对照验证，输出多算法一致性、STRING 验证率、新发现候选互作、综合评分、可视化图和 Cytoscape/Excel 文件。

## 目录结构

- `1.算法结果汇总/原始结果`：复制后的 6 个算法原始输出
- `1.算法结果汇总/标准化结果`：统一格式结果
- `2.STRING数据库/原始下载`：STRING API 或批量文件下载
- `2.STRING数据库/物种映射`：预留物种和 ID 映射辅助文件
- `2.STRING数据库/处理后数据`：处理后的 STRING links 文件，默认列为 `protein1, protein2, combined_score`
- `3.蛋白ID映射/映射表`：蛋白 ID 到 STRING Ensembl Protein ID 映射表
- `3.蛋白ID映射/映射日志`：映射统计日志
- `4.算法交叉验证`：算法重叠矩阵、一致性评分、Venn/UpSet 数据
- `5.STRING验证对照`：验证明细、验证率统计、新发现互作
- `6.综合评分`：综合评分和高置信互作
- `7.可视化`：重叠图、STRING 验证图、综合网络图、性能对比图
- `8.导出结果`：Excel、Cytoscape、HTML 报告
- `脚本`：全部 R 脚本和公共函数库

## 依赖包

```r
install.packages(c(
  "dplyr", "tidyr", "httr", "jsonlite", "ggplot2",
  "VennDiagram", "ggvenn", "pheatmap", "igraph", "openxlsx"
))
```

## 默认参数

- 物种：人类 `9606`
- STRING 版本：v12.0
- STRING 阈值：`400`
- 综合评分权重：`w1=0.4, w2=0.3, w3=0.3`
- 综合评分公式：`综合分数 = w1 * 支持算法数量归一化 + w2 * 归一化平均分数 + w3 * STRING验证加分`

## 执行顺序

在工作区根目录运行：

```powershell
Rscript "07_复合分析与整合/脚本/01_收集算法结果.R"
Rscript "07_复合分析与整合/脚本/02_标准化格式.R"
Rscript "07_复合分析与整合/脚本/03_下载STRING数据.R" --species=9606 --min_score=400 --mode=url
Rscript "07_复合分析与整合/脚本/04_蛋白ID映射.R" --species=9606 --mapping_file=你的蛋白到STRING映射表.tsv
Rscript "07_复合分析与整合/脚本/05_算法交叉验证.R"
Rscript "07_复合分析与整合/脚本/06_STRING验证对照.R" --min_score=400 --string_file="07_复合分析与整合/2.STRING数据库/处理后数据/string_links_filtered.tsv"
Rscript "07_复合分析与整合/脚本/07_综合评分.R" --w1=0.4 --w2=0.3 --w3=0.3 --high_conf=0.7
Rscript "07_复合分析与整合/脚本/08_可视化_重叠分析.R"
Rscript "07_复合分析与整合/脚本/09_可视化_STRING验证.R"
Rscript "07_复合分析与整合/脚本/10_可视化_综合网络.R" --top_n=200
Rscript "07_复合分析与整合/脚本/11_导出Excel报告.R"
Rscript "07_复合分析与整合/脚本/12_导出Cytoscape.R"
```

## 参数说明

所有脚本均使用 `commandArgs(trailingOnly = TRUE)` 解析 `--key=value` 参数。

- `--species=9606`：STRING 物种编号
- `--min_score=400`：STRING combined score 或 API score 最低阈值
- `--mode=url|api`：STRING 下载方式，`url` 为批量文件，`api` 为基于当前蛋白列表查询
- `--mapping_file=...`：两列表格，第一列为本项目蛋白 ID，第二列为 STRING ID，例如 `9606.ENSP000...`
- `--string_file=...`：STRING links 文件，支持列 `protein1, protein2, combined_score` 或 `stringId_A, stringId_B, score`
- `--w1/--w2/--w3`：综合评分权重
- `--high_conf=0.7`：高置信互作输出阈值
- `--top_n=200`：网络图最多绘制的边数

## STRING 验证逻辑

`06_STRING验证对照.R` 是核心脚本，主要步骤如下：

1. 读取标准化后的六算法互作表。
2. 读取蛋白到 STRING Ensembl Protein ID 的映射表。
3. 读取 STRING links，并按 `combined_score >= min_score` 过滤。
4. 将 bait/prey 映射为 STRING ID，构建无方向蛋白对 key。
5. 与 STRING 蛋白对匹配，标记 `mapped_to_string` 和 `string_validated`。
6. 输出按算法分层的验证率。
7. 输出按 1-6 个算法支持数分层的验证率。
8. 将已映射但未被 STRING 验证的多算法支持蛋白对导出为候选新发现互作。

## 输入文件要求

六个算法默认输入路径如下；若不存在，`01_收集算法结果.R` 会跳过并写入日志，后续可替换为真实文件。

- `01_CompPASS_互作评分/6.CompPASS评分结果/原始结果/comppass_results.tsv`
- `02_MiST_互作评分/6.MiST评分结果/preprocessed_NoC_MAT_MIST.txt`
- `03_HGSCore_辅助评分/6.HGSCore评分结果/hgscore_filtered.tsv`
- `04_CS_Score_辅助评分/6.CS_Score评分结果/cs_score_filtered_pairs.tsv`
- `05_CRAPome_污染过滤/6.CRAPome过滤结果/crapome_filtered.tsv`
- `06_PPIrank_网络补充排序/6.PPIrank排序结果/ppirank_filtered.tsv`

## 关键输出

- `1.算法结果汇总/标准化结果/all_algorithms_standardized.tsv`
- `4.算法交叉验证/重叠分析/algorithm_overlap_matrix.tsv`
- `4.算法交叉验证/一致性评分/pair_consistency_scores.tsv`
- `5.STRING验证对照/验证结果/algorithm_pairs_string_validation.tsv`
- `5.STRING验证对照/验证率统计/validation_rate_by_algorithm.tsv`
- `5.STRING验证对照/验证率统计/validation_rate_by_support_count.tsv`
- `5.STRING验证对照/新发现互作/candidate_novel_interactions.tsv`
- `6.综合评分/评分结果/integrated_interaction_scores.tsv`
- `6.综合评分/高置信互作/high_confidence_interactions.tsv`
- `8.导出结果/Excel报告/phase4_integration_report.xlsx`
- `8.导出结果/Cytoscape文件/high_confidence_interactions.sif`
