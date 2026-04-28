# CS_Score 辅助评分 R 工作流

本目录提供从 `proteinGroups.txt` 到 CS_Score 蛋白对评分结果与可视化的完整 R 工作流。

## 算法依据

论文 `A Novel Scoring Approach for Protein Co-Purification Data Reveals High Interaction Specificity` 的核心思想是：对 AP/MS purification 数据中的蛋白对统计共现次数，再通过受约束随机打乱构建背景分布，最后用 Z-score 表示实际共现相对随机背景的偏离程度。论文补充表说明只报告观察共现次数大于 1 的蛋白对，因此本实现默认仅对 `observed_cooccurrence >= 2` 的蛋白对输出评分。

本实现的受约束随机打乱保持每个 purification 的成员数量不变，在全部 prey 蛋白背景中随机分配成员，重复 `N` 次，计算：

```text
CS_Score Z = (observed_cooccurrence - mean(random_cooccurrence)) / sd(random_cooccurrence)
```

## 目录结构

```text
04_CS_Score_辅助评分/
├── 1.原始输入/
├── 2.样本信息与分组/
├── 3.数据清洗与预处理/
├── 4.CS_Score输入文件/
├── 5.质控/
├── 6.CS_Score评分结果/
├── 7.可视化/
└── 脚本/
```

## 运行顺序

在 `04_CS_Score_辅助评分` 目录下运行：

```powershell
Rscript ./脚本/01_数据清洗.R --input=./1.原始输入/proteinGroups.txt
Rscript ./脚本/02_格式转换_CS_Score.R --quant_type=auto --presence_threshold=0
Rscript ./脚本/03_计算CS_Score.R --n_perm=1000 --seed=1 --min_observed=2
Rscript ./脚本/04_结果筛选.R --min_cooccurrence=2 --min_z=2
Rscript ./脚本/05_可视化.R --top_n=30
```

## 输入

默认输入为 MaxQuant `proteinGroups.txt`：

- `1.原始输入/proteinGroups.txt`
- 可选样本注释：`../../00_统一数据准备/2.样本注释/sample_annotation_master.tsv`
- 可选统一清洗结果：`../../00_统一数据准备/3.清洗后数据/proteinGroups_cleaned.tsv`

定量列会自动识别以下前缀：`LFQ intensity `、`Intensity `、`iBAQ `、`MS/MS count `。也可以用 `--quant_type=` 指定。

## 输出

- `3.数据清洗与预处理/proteinGroups_cleaned.tsv`
- `4.CS_Score输入文件/cs_score_member_long.tsv`
- `4.CS_Score输入文件/cs_score_binary_matrix.tsv`
- `6.CS_Score评分结果/cs_score_all_pairs.tsv`
- `6.CS_Score评分结果/cs_score_filtered_pairs.tsv`
- `7.可视化/01_共现次数分布.png`
- `7.可视化/02_Zscore分布.png`
- `7.可视化/03_Top蛋白对.png`
- `7.可视化/04_Top蛋白对热图.png`

## 依赖

R 包：`dplyr`、`tidyr`、`ggplot2`。

