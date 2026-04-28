# PPIrank 网络补充排序

PPIrank 是 TAP/MS 数据的网络层补充排序方法，不直接替代主打分算法，而是把 `spectral counts`、`negative controls` 和重复性整合成可解释的优先级结果。

## 核心原理

PPIrank 定位为网络层补充排序，强调：
- spectral counts
- negative controls
- biological replicates
- 互作在重复中的稳定性

它不是直接读取宽表评分，而是要求先整理成样本级长表。

## 目录结构

- `1.原始输入`：原始 `proteinGroups.txt`、原始导出表、补充元数据
- `2.样本信息与分组`：bait、replicate、control、condition 的样本注释
- `3.数据清洗与预处理`：去除 `Reverse`、`Potential contaminant`、`Only identified by site` 等噪声
- `4.PPIrank输入文件`：转成 PPIrank 可直接读取的长表与附加注释表
- `5.质控`：重复性、缺失值、负对照覆盖度、单次谱图占比检查
- `6.PPIrank排序结果`：打分结果、阈值结果、候选互作列表
- `7.可视化`：网络图、排序分布、与参考库/其他算法的重叠图
- `8.文档与参考`：论文原文、官方页面、下载记录与方法说明
- `脚本`：按任务拆分的独立脚本

## 执行流程

### 流程顺序

1. 读入 `proteinGroups.txt` 和样本注释
2. 去除 `Reverse`、`Potential contaminant`、`Only identified by site`
3. 统一 bait 名、重复号、对照标记和 prey ID
4. 生成 PPIrank 输入长表，必要时补蛋白长度和对照汇总
5. 做质控，重点看重复性、负对照覆盖度和单次谱图计数
6. 按 PPIrank 结果阈值筛选候选互作
7. 输出网络图、Top prey 条形图、方法间重叠图

## 输入要求

### 关键方法约束

- PPIrank 依赖 `spectral counts`
- 需要 `negative controls`，并且至少 `3` 个生物学重复
- 论文明确提到：若某个互作在三个重复里都只有 `1` 个谱图计数，容易被视为噪声，应过滤
- 文献展示的核心图形是：bait 标绿的互作网络图、与已知互作的重叠比较图、不同方法的 pairwise overlap 图

### 最小输入字段

**必需字段：**
- `bait_name`
- `replicate_id`
- `prey_id`
- `spectral_count`
- `is_control`

**推荐字段：**
- `sample_id`
- `condition`
- `protein_length`
- `gene_symbol`
- `source_file`

### 本地约定输入

- 原始清洗前数据：`1.原始输入/proteinGroups.txt`
- 样本注释：`2.样本信息与分组/sample_annotation_template.tsv`
- 清洗后中间表：`3.数据清洗与预处理/清洗后蛋白组.tsv`
- PPIrank 输入长表：`4.PPIrank输入文件/PPIrank_输入_长表.tsv`
- PPIrank 负对照子表：`4.PPIrank输入文件/PPIrank_输入_负对照.tsv`
- PPIrank 样本汇总：`4.PPIrank输入文件/PPIrank_输入_样本汇总.tsv`
- PPIrank 最终排序：`6.PPIrank排序结果/PPIrank_排序结果.tsv`

## 结果解释

- **高分结果**：更适合作为网络补充候选
- **中分结果**：保留给人工复核或跨算法一致性检查
- **低分结果**：通常只在和已知互作重叠时保留

## 脚本清单

- `01_数据清洗脚本.py`
- `02_格式转换脚本.py`
- `03_输入文件确认脚本.py`
- `04_结果筛选说明.md`
- `05_可视化脚本.py`

## 文献依据

- PPIrank 原文：`PPIRank - an advanced method for ranking protein-protein interactions in TAP/MS data`
- 关键展示：
  - bait 标绿的网络图
  - 与 BioGRID 已知互作的重叠比较
  - 不同方法的 pairwise overlap
