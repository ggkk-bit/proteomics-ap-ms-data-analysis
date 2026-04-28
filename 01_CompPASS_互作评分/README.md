# CompPASS 互作评分

CompPASS 是用于大规模 AP-MS 数据集的互作评分方法，特别适合有较多不同 bait、能够估计 prey 在整个数据集背景出现频率的场景。

## 核心原理

CompPASS 强调：
- 某个 prey 在不同 bait 实验中的出现分布
- 同一 bait 下的重复支持
- 与背景型高频蛋白的区分

经典 CompPASS 依赖以下信息：
- bait 与 prey 的配对关系
- 每次 IP/run 中 prey 的丰度或谱数
- 同一 bait 的重复
- prey 在不同 bait 间的出现频率

## 目录结构

- `1.原始输入`：共享输入副本和原始 `proteinGroups.txt`
- `2.样本信息与分组`：样本注释、bait 映射、重复分组
- `3.数据清洗与预处理`：去除反向库、潜在污染物、仅位点鉴定蛋白后的中间结果
- `4.CompPASS输入文件`：最终用于 CompPASS 的长表输入
- `5.质控`：输入检查报告
- `6.CompPASS评分结果`：CompPASS 原始结果和筛选结果
- `7.可视化`：评分分布、Top prey、热图等图形输出
- `8.文档与参考`：输入格式、方法习惯、参考材料
- `脚本`：按任务拆分的脚本

## 执行流程

### 推荐执行顺序

1. 准备 `1.原始输入/proteinGroups.txt` 和 `archiveSummary.tsv`
2. 运行 `脚本/00_生成样本注释表_CompPASS.py`
3. 抽查并修正 `2.样本信息与分组/sample_annotation.tsv`
4. 运行 `脚本/01_数据清洗_CompPASS.py`
5. 运行 `脚本/02_格式转换_CompPASS.py`
6. 运行 `脚本/03_确认CompPASS输入文件.py`
7. 使用 R 包 `cRomppass` 完成评分
8. 参考 `脚本/04_结果筛选说明.md` 进行结果分层
9. 运行 `脚本/05_可视化_CompPASS.py`

## 输入要求

### CompPASS 实际输入格式

CompPASS 经典输入主表为一个长表，每行表示一个 `run-bait-prey` 组合，至少包含：

- `idRun`：一次 IP/run 的唯一标识
- `idBait`：该 run 对应的 bait
- `idPrey`：prey 蛋白标识
- `countPrey`：该 prey 在该 run 中的计数或丰度

本工作区的格式转换脚本会生成：
- `4.CompPASS输入文件/最终输入/comppass_input.tsv`
- `4.CompPASS输入文件/最终输入/bait_run_map.tsv`
- `4.CompPASS输入文件/最终输入/prey_reference.tsv`

### 样本注释要求

当前最关键的是 `2.样本信息与分组/sample_annotation.tsv`。至少要补齐：
- `sample_id`
- `ip_name`
- `bait_name`
- `replicate_group`
- `keep`

`脚本/00_生成样本注释表_CompPASS.py` 会优先从 `proteinGroups.txt` 中识别 `sample_id`，再尝试用 `archiveSummary.tsv` 自动回填 `raw_file_name`、`bait_name`、`plate`、`well`、`replicate_letter` 和 `replicate_group`。

## 完整执行命令

本说明对应目录：
```
C:\Users\dakey\Desktop\01_进行中\古老蛋白算法调研\多算法分析工作区\01_CompPASS_互作评分
```

建议统一在 `01_CompPASS_互作评分` 根目录执行命令。

### 运行前需要准备的文件

- `1.原始输入/proteinGroups.txt`
- `1.原始输入/archiveSummary.tsv`
- `8.文档与参考/03_源码或软件/cRomppass_master.zip`

脚本会生成或使用以下文件：
- `2.样本信息与分组/sample_annotation.tsv`
- `3.数据清洗与预处理/中间结果/proteinGroups.cleaned.tsv`
- `3.数据清洗与预处理/过滤记录/过滤统计.json`
- `4.CompPASS输入文件/最终输入/comppass_input.tsv`
- `5.质控/输入检查/CompPASS输入检查报告.txt`
- `6.CompPASS评分结果/原始结果/comppass_results.tsv`
- `7.可视化/图形输出/*`

### 进入项目根目录

```powershell
Set-Location 'C:\Users\dakey\Desktop\01_进行中\古老蛋白算法调研\多算法分析工作区\01_CompPASS_互作评分'
```

### 第 1 步：生成样本注释表

```powershell
python '.\脚本\00_生成样本注释表_CompPASS.py'
```

默认读取：
- `1.原始输入/proteinGroups.txt`
- `1.原始输入/archiveSummary.tsv`

默认输出：
- `2.样本信息与分组/sample_annotation.tsv`

作用：
- 从 `proteinGroups.txt` 识别样本列，例如 `Intensity b10494`
- 用 `archiveSummary.tsv` 自动回填 `bait_name`
- 自动回填 `raw_file_name`、`plate`、`well`、`replicate_letter`
- 自动生成 `ip_name` 和 `replicate_group`

可选参数示例：
```powershell
python '.\脚本\00_生成样本注释表_CompPASS.py' --input '.\1.原始输入\proteinGroups.txt' --archive '.\1.原始输入\archiveSummary.tsv' --output '.\2.样本信息与分组\sample_annotation.tsv'
```

### 第 2 步：人工抽查样本注释表

重点检查：
- `sample_id` 是否和 `proteinGroups.txt` 样本列后缀一致
- `bait_name` 是否正确
- `replicate_group` 是否符合实验设计
- `replicate_letter` 是否合理
- `keep` 是否都应保留

### 第 3 步：清洗 proteinGroups

```powershell
python '.\脚本\01_数据清洗_CompPASS.py'
```

默认输入：
- `1.原始输入/proteinGroups.txt`
- `2.样本信息与分组/sample_annotation.tsv`

默认输出：
- `3.数据清洗与预处理/中间结果/proteinGroups.cleaned.tsv`
- `3.数据清洗与预处理/过滤记录/过滤统计.json`

作用：
- 删除 `Reverse`
- 删除 `Potential contaminant`
- 删除 `Only identified by site`
- 删除没有 prey ID 的行
- 删除在所有保留样本中定量值都为 0 的行

### 第 4 步：转换为 CompPASS 输入长表

```powershell
python '.\脚本\02_格式转换_CompPASS.py'
```

默认输出：
- `4.CompPASS输入文件/最终输入/comppass_input.tsv`
- `4.CompPASS输入文件/最终输入/bait_run_map.tsv`
- `4.CompPASS输入文件/最终输入/prey_reference.tsv`

说明：
- `comppass_input.tsv` 是真正给官方 `cRomppass` 用的主输入来源
- 核心列为 `idRun`、`idBait`、`idPrey`、`countPrey`
- 当前脚本优先使用 `MS/MS Count`，如果没有，再退回 `Intensity` 或 `LFQ intensity`

### 第 5 步：检查输入文件

```powershell
python '.\脚本\03_确认CompPASS输入文件.py'
```

默认输出：
- `5.质控/输入检查/CompPASS输入检查报告.txt`

重点检查：
- 是否有核心列
- `countPrey` 是否全为正数值
- 一个 `idRun` 是否只对应一个 `idBait`
- 是否有重复的 `run-bait-prey` 记录
- bait 数量和重复数是否太少

### 第 6 步：用官方 cRomppass 做正式评分

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' '.\脚本\06_运行官方CompPASS_cRomppass.R'
```

默认读取：
- `4.CompPASS输入文件/最终输入/comppass_input.tsv`
- `4.CompPASS输入文件/最终输入/bait_run_map.tsv`
- `4.CompPASS输入文件/最终输入/prey_reference.tsv`
- `8.文档与参考/03_源码或软件/cRomppass_master.zip`

默认输出：
- `6.CompPASS评分结果/原始结果/comppass_official_input.tsv`
- `6.CompPASS评分结果/原始结果/comppass_results.tsv`

说明：
- 这里直接使用官方 `cRomppass`
- 该 R 脚本会把当前长表转换成官方实现需要的六列表结构
- 该 R 脚本依赖 `dplyr` 和 `magrittr`

可选参数示例：
```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' '.\脚本\06_运行官方CompPASS_cRomppass.R' --norm_factor=0.98
```

### 第 7 步：对评分结果做可视化

只有在这一步文件已经存在时才能运行：
- `6.CompPASS评分结果/原始结果/comppass_results.tsv`

运行命令：
```powershell
python '.\脚本\05_可视化_CompPASS.py'
```

默认读取：
- `6.CompPASS评分结果/原始结果/comppass_results.tsv`

默认输出：
- `7.可视化/图形输出/01_WD分数分布.png`
- `7.可视化/图形输出/02_WD_Z散点图.png`
- `7.可视化/图形输出/03_top_prey_*.png`
- `7.可视化/图形输出/04_top_bait_prey_wd_heatmap.png`
- `7.可视化/图形输出/05_ppi_global_network.png`
- `7.可视化/图形输出/06_ppi_core_module.png`

说明：
- `05_ppi_global_network.png`：高分 `bait-prey` 总览图。每个 bait 作为中心节点，周围是该 bait 的高置信 prey
- `06_ppi_core_module.png`：多 bait 的重点星形图。每个小面板只看一个 bait，中心是 bait，周围是该 bait 的 top prey

### 一键顺序执行命令

```powershell
Set-Location 'C:\Users\dakey\Desktop\01_进行中\古老蛋白算法调研\多算法分析工作区\01_CompPASS_互作评分'
python '.\脚本\00_生成样本注释表_CompPASS.py'
python '.\脚本\01_数据清洗_CompPASS.py'
python '.\脚本\02_格式转换_CompPASS.py'
python '.\脚本\03_确认CompPASS输入文件.py'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' '.\脚本\06_运行官方CompPASS_cRomppass.R'
python '.\脚本\05_可视化_CompPASS.py'
```

## 结果解释

### 关键输出文件

- 输入主表：`4.CompPASS输入文件/最终输入/comppass_input.tsv`
- 官方评分结果：`6.CompPASS评分结果/原始结果/comppass_results.tsv`
- 可视化输入：`6.CompPASS评分结果/原始结果/comppass_results.tsv`

### 筛选建议

CompPASS 没有一个跨数据集完全统一的硬阈值。实践中通常：
- 以 `scoreWD` 或 `NWD/WD` 作为主排序指标
- 结合 `scoreZ` 辅助判断 enrichment
- 再结合重复、一致性、污染过滤和生物学合理性做分层

如果后续要与 MiST、SAINTexpress、CRAPome 联合，建议把 CompPASS 结果作为主评分证据层之一，而不是单独做刚性截断。

详细筛选说明见：[脚本/04_结果筛选说明.md](脚本/04_结果筛选说明.md)

## 常见问题

### 1. 为什么 `05_可视化_CompPASS.py` 会报找不到 `comppass_results.tsv`？

因为它画的是正式评分后的结果，不是输入长表。必须先跑第 6 步。

### 2. 为什么 `02_格式转换_CompPASS.py` 以前很慢？

之前版本是逐行循环，现在已经改成更快的写法，正常情况下很快就能完成。

### 3. `archiveSummary.tsv` 是必须的吗？

对自动生成 `sample_annotation.tsv` 来说，强烈建议提供。如果没有它，也能手工做样本注释，但你需要自己填写 `bait_name`、重复和分组信息。

### 4. 哪个结果文件最重要？

如果你要继续做图和筛选，最重要的是：
- `6.CompPASS评分结果/原始结果/comppass_results.tsv`

## 脚本清单

- `00_生成样本注释表_CompPASS.py`
- `01_数据清洗_CompPASS.py`
- `02_格式转换_CompPASS.py`
- `03_确认CompPASS输入文件.py`
- `05_可视化_CompPASS.py`
- `06_运行官方CompPASS_cRomppass.R`

## 参考文档

- [脚本/03_输入确认说明.md](脚本/03_输入确认说明.md) - 输入确认逻辑和通过标准
- [脚本/04_结果筛选说明.md](脚本/04_结果筛选说明.md) - WD/Z 等分数的筛选思路和阈值习惯
- `8.文档与参考/` - 论文、源码、参考材料

## 注意事项

经典实现通常以谱数或肽段计数作为 `countPrey`。如果当前数据只有 `Intensity` 或 `LFQ intensity`，可以先作为近似输入进行方法测试，但应在方法部分明确说明这不是最标准的 CompPASS 计数来源。
