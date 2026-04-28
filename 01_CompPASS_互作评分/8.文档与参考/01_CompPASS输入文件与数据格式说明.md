# CompPASS 输入文件与数据格式说明

## 1. CompPASS 经典输入要求

根据 CompPASS 方法文献和 `cRomppass` 这类 R 包实现的接口，核心输入是一个数据框或制表文件，每行表示一个唯一的 bait-prey-run 观测，至少包含四列：

- `idRun`：一次 AP-MS 实验运行或一个 IP 样本的唯一标识
- `idBait`：该 run 对应的 bait 名称
- `idPrey`：被捕获 prey 的蛋白标识
- `countPrey`：prey 在该 run 中的计数值

在 `cRomppass` 这类实现里，这四列是最核心的必要字段。需要注意：`cRomppass` 是一个 R 包实现，不等于原始 CompPASS 软件本体。输出中常见字段包括：

- `AvePSM`
- `scoreZ`
- `scoreS`
- `scoreD`
- `Entropy`
- `scoreWD`

## 2. 本工作区的落地输入

### 2.1 原始输入

- `1.原始输入/proteinGroups.txt`
- `1.原始输入/archiveSummary.tsv`
- `2.样本信息与分组/sample_annotation.tsv`

### 2.2 样本注释表建议字段

- `sample_id`：必须能映射到 `proteinGroups.txt` 的样本列后缀
- `raw_file_name`：原始 raw 文件名，便于追踪
- `ip_name`：每个 IP/run 的唯一名称，将写入 `idRun`
- `bait_name`：bait 名称，将写入 `idBait`
- `replicate_group`：同一 bait 的重复组
- `replicate_letter`：A/B/C 等重复标签
- `keep`：是否纳入分析，建议使用 `1/0`、`TRUE/FALSE` 或 `yes/no`

如果 `proteinGroups.txt` 的样本列是 `b13358` 这类 BioPlex 风格样本 ID，而不是完整 raw 文件名，建议先使用 `archiveSummary.tsv` 自动生成样本注释表。当前工作区提供：

- `脚本/00_生成样本注释表_CompPASS.py`

它会按 `newName -> sample_id` 的方式匹配，例如 `b13358.raw -> b13358`，并自动解析 `info` 字段里的 bait、板号和孔位信息。

## 3. 由 proteinGroups 到 CompPASS 输入的映射

### 3.1 prey 标识

建议优先顺序：

1. `Majority protein IDs` 的首个条目
2. `Protein IDs` 的首个条目
3. 同时保留 `Gene names` 作为辅助注释

### 3.2 countPrey 来源

CompPASS 经典上更偏向：

- 谱数（spectral count）
- 肽段计数
- 其他代表 prey 丰度且适于 run 间比较的计数

本工作区脚本使用以下优先级自动选择样本列：

1. `MS/MS Count <sample>` 或 `MS/MS count <sample>`
2. `Intensity.<sample>` 或 `Intensity <sample>`
3. `LFQ intensity <sample>`

如果实际使用了 `Intensity` 或 `LFQ intensity`，请在方法说明中写明“用于 CompPASS 风格排序测试，不是经典谱数输入”。

## 4. 本工作区输出的 CompPASS 输入文件

### 4.1 主输入文件

`4.CompPASS输入文件/最终输入/comppass_input.tsv`

字段：

- `idRun`
- `idBait`
- `idPrey`
- `countPrey`
- `countType`
- `sourceColumn`
- `preyGene`

其中前四列是 CompPASS 实际评分所需核心列。

### 4.2 辅助文件

`4.CompPASS输入文件/最终输入/bait_run_map.tsv`

字段：

- `idRun`
- `idBait`
- `sample_id`
- `replicate_group`
- `replicate_letter`
- `raw_file_name`

`4.CompPASS输入文件/最终输入/prey_reference.tsv`

字段：

- `idPrey`
- `preyGene`
- `majority_protein_ids`

## 5. 输入构建前的必要过滤

至少应删除：

- `Reverse` / `反向匹配`
- `Potential contaminant` / `潜在污染物`
- `Only identified by site` / `仅位点鉴定蛋白`

建议追加删除：

- 没有可用 prey 标识的条目
- 在全部保留样本中定量值都为 0 或缺失的条目

追加过滤的依据是：这些条目无法为 CompPASS 提供可解释的 bait-prey-run 观测，保留会增加背景噪音。

## 6. 对结果稳定性的提醒

CompPASS 依赖“整个 bait 集合中的 prey 频率”。因此：

- bait 数量太少时，背景频率估计会不稳
- 全部 bait 高度相关时，特异性信息会变弱
- 没有重复时，结果可解释性会下降

建议在运行前先用 `脚本/03_确认CompPASS输入文件.py` 检查 bait 数量、重复数和输入完整性。

