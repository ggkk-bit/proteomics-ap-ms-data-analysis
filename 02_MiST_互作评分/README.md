# MiST 互作评分

MiST 是用于 AP-MS / Affinity Purification-MS 数据的互作优先级评分方法，核心依赖三个维度：

- **abundance**：候选互作蛋白的丰度
- **reproducibility**：跨重复的一致性
- **specificity**：相对于其他 bait 的特异性

MiST 适合把已经做过基础清洗的蛋白鉴定结果，整理成可解释的 bait-prey 互作列表，再输出高置信互作网络。

## 目录结构

- `1.原始输入`：原始 `proteinGroups.txt`、样本表、补充说明
- `2.样本信息与分组`：bait、重复、对照、组别定义
- `3.数据清洗与预处理`：基础过滤、字段统一、质量检查
- `4.MiST输入文件`：MiST 官方输入文件与可选控制文件
- `5.质控`：输入检查、重复相关性、缺失率、分布检查
- `6.MiST评分结果`：MiST 原始结果与筛选结果
- `7.可视化`：阈值、分布、网络、Top prey 图
- `8.文档与参考`：论文、官方实现、格式依据
- `脚本`：按任务拆分的脚本和说明

## 执行流程

### 推荐执行顺序

1. 先确认 `proteinGroups.txt` 里是否有可用的强度列、分子量列和样本列
2. 按 contaminant / reverse / site-only 规则清洗
3. 根据样本注释表生成 `data.txt` 和 `keys.txt`
4. 运行 MiST 预处理和评分
5. 按阈值与重复支持筛出高置信互作
6. 输出图表和最终表格

### 评分原则

- 先保留重复中稳定出现的 prey
- 再优先保留高特异性、非背景型互作
- 最后用 `MiST >= 0.75` 作为默认高置信起点

## 输入要求

### proteinGroups.txt 优先使用的字段

- `Protein IDs` / `Majority protein IDs`
- `Gene names` / `Protein names`
- `Unique peptides` / `Razor + unique peptides` / `Peptides`
- `Mol. weight [kDa]`
- `Intensity <sample>` 系列列

### 基础清洗规则

必须删除：
- `Reverse`
- `Potential contaminant`
- `Only identified by site`

建议额外删除：
- 全部保留样本中定量值都为 0 的蛋白
- 缺少稳定蛋白 ID 的蛋白
- 如果项目更严格，可以要求 `Unique peptides >= 1` 或 `>= 2`

### MiST 实际输入要求

**data.txt 推荐字段：**
- `id`
- `ms_uniprot_ac`
- `ms_num_unique_peptide`
- `ms_protein_mw`
- 可选：`ms_protein_name`、排序或丰度列

**keys.txt 最少字段：**
- 第 1 列：`id`
- 第 2 列：`bait_name`

**可选辅助文件：**
- `remove.txt`
- `collapse.txt`
- `specificity_exclusions.txt`

## 完整执行命令

本说明对应目录：
```
C:\Users\dakey\Desktop\01_进行中\古老蛋白算法调研\多算法分析工作区\02_MiST_互作评分
```

建议所有命令都在 `02_MiST_互作评分` 根目录执行。脚本默认使用相对路径。

### 运行前检查

至少需要这些文件：
- `1.原始输入/proteinGroups.txt`
- `1.原始输入/archiveSummary.tsv`
- `4.MiST输入文件/mist.yml`
- `8.文档与参考/03_源码或软件/mist_master/mist-master/main.R`

关键输出会写到：
- `2.样本信息与分组/sample_annotation_auto_R.tsv`
- `3.数据清洗与预处理/proteinGroups_已清洗.tsv`
- `4.MiST输入文件/data.txt`
- `4.MiST输入文件/keys.txt`
- `6.MiST评分结果/preprocessed.txt`
- `6.MiST评分结果/mist_run.log`

### 进入目录

```powershell
Set-Location 'C:\Users\dakey\Desktop\01_进行中\古老蛋白算法调研\多算法分析工作区\02_MiST_互作评分'
$R = 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe'
```

### 第 1 步：生成样本注释表

```powershell
& $R '.\脚本\00_生成样本注释表.R'
```

生成后人工检查：
- `sample_id` 是否对应 `proteinGroups.txt` 里的 `Intensity <sample>` 列
- `bait_name` 是否正确
- `keep` 是否只保留要参与 MiST 的样本

说明：这里识别到的是样品/run 数，不是唯一 bait 数。比如当前数据有 1612 个 `Intensity <sample>` 列，对应 1612 个样品/run；自动注释后唯一 `bait_name` 是 806 个。

### 第 2 步：清洗 proteinGroups

```powershell
& $R '.\脚本\01_数据清洗脚本.R'
```

默认会过滤：
- `Reverse`
- `Potential contaminant`
- `Only identified by site`
- 缺少蛋白 ID 的行

脚本会把空的污染标记列当作"未标记"，不会因为 `Reverse` 等列为空就误删全部数据。清洗后的超宽 `proteinGroups` 会按原始 TSV 行筛选写出，避免 8000 多列表格被 `write.table()` 写成空字段。

### 第 3 步：转换为 MiST 输入

```powershell
& $R '.\脚本\02_格式转换脚本.R'
```

默认生成：
- `4.MiST输入文件/data.txt`
- `4.MiST输入文件/keys.txt`
- `4.MiST输入文件/remove.txt`
- `4.MiST输入文件/collapse.txt`
- `4.MiST输入文件/specificity_exclusions.txt`

脚本会主动检查 `data.txt` 是否为 0 行。如果清洗后的表没有蛋白 ID、没有正强度值，或样本列匹配失败，会直接报错，不再生成空输入后让 MiST 主程序崩溃。

`Mol. weight [kDa]` 会自动从 kDa 转成 Da，因为 MiST 官方代码按 `Da / 110` 估算蛋白长度。

### 第 4 步：安装依赖

第一次运行，或更换 R 版本后执行一次：

```powershell
& $R '.\脚本\03_安装MiST依赖.R'
```

### 第 5 步：运行 MiST

普通运行：

```powershell
& $R '.\脚本\03_运行MiST.R'
```

带日志运行：

```powershell
& $R '.\脚本\03_运行MiST_日志版.R'
```

日志位置：
- `6.MiST评分结果/mist_run.log`
- `6.MiST评分结果/mist_status.tsv`

运行时可以另开一个 PowerShell 窗口查看状态：

```powershell
Set-Location 'C:\Users\dakey\Desktop\01_进行中\古老蛋白算法调研\多算法分析工作区\02_MiST_互作评分'
$R = 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe'
& $R '.\脚本\03_查看MiST状态.R'
```

连续监控，每 60 秒刷新一次：

```powershell
& $R '.\脚本\03_查看MiST状态.R' --watch=true --interval=60 --cpu_interval=5 --active_cpu_delta=1
```

状态脚本会检查：
- 是否存在 MiST 的 `Rscript` 进程
- `Rscript`/`Rterm` 的 CPU 是否还在增长
- `mist_run.log` 最后更新时间
- 当前推断阶段：读取输入、合并 keys、构建矩阵、QC、MiST 打分、完成或报错
- 最后一行会直接给出结论：正在正常运行、疑似卡住、已完成、已失败、或没有运行

默认用 5 秒 CPU 采样判断是否仍在计算。默认超过 10 分钟日志没有更新且 CPU 没有明显增长，会提示"疑似卡住"。可以改阈值：

```powershell
& $R '.\脚本\03_查看MiST状态.R' --stale_minutes=30 --tail=50 --cpu_interval=10 --active_cpu_delta=2
```

### 一键运行

```powershell
Set-Location 'C:\Users\dakey\Desktop\01_进行中\古老蛋白算法调研\多算法分析工作区\02_MiST_互作评分'
$R = 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe'

& $R '.\脚本\00_生成样本注释表.R'
& $R '.\脚本\01_数据清洗脚本.R'
& $R '.\脚本\02_格式转换脚本.R'
& $R '.\脚本\03_安装MiST依赖.R'
& $R '.\脚本\03_运行MiST_日志版.R'
```

如果样本注释表还没有人工检查，不建议直接用一键运行。

### 常用自定义命令

指定输入输出：

```powershell
& $R '.\脚本\01_数据清洗脚本.R' --input='.\1.原始输入\proteinGroups.txt' --output='.\3.数据清洗与预处理\proteinGroups_已清洗.tsv'
```

指定 MiST 配置：

```powershell
& $R '.\脚本\03_运行MiST.R' --config='.\4.MiST输入文件\mist.yml'
```

指定官方 MiST 主程序：

```powershell
& $R '.\脚本\03_运行MiST.R' --main='.\8.文档与参考\03_源码或软件\mist_master\mist-master\main.R'
```

### 卡住时诊断

只跑预处理：

```powershell
& $R '.\脚本\03_运行MiST_诊断版.R' --stage=preprocess
```

只跑 QC：

```powershell
& $R '.\脚本\03_运行MiST_诊断版.R' --stage=qc
```

只跑 MiST 打分：

```powershell
& $R '.\脚本\03_运行MiST_诊断版.R' --stage=mist
```

诊断日志在：`6.MiST评分结果/diagnostic/`

## 常见报错

- 找不到 `data.txt` 或 `keys.txt`：先运行 `02_格式转换脚本.R`
- 找不到 `mist.yml`：确认 `4.MiST输入文件/mist.yml` 存在，或用 `--config=` 指定
- 找不到 `main.R`：确认 MiST 源码已解压到 `8.文档与参考/03_源码或软件/mist_master/mist-master/`，或用 `--main=` 指定
- 样本匹配失败：检查 `sample_annotation_auto_R.tsv` 里的 `sample_id` 是否能对应 `Intensity <sample>` 列
- `输出 data.txt 行数: 0`：说明没有任何正强度 bait-prey 记录。优先检查 `3.数据清洗与预处理/proteinGroups_已清洗.tsv` 是否有蛋白 ID 和正强度
- MiST 报 `missing value where TRUE/FALSE needed`：常见原因是分子量/长度为 0 或 NA。当前格式转换脚本已自动把 `Mol. weight [kDa]` 转成 Da 并补齐无效值
- 不确定是否卡住：运行 `脚本/03_查看MiST状态.R`。如果进程存在、日志近期更新，一般是在正常跑；如果进程不存在且日志没有 `SCORING FINISHED`，多半是已经报错退出

## 输出说明

### 关键输出文件

- `4.MiST输入文件/data.txt`
- `4.MiST输入文件/keys.txt`
- `4.MiST输入文件/remove.txt`（可选）
- `4.MiST输入文件/collapse.txt`（可选）
- `4.MiST输入文件/specificity_exclusions.txt`（可选）
- `6.MiST评分结果/preprocessed.txt`
- `6.MiST评分结果/mist_run.log`

### 结果筛选

优先查看：
- `6.MiST评分结果/preprocessed.txt`
- `6.MiST评分结果/mist_run.log`

筛选建议：
- 高置信起点可先看 `MiST >= 0.75`
- 真正阈值应结合你项目里的正负集、重复数和背景复杂度再调整

详细筛选说明见：[脚本/04_结果筛选说明.md](脚本/04_结果筛选说明.md)

## 脚本清单

- `00_生成样本注释表.R` / `.py`
- `01_数据清洗脚本.R`
- `02_格式转换脚本.R`
- `03_安装MiST依赖.R`
- `03_查看MiST状态.R`
- `03_运行MiST.R`
- `03_运行MiST_日志版.R`
- `03_运行MiST_诊断版.R`
- `05_可视化脚本.R`

## 参考文档

- [脚本/03_输入文件确认说明.md](脚本/03_输入文件确认说明.md) - 输入验证标准
- [脚本/04_结果筛选说明.md](脚本/04_结果筛选说明.md) - 筛选阈值和解释
- `8.文档与参考/` - 论文、官方实现、源码

## 当前约定

- 默认输入模板来自 MaxQuant `proteinGroups.txt`
- 默认推荐阈值先用 `MiST >= 0.75`
- 如果数据集存在明确正负样本，可再用 ROC / PR 曲线复核阈值
