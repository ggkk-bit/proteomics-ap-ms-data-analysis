# HGSCore 辅助评分工作流

HGSCore 是基于超几何分布的蛋白互作评分算法，使用官方 SMAD R 包实现。

## 算法简介

HGSCore 通过超几何检验评估蛋白对的共现显著性，结合 NSAF 归一化，输出蛋白对的互作概率分数。

**参考文献**：
- Guruharsha et al., Cell 2011 - Drosophila 蛋白复合物网络
- Hart et al., BMC Bioinformatics 2007 - 超几何分布误差模型

## 目录结构

```
03_HGSCore_辅助评分/
├── 1.原始输入/              # proteinGroups.txt
├── 2.样本信息与分组/        # 样本注释（可选）
├── 3.数据清洗与预处理/      # 清洗后的数据
├── 4.HGSCore输入文件/       # HG() 函数输入格式
├── 5.质控/                  # 清洗日志
├── 6.HGSCore评分结果/       # HG 分数结果
├── 7.可视化/                # 图表
├── 8.文档与参考/            # SMAD 源码和文档
└── 脚本/                    # R 脚本
    ├── 01_数据清洗.R
    ├── 02_格式转换_HGSCore.R
    ├── 03_运行HGSCore.R
    ├── 04_结果筛选.R
    └── 05_可视化.R
```

## 使用方式

### 方式 1：使用统一准备层（推荐）

先运行统一准备层，再使用 HGSCore：

```powershell
# 1. 统一准备层
cd ../00_统一数据准备
Rscript 脚本/02_generate_sample_annotation.R
Rscript 脚本/01_clean_protein_groups.R

# 2. HGSCore 分析
cd ../03_HGSCore_辅助评分
Rscript 脚本/01_数据清洗.R
Rscript 脚本/02_格式转换_HGSCore.R
Rscript 脚本/03_运行HGSCore.R
Rscript 脚本/04_结果筛选.R --min_score=3
Rscript 脚本/05_可视化.R
```

### 方式 2：单独运行

从本地 proteinGroups.txt 开始：

```powershell
cd 03_HGSCore_辅助评分
Rscript 脚本/01_数据清洗.R --input="./1.原始输入/proteinGroups.txt"
Rscript 脚本/02_格式转换_HGSCore.R
Rscript 脚本/03_运行HGSCore.R
Rscript 脚本/04_结果筛选.R --min_score=3
Rscript 脚本/05_可视化.R
```

## 输入输出格式

### 输入格式（HG 函数）

| 列名 | 说明 | 示例 |
|------|------|------|
| idRun | 样本 ID | Sample_1 |
| idPrey | 蛋白 ID | P12345 |
| countPrey | 谱数或强度 | 1234.5 |
| lenPrey | 蛋白长度 | 350 |

### 输出格式

| 列名 | 说明 |
|------|------|
| Protein_A | 蛋白 A |
| Protein_B | 蛋白 B |
| HG | HG 分数（负对数 p 值）|

## 依赖包

```r
# 核心依赖
install.packages(c("dplyr", "tidyr", "magrittr", "devtools"))
install.packages("RcppAlgos")

# 可视化
install.packages("ggplot2")

# SMAD 包（从本地加载）
# 已包含在 8.文档与参考/03_源码或软件/SMAD-master/
```

## 参数说明

### 01_数据清洗.R
- `--unified_cleaned`: 统一准备层清洗结果路径
- `--input`: 本地 proteinGroups.txt 路径
- `--output`: 输出路径

### 02_格式转换_HGSCore.R
- `--cleaned`: 清洗后的数据路径
- `--annotation`: 样本注释路径（可选）
- `--output`: 输出路径
- `--quant_type`: 定量类型（auto/Intensity /LFQ intensity ）

### 03_运行HGSCore.R
- `--input`: HGSCore 输入文件
- `--output`: 结果输出路径
- `--smad_path`: SMAD 包路径

### 04_结果筛选.R
- `--input`: HGSCore 结果路径
- `--output`: 筛选后结果路径
- `--min_score`: 最小 HG 分数阈值（默认 3）

### 05_可视化.R
- `--input`: 筛选后结果路径
- `--output_dir`: 图表输出目录

## 注意事项

1. **SMAD 包加载**：脚本会自动从本地加载 SMAD 包，无需安装到系统
2. **蛋白长度**：必须有 `Sequence length` 列
3. **定量类型**：自动检测 Intensity/LFQ/iBAQ/MS/MS count
4. **筛选阈值**：HG >= 3 是常用阈值，可根据数据调整

## 与其他算法的关系

- **CompPASS/MiST**：主要评分算法，HGSCore 作为辅助验证
- **复合分析**：HGSCore 结果可与其他算法结果交叉验证
