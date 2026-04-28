# CRAPome 输入与筛选说明

## 1. CRAPome 需要什么输入
CRAPome 的上传分析通常使用长表，推荐四列：
- `Bait Name`
- `AP Name`
- `Prey Name`
- `Spectral Count`

其中：
- `Bait Name` 是 bait 名称
- `AP Name` 是亲和纯化样本名，同一 bait 的不同重复必须能区分
- `Prey Name` 是 prey 蛋白或基因名
- `Spectral Count` 是该 prey 在该 AP 中的谱数

## 2. 从 MaxQuant 到 CRAPome 的现实问题
`proteinGroups.txt` 适合做基础过滤和蛋白映射，但通常不能单独生成完整的 CRAPome 上传表。

实际转换时一般还需要：
- `evidence.txt` 或样本级谱数来源
- 样本注释表
- 必要时的蛋白 ID 映射表

如果目前只有 `proteinGroups.txt`，它更适合做：
- 反向/污染物/位点蛋白过滤
- 蛋白 ID 统一
- 和下游谱数表做 join 的锚点

## 3. 基础清洗规则
必须删除：
- `Reverse`
- `Potential contaminant`
- `Only identified by site`

建议额外检查：
- 缺失稳定 ID 的蛋白组
- 识别量过低的 run
- 同一 bait 的重复命名是否一致

## 4. 结果筛选怎么做
建议把结果分成三层：
- 高置信：本项目 bait 中多重复命中，且公共背景频率低
- 候选：分数中等，但重复性或富集趋势清楚
- 探索性：只在单次重复出现，或背景频率偏高

常用筛选维度：
- CRAPome 背景频率低
- bait 端谱数高于对照端
- 多个重复中一致出现
- 不与常见背景蛋白重叠
