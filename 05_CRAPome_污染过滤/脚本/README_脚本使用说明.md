# CRAPome 脚本使用说明

这份说明对应 `05_CRAPome_污染过滤\脚本` 目录下的脚本实际使用顺序。

当前这套工作区里，`03_格式转换_CRAPome输入.py` 默认建议直接使用 `02_MiST_互作评分` 已生成的样本级矩阵 `preprocessed_MAT.txt`。  
不要默认使用 `..\3.数据清洗与预处理\sample_counts_wide.tsv`，因为这个文件当前并没有自动生成。

## 0. 进入脚本目录

```powershell
Set-Location "C:\Users\dakey\Desktop\01_进行中\古老蛋白算法调研\下载的文献进行测试\多算法分析工作区\05_CRAPome_污染过滤\脚本"
```

## 1. 安装依赖

```powershell
python -m pip install pandas matplotlib seaborn
```

说明：

- `01`、`02`、`03` 脚本需要 `pandas`
- `05_可视化.py` 还需要 `matplotlib` 和 `seaborn`

## 2. 需要哪些输入文件

当前流程实际使用这些文件：

- `..\1.原始输入\proteinGroups.txt`
- `..\2.样本信息与分组\sample_annotation.tsv` 或者直接使用 MiST 那边的真实样本注释
- `..\..\02_MiST_互作评分\6.MiST评分结果\preprocessed_MAT.txt`
- `..\..\02_MiST_互作评分\2.样本信息与分组\sample_annotation_auto_R.tsv`

注意：

- `..\2.样本信息与分组\sample_annotation.tsv` 如果还是模板示例，就不要拿它做正式转换
- `..\3.数据清洗与预处理\proteinGroups_cleaned.txt` 只是清洗后的蛋白表，不等于样本级谱数矩阵
- `..\3.数据清洗与预处理\sample_counts_wide.tsv` 只有在你自己额外准备了样本级宽表时才存在

如果你只是想先复制一个本地模板文件名出来，可以运行：

```powershell
Copy-Item ..\2.样本信息与分组\sample_annotation_template.tsv ..\2.样本信息与分组\sample_annotation.tsv
```

## 3. 第一步：检查输入文件

这个脚本只做检查，不改数据。

如果你已经有正式样本注释：

```powershell
python .\01_输入文件确认与说明.py `
  --protein-groups ..\1.原始输入\proteinGroups.txt `
  --sample-annotation ..\2.样本信息与分组\sample_annotation.tsv `
  --output-report ..\8.文档与参考\01_CRAPome输入检查报告.txt
```

如果你暂时还没有正式样本注释，可以先不带 `--sample-annotation`：

```powershell
python .\01_输入文件确认与说明.py `
  --protein-groups ..\1.原始输入\proteinGroups.txt `
  --output-report ..\8.文档与参考\01_CRAPome输入检查报告.txt
```

## 4. 第二步：清洗 proteinGroups.txt

默认会去掉：

- `Reverse`
- `Potential contaminant`
- `Only identified by site`
- 空的 `Gene names`

命令：

```powershell
python .\02_数据清洗.py `
  --input ..\1.原始输入\proteinGroups.txt `
  --output ..\3.数据清洗与预处理\proteinGroups_cleaned.txt `
  --min-unique-peptides 1
```

如果你不想加 `Unique peptides` 阈值：

```powershell
python .\02_数据清洗.py `
  --input ..\1.原始输入\proteinGroups.txt `
  --output ..\3.数据清洗与预处理\proteinGroups_cleaned.txt
```

## 5. 第三步：转换成 CRAPome 上传四列表

这是当前工作区最应该执行的命令：

```powershell
python .\03_格式转换_CRAPome输入.py `
  --counts-wide ..\..\02_MiST_互作评分\6.MiST评分结果\preprocessed_MAT.txt `
  --sample-annotation ..\..\02_MiST_互作评分\2.样本信息与分组\sample_annotation_auto_R.tsv `
  --output ..\4.CRAPome输入文件\crapome_input.tsv
```

这一步会生成 CRAPome 要求的四列表：

- `Bait Name`
- `AP Name`
- `Prey Name`
- `Spectral Count`

输出文件：

- `..\4.CRAPome输入文件\crapome_input.tsv`

检查输出文件：

```powershell
Get-Item ..\4.CRAPome输入文件\crapome_input.tsv
Get-Content -TotalCount 5 ..\4.CRAPome输入文件\crapome_input.tsv
```

注意：

- PowerShell 里不要把 Markdown 代码块标记 ``` 一起复制执行
- 当前这份说明不再提供 `sample_counts_wide.tsv` 那条旧命令，原因是你当前工作区并没有这份文件

## 6. 第四步：在 CRAPome 网站上跑 Workflow 3

根据 CRAPome 官方流程，在线上传分析属于 `Workflow 3`，这一步需要注册用户登录后才能使用。  
官方页面：

- `https://www.reprint-apms.org/?q=chooseworkflow`
- `https://www.reprint-apms.org/?q=tutorial`

操作顺序：

1. 打开 `https://www.reprint-apms.org/?q=chooseworkflow`
2. 注册账号并登录，如果你还没有账号，先完成注册
3. 进入 `Workflow 3`，这是用于上传你自己的 AP-MS 数据进行在线分析的入口
4. 在本地准备好刚生成的文件：`..\4.CRAPome输入文件\crapome_input.tsv`
5. 在网站上传这个 `crapome_input.tsv`
6. 按网页提示选择或确认对照集、背景库和分析参数
7. 提交任务，等待网站返回分析结果页面
8. 下载结果表，保存到本地目录：`..\6.CRAPome过滤结果\`

建议：

- 下载后的结果文件统一命名为 `crapome_results.tsv`
- 如果网站返回多个表，主结果表放在 `..\6.CRAPome过滤结果\`，其他辅助文件也放同一目录

## 7. 第五步：导出本地筛选建议

这个脚本不会自动过滤，只是把推荐阈值写出来，方便你和网站结果一起看。

```powershell
python .\04_结果筛选说明.py > ..\8.文档与参考\02_CRAPome结果筛选建议.txt
```

输出文件：

- `..\8.文档与参考\02_CRAPome结果筛选建议.txt`

## 8. 第六步：对 CRAPome 结果做可视化

这一步要在你已经从网站下载结果之后再运行。

假设你把结果保存为：

- `..\6.CRAPome过滤结果\crapome_results.tsv`

运行：

```powershell
python .\05_可视化.py `
  --input ..\6.CRAPome过滤结果\crapome_results.tsv `
  --output-dir ..\7.可视化
```

输出目录：

- `..\7.可视化`

如果提示缺列，说明下载回来的结果表列名和脚本预期不一致，需要先统一列名或调整脚本。

## 9. 当前推荐的完整执行顺序

按当前工作区实际情况，推荐顺序如下：

```powershell
Set-Location "C:\Users\dakey\Desktop\01_进行中\古老蛋白算法调研\下载的文献进行测试\多算法分析工作区\05_CRAPome_污染过滤\脚本"
python -m pip install pandas matplotlib seaborn
python .\01_输入文件确认与说明.py --protein-groups ..\1.原始输入\proteinGroups.txt --output-report ..\8.文档与参考\01_CRAPome输入检查报告.txt
python .\02_数据清洗.py --input ..\1.原始输入\proteinGroups.txt --output ..\3.数据清洗与预处理\proteinGroups_cleaned.txt --min-unique-peptides 1
python .\03_格式转换_CRAPome输入.py --counts-wide ..\..\02_MiST_互作评分\6.MiST评分结果\preprocessed_MAT.txt --sample-annotation ..\..\02_MiST_互作评分\2.样本信息与分组\sample_annotation_auto_R.tsv --output ..\4.CRAPome输入文件\crapome_input.tsv
python .\04_结果筛选说明.py > ..\8.文档与参考\02_CRAPome结果筛选建议.txt
```

然后：

1. 登录 CRAPome 网站
2. 在 Workflow 3 上传 `..\4.CRAPome输入文件\crapome_input.tsv`
3. 下载结果到 `..\6.CRAPome过滤结果\crapome_results.tsv`
4. 再运行 `05_可视化.py`

## 10. 这套脚本当前不会自动完成的事情

- 不会自动从 `proteinGroups.txt` 还原出样本级谱数宽表
- 不会自动登录、上传、运行或下载 CRAPome 网站结果
- `04_结果筛选说明.py` 只输出推荐规则，不会自动执行过滤

如果后面需要，我可以继续把 `05_CRAPome_污染过滤\README.md` 也同步改成这一版说明。
