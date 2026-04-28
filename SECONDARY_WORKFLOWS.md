# Secondary AP-MS workflow runner

This project now has one reproducible entry point for the secondary AP-MS
analyses shown in these folders:

- `03_HGSCore_辅助评分`
- `04_CS_Score_辅助评分`
- `05_CRAPome_污染过滤`
- `06_PPIrank_网络补充排序`

Run from the project root:

```powershell
python run_secondary_workflows.py
```

Inputs:

- `proteinGroups.txt`
- `00_统一数据准备/2.样本注释/sample_annotation_master.tsv`

The shared cleaning step removes MaxQuant rows marked as:

- `Reverse`
- `Potential contaminant`
- `Only identified by site`

It also removes rows without a stable protein ID and rows that are zero across
all kept samples in the sample annotation table.

Outputs:

- `03_HGSCore_辅助评分/6.HGSCore评分结果/hgscore_candidates.tsv`
- `03_HGSCore_辅助评分/7.可视化结果/hgscore_score_distribution.png`
- `04_CS_Score_辅助评分/6.CS_Score评分结果/cs_score_candidates.tsv`
- `04_CS_Score_辅助评分/7.可视化结果/cs_score_distribution.png`
- `05_CRAPome_污染过滤/6.CRAPome过滤结果/crapome_background_frequency.tsv`
- `05_CRAPome_污染过滤/6.CRAPome过滤结果/crapome_filtered.tsv`
- `05_CRAPome_污染过滤/7.可视化结果/background_frequency_distribution.png`
- `06_PPIrank_网络补充排序/6.PPIrank排序结果/ppirank_ranked.tsv`
- `06_PPIrank_网络补充排序/7.可视化结果/ppirank_score_distribution.png`
- `06_PPIrank_网络补充排序/7.可视化结果/top_baits_by_candidate_count.png`

Important caveat:

The CRAPome stage currently uses project-level prey background frequency because
the repository does not include an external CRAPome reference table or the full
spectral-count upload workflow. This is suitable as a local contaminant filter,
but it should be replaced or supplemented with external CRAPome reference data
for publication-grade contaminant annotation.

Useful options:

```powershell
python run_secondary_workflows.py --top-n-per-bait 50
python run_secondary_workflows.py --contaminant-frequency 0.20
python run_secondary_workflows.py --protein-groups "path\to\proteinGroups.txt" --annotation "path\to\sample_annotation_master.tsv"
```

Each run writes `secondary_workflow_validation.json` locally with cleaning
counts, sample matching counts, output row counts, and assumptions. That file is
ignored by Git because it contains local absolute paths.
