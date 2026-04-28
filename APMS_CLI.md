# AP-MS reproducible CLI

This is the first maintainable command-line layer on top of the existing
algorithm work folders.

Run from the project root:

```powershell
python apms.py run --config config/project.yaml
```

The CLI creates one run folder per analysis:

```text
runs/<run_id>/
  config.yaml
  provenance.json
  report.html
  tables/
    standard_apms_table.tsv
    hgscore_candidates.tsv
    cs_score_candidates.tsv
    crapome_background_frequency.tsv
    crapome_filtered.tsv
    ppirank_ranked.tsv
  figures/
    hgscore_score_distribution.png
    cs_score_distribution.png
    background_frequency_distribution.png
    ppirank_score_distribution.png
    top_baits_by_candidate_count.png
```

The run folder is ignored by Git because it contains local data and large
outputs. The code, config, and docs are tracked.

Useful commands:

```powershell
python apms.py run --config config/project.yaml --run-id test_run
python apms.py run --config config/project.yaml --top-n-per-bait 50
python apms.py run --config config/project.yaml --contaminant-frequency 0.20
```

What this layer solves:

- one command produces a complete analysis package;
- each run preserves its config, Git commit, inputs, parameters, tables, figures,
  and HTML report;
- downstream app/front-end work can read from `runs/<run_id>/` instead of
  scraping many algorithm directories;
- future algorithms can be added by producing standard result tables and copying
  them into the run package.

The current implementation reuses the validated secondary workflow runner in
`run_secondary_workflows.py`. It does not yet replace the legacy algorithm
folders; it wraps them with a reproducible run structure.
