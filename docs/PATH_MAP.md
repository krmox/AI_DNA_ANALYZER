# Path map (repository regrouping)

Files were moved without changing their content. Older reports, manifests and
frozen result files still cite the **old** paths; use this table to find them.
Scripts keep working when run from the repository root.

| Old path | New path |
|---|---|
| `FINAL_RESEARCH_REPORT.md`, `FINAL_RESEARCH_SUMMARY.html`, `FINAL_SCIENTIFIC_AUDIT.md`, `FINAL_BENCHMARK_RESULTS.csv`, `METHOD_C_FINAL_REVALIDATION.md`, `DEEPMIND_EXTERNAL_VALIDATION.md`, `AI_DNA_ANALYZER_Scientific_Report.pdf` | `docs/reports/<same name>` |
| `History/` | `docs/history/` |
| `download_giab_chr21.sh`, `fetch_v1*.sh`, `prepare_v1*.sh`, `run_v1*_extract.sh` | `scripts/data/<same name>` |
| `logs_v13` … `logs_v19`, `logs_chr20`, `logs_hg005` | `results/logs/v13` … `v19`, `chr20`, `hg005` |
| `HG00{2,3,4,5}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi` | `datasets/giab_truth_index/<same name>` |

Path strings inside the moved shell scripts (`logs_*` → `results/logs/*`) and in
`experimental/chr20_validation/run_extract_chr20.sh` were updated; no Python file,
report or result file content was changed. Not moved on purpose: root `*.py`,
`checkpoints*/`, `native*/`, `VALIDATION_ROUND_2.md` (see README, "Repository layout").
