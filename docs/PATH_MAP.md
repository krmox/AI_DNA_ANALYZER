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

## Second regrouping (conservative; contents unchanged except where noted)

Files moved with plain `mv` (no content change) except for the lines listed under *Edits*.
Older reports, manifests and `results/` still cite the old paths.

| Old path (repository root) | New path |
|---|---|
| `bench_v15_safety.py`, `bench_v16_unseen.py`, `bench_v17_unseen.py`, `bench_v18_unseen.py`, `bench_v20_chr20.py`, `safety_layer.py`, `analyse_v15_strata.py`, `v19_timing.py`, `select_v14_regions.py`, `select_v15_regions.py`, `select_v16_region.py`, `select_v17_region.py`, `select_v19_regions.py`, `summarize_v14.py` … `summarize_v19.py`, `test_bench_v13.py` … `test_bench_v17.py`, `test_safety_layer.py` | `benchmarks/<same name>` |
| `test_binomial_baseline.py`, `test_quality_error.py`, `test_locus_evidence.py`, `test_providers.py`, `test_dataloader.py`, `test_generate_mock_bam.py`, `test_generate_mock_vcf.py`, `test_residual_binomial.py` | `tests/<same name>` |
| neural-era scripts: `train_{pb_residual,raw_pileup_llr,read_level,residual_binomial}.py`, `evaluate_{cheap_pipeline,checkpoint,llr_experiment,read_level,residual_experiment}.py`, `sanity_*.py`, `cache_*.py`, `model_raw_pileup_llr.py`, `model_read_level.py`, `read_level_controls.py`, `baseline_raw_pileup.py`, `diagnose_coverage_regime.py`, `dna_mutation_dataset.py`, `main.py`, `select_residual_lambda.py`, `summarise_label_distribution.py`, `summarize_{read_level,residual}_experiment.py`, `test_bimamba.py`, `test_read_level.py` | `deprecated/neural/<same name>` |
| `checkpoints_readlevel/`, `checkpoints_residual_binomial/`, `checkpoints_rawpileup_llr/`, `checkpoints_rawpileup_15x/` | `deprecated/checkpoints/{readlevel,residual_binomial,rawpileup_llr,rawpileup_15x}/` |

*Edits.* `benchmarks/*` and `deprecated/neural/*` scripts gained one line, `import _bootstrap`, so
`python3 <dir>/<script>.py` still finds the root modules. `benchmarks/bench_v20_chr20.py`: `parents[0]` →
`parents[1]` (repository root). `--checkpoint-dir` defaults of `train_read_level.py`,
`train_residual_binomial.py`, `train_raw_pileup_llr.py` now point to `deprecated/checkpoints/...`.
Not moved because kept code depends on them: see `deprecated/README.md` and `benchmarks/README.md`.
