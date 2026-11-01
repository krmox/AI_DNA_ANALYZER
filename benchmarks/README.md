# benchmarks/ — V1 regional benchmark rounds v15–v20 (historical, frozen)

Benchmark drivers, region selectors, summarisers and their tests for the V1 cascade rounds
**v15, v16, v17, v18, v20**, plus the shared `safety_layer`, `analyse_v15_strata` and `v19_timing`.
Recorded outputs are in `results/bench_v15 … bench_v20/` and `results/logs/` (frozen; not modified).

| Files | Purpose |
|---|---|
| `select_v14_regions.py … select_v19_regions.py` | pre-registered region selection |
| `bench_v15_safety.py`, `bench_v16/17/18_unseen.py`, `bench_v20_chr20.py` | benchmark drivers |
| `summarize_v14.py … summarize_v19.py` | summary tables from `results/` |
| `safety_layer.py`, `analyse_v15_strata.py`, `v19_timing.py` | helpers |
| `test_bench_v13.py … test_bench_v17.py`, `test_safety_layer.py` | unit tests |

**Not here:** `bench_v12_stage1/stage23`, `bench_v13_robustness`, `bench_v14_crosschrom`, `bench_v19_unseen`,
`extract_bench_v12`, `robustness_benchmark` stay at the repository root because frozen scripts under
`experimental/` and the `verify.sh` V1 tests import them there.

**Run from the repository root** (data/results paths are relative to it):
`python3 benchmarks/bench_v16_unseen.py …`, tests: `python3 -m pytest benchmarks`.
`_bootstrap.py` / `conftest.py` put the root on `sys.path`.
