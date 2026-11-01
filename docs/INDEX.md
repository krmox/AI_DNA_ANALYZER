# Documentation index

| I want to … | Go to |
|---|---|
| understand the project, versions, results | `README.md` |
| install and build | `README.md` → *Installation*; `docs/PLATFORM_SUPPORT.md`; `environment.yml`, `Dockerfile` |
| verify a checkout | `bash verify.sh full` (see `verify.sh --help`) |
| find the SNP caller | R2/V2.x/V2 engines: `native_r2/`, `native_v2x/`, `native_v2/`; V1 Python/C: root `*.py`, `native/` |
| find the indel caller | `nindels/` (`README.md`, `REPORT.md`, `EXPERIMENTS.md`) |
| find tests | root `test_*.py` (the `verify.sh` V1 set), `tests/`, `benchmarks/test_*`, `native_*/tests`, `nindels/tests` |
| reproduce benchmarks | `README.md` → *Reproducibility*; `benchmarks/README.md`; `scripts/data/` |
| scientific methods and manuscripts | `research/` (`CLAIM_EVIDENCE_MAP.csv`, `PAPER_MANUSCRIPT_V2.md`, `latex/`, `stage4/`, `final_validation/`) |
| recorded outputs | `results/` (V1), `research/` (R2/V2.x), `nindels/results*` |
| old reports and development logs | `docs/reports/`, `docs/history/`, `VALIDATION_ROUND_2.md` |
| abandoned neural work | `deprecated/` |
| frozen experiments | `experimental/` |
| where did file X go? | `docs/PATH_MAP.md` |
| build/test status, release readiness | `docs/BUILD_AND_TEST_REPORT.md`, `docs/RELEASE_READINESS.md`, `docs/PORTABILITY_AUDIT.md` |

## Component status

| Component | Status |
|---|---|
| `native_r2/` (R2), `nindels/` (Baseline/A/B/A+B) | current; R2 frozen |
| `native_v2x/`, `native_v2/`, V1 root modules, `native/` | stable / frozen releases (`v1.0.0`, `v2.0.0`) |
| `benchmarks/`, `tests/` | frozen V1 evaluation code |
| `experimental/` | historical experiments, frozen |
| `deprecated/`, `checkpoints*/` | abandoned neural line, unused |
| R3.1 | branch `r3.1` only (not on `main`) |
