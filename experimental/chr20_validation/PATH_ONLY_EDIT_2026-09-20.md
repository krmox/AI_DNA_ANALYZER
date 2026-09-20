# Path-only edits to frozen chr20 scripts (2026-09-20)

The SHA-256 values in `ANALYSIS_FREEZE.json` / `FROZEN_MANIFEST.json` describe the scripts **as run** (2026-09-20 morning). To remove machine-specific absolute paths and the dependency on the worktree `.claude/worktrees/opt-chr20`, the scripts below were edited afterwards. Only path resolution changed (absolute path literals replaced by a project root derived from the script location); no threshold, statistic, algorithm or output was touched. Pre-edit copies: `/home/mark/Documents/Projects/AI_DNA_ANALYZER_backup_pre_freeze_20260920/experimental_scripts_only/` and `/home/mark/Documents/Projects/AI_DNA_ANALYZER_backup_pre_freeze_20260920/bench_v20_chr20.py`.

| file | sha256 immediately before this edit | sha256 now |
|---|---|---|
| `bench_v20_chr20.py` | 74a0f3f249cf1a1e… | 06ef5d008725dd24… |
| `experimental/chr20_validation/run_chr20_pipeline.py` | b85d64b3fb11d8d0… | 12666af837a43c2b… |
| `experimental/chr20_validation/run_happy_chr20.sh` | 85c50e419fde2ed9… | f4306e3d2093b338… |
| `experimental/chr20_validation/run_extract_chr20.sh` | 295dbfcf79b6dd18… | 8e397f2d7ada2b03… |

Note: for `bench_v20_chr20.py` (74a0f3f249cf1a1e) and `run_chr20_pipeline.py` (b85d64b3fb11d8d0) the pre-edit hash equals the value in the freeze record. For `run_happy_chr20.sh` the freeze record lists d94143a977952e32, which differs from the pre-edit file (85c50e419fde2ed9), i.e. that script was already modified between freezing and this edit; the content of that change is not recoverable from the record and the hap.py outputs themselves (`happy/*/happy.summary.csv`) are unaffected by this path edit.

## Diff of the changed lines (old → new)
```diff
< M = "/home/mark/Documents/Projects/AI_DNA_ANALYZER"
> import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[0])  # project root (path-independent)
> M = _ROOT
< REPO = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/.claude/worktrees/opt-chr20"
> import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[2])  # project root (path-independent)
> REPO = _ROOT
< HERE=/home/mark/Documents/Projects/AI_DNA_ANALYZER/.claude/worktrees/opt-chr20/experimental/chr20_validation
< D=/home/mark/Documents/Projects/AI_DNA_ANALYZER/data/giab_hg002_chr20
> ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../.." && pwd)"  # project root (path-independent)
> HERE=$ROOT/experimental/chr20_validation
> D=$ROOT/data/giab_hg002_chr20
< M=/home/mark/Documents/Projects/AI_DNA_ANALYZER
> ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/../.." && pwd)"  # project root (path-independent)
> M=$ROOT
```
