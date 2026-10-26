# Varion R2 — implementation lock (frozen 2026-10-26, before any BLIND run)

| Item | Value |
|---|---|
| Base | R1.1 commit `659df82` (unchanged, `native_v2x/`); R2 lives in untracked `native_r2/` (no commit made) |
| Source | `native_r2/{src,include,tail,Makefile}` per-file sha256: `/mnt/archive/AI_DNA_ANALYZER_r2/gates/source_sha256.txt` (21 files); concatenated sha256 `b4c93126a46abadae246abd897d9644af05676bac942b0b7b97a506c94daaabe` |
| Binary | `/mnt/archive/AI_DNA_ANALYZER_r2/build/dnav2` sha256 `6c675174854917ca3bbbfc1192c9205ac0043a696bd37d6ff2b2c08c3a3304be` (g++ 15.3.1, `-O3 -std=c++20 -ffp-contract=off -fno-fast-math -fno-associative-math -frounding-math`, same as R1.1) |
| Model | `research/stage4/frozen_R2.model` (= arm A3R, C=0.3), sha256 `649215d52e57cf00d9f69b8631d7c0d9a1848db2550fc593b7d464c1b4cab6ac`; threshold z >= -0.1 fixed by blocked CV on TRAIN |
| Command | `dnav2 run ... --mode cascade --model frozen_R2.model --r2 --threads 8` (calibration 150 blocks x 64 windows; `--ood auto`; tail strength 1; no researcher regions) |
| Training data | TRAIN dumps of HG002 chr20, HG003 chr8, HG004 chr20/21 (see PREREG_R2.md); fit script `scripts/fit_arms.py` |
| Fixed choices | gate thresholds (unused by A3R), OOD envelope eps [0.0005,0.02], MAPQ p90 [40,80], BQ p50 [15,60], depth [4,800] — set before any BLIND run |

## Gates passed on this binary
* No flags: whole VCF (incl. header) byte-identical to R1.1 (HG004 chr21 30x, 45,723 records).
* `--r2` with R1.model: VCF body byte-identical to R1.1 (only the `##varion_r2` header line is added).
* Unit tests `make test-tail`: 40,186 checks, 0 failures (quantiles, OOD envelope, pattern layer vs brute force + fuzz, eps-aware binomial exactness, model tail/gate/cap).
* Calibration and decisions independent of thread count (A3R, 1 vs 8 threads: byte-identical VCF + feature dump).
* Output independent of the truth argument (a nonexistent truth path gives the identical VCF).
* AddressSanitizer + UBSan(trap) build: same VCF, 0 reports on a 80k-window slice with the A3R model. TSan and the UBSan runtime are not installed on this machine (no root): races are excluded by design (calibration uses per-thread accumulators merged by integer sums; the model and tail parameters are read-only after load) and by the thread-count determinism test, not by a race detector.
* Strict warnings (`-Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion -Wformat=2 -Wnull-dereference -Wdouble-promotion`): 0 warnings in new code (`tail/`); inherited R1.1 code keeps its original warning level.
* `--ood force` produces the conservative Binomial caller (records the reason in the header).

## VALID result (real binary + hap.py, HG004 chr22, SNP PASS) — not blind
R2 vs R1.1 F1: 10x 0.91574/0.91486, 15x 0.96482/0.96303, 20x 0.97985/0.97698, 30x 0.98807/0.98440, 68x 0.99156/0.98749. bcftools (Stage 3): 30x 0.9854, 68x 0.9888.
