# Claims audit (2026-10-29)

Each reported value was checked against the underlying artifacts where they exist in this checkout.

## Reproduced from artifacts
- **HG007 table** (chr17/19 × 10/20/30× × R2/R1.1/bcftools, plus native): all 24 means reproduce from `hg007/RESULTS.jsonl` to 4 d.p. (seeds n=3; native n=1).
- **HG001 table**: 17 of 18 reproduce exactly. HG001 chr19 30× R2 mean is 0.99325 in the raw rows (displays 0.9932 or 0.9933 by rounding); no substantive difference.
- Tests: `make test-tail` 40,186 checks / 0 failures and `make test-frame` pass (re-run this campaign). Full `make test` blocked: `torch` not installed. NOT re-run in this phase.
- R2 locked identity: 21 sources / binary / model hashes verified earlier in this campaign; full hashes sit in `research/stage4/R2_LOCK.md`.

## Claims to weaken or relabel
1. **"Non-inferior to bcftools" (P2).** Only 10/20/30× were tested on HG001/HG007, not the preregistered 10/15/20/30/50×. HG007 has no CIs at all. Report P2 as *partly assessed* (HG001) and *unassessed* (HG007), not "met".
2. **"Independent" labels.** HG003 chr8 was called independent for V2.x, but it is R2 TRAIN. HG004 chr22 is R2 VALID. HG005 was used to debug. HG006 chr16 and now HG001/HG007 are exposed. No sample is untouched.
3. **HG006 chr16 blind-new.** CIs "NOT yet computed"; seeds 102–103 skipped at 10/15/20/50×. It can't count as a blind result for any future claim.
4. **Baseline comparability.** Stage 3 and HG006 used bcftools 1.23.1 (older docs mention 1.23.1); this session's HG001/HG007 runs used bcftools 1.21 (conda). Same command line, different versions. State this wherever the tables are pooled.
5. **DeepVariant comparisons.** Existing ones are HG002 only (chr20 300×; chr21 15×), via containers that are not available now (no docker/podman on this machine). The chr21 15× table mixes the locus evaluator (Varion) with hap.py (baselines), so it is not like-for-like. No DeepVariant result exists for any held-out sample.
6. **Runtime claims.** "~5–6× less wall time than bcftools" is hardware/workload specific (8 threads, 32×8 bcftools shards); chr19 timings overlapped with a download. Single runs, no repeats.
7. **Confidence intervals** are 2-Mb block bootstraps within one contig. They are not between-genome intervals; seeds are subsamples of one BAM.
8. **Effect sizes.** +0.0011 is an absolute F1 difference (0.11 points), not a 0.1% relative gain. The preregistered +0.003 was not met on HG001, HG007 or HG006.
9. **HG005 numbers** (post-fix) are record-only; raw artefacts were lost (README.md:166).
10. **Scope.** SNPs only; GIAB high-confidence regions; short-read Illumina; GRCh38; ≤2 contigs per new sample; no whole-genome run for any sample.

## Supported (with the stated scope)
- On HG001 and HG007 chr17/chr19, R2 ≥ R1.1 at all tested depths by +0.0005…+0.0027 F1, and R2 > bcftools by ≈+0.002…+0.035 depending on depth; the largest margins are at 10×.
- R2 reaches F1 ≥ 0.990 at 30× on all four sample×contig cells.
- The +0.003 R2 improvement target was *not* met anywhere it was tested.

## Missing evidence (not estimated)
- R2 BLIND chr18 (HG002/3/4) and HG005 results: no artifacts found.
- Strelka2, FreeBayes: never run.
- ASan/UBSan: recorded in R2_LOCK.md; TSan, valgrind: not available/not run.
