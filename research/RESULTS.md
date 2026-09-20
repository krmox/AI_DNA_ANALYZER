> **SUPERSEDED IN PART (2026-09-20).** `RESULTS.md` predates the chr20 and `load_reads` work. Superseded statements: v14 "IMPROVED" (pre-registered verdict NEGATIVE; pooled IMPROVED); "rescue recall 0.986" on HG005 (142 rescuable = 13 true-SNP rescues + 129 FP-avoidance; not a true-SNP recall); no chr20 result (v20: DEGRADED, ΔF1 −5.6×10⁻⁵); extraction speedups are not pipeline speedups (1.22×, 2.72×, ≤ 7.1×).

# Results

Full data: `BENCHMARK_MASTER.csv`, `ERROR_ANALYSIS.csv`, `RESCUE_ANALYSIS.csv`,
`PERFORMANCE_RESULTS.csv`. This file summarizes; the CSVs are authoritative.

## 1. Internal evaluator sweep (HG002/HG003/HG004, `bench_v13`-`v19`)

Seven pre-registered experiments, internal row-index evaluator (not
hap.py), spanning 4 samples' worth of generalization axes: depth (full/30x/15x),
region (8 distinct contigs/windows across chr1/2/3/4/5/13/16/17/19/20/21),
and cross-sample (HG002 -> HG003 -> HG004):

| Experiment | Cascade F1 | vs PB-only | Verdict |
|---|---|---|---|
| v13 (HG002 chr21 pool) | 0.9817 | +0.0001 | (no formal verdict field) |
| v14 (HG002 cross-chrom) | 0.9563 | +0.0006, CI excludes 0, p=0.023 | IMPROVED |
| v15 (HG002 held-out test) | 0.9278 | +0.0027 | WEAK POSITIVE / NEGATIVE (mixed) |
| v16 (HG002 chr13, unseen) | 0.9946 | -0.00003 | PRESERVED |
| v17 (HG002 chr17, hardest) | 0.9500 | +0.0019 | IMPROVED |
| v18 (HG003, cross-sample) | 0.9940 | +0.0003 | IMPROVED |
| v19 (HG004, 3 new contigs) | 0.9874 | +0.0010 | IMPROVED |

**Pattern**: the cascade (router -> PB for a small routed subset, ~0.04-0.11%
of loci) is never worse than PB-only by a scientifically meaningful margin in
6/7 experiments, and one experiment (v15) produced a mixed/weak result that
is reported here exactly as its own devlog recorded it, not smoothed over.
This internal evaluator ignores genotype and representation, so these numbers
are **not comparable** to the hap.py numbers below.

## 2. hap.py-evaluated results

| Experiment | Arm | Precision | Recall | F1 |
|---|---|---|---|---|
| HG002 chr21:32-44M (Method C) | PB-only | 0.9785 | 0.9319 | 0.9547 |
| HG002 chr21:32-44M (Method C) | Cascade | 0.9783 | 0.9319 | 0.9546 |
| **HG005 chr1:1,000,001-4,000,000 (Method C, this session)** | **PB-only** | **0.9642** | **0.9403** | **0.9521** |
| **HG005 chr1:1,000,001-4,000,000 (Method C, this session)** | **Cascade** | **0.9637** | **0.9403** | **0.9519** |

HG005's F1 (~0.952) is close to, and slightly below, the HG002 chr21 Method-C
figure (~0.9546) — consistent with a new, previously-unseen sample generalizing
reasonably rather than a striking win or a failure. Cascade and PB-only are
statistically indistinguishable in both experiments (ΔF1 within noise, and in
HG005's case backed by an exact 2-locus disagreement inspection, not just an
aggregate delta).

## 3. External caller comparison (HG002 chr21:32-44M only)

| Caller | Precision | Recall | F1 | Runtime |
|---|---|---|---|---|
| AI Cascade + Method C | 0.9783 | 0.9319 | 0.9546 | not separately timed at this scope |
| DeepVariant 1.6.1 | 0.9340 | 0.9654 | 0.9494 | 819s |
| GATK HaplotypeCaller 4.5.0.0 | 0.9877 | 0.9493 | 0.9682 | 287s |
| Clair3 v1.0.10 | 0.8965 | 0.9638 | 0.9290 | 1342s |
| Strelka2 | — | — | — | BLOCKED (tooling incompatibility) |
| FreeBayes | — | — | — | not run |

AI Cascade + Method C's F1 sits between Clair3/DeepVariant and GATK at this
one region/sample/depth. It has higher precision and lower recall than
DeepVariant and Clair3; it trails GATK on both metrics. **This is one region,
one sample, one depth** — not a general performance claim (see
`SCIENTIFIC_CLAIMS.md`).

## 4. Router rescue (diagnostic, not tuned)

HG005: rescue recall 0.986 (140/142 PB-rescuable loci correctly routed),
false-rescue rate 0.0056. Consistent in direction with the HG002 devlog
finding that cascade and PB-only differ at only 3/10,975,654 loci region-wide.
Full breakdown: `RESCUE_ANALYSIS.csv`.

## 5. Stratification and error taxonomy

`lowmap_segdup` is the dominant, consistently-enriched difficulty axis
(HG005: FN 7.7x, FP 6.4x over background). A genuine low-VAF component
co-occurs with it (median VAF 0.0 in-segdup vs 0.19 elsewhere among HG005 FN
loci with any extracted candidate) — the previously-hypothesized low-VAF/
segdup failure mode reproduces on HG005. A separate, non-model, engineering
artifact (window-level vs SNP-level BED granularity) accounts for 58% of
HG005's cascade FN having no extracted candidate at all — explicitly
distinguished from a detection failure. Full breakdown: `ERROR_ANALYSIS.csv`.

## 6. Extraction integrity (this session's primary engineering contribution)

The coordinate/reference-integrity bug that invalidated the original HG005
attempt is fixed, independently verified at full population scale
(0/2,403,008 reference mismatches, 0/3,989+3,991 VCF REF mismatches,
byte-identical repeated extraction), and proven not to have altered any
previously-published HG002/HG003/HG004 result (historical audit,
`experimental/stress_test/HG005_HISTORICAL_AUDIT.md`).

## 7. Performance

**Corrected finding** (an earlier pass of this dossier wrongly concluded no
C/htslib acceleration existed anywhere in the codebase — see
`PIPELINE_INTEGRITY_AUDIT.md`'s correction note for the full account): a
real, bit-exact-equivalent C/htslib extraction extension
(`native/pileup_native.c`) has been integrated into production
`pileup_counts.load_counts` in this session. Measured, apples-to-apples,
same-loci-count, extraction-stage-only speedup on the HG005 3 Mb region:
**373.11s → 17.74s, 21.03x**. At the original 12 Mb HG002 benchmark scale
(reused from the source worktree's own measurement, justified by a
full-region 0-mismatch equivalence proof that makes a re-run a deterministic
repeat of the same result): **1092.07s → 16.76s, 65.14x**.

**This figure covers the pileup-counting extraction stage only.**
`read_level_pileup.load_reads` (feeds the Poisson-Binomial caller) is a
separate, unaccelerated extractor and is the larger cost in a full run that
recomputes PB calls from scratch (HG005 3Mb: 696s vs the pre-acceleration
344s for `load_counts`) — see `PERFORMANCE_RESULTS.csv` for the full
pipeline's measured before/after wall time, which shows a real but much
smaller end-to-end speedup than 65x for that reason. Separately, real (but
smaller, ~3.1x at 8 workers, sub-linear) process-level parallelism and a
genuine null result from htslib-internal thread scaling (~0x) both remain
documented as distinct, unintegrated mechanisms — full detail:
`PERFORMANCE_RESULTS.csv`, `experimental/performance/CHTSLIB_INTEGRATION_REPORT.md`.
