# Varion research history (reconstructed 2026-10-29, Phase 1)

Method: repo + Git log (42 commits, 2026-08-04 → 2026-10-26, HEAD `3cdfe5a`), the devlogs in `History/`, `research/*.md`, `research/stage3|stage4`, `BENCHMARK_MASTER.csv`, and the result files of the HG007/HG001 campaigns. HDD artifact existence was spot-checked, not exhaustively. Anything not confirmed is marked NOT_VERIFIED. Machine-readable version: `experiment_ledger.csv`.

## Timeline
| Period | Work | Data | Status |
|---|---|---|---|
| 2026-08-04..12 | BiMamba neural prototype, mock BAM/VCF, loss experiments, read-level and PB-residual nulls | mock + HG002 | development; mostly negative/NULL |
| 2026-08-12..19 | Frozen cheap-router → Poisson-binomial cascade, v13–v19 benchmarks | HG002 (chr21, 20, 19, 4, 1, 16, 15, 7, 14, 13, 17), HG003 chr13, HG004 chr2/3/5 | development + internal validation; internal locus evaluator, not hap.py |
| 2026-09 | Method C genotype layer; HG005 stress test (first result F1 0.0085 invalidated, extraction bug fixed); chr20 300× benchmark vs DeepVariant 1.10.0 / GATK 4.6.2.0 (containers) | HG002 chr20, HG005 chr1 3 Mb | debugging + same-sample replication |
| 2026-09-26 | V2 native engine, V2.x A+C fitted on HG002 chr20 even blocks; preregistered final eval on HG003 chr8 (F1 0.989012) | HG002, HG003 | tuning; "independent" only relative to V2.x |
| 2026-10 | Stage 3: R1.1 vs bcftools on chr22 of HG002/3/4 | trio | exposed sample, new contig |
| 2026-10-26 | R2: PREREG_R2 (TRAIN HG002 chr20, HG003 chr8, HG004 chr20/21; VALID HG004 chr22), model A3R chosen, frozen; HG006 chr16 blind run | trio, HG006 | tuning + one blind contig |
| 2026-10-28..29 | This session: cleanup, HG007 chr17/19, HG001 chr17/19 | HG007, HG001 | independent of development; now exposed |

## Sample exposure (the key table)
| Sample | Role in development/evaluation | Untouched? |
|---|---|---|
| HG002 | V1 constants, V2.x fit (chr20), R2 TRAIN (chr20), many benchmarks | NO |
| HG003 | v18, final eval chr8 (for V2.x), then **R2 TRAIN chr8** | NO |
| HG004 | v19, R2 TRAIN (chr20/21), **R2 VALID (chr22)** | NO |
| HG005 | stress test + extraction-bug debugging (code was fixed because of it); planned R2 BLIND never run | NO (debugged on) |
| HG006 | R2 blind-new chr16 run (single contig, one reduced schedule) | NO (exposed) |
| HG007 | this session's chr17/chr19 campaign (never in development) | NO under the new rule (evaluated) |
| HG001 | this session's chr17/chr19 campaign (never in development) | NO under the new rule (evaluated) |
| HG008/HG009 | cancer tumor/normal; no germline v4.2.1 GRCh38 benchmark found on the GIAB FTP (checked 2026-10-29) | n/a — unusable |

## Grouping
- **A development/debug:** v13–v19; HG005 stress test and fix; mock-data work.
- **B tuning/selection:** V2.x on HG002 chr20; R2 TRAIN and VALID; depth_robust stage 1; recall_recovery.
- **C internal validation:** chr20 300× head-to-head (same sample); HG003 chr8 final eval; Stage 3 chr22; HG004 chr22 reproduction.
- **D independent evaluation (relative to development):** HG006 chr16 (exposed), HG007 chr17/19, HG001 chr17/19.
- **E incomplete/invalid:** HG005 pre-fix result; HG005 raw artefacts lost; R2 BLIND chr18 (HG002/3/4) and HG005 never run; Strelka2/FreeBayes not run; harness-bug rows archived (HG007 empty-INDEL parse; HG001 index name, missing `QUERY.TP`, commented-out line).

## What is a new sample vs. a new chromosome vs. a rerun
HG007 and HG001 are new samples. HG003 chr8 and HG006 chr16 are new contigs. Seeds 101–103 are re-subsamples of the same BAM, not replicates.
