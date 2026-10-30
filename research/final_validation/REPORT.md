# Frozen Varion R2 — cleanup + fresh-data validation (2026-10-28)

## Frozen model identity (unchanged)
Commit 3cdfe5a (clean tree). All 21 files match `gates/source_sha256.txt`. Original binary sha256 `6c675174…3304be` matches R2_LOCK.md and was never modified or rebuilt in place. `frozen_R2.model` sha256 `649215d5…cab6ac` matches. No threshold/parameter/source change.
Environment: micromamba env `vt` (htslib/samtools/bcftools 1.21, pysam) and `hp` (hap.py 0.3.15, py2.7) under `~/varion_mamba` (user space, no root).
Original binary runs unmodified with `LD_LIBRARY_PATH=~/varion_mamba/envs/vt/lib`.

## Cleanup
Permanently deleted (quarantine contents matched the report: 41 files, 33,699,792,585 B): 11 HG004 downsample BAMs (stage2 f053/f010 chr20/21, stage3 chr22 30x, r2 *_s301), 5 depth_robust BAMs, benchmark `withRG.bam`, benchmark `.fna.gz` (hash-identical to retained `.fna`), 5 verified-duplicate files, the incomplete chr18 `.part`. Space: free 14.05 GB -> 47.75 GB (+33.70 GB). Then 10 GB used for fresh HG007 data (35 GB free now).
Kept: frozen source/model/locks, scripts, reports, TRAIN/VALID feature dumps (`r2/dump2`), natives (HG004 chr20/21/22, HG006 chr16, HG003 chr8 300x, HG002 chr20 300x), reference + truth, `.npz` caches, `container_storage.img` (64 GB, hap.py image; no longer needed for hap.py — candidate for you to remove), depth30 and v2/v2x directories (unclear role, skipped). Full file lists: `project_cleanup.log`, `project_cleanup_candidates.tsv`.

## Fresh data (official GIAB FTP)
HG007 (NA24695), HiSeq 2x148 novoalign 100x GRCh38 (`…/NA24695_Mother_HiSeq100x/…/HG007.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.100x.bam`), streamed per contig with `fetch_hg007.sh`: chr19 (4.24 GB, sha256 7570bc44…), chr17 (6.30 GB, sha256 6925d373…); `samtools quickcheck` OK. Truth `release/ChineseTrio/HG007_NA24695_mother/NISTv4.2.1/GRCh38/HG007_GRCh38_1_22_v4.2.1_benchmark.{vcf.gz,bed}` (no `_noinconsistent.bed` is published for HG007; same `.bed` convention as HG005/HG006). Reference: GRCh38 no_alt+hs38d1 analysis set (same as all earlier runs). HG007 and chr17/chr19-of-HG007 were never used in R2 development (TRAIN: HG002 chr20, HG003 chr8, HG004 chr20/21; VALID: HG004 chr22; HG006 chr16 = already used, NOT blind any more).

## Experiments completed (SNP PASS, hap.py 0.3.15, confident BED, one contig each)
Matrix: 2 contigs × {10x,20x,30x}×seeds{101,102,103} + native(~100x) × {R2 frozen, R1.1 frozen, bcftools baseline (same flags as research/depth30_baseline)} = 24 conditions-callers per contig, 0 failures (4 early chr17 rows hit a parser bug in MY driver — empty INDEL precision — archived in `hg007/RESULTS_failed_parse_archived.jsonl`, fixed, and re-run).
Full table, per-run counts, wall/CPU/RSS: `hg007/SUMMARY.md`; raw rows incl. commands, BAM sha256, samtools -s flags: `hg007/RESULTS.jsonl`; per-run hap.py outputs: `hg007/runs/`.

| contig | depth | R2 F1 | R1.1 F1 | bcftools F1 | R2−bcftools (min over seeds) |
|---|---|---|---|---|---|
| chr17 | 10x | 0.94242 | 0.94083 | 0.90828 | +0.0341 (+0.0338) |
| chr17 | 20x | 0.98578 | 0.98305 | 0.97814 | +0.0076 (+0.0073) |
| chr17 | 30x | 0.99073 | 0.98847 | 0.98654 | +0.0042 (+0.0037) |
| chr17 | ~100x | 0.99275 | 0.99187 | 0.98897 | +0.0038 |
| chr19 | 10x | 0.94051 | 0.93978 | 0.90403 | +0.0365 (+0.0357) |
| chr19 | 20x | 0.98782 | 0.98628 | 0.98248 | +0.0053 (+0.0052) |
| chr19 | 30x | 0.99381 | 0.99287 | 0.99263 | +0.0012 (+0.0010) |
| chr19 | ~100x | 0.99603 | 0.99603 | 0.99557 | +0.0005 |

Cost (8 threads): R2 ≈ R1.1 in wall/CPU/RSS (e.g. chr17 30x: 16.9 s / 118 s CPU / 233 MB vs 16.9 s / 117 s / 234 MB); bcftools 5–6× the wall and ~6× the CPU at ~130 MB.
Reproduction check (VALID, seen data): HG004 chr22 native, R2 in the new environment: SNP F1 0.991563, TP 37,956 / FN 374 / FP 272 — equals the locked 0.99156.

## Tests
| category | verdict |
|---|---|
| Source/binary/model identity | PASSED |
| `make test-tail` (R2 tail unit tests, rebuilt out-of-tree with htslib 1.21) | PASSED: 40,186 checks, 0 failures (= R2_LOCK) |
| `make test-frame` | PASSED |
| Rebuilt binary vs original (HG007 chr19 native, same command) | PASSED: byte-identical VCF (sha256 ef6f069b…). Rebuilt binary hash differs (`b04b6984…` vs original `6c675174…`; different compiler/headers), original preserved |
| `make test` (Python kernel/oracle tests) | BLOCKED: needs `torch` (legacy v1 Python stack, `quality_error_model.py`) — not installed. Minimum action: an env with torch (large, not attempted) |
| ASan/TSan/valgrind gates | not rerun |
| Accuracy on fresh independent sample/chromosomes | COMPLETED (point estimates) |

Against the PREREG endpoints (point estimates only; the preregistered 2-Mb block-bootstrap CIs were NOT computed, so none of these is a formal pass): P1 (+0.003 at 30x vs R1.1) **not met** (chr17 +0.0023, chr19 +0.0009); P2 (non-inferior to bcftools at 10/20/30x) met in every condition by point estimate; P3 (F1≥0.99 at ≥30x) reached on both contigs by R2, by bcftools only on chr19; P4 (no low-depth regression) met (R2 ≥ R1.1 at 10x).

## Limitations
One independent sample, two chromosomes, SNP-only scoring (INDEL not assessed, R2 emits SNPs), same lab/platform/aligner family as training, seeds share one BAM (ranges are not CIs), 8-thread timings on a desktop with a background HDD copy during early chr19 runs. This does not establish general superiority; the R2 gain over R1.1 is small (+0.001–0.003) and over bcftools shrinks to ≤0.005 at ≥30x.

## Defects
Confirmed: none in the frozen caller. Environment/portability notes (not algorithm defects): Makefile hard-codes `/usr/lib64` for libhts and takes headers from pysam (pysam's bundled headers are newer than libhts 1.21 → link error `tbx_itr_next1`); workaround was `HTS_INC`/`HTS_LIB` overrides. Scripts hard-code `/mnt/archive` and `/home/mark/...`.

## Remaining blockers / next steps
1. `make test` needs torch.  2. Optional: more fresh contigs/samples (e.g. HG007 chr12, HG005/HG006 other chromosomes) and the preregistered block-bootstrap CIs.  3. You may delete `container_storage.img` (64 GB) and the large natives you no longer need.
