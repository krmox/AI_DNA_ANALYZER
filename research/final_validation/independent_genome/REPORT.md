# Frozen Varion R2 on an independent genome: HG001 (NA12878), chr17 + chr19

Date: 2026-10-28/29. Details: ENVIRONMENT.md (hashes, versions, commands), SUMMARY.md (every run), ENDPOINTS.md, BOOTSTRAP.tsv, RESULTS.jsonl, FAILED_RUNS.md, logs/, tests/.

## 1. Sample selection (HG009 was rejected)
HG009 is **not usable** for this protocol. Under GIAB `data_somatic/` it is a tumor/normal cancer cell-line pair (HG008/HG009, "Cancer Genome in a Bottle"); no germline v4.2.1 GRCh38 benchmark of the HG001–HG007 type exists for it, so SNP-PASS scoring against a germline truth set would be invalid. Substitute: **HG001/NA12878**, because (a) project records (History/18_DEVLOG.md "never touched (HG001, …)", VALIDATION_ROUND_2.md) show it was never used for R1.1/R2 development, tuning, training, validation or blind runs; (b) GIAB provides an Illumina HiSeq 2x148 300x novoalign GRCh38 BAM (same alignment recipe as HG002–HG007) and a v4.2.1 GRCh38 truth VCF + confident BED; (c) it is a different genome from the trio-related samples used so far (HG002–HG004 AJ trio; HG005–HG007 Chinese trio). Caveat: same lab, platform family and aligner as the development data.
Contigs: chr17 and chr19 (not R2 TRAIN/VALID/BLIND contigs of any sample; chr18, chr20, chr21, chr22, chr8 avoided). chr19 and chr17 were chosen for adequate confident-region size at tractable download volume (300x streams: 10.3 and 15.2 GB); chr17 is segdup-rich. Both were also evaluated on HG007 earlier (separate report). Regions: contig ∩ HG001 v4.2.1 confident BED, identical for every caller.

## 2. Methods (all frozen; nothing tuned)
Frozen R2 identity re-verified before and after (tests/identity_after_campaign.txt): commit 3cdfe5a, 21/21 source hashes, R2 binary sha256 6c675174…, model 649215d5…, R1.1 binary e83d87fd…, R1.model 63816707…. Original binary used unmodified (needs libhts.so.3 → micromamba htslib 1.21 via LD_LIBRARY_PATH). Callers on identical BAMs/BED/reference: R2 (dnav2 --r2 + frozen_R2.model), R1.1 (frozen binary+model), bcftools (stage3/depth30 baseline flags, 32 shards). Downsampling: samtools view -s SEED.FFFF by BED-mean depth (native HG001 BED-mean depth 267.3x chr19; see RESULTS.jsonl for chr17), seeds 101–103, achieved depth recorded. Scoring: hap.py 0.3.15 (xcmp), SNP Filter=PASS primary; INDELs reported separately in SUMMARY.md (Varion emits SNPs only). Uncertainty: PREREG 2-Mb block bootstrap (4,000 resamples, seed 20261027, percentile 95%), paired. Aggregation over contigs×seeds is **not** numerically specified in the PREREG: documented in analyze.py (blocks resampled per contig, one draw shared by the contig's seeds/callers, ΔF1 averaged over contig×seed cells). Seeds share one BAM and blocks are spatially correlated, so intervals reflect within-contig spatial variability only, not between-genome variability.
Matrix executed: 2 contigs × {10,20,30}x × 3 seeds × 3 callers = 54 runs, all exit 0, hap.py OK, 0 metric inconsistencies. Not run: native 300x (impractical), 15x/50x PREREG depths.

## 3. Results (SNP PASS F1, mean of 3 seeds; full per-run table in SUMMARY.md)
| contig | depth | R2 | R1.1 | bcftools | R2−R1.1 | R2−bcftools |
|---|---|---|---|---|---|---|
| chr17 | 10x | 0.94634 | 0.94539 | 0.91324 | +0.00095 | +0.03310 |
| chr17 | 20x | 0.98800 | 0.98624 | 0.98156 | +0.00175 | +0.00643 |
| chr17 | 30x | 0.99217 | 0.99054 | 0.98886 | +0.00163 | +0.00331 |
| chr19 | 10x | 0.94150 | 0.94103 | 0.90468 | +0.00047 | +0.03681 |
| chr19 | 20x | 0.98787 | 0.98697 | 0.98262 | +0.00090 | +0.00524 |
| chr19 | 30x | 0.99325 | 0.99271 | 0.99198 | +0.00054 | +0.00127 |

### Preregistered endpoints (thresholds unchanged; see ENDPOINTS.md)
* **P1** (30x, ΔF1 R2−R1.1 ≥ +0.003 and CI lower > 0): ΔF1 = +0.00109, CI [+0.00042, +0.00207] → **NOT MET**. The gain is statistically positive within these contigs but about one third of the preregistered target.
* **P2** (R2−bcftools lower bound ≥ −0.0005) at the tested depths: 10x +0.0350 [+0.0331, +0.0368], 20x +0.0058 [+0.0049, +0.0072], 30x +0.0023 [+0.0011, +0.0040] → **MET** at 10/20/30x (15x and 50x not tested, so P2 is only partially assessed).
* **P3** (R2 ≥ 0.990 at 30x): reached on both contigs (0.9922, 0.9933); bcftools also reaches it on chr19 (ceiling claim allowed) but not on chr17 (0.9889), so the target is not claimable as a data ceiling there.
* **P4** (no 10x regression vs R1.1): **MET** on both contigs (+0.0009, +0.0005); 15x not tested.
* **Guardrails vs R1.1** (indicative; runs were sequential single runs, not 5 interleaved reps): CPU within ±1.5% and peak RSS within ±1.4% everywhere; wall within +10% / 2 s except chr17 20x (+16.6%, +2.2 s, with equal CPU: a single noisy-wall comparison, not a measured slowdown).

## 4. Interpretation
* **R2 performs well on this independent sample:** F1 0.941–0.943 at 10x, 0.988 at 20x, 0.992–0.993 at 30x on both chromosomes, with zero failed runs.
* **Generalization across the tested chromosomes and depths:** consistent sign and magnitude of R2's small advantage over R1.1 and its larger advantage over bcftools on two chromosomes and three depths, and consistent with the earlier HG007 results. This is two chromosomes of one extra genome; it does not establish genome-wide or cross-platform generalization.
* **Runtime/memory:** R2 costs the same as R1.1 (CPU ±1.5%, RSS ±1.4%). Versus bcftools, R2 used ~5–6x less wall time and ~6x less CPU at ~2x the peak RSS (200–240 MB vs 100–135 MB), at 8 threads; bcftools is the 32-shard parallel baseline, not a tuned fast configuration.
* **Superiority over a baseline:** R2 is better than bcftools in every condition, with intervals excluding zero, but the margin shrinks to +0.001 (chr19) to +0.003 (chr17) at 30x; at 30x on chr19 both are near the data ceiling. R2's gain over R1.1 is real but small (+0.0005 to +0.0018) and does not reach the preregistered +0.003.
* **Unanswered:** other genomes/platforms/aligners/library types; other chromosomes (incl. chrX/Y, large chr1-5); INDELs (not called); depths below 10x or above 30x and 15x/50x endpoints; variance between genomes (the CIs are within-contig); a fully interleaved timing study; the Python test suite (`make test`, needs torch).
* No claim of universal superiority, clinical validity or broad generalization is made.

## 5. Tests and reproducibility
Unit tests rebuilt out-of-tree (original binary untouched): `make test-tail` 40,186 checks, 0 failures; `make test-frame` passed; `make test` (Python kernels) **blocked** (torch missing). Outputs in this directory: run_campaign.py, matrix.sh, fetch_hg001.sh, analyze.py, logs/. Reproduce: `fetch_hg001.sh chr19 chr17` → `matrix.sh chr19; matrix.sh chr17` → `python analyze.py` (needs the micromamba envs `vt`, `hp`).
Storage: HG001 natives (25.5 GB) deleted after completion; they are recoverable with fetch_hg001.sh and their sha256 are in logs/fetch_sha256.log. Truth files and per-run hap.py summaries/VCFs retained (/mnt/hdd/AI_DNA_ANALYZER_hg001/truth, runs/). Nothing committed.
