# HG005 Max-Feasible Stress Test — chr1:1,000,001–4,000,000 (large regional/chromosome-scale HG005 stress test)

**Status: ACCURACY RESULTS INVALIDATED BY A NEWLY-DISCOVERED EXTRACTION-PIPELINE DATA-INTEGRITY BUG.**
This is not a whole-genome validation, and (per the finding below) it is not a trustworthy
regional accuracy validation either. The frozen production model (thresholds, router,
Method C) was **not modified, retrained, or tuned** at any point in this exercise.

## Executive Summary

- **Goal**: largest scientifically defensible HG005 (NA24631) Illumina/GRCh38 regional
  benchmark of the frozen AI_DNA_ANALYZER cascade that fits a ~12 GB disk budget.
- **Dataset found and verified**: GIAB NHGRI Illumina 300x novoalign BAM for HG005 on
  GRCh38 (indexed, HTTP range-request capable), matching NIST v4.2.1 truth VCF + BED.
- **Disk-feasible scope**: up to ~220–430 Mb of contiguous chr1 sequence at a 30x-equivalent
  downsample fit the disk budget with margin.
- **Actual binding constraint turned out to be CPU throughput, not disk**: the existing
  Python/pysam feature-extraction path runs at ~16 CPU-minutes/Mb. A 220 Mb region would
  take on the order of 50+ CPU-hours; even 24 Mb, parallelized across 8 threads, took ~44
  minutes wall-clock.
- **Critical finding**: while chasing why the 24 Mb hap.py run scored a near-zero F1
  (~0.0085), independent verification against `pysam.FastaFile`/`samtools faidx` showed
  that **70–78% of loci in every 3 Mb chunk had the wrong reference base recorded** by
  `extract_bench_v12.py`/`pileup_counts.py` for this BAM/region — reproduced in both the
  8-way parallel run and a clean single-process serial re-run of the same 3 Mb window.
  This is an **upstream data-integrity bug**, not a router/threshold/Method-C problem
  (those stages are pure functions of the LLR/count arrays this bug corrupts). Because the
  frozen model is contractually off-limits for changes in this exercise, and root-causing
  the bug fully was out of scope for a no-tuning stress test, **no accuracy, genotype, or
  rescue number produced by this run can be trusted**, and none is reported as valid.
- **What IS trustworthy**: dataset provenance, disk/compute-budget calculations, the
  pilot/calibration methodology, the discovery and reproduction of the reference-base bug,
  and the hap.py/Docker/tooling setup (all of which worked correctly).

## 1. Dataset Provenance

- Sample: **HG005 / NA24631** (Chinese trio son), GIAB.
- Alignment: `HG005.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam`
  (NHGRI Illumina 300x, novoalign, GRCh38 + hs38d1 decoy minus ALTs), at
  `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/`.
  Confirmed via `curl -I`: `Accept-Ranges: bytes`, `Content-Length: 668998775517`
  (~669 GB total file). Matching `.bai` present (15.1 MB). Indexed HTTP-range regional
  extraction via `samtools view -b <URL> <region>` confirmed working — the whole file was
  never downloaded.
- Truth: `HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` (+`.tbi`) and
  `HG005_GRCh38_1_22_v4.2.1_benchmark.bed` (NOTE: HG005's high-confidence BED is named
  `_benchmark.bed`, not `_benchmark_noinconsistent.bed` as for HG002/3/4 — verified via
  directory listing before use) at
  `https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/`.
- Reference: `data/reference/chr1_full.fa` (pre-existing local Ensembl-style GRCh38 chr1,
  contig name `"1"`). BAM/VCF use UCSC-style `"chr1"`; both prefix conventions are
  normalised transparently by the existing `providers.py:_normalise_contig` — verified
  correct for FASTA/BAM/VCF contig resolution, not a new addition.

## 2. Region Selection (mechanical, independent of accuracy)

- Chosen mechanically: start at the round offset chr1:1,000,001 (clear of the chr1
  telomeric gap), contiguous, extended purely by disk-budget arithmetic — never adjusted
  after seeing any accuracy result.
- Downsample: `samtools view -s 42.09` (seed 42, keep 9%), targeting ~30x from a measured
  native depth of ~300–340x (calibrated via pilot: 1 Mb native = 197.7 MB; 10 Mb at frac
  0.09 = 187.2 MB → 18.7 MB/Mb, matching the historical 30x-tag ratio used for HG002/3/4
  in this repo, e.g. `cache/bench_v19/chr2_ordinary_30x.npz` = 44.7 MB/3 Mb).

## 3. Disk Budget

| Quantity | Value |
|---|---|
| Free disk at start of this task | 12 GB (13 GB per `df`) |
| Required safety margin | ≥4.5 GB (per instruction) |
| Measured BAM footprint | 18.7–19.1 MB/Mb (30x-downsampled) |
| Measured feature-cache (`--features none`) footprint | ~15 MB/Mb (historical) |
| Disk-feasible contiguous length (peak, BAM+cache coexisting) | ~220 Mb |
| Disk-feasible length with staged download+cache+delete-raw-BAM | up to ~430 Mb (theoretical; not needed, CPU was binding) |
| **Actual region attempted** | 24 Mb (chr1:1,000,001–25,000,000), reduced from disk-feasible 220 Mb because of the CPU-throughput finding below |
| **Region used for the final (still-invalidated) accuracy attempt** | 3 Mb (chr1:1,000,001–4,000,000), reduced further to make a clean serial re-run tractable |
| 24 Mb region BAM downloaded | 458.5 MB |
| Truth VCF slice / highconf BED slice | 1.2 MB / 204 KB |
| hap.py Docker image (`quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0`, same digest as prior benchmarks: `sha256:d63b963a6cb0...`) | 1.09 GB |
| Free disk at end | 11 GB |
| Temporary intermediate deleted this session (disposable, created by this task only) | two calibration pilot BAMs (1 Mb, 10 Mb); the 24 Mb raw region BAM's per-chunk raw sub-BAMs were never separately materialized (extraction read directly from the one 458 MB region BAM) |

No existing project data, frozen source files, or prior devlog artifacts were touched or
deleted.

## 4. Frozen Configuration (verified unchanged)

```
FROZEN_BINOMIAL_THRESHOLD = 7.0
FROZEN_PB_THRESHOLD       = 10.5
FROZEN_ROUTER_CUTOFF      = 5.411872376933351
```
(from `cascade.py`, byte-identical to the values specified for this exercise.) Two-stage
routing formula reused verbatim from `robustness_benchmark.py::route_mask`/`evaluate_arms`
(binomial+PB only — no Mamba/Stage-2, matching `--features none` and the precedent set by
`experimental/unified_happy/build_vcfs.py`). Method C genotyping formula reused verbatim
from `experimental/genotype_layer/method_c_regression.py::optimized_method_c` (fixed
EPS=0.01 binomial posterior argmax over {0/0,0/1,1/1}; GQ capped at 99; 0/0→0/1 forced-call
safety rule with GQ capped at 5). No formula, threshold, or decision rule was edited.

## 5. Compute-Throughput Finding (the actual binding constraint)

The existing `extract_bench_v12.py` (Python/pysam, per-locus binomial+PB feature
extraction) runs at **~16 CPU-minutes per Mb** single-threaded (measured: a 3 Mb chunk
took 47.6 min in a prior v19 run under 9-way contention, and 22.4 min in this session's
clean single-process run of chr1:1,000,001–4,000,000 — still ~7.5 min/500 kb sub-window).
Extrapolated to the disk-feasible 220 Mb region, this is **~59 hours single-threaded** or
roughly **~7–8 hours even with full 8-way parallelism** — not completable in this session.
8-way parallel extraction of the 24 Mb region did complete in ~44 minutes wall-clock, which
is what enabled testing at 24 Mb at all. **Disk was never the limiting factor; CPU
throughput of the existing pysam-based extractor was.**

## 6. Critical Finding: Reference-Base Corruption in Full-Chunk Extraction

### What was observed
`ai_cascade_C`/`ai_pb_only_C` hap.py runs against the 24 Mb region's HG005 truth scored:

```
SNP  TRUTH.TOTAL=28435  TRUTH.TP=203  QUERY.FP=18902  Recall=0.0071  Precision=0.0106  F1=0.0085
```

This is a catastrophic, near-chance result — inconsistent with every prior HG002/3/4
benchmark in this repo (which used the identical routing/Method-C code and scored well).

### Root-cause investigation (not a fix — investigation only, per "do not tune")
1. Spot-checked individual query VCF records against the raw truth VCF: the first several
   records (chr1:1000079, 1000112, 1000291, 1000453, ...) matched **exactly** (POS, REF,
   ALT, and hap.py's own per-record `BD=TP:TP` annotation).
2. Direct set-intersection of all query vs. truth POS values: only 254/28,149 query
   positions (0.9%) appear anywhere in the truth VCF at all — this is a property of the
   raw files, not a hap.py engine artifact.
3. For non-matching query calls, compared the VCF's `REF` column against
   `samtools faidx`/`pysam.FastaFile.fetch` on the exact same reference file at the exact
   same coordinate:
   - chr1:1,001,742 — query REF=`C`, actual genome=`G` (**wrong**)
   - chr1:1,004,113 — query REF=`A`, actual genome=`G` (**wrong**)
   - chr1:1,002,321 — query REF=`A`, actual genome=`A` (correct)
4. Isolated a 500 bp ad hoc extraction (`--region 1001500 1002000`, smaller than one
   `CHUNK_BP`=500,000 chunk) and got the **correct** reference base (G) at chr1:1,001,741.
5. Re-ran chr1:1,000,001–4,000,000 **serially, single-process, no contention**, and sampled
   2,000 loci against `pysam.FastaFile`: **74.85% mismatch rate** — i.e. the corruption is
   **not** a parallelism/race-condition artifact (it reproduces identically with zero
   concurrency). It reproduces specifically when the extractor processes full,
   `CHUNK_BP`-sized (500,000 bp) chunks starting at round genome offsets, as every real
   invocation of this script does; it did *not* reproduce in the one ad hoc sub-chunk-sized
   test. All 7 other 24 Mb parallel chunks were independently sampled and showed the same
   70–78% mismatch rate.
6. Traced the code path: `pileup_counts.py::_count_locus`/`window_counts` derive
   `reference_index` **strictly** from the `reference` string slice passed in by the
   caller, indexed by `offset = column.reference_pos - window.start` — this part of the
   code is correct by inspection. The bug therefore lies further upstream, in how
   `providers.py`/`load_counts`/`load_reads` fetch and hand off the reference slice for a
   full 500,000 bp window on **this** BAM/contig, not in the reference-index bookkeeping
   itself. Full root-cause (the exact line) was not identified — it was out of scope to
   chase further under a no-tuning, no-code-change stress test, and the investigation
   already consumed the bulk of this session's compute/time budget.

### Why this matters
- Every downstream number computed from `counts` (binomial LLR, PB LLR, alt-allele choice,
  Method C genotype, and therefore all hap.py accuracy metrics) is built on the
  `reference_index` column. With ~75% of reference bases wrong in this run, **all
  accuracy/genotype/rescue numbers produced by this session's pipeline are invalid** and
  are **not reported** below as if they were real findings.
- This is **not** evidence that the frozen thresholds/router/Method-C generalize poorly to
  HG005 — those stages never saw wrong data by their own doing; they are pure functions of
  whatever `counts`/LLR arrays they are given, and the arrays here are wrong at the input
  layer.
- This is also **not proven to be HG005-specific or chr1-specific** — it may be a latent
  defect in `pileup_counts.py`/`providers.py`'s full-chunk reference handling that simply
  never manifested on the specific BAMs/regions previously exercised (chr21, chr2, chr3,
  chr5, chr13, chr17 3 Mb windows on HG002/3/4's ~50–70x novoalign BAMs). It could affect
  any future large-region extraction with this code, on any sample.
- Per the explicit instructions for this exercise, **this bug was not fixed and the frozen
  pipeline was not touched**. It is reported here as the primary finding of the stress
  test.

## 7. Accuracy / Genotyping / Stratification / Rescue / Error Taxonomy

**Not reported.** Every metric that depends on `counts`/reference-index (detection, allele,
genotype accuracy, hap.py precision/recall/F1, stratified breakdowns, rescue recall, false
rescue rate, error taxonomy) was computed from corrupted input and is withheld rather than
presented as if valid. The routing-fraction figure below is likewise unreliable because
`binomial_llr`/`pb_llr` are themselves functions of the corrupted counts, and is reported
only as an artifact of the run, not a validated router-behavior finding:

```
PB-routed (24 Mb run, INVALID INPUT — reported for completeness only): 8,601 / 20,799,919 loci (0.041%)
```

## 8. External Callers, OLD-vs-NEW Equivalence, Statistics

- **External callers (DeepVariant/GATK/Clair3)**: **BLOCKED**. Given the reference-base
  bug invalidates the AI_DNA_ANALYZER side of any comparison, and given the remaining
  session budget, no external caller was run. This is a scope decision made after
  discovering the bug, not a resource/tooling failure — the hap.py tooling itself (Docker,
  image, invocation) worked correctly and is reusable once the extraction bug is fixed.
- **OLD-vs-NEW (pysam vs C/htslib) equivalence**: not exercised in this session. The
  previously-documented 12 Mb bit-exact equivalence result (cited in prior devlogs) is a
  separate, already-completed check on HG002/3/4 data and is unaffected by this finding,
  but it does not by itself validate the reference-index path on new BAMs/regions — indeed
  this stress test shows that path can silently fail on new inputs.
- **Confidence intervals / bootstrap**: not computed — no valid accuracy numbers exist to
  bound.

## 9. Failure Classification

| Item | Classification |
|---|---|
| Whole-genome / large (≥100 Mb) HG005 extraction | **RESOURCE LIMIT** (CPU throughput of existing extractor, not disk) |
| 24 Mb / 3 Mb hap.py accuracy result | **ENGINEERING FAILURE** (reference-base corruption in extraction, root cause not fully isolated) |
| External caller comparison | **SCOPE-LIMITED, not attempted** (would inherit the same invalid input) |
| Frozen thresholds / Method C / router code | **UNCHANGED, UNIMPLICATED** — verified via SHA-visible source inspection, not modified |

## 10. Limitations (explicit)

- No trustworthy accuracy, genotype, or rescue numbers for HG005 exist from this session.
- The extraction bug's exact root cause (which function, which line) is not identified.
- Whether the bug is HG005-specific, chr1-specific, or a general latent defect is unknown.
- Whether previously-published chr21/chr2/chr3/chr5/chr13/chr17 HG002/3/4 results are
  themselves affected by a milder/undetected version of this bug was **not** re-checked in
  this session (out of scope; those results predate this discovery and were not
  re-verified against `pysam.FastaFile` here).
- Disk-budget and compute-throughput figures in Sections 3 and 5 ARE trustworthy
  (measured directly, independent of the reference-index bug).

## 11. Claims Supported by the Data

- The GIAB HG005 NHGRI Illumina 300x BAM on GRCh38 supports indexed, HTTP-range regional
  extraction without downloading the full 669 GB file.
- A contiguous chr1 region of ~220–430 Mb would fit the measured disk budget with the
  required safety margin; the realized limiting factor for this session was CPU throughput
  of the existing Python/pysam extractor (~16 CPU-min/Mb), not disk.
- A previously-latent reference-base assignment defect in `pileup_counts.py`/`providers.py`
  (or their interaction) was discovered, reproduced in both parallel and serial execution,
  and localized to full-`CHUNK_BP`-sized extraction windows on the HG005 chr1 BAM; it was
  not present in one small sub-chunk-sized ad hoc test.
- The hap.py 0.3.15 biocontainer tooling, Docker invocation conventions, and Method
  C/cascade/router code paths used in prior benchmarks all executed successfully against
  this new sample/region once given (corrupted) input.

## 12. Claims NOT Supported by the Data

- Do **not** claim the frozen AI_DNA_ANALYZER cascade "generalizes" or "fails to
  generalize" to HG005 — no valid accuracy measurement was obtained.
- Do **not** cite the 0.0085 F1 / 0.7% recall figures as a real performance number for
  anything.
- Do **not** claim whole-genome validation of any kind.
- Do **not** claim the router/rescue behavior was validated on HG005 — the routing
  percentage reported is an artifact of corrupted LLR inputs.
- Do **not** claim this stress test proves or disproves the correctness of prior
  HG002/3/4 chr21/chr2/chr3/chr5/chr13/chr17 results; they were not re-audited here.

## 13. Reproducibility

- Git: branch `main`, HEAD `a3d5761` at start, untouched throughout (verified via
  `git status --short` — only new untracked files under `data/giab_hg005_stress/`,
  `cache/hg005_stress/`, `logs_hg005/`, `experimental/stress_test/`, plus two stray
  `.bai`/`.tbi` index files accidentally written to the repo root by an early `samtools`
  invocation — harmless, not deleted per the no-delete instruction).
- Python 3.14.3, pysam 0.24.0, numpy 2.4.4, scipy 1.17.1, samtools/bcftools 1.23.1.
- hap.py: `quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0`,
  digest `sha256:d63b963a6cb01b4830393b22369e7b91d298e4156dde353739e74e4cfa4f96d0`
  (identical to prior unified_happy benchmark).
- Downsample seed: 42, fraction 0.09 (`samtools view -s 42.09`).
- Exact commands: see `experimental/stress_test/run_hg005_pipeline.py` (merge + frozen
  cascade + Method C) and the shell invocations recorded in this document (Sections 2–4,
  6). Per-chunk extraction logs: `logs_hg005/chr1_*.log`
  (parallel) and `logs_hg005/chr1_1000000_4000000.serial.log` (clean serial re-run).
- Region/truth/BED files: `data/giab_hg005_stress/` (BAM deleted is not applicable — the
  single 458.5 MB region BAM was kept; no chunk-level raw BAMs existed to delete).
- Manifest: `HG005_MANIFEST.json` in this directory.
