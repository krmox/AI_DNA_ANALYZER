# Does the frozen binomial → router → PB cascade hold on a representative unseen region? — chr13:70–73 Mb benchmark

**Date:** 2026-08-17

**Status of sections:** §1–§10 are the pre-registration. They were written and
saved **before any caller, router or metric was run on any chr13 locus**, and
before any chr13 BAM slice had finished downloading. Everything from §11 onwards
was written after the corresponding measurement.

Devlogs 3–15 are not edited by this experiment. `cascade.py`, `cheap_router.py`,
`binomial_baseline.py`, `quality_error_model.py`, `robustness_benchmark.py`,
`extract_bench_v12.py`, `pileup_counts.py` and `read_level_pileup.py` are not
edited by this experiment; §17 verifies that mechanically. Nothing is committed
or pushed.

---

## 1. Research question

The production architecture is frozen:

```
ALL LOCI → binomial LLR (ε = 0.01)
             ├── |binomial_LLR − 7.0| >  5.411872376933351  → cheap answer
             └── |binomial_LLR − 7.0| ≤  5.411872376933351  → PB → final call
```

Devlogs 11–13 validated it on chr21 (19.4 M loci, ΔF1 = 0). Devlog 14 validated
it across chr20/chr19/chr4/chr1 (10.37 M loci) and found exactly one failure
mode, in segmental duplication. Devlog 15 tested a candidate ε-sensitivity
safety layer on chr16/chr15/chr7/chr14; that layer is **not** part of this
benchmark, and neither is the PB → Mamba branch.

This experiment is not an attempt to improve anything. It asks one practical
question on one region no experiment in this project has ever touched:

1. Does the frozen cascade reproduce PB-only accuracy there?
2. How much PB computation does it remove?
3. What is the real wall-clock speedup?
4. Where and why does it disagree with PB?
5. Are the disagreements concentrated in particular depth / VAF / base-quality /
   mapping-quality / genomic-difficulty regimes?
6. Does the router beat matched-budget random and reverse-confidence routing?
7. What are the operating boundaries on this region?

A deliberate difference from devlogs 14 and 15: those chose regions at the
*extremes* of the annotation (maximum GC, minimum GC, maximum segdup, or an
explicitly easy neutral window). This one chooses a **representative** region —
median-difficulty contig, contig-typical difficulty window — because the
practical question "what does this architecture do on ordinary unseen sequence"
has never been asked with a rule that selects for typicality.

## 2. Hypotheses

* **H1 (accuracy transfers).** ΔF1 (cascade − PB-only) is indistinguishable
  from 0 at every depth and pooled, by the §7 decision rule. Rationale: the
  router's input is a binomial LLR over read counts and a fixed ε; it carries no
  chromosome-specific information, and devlogs 11–15 found the only departures
  in segmental duplication, which is 1.9% of this region's high-confidence
  bases.
* **H2 (compute).** The routed fraction lands in the range devlogs 13–14
  observed on ordinary-structure regions, roughly 3×10⁻⁵ – 3×10⁻³ of loci, and
  rises as depth falls (fewer reads → more loci in the uncertain band). PB
  compute fraction is therefore ≲ 0.5% and the projected end-to-end speedup is
  ≳ 100×.
* **H3 (disagreements are rare and mostly favourable).** Disagreement rate
  ≲ 10⁻⁵. Devlog 14 found departures net-favourable (26 avoided PB false
  positives vs 5 lost true SNP). Prediction: on a region with 1.9% segdup, the
  known failure mode (fixed-ε binomial confidently wrong at low VAF in
  duplicated sequence) is present but rare, so the absolute count of lost true
  SNP is 0–2 and the pooled ΔF1 sign is not predictable in advance. The sign is
  explicitly **not** part of any hypothesis: only |ΔF1| matters.
* **H4 (controls).** At matched PB budget, random routing and
  reverse-confidence routing both do worse than the frozen router — or the
  routed set is too small on this region for any of the three to differ, in
  which case the control is reported as underpowered rather than as support.

H4's second branch is written down now because devlog 14's routed sets were
50–1,300 loci per cell and its controls were frequently uninformative. If that
repeats here, it is a limitation, not a finding.

## 3. Frozen configuration — not tuned, not refit, not swept for production

Identical to devlog 13 §3.1, devlog 14 §3.1 and devlog 15 §3. Nothing below is
re-derived, re-fit or re-selected in this session.

| Constant | Value | Source |
|---|---|---|
| `FROZEN_BINOMIAL_THRESHOLD` | 7.0 | imported from `cascade.py` |
| `FROZEN_PB_THRESHOLD` | 10.5 | imported from `cascade.py` |
| `FROZEN_ROUTER_CUTOFF` | 5.411872376933351 | imported from `cascade.py` |
| binomial ε | 0.01, `error_floor` 1e-3 | `BinomialVariantCaller` defaults |
| PB ε floor / ceiling | 1e-4 / 0.25 | `quality_error_model` |
| `MAX_READS` | 48 | `read_level_pileup` |
| PB implementation | the repaired one (devlog 13 §9) | `quality_error_model.poisson_binomial_llr` |
| extraction chunk | 500,000 bp | `extract_bench_v12.CHUNK_BP` |
| window length | 64 | `SEQ_LEN` |

Baseline integrity, checked at the start of this session, at `HEAD = 1dfa759`,
tree clean:

```
git status --porcelain            -> (empty)
cascade.py              sha256 b0ee9f4b24fe06dc…  IDENTICAL-TO-HEAD
cheap_router.py         sha256 5e13cca000f4565a…  IDENTICAL-TO-HEAD
binomial_baseline.py    sha256 f0408ecae0929aef…  IDENTICAL-TO-HEAD
quality_error_model.py  sha256 a69cb49591cbac13…  IDENTICAL-TO-HEAD
robustness_benchmark.py sha256 6fd5c4df39806e05…  IDENTICAL-TO-HEAD
bench_v13_robustness.py sha256 839f105a45c3e924…  IDENTICAL-TO-HEAD
bench_v12_stage1.py     sha256 745d7d694653ebd9…  IDENTICAL-TO-HEAD
extract_bench_v12.py    sha256 8ec91cbd646cfdeb…  IDENTICAL-TO-HEAD
pileup_counts.py        sha256 70b141854cb45fa3…  IDENTICAL-TO-HEAD
read_level_pileup.py    sha256 2a63c448eb271538…  IDENTICAL-TO-HEAD
```

The routing rule is used verbatim through `robustness_benchmark.route_mask` with
`FROZEN_ROUTER_CUTOFF`; the three constants are only ever read.

New code written for this experiment is benchmark infrastructure only:
`select_v16_region.py`, `fetch_v16_reference.sh`, `fetch_v16_data.sh`,
`prepare_v16_depths.sh`, `run_v16_extract.sh`, `bench_v16_unseen.py`,
`summarize_v16.py`, `test_bench_v16.py`. None of it is imported by production
code paths.

## 4. New-region selection and provenance

### 4.1 The rules, fixed before execution

`select_v16_region.py` reads only the GIAB high-confidence BED, the GIAB v3.1
stratification BEDs and the reference FASTA. It cannot see a caller score, a
router decision or a metric.

* **Rule C16 (contig).** Among GRCh38 autosomes untouched by any previous
  devlog — excluding `chr21` (devlogs 3–13), `chr20`, `chr19`, `chr4`, `chr1`
  (devlog 14) and `chr16`, `chr15`, `chr7`, `chr14` (devlog 15) — rank the 13
  survivors by `segdup_hc_fraction = |HC ∩ lowmap_segdup| / |HC|`, descending,
  and take the **median rank** (index 6 of 13). Ties by smaller contig number.
* **Rule R16 (window).** Among 3 Mb windows at whole-Mb offsets with
  `hc_fraction ≥ 0.50` and `n_fraction == 0`, take the window whose
  `alldifficult_hc_fraction` is closest to the **contig-wide**
  `alldifficult_hc_fraction`. Ties by smaller coordinate.

Rule C16 was run to completion **before the chr13 FASTA existed**; its inputs
are two BED files. Rule R16 ran on the FASTA immediately after, still before any
BAM slice or caller. 3 Mb is fixed for statistical power, from devlog 14's
measured ~1,000 SNP/Mb and ΔF1 CI half-width ~0.0025 per 1 Mb cell → ~0.0015
expected here, inside the §7 power criterion of ≤ 0.01. It is not chosen from
any property of chr13.

### 4.2 Rule C16 output

Descending `segdup_hc_fraction` over the 13 candidates:

```
chr17 .08609  chr22 .08265  chr9 .07121  chr10 .06456  chr11 .06367
chr2  .05125  chr13 .04867  <- median rank 6/13, SELECTED
chr5  .04842  chr6 .04667  chr18 .04321  chr12 .04309  chr3 .04151  chr8 .03692
```

chr13: length 114,364,328; HC bases 92,390,841 (80.8%); segdup∩HC 4.87%;
alldifficult∩HC 16.92%.

### 4.3 The selected region

| Field | Value |
|---|---|
| Region (GRCh38, 0-based half-open) | **chr13:70,000,000–73,000,000** |
| Size | 3,000,000 bp |
| Sample | GIAB HG002 (NA24385 son) |
| Reference build | GRCh38, Ensembl release-110 chromosome FASTA |
| Platform | Illumina 2×250, novoalign (`HG002.GRCh38.2x250.bam`) |
| Truth | GIAB v4.2.1 benchmark VCF |
| GC | 0.3525 |
| N fraction | 0.0 |
| HC coverage | 0.9783 |
| alldifficult (bases / ∩HC) | 0.1781 / 0.1707 (contig-wide ∩HC: 0.1692) |
| lowmap_segdup (bases / ∩HC) | 0.0208 / 0.0190 |
| tandemrepeats (bases / ∩HC) | 0.0553 / 0.0486 |
| Expected SNP (from devlog 14's ~1,000/Mb) | ~3,000 per depth cell |

Difficulty gap to the contig-wide value: 0.0015 — i.e. this window's difficulty
is within 0.9% relative of typical chr13. That is the whole point of Rule R16.

### 4.4 Depth realisations

Native (~full), 30x and 15x, by `samtools view -s 42.FRACTION` with FRACTION
derived from the region's **measured** native mean depth. These are downsamples
of one library: three depth regimes of one region, never three independent
regions. 3 cells total.

### 4.5 Provenance

* Reads: `samtools view -b <GIAB URL> chr13:70000001-73000000` — same URL and
  recipe as devlogs 12–15 (`fetch_v16_data.sh`).
* Truth: `bcftools view -r chr13:70000001-73000000` on the v4.2.1 VCF.
* HC BED: `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`, filtered to
  the span.
* Reference: Ensembl release-110 `Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz`
  (`fetch_v16_reference.sh`), same source/release as every existing
  `data/reference/chr*_full.fa`.
* Difficulty annotation: GIAB genome-stratifications v3.1 GRCh38 —
  `lowmap_segdup`, `alldifficult`, `tandemrepeats`, the same three files
  devlogs 14–15 used.
* Retrieved 2026-08-17.

Nothing is synthesised; no artificial difficulty is injected; no read is
modified. Cross-sample and cross-platform transfer remain untested and no claim
will be made about them.

## 5. Independence audit

Every region this project has ever used lies on chr21, chr20, chr19, chr4, chr1,
chr16, chr15, chr7 or chr14. chr13 is none of them, so coordinate overlap with
any previous training, validation, tuning or test region is **impossible**, not
merely unlikely. Rule C16 could not have selected a used contig by construction.

`bench_v16_unseen.verify_caches` asserts mechanically, per cell, before any
metric is computed: the cached contig and coordinates equal the pre-registered
ones; the contig is not in `PREVIOUSLY_USED_CONTIGS`; the locus count is
window-aligned and equals an independent reconstruction of the window tiling;
no LLR is non-finite; and every reconstructed genomic coordinate lies inside
[70,000,000, 73,000,000). A failure aborts the benchmark rather than
downgrading it.

## 6. Experimental arms and metrics

Three arms on the identical scoring frame (`label ∈ {0, SNP}`; indel and no-call
loci excluded, as in every devlog since 9), identical loci and identical truth:

* **A. Binomial-only** — `binomial_llr ≥ 7.0`.
* **B. PB-only** — `pb_llr ≥ 10.5`. The reference arm.
* **C. Frozen router → PB** — PB where `|binomial_llr − 7.0| ≤ 5.411872376933351`,
  binomial elsewhere. **The architecture under evaluation.**

Primary comparison **C vs B**.

Recorded per arm and per cell: F1, precision, recall, TP, FP, FN, TN, accuracy,
balanced accuracy; ΔF1 vs PB; absolute loci scored and absolute loci routed to
PB (not only fractions); PB compute fraction; measured throughput and projected
cascade vs PB-only wall-clock and speedup; disagreement count and rate vs PB.
Disagreement count is an engineering statement about *reproducing PB*; it is
never used as an accuracy proxy, and vice versa.

**Primary endpoint: ΔF1 = F1(C) − F1(B)**, per cell and pooled, with a paired
bootstrap 95% CI (`evaluate_quality_error.paired_bootstrap`, 10,000 resamples,
module-default seed) — the same estimator, seed and code path as devlogs 11–15.

Stratified (§13), reusing existing helpers only, no new annotations invented:
depth (coarse and fine bins), VAF, mean base quality, mean mapping quality,
GIAB v3.1 `lowmap_segdup` / `alldifficult` / `tandemrepeats` inside vs outside,
and the existing Stage-1 `strata_masks`.

Controls (§14), `robustness_benchmark.controls` unchanged: matched-budget random
routing (5 seeds) and reverse-confidence routing.

Operating curve (§15): `cutoff_sweep` over the existing cutoff grid, diagnostic
only. It costs no new extraction — all LLRs are already cached — so it is cheap.
**No row of it may be promoted to a production threshold.**

## 7. Pre-registered decision rule

Quoted verbatim from devlog 13 §3.4 / devlog 14 §3.4 so it cannot drift. Per
cell and for each pooled result:

* **PRESERVED** — the paired-bootstrap 95% CI for ΔF1 contains 0 **and**
  |ΔF1| < 0.001.
* **DEGRADED** — the CI excludes 0 with ΔF1 < 0, or |ΔF1| ≥ 0.001.
* **IMPROVED** — the CI excludes 0 with ΔF1 > 0.
* **UNDERPOWERED** — < 100 SNP, or CI half-width > 0.01. Reported as
  underpowered, never rounded to "preserved".

Final verdict vocabulary, mapped now:

* **STRONG POSITIVE** — every powered cell and the pooled result PRESERVED (or
  IMPROVED), PB compute fraction ≤ 1%, projected speedup ≥ 50×, and no new
  failure mode.
* **POSITIVE** — the above but with a material limitation: some cells
  UNDERPOWERED, or the controls uninformative, or the compute reduction between
  10× and 50×.
* **MIXED / PARTIAL GENERALIZATION** — pooled result PRESERVED/IMPROVED but at
  least one powered stratum DEGRADED, or a new failure regime appears.
* **NEGATIVE** — any powered *cell* DEGRADED, or the projected speedup < 10×.
* **INCONCLUSIVE** — too few powered cells to distinguish the above.

Compute thresholds are stated as absolute numbers rather than "within an order
of magnitude of chr21" so that they cannot be re-anchored after the fact.

### 7.1 What would falsify H1

Any powered cell whose ΔF1 CI excludes 0 on the negative side; or a disagreement
rate above 10⁻⁴ (chr21 saw ~6×10⁻⁷; devlog 14's ordinary regions ~10⁻⁶).

## 8. Stopping criteria

Extraction runs **once** per cell. Metrics are computed **once**. No region,
depth, threshold, cutoff or evaluation criterion is added, removed or changed
after any result is read, except to fix an outright software crash — any such
re-run is disclosed in §17. If the frozen cutoff fails here, the failure is
reported and the cutoff stays at 5.411872376933351.

The experiment stops when: 3 cells are extracted and scored; all three arms and
all metrics of §6 are computed; controls are run or documented as underpowered;
the disagreement analysis is complete; compute accounting is verified; the test
suite passes; and §17 records reproducibility. No further cycle is run
regardless of outcome.

## 9. Runtime plan

Reused without regeneration: the GIAB v3.1 stratification BEDs, the HC BED, the
v4.2.1 truth VCF index, the extractor, every metric helper, the bootstrap, the
controls and the sweep. Regenerated: nothing that already exists.

New compute is exactly: one reference FASTA download (done, 82 s), one 3 Mb BAM
slice download, two `samtools view -s` downsamples, three extractions (the only
expensive step — each computes the binomial and the repaired PB over every
locus, which is also what supplies the *measured* PB cost the speedup claim
rests on), and one scoring pass over cached arrays. Expected wall-clock: well
under two hours, dominated by extraction. No agent fan-out; three concurrent
extraction processes at most, on 8 cores.

## 10. What this experiment cannot show

Fixed by design, stated before results exist: one sample (HG002), one platform
(Illumina 2×250, novoalign), one aligner, one truth set (v4.2.1), one reference
build, SNP only (indels are excluded from the scoring frame), and one 3 Mb
window on one chromosome. A positive result licenses no claim about other
samples, platforms, variant classes or the whole genome. §18 restates this
against whatever is actually found.

---

*(§11 onward is written after the corresponding measurement.)*

## 11. Data extraction and QC

| Step | Result |
|---|---|
| Reference FASTA | `data/reference/chr13_full.fa`, 114,364,328 bp, contig name `13`, 82 s |
| BAM slice | `chr13_representative_full.bam`, 119.4 MB, **851,891** mapped reads in region, index OK |
| Native mean depth (measured) | **69.63×** |
| 30x realisation | `samtools view -s 42.4309` → achieved **29.93×** |
| 15x realisation | `samtools view -s 42.2154` → achieved **14.99×** |
| Truth VCF slice | 5,743 records, **4,889 SNV**, 0 records outside the region |
| REF allele consistency | **0 mismatches** between VCF `REF` and the FASTA at all 5,743 records |
| High-confidence BED | 391 intervals, 2,934,993 bp in region = **97.83%** (equals the selection statistic) |
| Reference N fraction / GC | 0.0 / 0.35249 (equals the selection statistic) |
| Smoke test | 64 kb extraction, 61,760 loci, 183 SNP, 47 s, exit 0 — run before the full extraction |
| Extraction | 3 cells × 2,910,592 loci × 4,758 SNP, identical locus counts; 15x 15.7 min, 30x 23.6 min, full 39.7 min (three processes concurrently) |
| Partial-output guard | extractions wrote to `*.partial.npz` and were moved into place only on exit 0; no partial artefact survived |
| Orphan processes | none after the run |

The three cells are the same loci at three depths, which is what makes the depth
comparison paired.

## 12. Primary results

Loci scored per cell: **2,908,917** (1,675 loci of the 2,910,592 cached carry an
indel or no-call label and are outside the scoring frame, as in every devlog
since 9). SNP per cell: **4,758**. Pooled: **8,726,751** loci, **14,274** SNP.

| Cell | loci | SNP | F1 PB-only | F1 cascade | ΔF1 | ΔF1 95% CI | p | routed→PB | disagree | verdict |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---|
| full (69.6×) | 2,908,917 | 4,758 | 0.997380 | 0.997275 | −0.000105 | [−0.000511, +0.000212] | 0.79 | 0.0034% (99) | 3 | **PRESERVED** |
| 30x | 2,908,917 | 4,758 | 0.997272 | 0.997272 | **0.000000** | [0, 0] | — | 0.0131% (382) | 0 | **PRESERVED** |
| 15x | 2,908,917 | 4,758 | 0.989095 | 0.989095 | **0.000000** | [0, 0] | — | 0.1011% (2,941) | 0 | **PRESERVED** |
| **pooled** | 8,726,751 | 14,274 | 0.994600 | 0.994565 | −0.0000349 | [−0.000171, +0.0000703] | 0.79 | 0.0392% (3,422) | 3 | **PRESERVED** |

All three arms, pooled:

| arm | F1 | precision | recall | TP | FP | FN | TN | accuracy | balanced acc |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| binomial-only | 0.989686 | 0.991216 | 0.988160 | 14,105 | 125 | 169 | 8,712,352 | 0.999966 | 0.994073 |
| PB-only | 0.994600 | 0.995647 | 0.993555 | 14,182 | 62 | 92 | 8,712,415 | 0.999982 | 0.996774 |
| cascade | 0.994565 | 0.995577 | 0.993555 | 14,182 | **63** | 92 | 8,712,414 | 0.999982 | 0.996774 |

The cascade's entire pooled departure from PB is **one extra false positive** out
of 8.7 M loci. TP, FN and recall are identical to PB in every cell.

The cheap caller alone is *not* adequate, which is what makes the result
non-trivial: binomial-only vs PB-only is ΔF1 = −0.00094 (full, CI
[−0.00190, 0], p = 0.066), −0.00272 (30x, CI [−0.00385, −0.00168], p < 10⁻⁴)
and −0.01122 (15x, CI [−0.01356, −0.00904], p < 10⁻⁴). At 15x the binomial
loses 74 true SNP that PB finds; the cascade recovers **all 74** while running
PB on 0.101% of loci.

### 12.1 Compute and wall-clock

Throughput measured in a quiet single process on this experiment's own BAMs,
after every extraction had exited (128 kb probe, 123,520 loci):

| depth | binomial loci/s | PB loci/s | PB : binomial cost ratio |
|---|---:|---:|---:|
| full | 2,879,404 | 8,718 | 330× |
| 30x | 3,183,891 | 9,611 | 331× |
| 15x | 2,997,106 | 11,022 | 272× |

| Cell | loci | PB loci | PB compute fraction | PB-only s | cascade s | speedup |
|---|---:|---:|---:|---:|---:|---:|
| full | 2,908,917 | 99 | 0.00340% | 333.7 | 1.02 | **326.6×** |
| 30x | 2,908,917 | 382 | 0.01313% | 302.7 | 0.95 | **317.5×** |
| 15x | 2,908,917 | 2,941 | 0.10111% | 263.9 | 1.24 | **213.3×** |

PB computation eliminated: **99.997% / 99.987% / 99.899%** of it. The speedups
are projections from measured per-locus rates applied to the real routed counts,
not stopwatch readings of two full pipelines — that is the same accounting
devlogs 11–15 used, and it is stated as a projection rather than a measurement.
The pileup/read-tensor I/O that both arms share is excluded from both sides; the
measured extraction cost of that shared stage was 195–646 s (counts) and
422–1,358 s (read tensors) per cell, so a whole-pipeline speedup would be
smaller than the caller-stage speedup reported here.

## 13. Stratified benchmark (same region)

Full tables in `results/bench_v16/summary.md`. What they show:

* **Depth (fine bins).** ΔF1 = 0 in every bin of the 30x and 15x cells and in
  four of five bins at native depth. The two non-zero bins are 48–59× (ΔF1
  −0.001042, 2 disagreements, 948 SNP) and 60–79× (+0.000178, 1 disagreement) —
  both driven by single events, both at native depth. Routed fraction falls
  monotonically with depth within every cell (15x: 0.645% at 1–9×, 0.025% at
  10–29×; full: 0.366% at 10–29×, 0.000% above 80×), exactly as H2 predicted.
* **VAF.** Every disagreement in the experiment sits at VAF ≤ 0.19. The
  vaf 0.15–0.25 bin is where routing concentrates (73.6% of that bin routed at
  full depth, 92.3% at 30x, 100.0% at 15x) and where 2 of 3 disagreements land.
  Above VAF 0.25 the cascade equals PB everywhere at every depth.
* **Base quality.** All 3 disagreements sit at mean base quality < 30, in bins
  holding 0.7% of loci. Above Q30 the cascade is identical to PB in all three
  cells.
* **Mapping quality.** Uninformative on this region: 99.98% of loci have mean
  MAPQ ≥ 60 (novoalign emits 70). The MQ40–60 bin holds 255–475 loci and 1–4
  SNP. Nothing about mapping-quality dependence can be concluded here, and
  nothing is.
* **Genomic difficulty (GIAB v3.1).** In `lowmap_segdup` (53,298 loci, 172 SNP)
  the cascade equals PB at all three depths, F1 1.000000 at full and 30x. In
  `alldifficult` (494,289 loci, 853 SNP) ΔF1 is −0.000576 at full (1
  disagreement) and 0 at 30x/15x. In `tandemrepeats` ΔF1 = 0 everywhere. Routed
  fraction is 2–3× higher inside the difficult annotations than outside, i.e.
  the router spends its budget where the sequence is harder.
* **Stage-1 hard strata.** `hard_pb_uncertain` (loci within 13.5 of the PB
  threshold) routes 20.5% / 6.1% / 3.6% of its loci and contains all 3
  disagreements; ΔF1 there is 0 at 30x and 15x. `hard_low_vaf` ΔF1 −0.005061 at
  full (71 SNP — thin), 0 at 30x/15x. `hard_low_depth` ΔF1 = 0 in both cells
  where it exists.

## 14. Controls (matched PB budget)

| Cell | PB budget | frozen router ΔF1 | random routing ΔF1 (mean ± sd, 5 seeds) | reverse-confidence ΔF1 |
|---|---:|---:|---:|---:|
| full | 99 | **−0.000105** | −0.000940 ± 0.000000 | −0.000940 |
| 30x | 382 | **0.000000** | −0.002716 ± 0.000000 | −0.002716 |
| 15x | 2,941 | **0.000000** | −0.011196 ± 0.000043 | −0.011218 |

Spending the identical PB budget at random, or on the loci the binomial is
*most* confident about, recovers essentially none of the gap to PB: at 15x the
frozen router closes it entirely (ΔF1 = 0) while random routing sits at
−0.0112, i.e. no better than binomial-only (−0.01122). The router's advantage is
therefore the confidence ordering, not the budget. The random control's zero
variance at full and 30x is expected rather than suspicious: with 99–382 loci
drawn at random from 2.9 M, no seed ever hits a locus where PB and the binomial
disagree, so every seed reproduces binomial-only exactly.

This control is informative here, unlike in devlog 14, because the gap being
closed (binomial-only vs PB) is itself statistically significant at 30x and 15x.

## 15. Operating curve (diagnostic only — the cutoff stays frozen)

| cell | frozen point | first cutoff reaching ΔF1 = 0 | cost of going there |
|---|---|---|---|
| full | 5.412 → 0.0034% routed, ΔF1 −0.000105, 3 disagreements | 6.0 (0.0038% routed) | 326.6× → 326.2× |
| 30x | 5.412 → 0.0131% routed, ΔF1 0, 0 disagreements | 3.5 (0.0046%) | already there |
| 15x | 5.412 → 0.1011% routed, ΔF1 0, 0 disagreements | **5.412** | already there |

At 15x the frozen cutoff is the *smallest* swept cutoff that reaches ΔF1 = 0:
5.000 gives −0.000107 and 4.000 gives −0.000109. At 30x it is comfortably past
the knee (3.5 would have sufficed). At native depth the knee is just above it
(6.0). The cutoff was selected on chr21 in devlog 11 and it lands on the knee of
a chr13 curve it has never seen — evidence that 5.412 is a property of the
statistic rather than of chr21, which is what devlog 14 §H4 asked and could not
answer cleanly. **No production threshold is changed by this observation.**

## 16. Disagreement analysis and failure modes

**3 disagreements in 8,726,751 scored loci — rate 3.44×10⁻⁷**, all in the
native-depth cell, none at 30x or 15x. Pre-registered falsification threshold
was 10⁻⁴; the observed rate is 290× below it.

| pos (chr13) | truth | routed | effect | depth | VAF | baseQ | MAPQ | binom LLR | PB LLR | annotation |
|---:|---|---|---|---:|---:|---:|---:|---:|---:|---|
| 71,748,201 | non-SNP | no | `router_avoids_fp` | 65 | 0.113 | 27.6 | 70 | −2.09 | +11.07 | — |
| 72,736,919 | non-SNP | no | `router_fp` | 58 | 0.167 | 24.9 | 70 | +12.71 | +2.28 | — |
| 72,944,485 | non-SNP | no | `router_fp` | 57 | 0.188 | 26.5 | 70 | +18.39 | +8.32 | alldifficult |

* **Zero disagreements at truth-SNP loci** (0 of 14,274). The cascade lost no
  true variant that PB-only found, at any depth. Devlog 14's `router_fn`
  failure mode — 5 true SNP lost in chr1 segdup at VAF 0.11–0.13 — **did not
  recur here**. That is consistent with devlog 14's own localisation of it to
  segmental duplication: this region is 1.9% segdup against chr1_segdup's 100%,
  and inside this region's segdup the cascade and PB agree everywhere.
* **Net effect: −1 false positive.** Two cascade FPs against one PB FP avoided.
  Devlog 14's pooled departure was net *favourable*; here it is net
  unfavourable by one event. Both are single-digit counts on different regions
  and neither supports a claim about the sign.
* **No new failure mode.** All three events sit in the regime devlogs 13–15
  already characterised: fixed-ε binomial confidently wrong at low VAF
  (0.11–0.19) and low mean base quality (24.9–27.6), with the router unable to
  see them because their router margins (5.7–11.4) are all *outside* the cutoff
  — none of the three was routed. This is devlog 14 §17.3's sentence reproduced
  on a fourth kind of sequence, and it confirms rather than contradicts it.
* **Boundary.** Pooled disagreement rate by base quality: 1.8×10⁻⁴ at Q20–25,
  3.4×10⁻⁵ at Q25–30, **0** above Q30. By VAF: 7.4×10⁻⁴ at 0.15–0.25,
  2.4×10⁻⁶ at 0–0.15, **0** above 0.25. The operating boundary of the
  architecture on this region is therefore: **VAF ≳ 0.25 or mean base quality
  ≳ Q30 → the cascade is byte-identical to PB.**

## 17. Statistical interpretation, reproducibility and QC

* The 30x and 15x cells are **exact equality**, not statistical equivalence: the
  cascade and PB produce literally identical call vectors over 2,908,917 loci,
  which is why their bootstrap CI is the degenerate [0, 0]. No inference is
  needed for those two cells and none is claimed from the degenerate interval.
* The native-depth cell and the pool are **statistically indistinguishable and
  practically equivalent**: ΔF1 −0.000105 (CI [−0.00051, +0.00021]) and
  −0.0000349 (CI [−0.00017, +0.00007]). Both intervals exclude any effect larger
  than ±0.0005 in either direction, so this is a bounded null rather than
  "p > 0.05 therefore equal": the effect size that the data can still hide is
  smaller than the pre-registered practical threshold of 0.001, and both
  contain 0.
* Effect direction is honestly negative at native depth (one extra FP) and the
  §7 rule calls it PRESERVED because |ΔF1| < 0.001 and the CI contains 0. It is
  not called "identical", which would be false for that cell.
* Power: every cell carries 4,758 SNP against the ≥ 100 requirement, and CI
  half-widths are 0.00036 / 0 / 0 against the ≤ 0.01 requirement. No cell is
  UNDERPOWERED. The 3 Mb window sizing in §4.1 was, if anything, conservative.
* **Not powered** to say anything about: mapping quality (99.98% of loci at
  MAPQ ≥ 60), segmental duplication beyond 172 SNP, or the *sign* of a
  single-event ΔF1.

QC actually run after execution:

```
independence audit                    all_ok = True over 3 cells
  contig chr13, region [70000000, 73000000] in all 3 cells
  loci_cached == loci_from_window_tiling == 2,910,592 in all 3 cells
  window-aligned, 0 non-finite binomial LLR, 0 non-finite PB LLR
  coordinate overlap with the 13 historical regions: [] (none)
non-finite scan of both result artefacts      -> none
arm-consistency re-derivation from the caches -> routed counts, disagreement
  counts, TP and loci_scored all reproduce the report exactly, in all 3 cells
pytest (12 suites, incl. the new test_bench_v16) -> 343 passed
git status --porcelain (tracked files)        -> no tracked file modified
orphan processes                              -> none
```

Frozen-file integrity re-verified after the run: `cascade.py`,
`cheap_router.py`, `binomial_baseline.py`, `quality_error_model.py`,
`robustness_benchmark.py`, `bench_v13_robustness.py`, `bench_v12_stage1.py`,
`extract_bench_v12.py`, `pileup_counts.py`, `read_level_pileup.py` all still
IDENTICAL-TO-HEAD (`1dfa759`). Nothing was committed or pushed.

### 17.1 Deviations from the pre-registration — disclosed

1. **Mapping-quality bin edges were wrong and were fixed after seeing the
   output.** `mapq_stratified` was written with edges topping out at 61,
   mirroring the base-quality stratifier. This library's MAPQ is 70, so those
   bins covered 0.02% of loci. The edges were corrected to `(0,20) (20,40)
   (40,60) (60,71)` and the axis recomputed into
   `results/bench_v16/mapq_stratified.json`; `benchmark_results.json` was left
   untouched and its `by_mapq` block is superseded by that file. This is a
   defect in a diagnostic axis: no arm, endpoint, threshold, control or verdict
   depends on it, and the corrected axis changes no conclusion (it is
   degenerate either way). Disclosed rather than silently re-run.
2. **A `KeyError` in `summarize_v16.py`** (the base-quality blocks store two F1
   values but no pre-computed delta) was fixed and the markdown re-rendered.
   Pure formatting; the JSON artefacts were not regenerated.

Nothing else was changed after any number was read. No region, depth, cutoff,
threshold or evaluation criterion was altered at any point.

Reproduction, in order: `select_v16_region.py --contig-only` → `select_v16_region.py`
→ `fetch_v16_reference.sh` → `fetch_v16_data.sh` → `prepare_v16_depths.sh` →
`run_v16_extract.sh` → `bench_v16_unseen.py` → `summarize_v16.py`. Determinism:
downsampling seed 42, bootstrap seed 20260812 (module default), control seeds
0–4, fixed cutoff grid. Machine: 8 cores, 31 GB RAM, Linux 6.19.12, Python
3.14, numpy 2.4.4, scipy 1.17.1, pysam 0.24.0, samtools/bcftools 1.23.1.
End-to-end wall clock for the whole experiment: ~54 min from first download to
final artefact (extraction 39.7 min of it), plus ~2.5 min of region selection
and reference download.

## 18. Limitations

Everything in §10 still holds, and is now specific:

* **One sample, one platform, one aligner, one truth set, one reference build.**
  HG002 / Illumina 2×250 / novoalign / GIAB v4.2.1 / GRCh38. Cross-sample and
  cross-platform transfer remain **untested** after six devlogs.
* **SNP only.** Indel and no-call loci are outside the scoring frame; 854 of the
  region's 5,743 truth records are non-SNV and are not scored by any arm.
* **One 3 Mb window.** Representative *of chr13 by an annotation rule*, not of
  the genome. It contains 1.9% segdup, so it tests the known failure mode only
  weakly — 172 SNP inside `lowmap_segdup`, where the cascade happened to be
  perfect. A region like devlog 14's chr1_segdup would stress it far harder, and
  this experiment does not repeat that stress.
* **The depth cells are downsamples of one library**, not independent
  sequencing runs; their agreement is not three independent confirmations.
* **The speedups are projections** from measured per-locus caller throughput,
  and exclude the pileup/read-tensor I/O both arms share.
* **Three disagreements** is a small sample for characterising a failure mode.
  The regime description in §16 is consistent with devlogs 13–15 but rests here
  on 3 events.
* **Mapping quality is untestable on this region** (§13).

## 19. Final verdict

**STRONG POSITIVE**, scoped to HG002 / Illumina 2×250 novoalign / GRCh38 /
GIAB v4.2.1 SNP calling on chr13:70–73 Mb.

Against the §7 criteria, checked in order: every powered cell PRESERVED (3/3,
none underpowered); pooled PRESERVED; PB compute fraction 0.0034%–0.1011%, all
≤ 1%; projected speedup 213×–327×, all ≥ 50×; no new failure mode — the three
disagreements fall inside the regime devlogs 13–15 already documented, and the
one failure devlog 14 warned about (true SNP lost in segdup) did not occur.

Concretely: on 8.73 M unseen loci the frozen cascade reproduced PB-only's calls
exactly at 30x and 15x, differed from PB at 3 loci out of 8.73 M at native
depth, cost one net false positive in total, lost no true variant, and did it
while running the expensive caller on 0.039% of loci.

This is **not** a claim that the architecture is proven. It is a claim that it
survived a genuinely unseen, mechanically selected, representative region under
a pre-registered decision rule — the sixth region family it has survived, and
still on one sample and one platform.

Devlog 14's §17.3 characterisation of the departure regime is **confirmed**, not
contradicted. No previous claim is contradicted by this experiment.

## 20. Future work (recorded, not implemented)

1. **Cross-sample and cross-platform transfer** — HG003/HG004, or PacBio HiFi /
   ONT. This is now the largest untested axis by a wide margin; every accuracy
   claim this project has made is about one sample on one instrument.
2. **A segdup-heavy unseen region at native depth**, to test the `router_fn`
   failure mode directly rather than incidentally. This region was too easy on
   that axis by construction.
3. **Indels.** Excluded from every devlog's scoring frame so far.
4. **Whole-pipeline wall-clock**, measuring the shared pileup/read-tensor stage
   end to end, so the speedup can be quoted for the pipeline rather than the
   caller stage.
5. **Observed but not acted on:** the frozen cutoff sits exactly on the 15x knee
   of a curve from a chromosome it never saw (§15). That is a testable claim
   about the cutoff being a property of the statistic; it deserves its own
   experiment and must not be inferred from this one.
6. The ε-sensitivity safety layer of devlog 15 was not exercised here. Whether
   it would have changed any of the 3 disagreements is unknown and was
   deliberately not measured, to keep this benchmark about the frozen
   production architecture.
