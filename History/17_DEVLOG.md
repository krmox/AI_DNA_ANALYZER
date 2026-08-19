# Does the frozen binomial → router → PB cascade hold on the hardest unseen contig? — chr17:18–21 Mb segmental-duplication benchmark

**Date:** 2026-08-19

**Status of sections:** §1–§10 are the pre-registration. They were written and
saved **before any caller, router or metric was run on any chr17 locus**, and
before any chr17 BAM slice had been downloaded. Everything from §11 onwards was
written after the corresponding measurement.

Devlogs 3–16 are not edited by this experiment. `cascade.py`, `cheap_router.py`,
`binomial_baseline.py`, `quality_error_model.py`, `robustness_benchmark.py`,
`extract_bench_v12.py`, `pileup_counts.py`, `read_level_pileup.py`,
`bench_v12_stage1.py`, `bench_v13_robustness.py` and `bench_v14_crosschrom.py`
are not edited by this experiment; §17 verifies that mechanically. Nothing is
committed or pushed.

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
benchmark, and neither is the PB → Mamba branch. Devlog 16 validated a
*median-difficulty* window on chr13.

This experiment is not an attempt to improve anything. It deliberately targets
the one regime devlog 14 flagged, on a contig no experiment in this project has
ever touched, and asks:

1. Does the frozen cascade reproduce PB-only accuracy on a segdup-saturated window?
2. How much PB computation does it remove there?
3. What is the real wall-clock speedup?
4. Where and why does it disagree with PB?
5. Are the disagreements concentrated in particular depth / VAF / base-quality /
   mapping-quality / genomic-difficulty regimes?
6. Does the router beat matched-budget random and reverse-confidence routing?
7. Does the devlog 14 segdup failure mode reappear, and at what rate?

## 2. Hypotheses

* **H1 (primary).** On chr17:18–21 Mb the cascade's SNP F1 is statistically
  indistinguishable from PB-only: the paired-bootstrap 95% CI for ΔF1 contains 0
  and |ΔF1| < 0.001, at native, 30× and 15× depth.
* **H2 (compute).** The frozen router sends ≤ 1% of loci to PB.
* **H3 (router signal).** At matched PB budget the frozen router beats random
  routing and reverse-confidence routing.
* **H4 (failure mode).** The devlog 14 segdup failure mode reappears here at an
  elevated but still ≤ 10⁻⁴ disagreement rate.

H4 is stated as a *prediction*, not a success criterion; §7 governs the verdict.

## 3. Frozen configuration — not tuned, not refit, not swept for production

Recorded in full in `results/bench_v17/freeze_audit.txt` before any chr17 data
was touched.

| item | value |
|---|---|
| git HEAD | `426c209d4d86226fa0c3227be61e32eb5aec0910` |
| working tree | clean except the new `results/bench_v17/` |
| binomial threshold | 7.0 |
| binomial error rate ε | 0.01 |
| PB threshold | 10.5 |
| frozen router cutoff | 5.411872376933351 |
| reference | Ensembl release-110 GRCh38, per-chromosome FASTA |
| truth | GIAB HG002 v4.2.1 benchmark VCF, SNP only |
| high-confidence BED | v4.2.1 `_noinconsistent` BED |
| stratifications | GIAB v3.1 (`lowmap_segdup`, `alldifficult`, `tandemrepeats`) |
| python / numpy / scipy / pysam | 3.14.3 / 2.4.4 / 1.17.1 / 0.24.0 |
| samtools / bcftools | 1.23.1 / 1.23.1 (htslib 1.23.1) |

SHA-256 of every frozen implementation file is in
`results/bench_v17/freeze_audit.txt`. §17 re-hashes them after the benchmark and
requires an exact match.

**Bug rule, fixed in advance.** If a defect is found in any frozen file during
this experiment it is **not** fixed here. It is recorded in §17 and every cell it
touches is classified `INVALID` and excluded from the verdict rather than
re-scored. Only an outright crash in *non-frozen* v17 glue code may be repaired,
and any such repair is disclosed in §17.

## 4. New-region selection and provenance

### 4.1 The rules, fixed before execution

Implemented in `select_v17_region.py`, whose docstring is the authoritative
statement. Both rules read only the HC BED, the stratification BEDs, and (only
after the contig is already committed) the reference FASTA. No caller score,
router decision or metric is an input.

* **Rule C17 (contig).** Among the GRCh38 autosomes no previous experiment has
  touched — every autosome except chr21, chr20, chr19, chr4, chr1, chr16, chr15,
  chr7, chr14, chr13 — rank by `segdup_hc_fraction = |HC ∩ lowmap_segdup| / |HC|`
  descending and take **rank 0**, the hardest remaining contig. Ties → smaller
  contig number.
* **Rule R17 (window).** Among all 3 Mb windows at whole-Mb offsets on that
  contig with `hc_fraction ≥ 0.50` and `n_fraction == 0`, take the window
  **maximising** `lowmap_segdup_hc_fraction`. Ties → smaller coordinate.

The `hc_fraction ≥ 0.50` floor and the zero-N requirement exist only so the
window carries scorable truth and callable sequence; both are verbatim from
devlog 15 Rule H / devlog 16 Rule R16. The 3 Mb width is copied from devlog 16
for statistical power and is not a property of this data.

Rule C17 is the deliberate opposite of devlog 16's median rule: this experiment
is a hard-data stress test, so it takes the worst remaining contig rather than a
representative one.

### 4.2 Rule C17 output

Candidates ranked by `segdup_hc_fraction` (top 5): chr17 0.0861, chr22 0.0827,
chr9 0.0712, chr10 0.0646, chr11 0.0637. **Rank 0 = chr17.**

### 4.3 The selected region

**chr17:18,000,000–21,000,000** (0-based half-open), 3,000,000 bp.

| property | value |
|---|---|
| HC fraction | 0.8510 |
| N fraction | 0.0 |
| GC fraction | 0.4628 |
| `lowmap_segdup` fraction (of window) | 0.4348 |
| `lowmap_segdup` fraction **of HC bases** | **0.3662** |
| `alldifficult` fraction of HC bases | **0.4514** |
| `tandemrepeats` fraction of HC bases | 0.0469 |
| eligible windows considered | 70 of 81 |

For comparison, devlog 16's chr13 window carried `lowmap_segdup_hc_fraction`
≈ 0.01–0.02 and `alldifficult_hc_fraction` ≈ 0.16. This window is roughly an
order of magnitude more segdup-saturated: 37% of its scorable bases sit in
low-mappability segmental duplication, and 45% are annotated difficult. This is
17p11.2 — the Smith–Magenis / *KRT*-proximal duplication block — one of the
hardest short-read regions in the genome that still carries GIAB truth.

Eleven windows were excluded by the pre-registered filters, all of them
centromere-proximal windows with non-zero N and `hc_fraction` far below 0.50
(e.g. start 24 Mb: `segdup_hc = 1.00` but `hc_fraction = 0.007`). Those windows
are harder still, but unscorable; the exclusion rule was fixed before the
statistics were computed.

### 4.4 Depth realisations

The native slice is downsampled to nominal 30× and 15× with
`samtools view -s 42.FRACTION`, the fraction derived from the region's
**measured** native mean depth (`samtools depth -a`), recipe verbatim from
`prepare_v14/v15/v16_depths.sh`. These are three depth regimes of **one** region
from **one** library, never three independent regions; §17 and §18 treat them as
correlated.

### 4.5 Provenance

| item | source |
|---|---|
| BAM | GIAB `AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam`, sliced by `samtools view -b <url> chr17:18000001-21000000` |
| truth VCF | GIAB `NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`, sliced by `bcftools view -r` |
| HC BED | `data/giab_hg002_v14/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`, awk-restricted to the window |
| reference | Ensembl release-110 `Homo_sapiens.GRCh38.dna.chromosome.17.fa.gz` |
| stratifications | GIAB v3.1, already present in `data/strat_v31/` |

Scripts: `fetch_v17_reference.sh`, `fetch_v17_data.sh`, `prepare_v17_depths.sh`,
`run_v17_extract.sh`, `bench_v17_unseen.py`, `summarize_v17.py`.

## 5. Independence audit

Contig-level exclusion already makes overlap impossible — Rule C17's candidate
set excludes every previously used contig — but the audit checks coordinates
mechanically anyway. `bench_v17_unseen.verify_caches` asserts, per cell, that the
cached contig is chr17, that the cached region is exactly [18000000, 21000000),
that every locus position lies inside it, that the contig is not in
`PREVIOUSLY_USED_CONTIGS`, and that the window has zero coordinate overlap with
the explicit list of every region any previous devlog extracted. It also checks
window alignment, cache/tiling locus-count agreement and finiteness of both LLR
arrays. The artifact is `results/bench_v17/independence_audit.json`; a single
failing check aborts scoring.

## 6. Experimental arms and metrics

Three arms per cell, all on the identical scoring frame (callable loci, SNP vs
non-SNP, indels excluded):

* **A — binomial only:** `binomial_LLR ≥ 7.0`.
* **B — PB only:** `pb_LLR ≥ 10.5` on every locus.
* **C — frozen cascade:** binomial answer where
  `|binomial_LLR − 7.0| > 5.411872376933351`, otherwise the PB answer.

Reported per arm: TP, FP, FN, TN, precision, recall, F1, accuracy, balanced
accuracy. Reported for C: ΔF1 vs B with a paired bootstrap 95% CI, PB routing
count and fraction, PB compute fraction, measured binomial and PB throughput,
projected wall-clock and speedup, cascade-vs-PB disagreement count and rate.

Stratifications: depth (coarse and fine), base quality, VAF, mapping quality,
and the three GIAB structural annotations (`lowmap_segdup`, `alldifficult`,
`tandemrepeats`) — the last being the point of this region.

Controls at matched PB budget: random routing and reverse-confidence routing.

Every cascade/PB disagreement is dumped individually with position, depth, VAF,
mean base quality, mean mapping quality, annotations, truth label and all three
arm calls.

## 7. Pre-registered decision rule

Quoted verbatim from devlog 13 §3.4 / devlog 14 §3.4 / devlog 16 §7 so it cannot
drift. Per cell and for each pooled result:

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
* **POSITIVE / PRESERVED** — the above but with a material limitation: some cells
  UNDERPOWERED, or the controls uninformative, or the compute reduction between
  10× and 50×, or small localized discrepancies.
* **MIXED / PARTIAL GENERALIZATION** — pooled result PRESERVED/IMPROVED but at
  least one powered stratum DEGRADED, or a new failure regime appears.
* **NEGATIVE** — any powered *cell* DEGRADED, or the projected speedup < 10×.
* **INCONCLUSIVE** — too few powered cells to distinguish the above.

### 7.1 What would falsify H1

Any powered cell whose ΔF1 CI excludes 0 on the negative side; or a disagreement
rate above 10⁻⁴ (chr21 saw ~6×10⁻⁷; devlog 14's ordinary regions ~10⁻⁶; devlog
16's chr13 window ~10⁻⁶).

## 8. Stopping criteria

Extraction runs **once** per cell. Metrics are computed **once**. No region,
depth, threshold, cutoff or evaluation criterion is added, removed or changed
after any result is read, except to fix an outright crash in non-frozen v17 glue
— any such re-run is disclosed in §17. **If this region performs poorly it is
reported as a negative result and no second region is drawn.** If the frozen
cutoff fails here, the failure is reported and the cutoff stays at
5.411872376933351.

The experiment stops when: 3 cells are extracted and scored; all three arms and
all metrics of §6 are computed; controls are run or documented as underpowered;
the disagreement analysis is complete; compute accounting is verified; the test
suite passes; and §17 records reproducibility. No further cycle is run regardless
of outcome.

## 9. Runtime plan

Reused without regeneration: the GIAB v3.1 stratification BEDs, the HC BED, the
v4.2.1 truth VCF index, the extractor, every metric helper, the bootstrap, the
controls and the sweep. Regenerated: nothing that already exists.

New compute is exactly: one reference FASTA download (done, 16 s), one 3 Mb BAM
slice download, two `samtools view -s` downsamples, three extractions (the only
expensive step), and one scoring pass over cached arrays. Three concurrent
extraction processes at most.

## 10. What this experiment cannot show

Fixed by design, stated before results exist: one sample (HG002), one platform
(Illumina 2×250, novoalign), one aligner, one truth set (v4.2.1), one reference
build, SNP only (indels excluded from the scoring frame), and one 3 Mb window on
one chromosome. The three depth cells are downsamples of one library and are
correlated, not independent replicates. GIAB truth is *itself* least reliable in
exactly the segdup regime this window targets, so the HC BED excludes the worst
of it — meaning this benchmark measures the cascade on the **scorable** hard
regions, not on the unscorable ones. A positive result licenses no claim about
other samples, platforms, variant classes or the whole genome. §18 restates this
against whatever is actually found.

---

*(§11 onward is written after the corresponding measurement.)*

## 11. Data extraction and QC

Everything below was measured after §1–§10 were saved.

| step | result |
|---|---|
| reference | Ensembl release-110 chr17, 84.6 MB, 16 s |
| BAM slice | `chr17:18000001-21000000`, 123.7 MB, 49 s |
| truth VCF slice | 3,679 records, of which **3,135 SNP** |
| HC BED slice | 1,381 intervals, 2,562,660 bp |
| REF-base concordance | **3,134 / 3,134** single-base VCF REF alleles match the Ensembl FASTA (0 mismatches) |
| native mean depth | 65.13× (`samtools depth -a` over the window) |
| 30× realisation | `-s 42.4606` → achieved 29.95× |
| 15× realisation | `-s 42.2303` → achieved 14.99× |
| extraction | 3 cells, all `rc=0`, 14–37 min each, run concurrently |
| cached loci | 2,467,776 per cell; 2,466,725 enter the scoring frame; **2,828 SNP** |

The scored SNP count (2,828) is below the VCF's 3,135 because the scoring frame
keeps only callable SNP/non-SNP loci inside the HC BED and drops indels.

**Difficulty is real, not just annotated.** PB-only F1 on this window is 0.9507
(native depth) against 0.99+ on chr21 and chr13. Inside `lowmap_segdup` PB-only
F1 falls to 0.8837; outside it, 0.9903 on the same cell. This benchmark is
roughly an order of magnitude harder than any previous one in this project.

## 12. Primary results

Frozen configuration, three arms, one region, three depths.

| cell | loci | SNP | F1 binomial | F1 PB | F1 cascade | ΔF1 vs PB | 95% CI | verdict |
|---|---:|---:|---:|---:|---:|---:|---|---|
| full (69.5×) | 2,466,725 | 2,828 | 0.952565 | 0.950682 | **0.954506** | **+0.003823** | [+0.002267, +0.005506] | **IMPROVED** |
| 30× (32.0×) | 2,466,725 | 2,828 | 0.945903 | 0.952078 | **0.953410** | **+0.001332** | [+0.000316, +0.002483] | **IMPROVED** |
| 15× (16.0×) | 2,466,725 | 2,828 | 0.927256 | 0.941345 | **0.941850** | **+0.000505** | [+0.000000, +0.001168] | **PRESERVED** |
| **pooled** | 7,400,175 | 8,484 | 0.941988 | 0.948070 | **0.949962** | **+0.001892** | [+0.001250, +0.002615] | **IMPROVED** |

Full confusion counts:

| cell | arm | P | R | TP | FP | FN | accuracy | balanced acc |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| full | binomial | 0.946578 | 0.958628 | 2711 | 153 | 117 | 0.99989054 | 0.979283 |
| full | PB | 0.940484 | 0.961103 | 2718 | 172 | 110 | 0.99988568 | 0.980517 |
| full | cascade | 0.948342 | 0.960750 | 2717 | 148 | 111 | 0.99989500 | 0.980345 |
| 30× | binomial | 0.942747 | 0.949081 | 2684 | 163 | 144 | 0.99987554 | 0.974507 |
| 30× | PB | 0.952246 | 0.951909 | 2692 | 135 | 136 | 0.99989014 | 0.975927 |
| 30× | cascade | 0.955272 | 0.951556 | 2691 | 126 | 137 | 0.99989338 | 0.975752 |
| 15× | binomial | 0.937477 | 0.917256 | 2594 | 173 | 234 | 0.99983500 | 0.958593 |
| 15× | PB | 0.952243 | 0.930693 | 2632 | 132 | 196 | 0.99986703 | 0.965320 |
| 15× | cascade | 0.953278 | 0.930693 | 2632 | 129 | 196 | 0.99986825 | 0.965320 |

No cell is DEGRADED. Two of three are IMPROVED — the cascade is *more accurate*
than PB-only here, which no previous devlog observed at cell level. §16 explains
why, and §17 argues why this must not be read as the cascade being a better
caller in general.

Balanced accuracy is the one metric where PB-only edges ahead at native depth
(0.980517 vs 0.980345): PB recovers 1 more true SNP than the cascade while
emitting 24 more false positives. F1 and accuracy both favour the cascade; the
recall difference is 1 SNP.

### 12.1 Compute and wall-clock

Throughput measured in a quiet single process on this experiment's own BAMs,
after every extraction had exited.

| cell | routed→PB | PB compute fraction | binomial loci/s | PB loci/s | PB-only s | cascade s | **speedup** |
|---|---:|---:|---:|---:|---:|---:|---:|
| full | 299 | **0.01212%** | 2,829,942 | 8,083 | 305.16 | 0.909 | **335.8×** |
| 30× | 797 | **0.03231%** | 3,092,361 | 9,058 | 272.32 | 0.886 | **307.5×** |
| 15× | 3,213 | **0.13025%** | 2,794,387 | 9,682 | 254.76 | 1.215 | **209.8×** |

H2 (≤ 1% routed) holds with three orders of magnitude of headroom even on the
hardest region tested. Routing load rises monotonically as depth falls — 0.012%
→ 0.032% → 0.130% — which is the expected behaviour of a margin-based router
when evidence thins, not a new effect.

## 13. Stratified benchmark

**Structure** — the point of this region:

| cell | stratum | loci | SNP | F1 PB | F1 cascade | ΔF1 | disagree |
|---|---|---:|---:|---:|---:|---:|---:|
| full | in `lowmap_segdup` | 875,964 | 1,048 | 0.883655 | 0.892381 | +0.008726 | 23 |
| full | not in `lowmap_segdup` | 1,590,761 | 1,780 | 0.990264 | 0.990815 | +0.000551 | 2 |
| full | in `alldifficult` | 1,091,026 | 1,321 | 0.901713 | 0.909091 | +0.007378 | 24 |
| full | not in `alldifficult` | 1,375,699 | 1,507 | 0.994063 | 0.994391 | +0.000328 | 1 |
| full | in `tandemrepeats` | 118,099 | 181 | 0.925450 | 0.932642 | +0.007193 | 3 |
| 30× | in `lowmap_segdup` | 875,964 | 1,048 | 0.881814 | 0.885119 | +0.003305 | 10 |
| 30× | not in `lowmap_segdup` | 1,590,761 | 1,780 | 0.992741 | 0.992741 | +0.000000 | 0 |
| 15× | in `lowmap_segdup` | 875,964 | 1,048 | 0.865185 | 0.866469 | +0.001284 | 3 |
| 15× | not in `lowmap_segdup` | 1,590,761 | 1,780 | 0.984581 | 0.984581 | +0.000000 | 0 |

Every disagreement in the whole experiment except one sits inside
`alldifficult`, and 36 of 38 inside `lowmap_segdup` — 36.6% of the bases carry
95% of the disagreements. That is exactly the concentration devlog 14 predicted.

**VAF** — the low-VAF bin is where everything happens:

| cell | VAF bin | loci | SNP | F1 PB | F1 cascade | ΔF1 | routed | disagree |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| full | 0.00–0.15 | 252,230 | 12 | 0.115385 | 0.148148 | +0.032764 | 0.03% | 25 |
| full | 0.15–0.25 | 231 | 15 | 0.291262 | 0.291262 | 0.000000 | 73.59% | 0 |
| full | 0.25–0.35 | 88 | 47 | 0.783333 | 0.783333 | 0.000000 | 18.18% | 0 |
| full | 0.35–0.45 | 314 | 294 | 0.967105 | 0.967105 | 0.000000 | 0.00% | 0 |
| full | ≥ 0.45 | 2,398 | 2,365 | ~0.99 | ~0.99 | 0.000000 | ≤1.5% | 0 |
| 15× | 0.15–0.25 | 2,617 | 75 | 0.635838 | 0.635838 | 0.000000 | **99.92%** | 0 |
| 15× | 0.25–0.35 | 371 | 130 | 0.848485 | 0.848485 | 0.000000 | 73.05% | 0 |

The router's behaviour in the genuinely ambiguous VAF band is textbook: at 15×
it routes **99.9%** of the 0.15–0.25 bin to PB and reproduces PB exactly there.
Every disagreement lives in the VAF < 0.15 bin, where the truth is
overwhelmingly non-SNP.

**Base quality** and **mapping quality**: ΔF1 ≥ 0 in every bin at every depth;
the lowest-MAPQ bin (20–40) has PB-only F1 0.68–0.80 and the cascade matches or
beats it (+0.0044 / +0.0046 / 0.0000). No stratum degrades.

**Depth (fine bins)**: ΔF1 ≥ 0 in all 15 populated bins. The largest gain is at
80×+ (+0.0079), where excess depth in collapsed duplications drives PB false
positives.

**Stage-1 hard strata**: `hard_low_depth` ΔF1 = 0.000000 at all three depths;
`hard_low_vaf` +0.0393 / +0.0137 / +0.0045; `hard_pb_uncertain` +0.0153 /
+0.0028 / +0.0030; `hard_low_quality` +0.0097 / 0.0000 / 0.0000. No hard stratum
loses accuracy.

## 14. Controls (matched PB budget)

| cell | frozen router F1 | random routing (5 seeds) | reverse-confidence F1 |
|---|---:|---:|---:|
| full | **0.954506** | 0.952565 ± 0.000000 | 0.952565 |
| 30× | **0.953410** | 0.945903 ± 0.000000 | 0.945903 |
| 15× | **0.941850** | 0.927256 ± 0.000000 | 0.927256 |

The frozen router beats both controls at every depth. Both controls collapse
exactly onto the binomial-only F1 with zero variance across seeds, and the reason
is mechanical rather than surprising: with 299–3,213 loci to spend out of 2.47 M,
a random or margin-inverted draw essentially never lands on one of the ≲ 40 loci
where PB and the binomial actually differ, so the PB budget is spent entirely on
loci where PB would have agreed anyway. **H3 holds**, and the control confirms
the frozen router's selection carries essentially all of the available signal:
its 299 chosen loci change the answer, an arbitrary 299 change nothing.

## 15. Operating curve (diagnostic only — the cutoff stays frozen)

**Secondary / exploratory.** This sweep does not and cannot alter §18's verdict;
the production cutoff remains 5.411872376933351.

Native-depth cutoff sweep: ΔF1 vs PB is positive across the entire swept range
0.5 → 20.0 (+0.0031 to +0.0050), and the frozen cutoff (5.4119) sits at the
sweep's **maximum**, ΔF1 = +0.003823, at a PB compute fraction of 0.0121% and
336× speedup. Raising the cutoff to 20 would cut disagreements from 25 to 3 but
also cut ΔF1 to +0.0005 and speedup to 34×. The frozen point was chosen on chr21
in devlog 13 and had never seen chr17; that it lands on this region's optimum is
a favourable coincidence, not evidence that the cutoff was tuned here.

## 16. Disagreement analysis and failure modes

38 cascade-vs-PB disagreements over 7,400,175 scored loci — pooled rate
**5.14 × 10⁻⁶**, well inside the pre-registered 10⁻⁴ falsification bound (§7.1)
though ~5–8× above the ordinary-region rate (~10⁻⁶) as expected on this data.

Every one of the 38 is the same structural event: a locus the router kept
**out** of PB (`routed_to_pb = false`) because the binomial was confident, where
PB-only would have called a SNP.

Biological classification, all 38 inspected individually:

| classification | count |
|---|---:|
| **cascade correct / PB wrong** | **36** |
| **PB correct / cascade wrong** | **2** |
| both wrong | 0 |
| unresolved | 0 |

| category | count |
|---|---:|
| PB false positive suppressed by the cascade | 36 |
| true SNP lost by the cascade (false negative) | 2 |
| non-SNP locus | 36 |
| true SNP locus | 2 |

The 36 suppressed events are one coherent biological population: depth 24–106,
**VAF 0.071–0.130**, mean base quality 24.8–38.9, MAPQ 36.6–70.0, and 35 of 36
annotated `lowmap_segdup`. That is the signature of **paralogous sequence
variants** — reads from the duplicate copy of a segmental duplication mismapping
into this one at a fixed low allele fraction. PB, which models per-read error,
reads the systematic low-VAF pile-up as evidence of a real variant; the binomial
LLR, comparing to a flat ε = 0.01, is confidently negative and the router never
hands the locus over. In this regime the cheap caller's cruder model is the more
robust one, and the cascade inherits that.

The 2 losses (the failure mode that matters):

| cell | position | truth | binom LLR | PB LLR | depth | VAF | baseQ | MAPQ | annotations |
|---|---:|---|---:|---:|---:|---:|---:|---:|---|
| full | chr17:18,640,872 | **SNP** | +1.57 | 11.78 | 34 | 0.129 | 36.7 | **45.2** | alldifficult, lowmap_segdup |
| 30× | chr17:18,461,522 | **SNP** | −0.70 | 12.27 | 35 | 0.115 | 38.8 | 65.4 | alldifficult, lowmap_segdup |

Both are real GIAB SNPs presenting at VAF 0.115–0.129 inside segmental
duplication — i.e. true variants whose read support is indistinguishable, at the
binomial's resolution, from the 36 paralogous artefacts above. Both have binomial
LLR near 0 but outside the ±5.41 router band on the confident-negative side, so
neither reaches PB. This is **the devlog 14 segdup failure mode, unchanged in
character** — H4's prediction — occurring here at 2 events in 2.47 M loci per
cell (8.1 × 10⁻⁷). It is not new.

Failure-boundary concentration: by VAF, all 38 disagreements fall in the
0.00–0.15 bin (rate 8.75 × 10⁻⁵ there, 0 everywhere else); by depth the rate
rises with depth (1.0 × 10⁻⁶ at 10–29× to 1.4 × 10⁻⁵ at 80×+), consistent with
collapsed-duplication pile-ups being the driver rather than evidence scarcity.

**No genuinely new failure mode was found.**

## 17. Statistical interpretation, reproducibility and QC

Post-run verification (`results/bench_v17/post_run_verification.txt`):

* SHA-256 of all 13 frozen implementation files is **byte-identical** to the
  pre-run audit. Nothing frozen was edited.
* Independence audit `all_ok = true` over 3 cells; contig chr17, region exactly
  [18000000, 21000000), zero coordinate overlap with all 14 historical regions,
  and chr17 is not in `PREVIOUSLY_USED_CONTIGS`.
* All four result JSONs parse; **no non-finite metric anywhere** in the report.
* Arm consistency re-derived independently from the raw `.npz` caches: the
  cascade call equals the frozen rule at every locus in all three cells; both LLR
  arrays finite everywhere.
* No orphan processes; all three extractions exited `rc=0` and the atomic
  `.partial.npz` staging left no truncated cache.
* Test suite: **273 passed** before the benchmark and **273 passed** after,
  across `test_cascade`, `test_cheap_router`, `test_bench_v12/13/14/15/16/17`,
  `test_robustness_benchmark`, `test_binomial_baseline`, `test_safety_layer`.
  `test_bench_v17.py` (17 tests) additionally verifies that the driver never
  assigns a frozen constant, that its region equals what `select_v17_region.py`
  independently chose, and that the independence audit actually fails on each of
  five kinds of corruption.

**On the IMPROVED verdicts.** The cascade beating PB by ΔF1 = +0.0019 pooled is
statistically solid (CI excludes 0) but must not be over-read. It is not evidence
that skipping PB improves calling. It is one specific interaction: in
segmental duplication, PB's per-read model is *miscalibrated upward* on
paralogous pile-ups, and the frozen router happens to keep those loci away from
it. Change the region, the aligner, or PB's calibration and the sign could
reverse. The defensible claim is the pre-registered one — **accuracy is not
lost** — with the improvement recorded as an observation about PB's segdup
behaviour, not as a property of the cascade.

**Deviations from the pre-registration:** none. No threshold, cutoff, region,
depth or criterion was changed after any result was read; no bug was found in any
frozen file; the experiment ran exactly once.

## 18. Limitations

As pre-registered in §10, plus what the results now make specific:

* One sample (HG002), one platform (Illumina 2×250, novoalign), one aligner, one
  truth set (v4.2.1), one reference build, SNP only, one 3 Mb window.
* The three depth cells are `samtools view -s` downsamples of one library. They
  are correlated; the pooled CI treats 8,484 SNP as if independent and is
  therefore optimistic. The per-cell CIs are the honest ones.
* **GIAB truth is weakest exactly where this experiment looked.** The HC BED
  excludes the worst of 17p11.2, so 14.9% of the window is unscorable. The
  cascade was measured on the *scorable* hard regions; the hardest segdup —
  including the windows Rule R17's `hc_fraction ≥ 0.50` filter excluded, one of
  which is 100% segdup — remains untested because no trustworthy truth exists
  there.
* 2,828 SNP per cell is enough for the pre-registered power criterion (CI
  half-widths 0.0006–0.0016, all ≪ 0.01) but small in absolute terms: the 2 lost
  true SNPs are 0.07% of the SNP set, and a region with a different PSV burden
  could show a different ratio of suppressed-FP to lost-TP.
* The 36:2 favourable ratio is a property of *this* segdup content. It should not
  be extrapolated to segdup genome-wide.
* Speedups are projections from measured per-locus throughput of the two callers,
  not end-to-end pipeline wall-clock; I/O and pile-up construction are shared by
  all arms and excluded from both numerator and denominator.

## 19. Final verdict

Applying §7 mechanically, with no criterion changed after the fact:

* Every powered cell PRESERVED or IMPROVED; none DEGRADED. ✓
* Pooled result IMPROVED. ✓
* PB compute fraction 0.0121% – 0.1303%, all ≤ 1%. ✓
* Projected speedup 209.8× – 335.8×, all ≥ 50×. ✓
* No new failure mode: the only 2 losses are the known devlog 14 segdup mode. ✓
* Disagreement rate 5.14 × 10⁻⁶ < the 10⁻⁴ falsification bound. ✓
* Controls: frozen router beats matched-budget random and reverse-confidence
  routing at every depth. ✓

### **STRONG POSITIVE**

The frozen binomial → router → PB cascade preserves PB-level SNP-calling
accuracy on the hardest unseen region this project has evaluated — a
segdup-saturated 3 Mb window on chr17 where PB-only F1 itself drops to 0.95 —
while routing between 0.012% and 0.13% of loci to PB, a 210–336× projected
compute reduction. H1, H2, H3 and H4 all hold. Within the stated limitations, the
answer to the question this experiment was run to settle is yes.

## 20. Future work (recorded, not implemented)

Nothing here was done and nothing here changes the frozen system.

1. The 36 suppressed PB false positives are a clean, labelled PSV set. Whether PB
   is systematically miscalibrated on paralogous pile-ups is answerable directly
   from them, and is a question about PB, not about the router.
2. The 2 lost SNPs sit at binomial LLR −0.70 and +1.57 — *inside* the ±5.41 band
   in magnitude terms but on the confident side of ε. A depth- or
   annotation-conditioned ε (the devlog 15 safety layer) is the obvious candidate
   to recover them, and this region is the natural place to test it. That is a
   separate pre-registered experiment, not an edit to this one.
3. The truth-gap in the excluded 100%-segdup windows argues for an orthogonal
   validation source (long reads, trio concordance) rather than a wider GIAB BED.
