# Does the frozen cheap-router → PB cascade leave chr21? — cross-chromosome robustness benchmark

**Date:** 2026-08-17

**Status of sections:** §1–§7 are the pre-registration. They were written and
saved **before any caller, router or metric was run on any chr20/chr19/chr4/chr1
locus**. Everything from §8 onwards was written after the corresponding
measurement. The boundary is marked in the text and is checkable in the file's
git history (the file is uncommitted, so "checkable" here means: §1–§7 were not
edited afterwards, and any edit that did happen is disclosed in §20).

---

## 1. Research question

Every accuracy claim this project has made about the frozen cascade

```
binomial LLR → route if |binomial_LLR − 7.0| ≤ 5.411872376933351 → PB → call
```

rests on **one chromosome**: chr21, HG002, Illumina 2×250. Devlog 11 validated
it on 19.0 Mb of chr21. Devlog 12 benchmarked a newly extracted chr21 region.
Devlog 13 repaired the PB depth-clipping bug and re-validated across 24.6 M
chr21 loci at 14x–127x. Devlog 13 §23 closes by listing
"cross-chromosome / cross-sample / cross-platform generalisation" as **NOT
TESTED**, and §24 puts "any chromosome other than 21" first on the untested
list.

This experiment tests exactly that, and nothing else:

1. Does the frozen router preserve PB-level accuracy on chromosomes it has never
   seen?
2. Does it still eliminate almost all PB compute there?
3. Does that hold across genuinely different genomic structure — GC-rich and
   Alu-dense sequence, AT-rich gene-poor sequence, and segmental-duplication /
   low-mappability sequence — rather than only in the ordinary euchromatin the
   chr21 benchmarks were dominated by?
4. Where does it break?

The purpose is to find the boundary, not to pass. A negative result here is
worth more than another chr21 confirmation, because the architecture is already
saturated on chr21 (devlog 13 §25: "six regions and 19.4 M loci at ΔF1 = 0 is a
saturated question").

## 2. Hypotheses

* **H1 (transfer).** The frozen cascade's ΔF1 versus PB-only is
  indistinguishable from 0 on all three ordinary-structure, non-chr21 regions,
  at all three depths. Rationale: the router's input is a binomial LLR, which is
  a function of read counts and a fixed ε — it carries no chromosome-specific
  information, so there is no obvious mechanism by which chromosome identity
  alone could break it.
* **H2 (compute).** The routed fraction stays within roughly an order of
  magnitude of the chr21 range (0.0013%–0.17% of loci), and rises with GC
  content and repeat density, because both raise the local error rate and so
  push more loci into the uncertain band.
* **H3 (the segdup region is where it breaks).** The chr1 segmental-duplication
  window is the pre-registered failure candidate. In segdup sequence, reads from
  a paralogous copy mismap into the locus and produce *systematically* wrong
  non-reference bases at high base quality and high mapping quality. Both
  callers are wrong there, but they are wrong *differently*: the binomial's
  fixed ε = 0.01 is pessimistic and suppresses such evidence, while PB trusts
  per-read quality and will call it. Prediction: the routed fraction rises
  sharply, the disagreement rate rises by orders of magnitude, and — unlike
  every chr21 region — some disagreements change a true call.
* **H4 (the cutoff is regime-dependent, or it is not).** If the frozen cutoff
  5.412 is a property of the *statistic* rather than of chr21, the
  accuracy-vs-compute knee should sit near it on every new region. If the knee
  moves substantially on the GC-rich or segdup regions, the cutoff is a chr21
  artefact. This is diagnostic and **cannot** change the cutoff in this
  experiment either way.

H3 is a prediction of failure. If it is wrong, that is recorded as such.

## 3. Pre-registration

### 3.1 Frozen configuration — not tuned, not refit, not swept for production

Identical to devlog 13 §3.1. Nothing below is re-derived, re-fit or re-selected.

| Constant | Value | Source |
|---|---|---|
| `FROZEN_BINOMIAL_THRESHOLD` | 7.0 | imported from `cascade.py` |
| `FROZEN_PB_THRESHOLD` | 10.5 | imported from `cascade.py` |
| `FROZEN_ROUTER_CUTOFF` | 5.411872376933351 | imported from `cascade.py` |
| binomial ε | 0.01 (`BinomialVariantCaller(error_rate=0.01)`) | unchanged |
| PB ε floor / ceiling | 1e-4 / 0.25 | `quality_error_model` |
| `MAX_READS` | 48 | `read_level_pileup` |
| PB implementation | the **repaired** one (devlog 13 §9) | `quality_error_model.poisson_binomial_llr` |

Routing rule, verbatim: `|binomial_llr − 7.0| ≤ 5.411872376933351` → PB,
otherwise the binomial answer. `cascade.py` and `cheap_router.py` are **not
edited in this session**; verified by `git diff` in §19.

### 3.2 Primary endpoint

**ΔF1 = F1(frozen router → PB) − F1(PB-only)**, per region/depth and pooled,
with a paired bootstrap 95% CI (`evaluate_quality_error.paired_bootstrap`,
10,000 resamples, module-default seed). Same estimator, same seed, same code
path as devlogs 11–13.

### 3.3 Secondary / engineering endpoints

Identical list to devlog 13 §3.3: disagreement count and rate (an engineering
statement about *reproducing PB*, never used as an accuracy proxy or vice
versa); precision, recall, TP/FP/FN/TN, accuracy, balanced accuracy; routed
fraction; measured caller throughput and the projected cascade vs PB-only
wall-clock; the diagnostic accuracy-vs-compute curve; depth-, quality- and
VAF-stratified metrics.

New to this devlog, because the regions are new:

* **Structure-stratified metrics.** Every locus is additionally labelled by
  three external GIAB v3.1 annotations — `lowmap_segdup`, `tandemrepeats`,
  `alldifficult` — and all arms are scored inside and outside each. These
  annotations are public BED files fixed before the experiment; they are not
  derived from any measurement made here.
* **Per-chromosome pooling**, so a single chromosome cannot be hidden inside a
  pooled average.

### 3.4 Decision rule, fixed in advance

**Unchanged from devlog 13 §3.4**, quoted so it cannot drift:

Per region/depth cell with ≥ 100 SNP, and for each pooled result:

* **PRESERVED** — the paired-bootstrap 95% CI for ΔF1 contains 0 **and**
  |ΔF1| < 0.001.
* **DEGRADED** — the CI excludes 0 with ΔF1 < 0, or |ΔF1| ≥ 0.001.
* **IMPROVED** — the CI excludes 0 with ΔF1 > 0.
* **UNDERPOWERED** — cell has < 100 SNP, or the CI half-width exceeds 0.01.
  Reported as underpowered, never rounded to "preserved".

The overall verdict is CONFIRMED only if every powered cell and every pooled
result is PRESERVED. A single powered DEGRADED cell makes the verdict NEGATIVE
or MIXED and is reported as such.

Verdict vocabulary for §21: STRONG POSITIVE / POSITIVE / NULL–INCONCLUSIVE /
NEGATIVE / STRONG NEGATIVE. Mapping, fixed now:

* **POSITIVE** — every powered cell PRESERVED, pooled PRESERVED, and the
  compute reduction is within an order of magnitude of chr21's.
* **STRONG POSITIVE** — the above *and* at least one powered cell in the
  deliberately hard (chr1 segdup) region is PRESERVED.
* **NULL / INCONCLUSIVE** — no powered cell DEGRADED but too few cells are
  powered to support a transfer claim.
* **NEGATIVE** — at least one powered cell DEGRADED.
* **STRONG NEGATIVE** — the pooled result DEGRADED, or the cascade loses true
  SNP that PB-only finds in more than one region.

### 3.5 Stopping criteria

Extraction runs **once** per pre-registered region/depth cell. Metrics are
computed **once**. No region is added, removed, re-extracted or re-scored after
any result is read, except to fix an outright software crash — any such re-run
is disclosed in §20. No threshold, cutoff or parameter is changed on the basis
of any number produced today. If the frozen cutoff fails, the failure is
reported and the cutoff stays at 5.411872376933351.

The experiment stops when: all 12 cells are extracted and scored; all three arms
and all metrics of §3.2–3.3 are computed; controls are run or documented as
underpowered; the disagreement and failure-boundary analyses are complete;
compute accounting is verified; the test suite passes; and §19 records
reproducibility. No further optimisation cycle is run regardless of outcome.

### 3.6 What would falsify H1

Any powered cell where the ΔF1 CI excludes 0 on the negative side, or where the
disagreement rate exceeds ~10⁻⁴ (chr21 saw ~6×10⁻⁷ over 24.6 M loci). Both are
recorded whatever they show.

## 4. Comparison arms

Three arms on the identical scoring frame (`label ∈ {0, SNP}`; indel and
no-call loci excluded, as in every devlog since 9):

* **A. Binomial-only** — `binomial_llr ≥ 7.0`. Establishes what the cheap caller
  costs in accuracy on its own.
* **B. PB-only** — `pb_llr ≥ 10.5`. The reference arm.
* **C. Frozen router → PB** — PB where `|binomial_llr − 7.0| ≤ 5.412`, binomial
  elsewhere.

Primary comparison **C vs B**.

### 4.1 Controls

Run per cell, reusing `robustness_benchmark.controls` unchanged:

* **Matched-coverage random routing** — route a uniformly random subset of loci
  of exactly the same size as the frozen router's routed set, 5 seeds. Tests
  whether the router's benefit is structure or merely PB budget.
* **Reverse-confidence routing** — route the loci *furthest* from the binomial
  threshold, same budget. Tests whether the confidence ordering is the active
  ingredient.

Where a cell's routed set or SNP count is too small for these to separate, that
is reported as underpowered rather than dressed up.

## 5. Dataset, provenance and region selection

### 5.1 Independence from every previous experiment

The frozen cutoff was selected on chr21:32–40 Mb. Every region ever used by this
project — devlogs 3 through 13 — lies on **chr21**. The four regions below lie
on chr20, chr19, chr4 and chr1. Coordinate overlap with any previous training,
validation, tuning or test region is therefore **impossible**, not merely
unlikely: they are different chromosomes. §7 records the mechanical audit.

The sample (HG002) and platform (Illumina 2×250, novoalign) are deliberately
held constant, so that "chromosome and genomic structure" is the only variable
that changes. Cross-sample and cross-platform transfer remain untested and no
claim will be made about them.

### 5.2 Region selection rule — mechanical and results-blind

Implemented in `select_v14_regions.py`, which reads only the GIAB
high-confidence BED, the GIAB v3.1 stratification BEDs and the reference FASTA.
It cannot see a caller score, and it was run to completion before any caller
was run. Rules are quoted in that file's docstring; in summary, over 1 Mb
windows at whole-Mb offsets:

* **chr20 (neutral)** — closest to the chromosome midpoint among windows with
  HC ≥ 0.90, difficult ≤ 0.15, no N.
* **chr19 (GC-rich)** — highest GC among windows with HC ≥ 0.90 and no N.
* **chr4 (AT-rich)** — lowest GC among windows with HC ≥ 0.90, difficult ≤ 0.15,
  no N.
* **chr1 (hard)** — maximum fraction of high-confidence bases inside
  `GRCh38_alllowmapandsegdupregions`, among windows with HC ≥ 0.50 and no N.

**Two rule revisions are disclosed here**, both made before any caller ran, with
only annotation statistics in existence at the time:

1. The original rule was "closest to the chromosome midpoint" for all three
   ordinary contigs. It placed the chr19 window at 28–29 Mb, which is
   pericentromeric at GC 43.0% — the opposite of the GC-rich regime chr19 was
   selected to represent. Replaced by explicit GC-extremum rules for chr19 and
   chr4, keeping the midpoint rule for chr20 as the neutral control.
2. The GC-rich rule initially inherited the `difficult ≤ 0.15` cap, which
   returned a 45.3% GC window — *lower* GC than chr20's neutral window, because
   GC-rich human sequence is Alu-dense by construction and is therefore
   annotated difficult. The cap was dropped for that arm alone.

Neither revision could have been informed by a result, because no result
existed. Both are recorded because pre-registration hygiene is worth more than
a tidy narrative.

### 5.3 The selected regions

| Key | Region (GRCh38) | Regime | GC | HC cov. | difficult | segdup∩HC | TR∩HC |
|---|---|---|---:|---:|---:|---:|---:|
| `chr20_neutral` | chr20:33,000,000–34,000,000 | ordinary euchromatin, neutral control | 0.477 | 0.961 | 0.130 | 0.013 | 0.062 |
| `chr19_gcrich` | chr19:1,000,000–2,000,000 | GC-rich, Alu/repeat-dense subtelomeric | 0.591 | 0.909 | 0.487 | 0.018 | 0.083 |
| `chr4_atrich` | chr4:101,000,000–102,000,000 | AT-rich, gene-poor isochore | 0.363 | 0.974 | 0.150 | 0.024 | 0.046 |
| `chr1_segdup` | chr1:120,000,000–121,000,000 | pericentromeric segmental duplication / low mappability | 0.406 | 0.722 | 0.998 | 1.000 | 0.020 |

GC spans 0.363–0.591 — a genuine isochore contrast, and wider than anything
chr21 offered. `chr1_segdup` has **100%** of its high-confidence bases inside
the low-mappability/segdup annotation, versus 1.3–2.4% for the ordinary
regions: this is a categorically different mapping regime, not a gradient.

### 5.4 Depth realisations

Each region is extracted at three depths: **full** (native ~70x), **30x** and
**15x**, the latter two by `samtools view -s` with a fixed seed and a fraction
derived from the region's measured native depth. As in devlog 13 §5, depth
realisations of one region are downsamples of **one library** and are treated as
three regimes of one region, never as three independent regions.

12 cells total: 4 regions × 3 depths.

### 5.5 Provenance

* Sample: GIAB **HG002** (NA24385 son), GRCh38.
* Reads: `HG002.GRCh38.2x250.bam`, NIST Illumina 2×250 novoalign, sliced
  remotely per region with `samtools view -b <URL> <contig>:<start+1>-<stop>`.
  Same URL and same recipe as `data/giab_hg002_real_30M/PROVENANCE.txt`.
* Truth: GIAB **v4.2.1** benchmark VCF
  (`HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`), region-sliced with `bcftools
  view -r`. Same release used by every previous devlog.
* High-confidence BED: `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`,
  filtered to the region's contig and span.
* Reference: Ensembl release-110 GRCh38 per-chromosome FASTA (contigs `20`,
  `19`, `4`, `1`), the same source and release as the existing
  `data/reference/chr21_full.fa`.
* Difficulty annotation: GIAB genome-stratifications **v3.1**, GRCh38 —
  `Union/GRCh38_alllowmapandsegdupregions.bed.gz`,
  `Union/GRCh38_alldifficultregions.bed.gz`,
  `LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz`.
* Retrieved 2026-08-17. Exact commands: `fetch_v14_data.sh`.

Nothing is synthesised. No artificial difficulty is injected. No read is
modified.

### 5.6 Regimes: what is and is not available

| Requested regime | Available? | What is used |
|---|---|---|
| ordinary depth | yes | 30x realisations of all four regions |
| low depth | yes | 15x realisations; plus the `depth ≤ 9` stratum everywhere |
| high depth | yes | native ~70x realisations (PB repaired in devlog 13, so this is now measurable) |
| low VAF | as a stratum | `0 < VAF ≤ 0.35` |
| low base quality | as a stratum | mean base quality < 30 |
| high sequencing-error burden | **proxy only** | no elevated-error BAM exists for HG002 at this platform; the GC-rich and segdup regions are the strongest *naturally occurring* high-error proxies, and are labelled as proxies. Unlike devlog 13, at least the proxy is now a genuine mapping-error regime rather than only a quality stratum. |
| difficult / local genomic regions | **yes, for the first time** | `chr1_segdup`, 100% of HC bases in segdup/low-mappability; plus the GIAB difficulty strata applied to every region |
| mixed-complexity regions | yes | `chr19_gcrich` (48.7% difficult) and whole-region pools |
| different chromosomes | **yes, for the first time** | chr20, chr19, chr4, chr1 |
| different sample | no | HG002 only — held constant deliberately |
| different platform | no | Illumina 2×250 only — held constant deliberately |
| indels / MNP / SV | no | SNP-only label set, as since devlog 9 |

## 6. Bug status

The PB depth-clipping defect (devlog 13 §8–§12) was diagnosed and fixed in the
current working tree, with 12 regression tests. **This experiment does not
re-open it.** §8 re-verifies that the fix is still present and still passing
before any high-depth data is used, because this benchmark's native-depth cells
are ~70x and would otherwise be invalid. If that verification fails, the
high-depth cells are dropped rather than reported.

## 7. Execution and leakage audit plan

1. Verify the PB fix is present and the existing suite passes (§8).
2. Extract 12 cells once each with `extract_bench_v12.py --contig <c>`.
3. Audit independence mechanically: assert every extracted region's contig ∉
   {chr21} and that no coordinate span intersects any previously used region
   (§9).
4. Score all cells, all three arms, all strata; run controls.
5. Accuracy-vs-compute curve (diagnostic only).
6. Disagreement dump and failure-boundary analysis.
7. Full test suite, artefact validation, git state.

### 7.1 Artefacts

| What | Path |
|---|---|
| Region selection | `select_v14_regions.py`, `results/bench_v14/region_selection.json` |
| Data acquisition | `fetch_v14_data.sh`, `prepare_v14_depths.sh` |
| Benchmark driver | `bench_v14_crosschrom.py` |
| Driver unit tests | `test_bench_v14.py` |
| Extraction caches | `cache/bench_v14/*.npz` |
| Primary results | `results/bench_v14/robustness_results.json` |
| Disagreement dump | `results/bench_v14/disagreements.json` |
| Leakage audit | `results/bench_v14/independence_audit.json` |

No file under `results/robustness_v11/`, `results/bench_v12/`,
`results/bench_v13/` or `History/` (other than this new file) is modified.

---

*Everything from here down was written after the corresponding measurement was
made. Sections 1–7 above were not edited afterwards, except as disclosed in §20.*

---

## 8. Execution plan and pre-flight validation (interrupted-run recovery)

### 8.1 What was interrupted, and what state it left

The first launch of `run_v14_extract.sh` was stopped by hand ~3 minutes in,
during batch 1 of 3 (the four native-depth cells). State audit at the moment of
the stop:

| Check | Result |
|---|---|
| Extractor processes alive afterwards | **none** (`pkill` confirmed; `rc=143` recorded per cell in `logs_v14/extract_driver.log`) |
| Files in `cache/bench_v14/` | **none** — zero bytes written |
| Partially written `.npz` | **impossible by construction**: `extract_bench_v12.main` accumulates in memory and calls `np.savez_compressed` exactly once, after the final chunk. A killed run leaves no output file, not a truncated one. |
| Historical artefacts touched | none (`results/robustness_v11`, `results/bench_v12`, `results/bench_v13`, `History/3–13`) |
| Downloaded inputs | complete and intact — see §8.2 |

So there is nothing to salvage and nothing to repair: the extraction is
**restartable from scratch**, and no partial artefact can silently contaminate
a result. Nothing was deleted to reach this state.

### 8.2 Pre-flight validation of the inputs (all passed before relaunch)

| Check | Result |
|---|---|
| 12 BAM cells present, indexed, `samtools quickcheck` | **12/12 OK** |
| Read counts | native 194k–300k; 30x 122k–125k; 15x 61k–63k reads per Mb |
| Achieved mean depth (`samtools depth -a`) | 14.9–15.1x, 29.8–30.2x, native 47.6–73.7x |
| Truth VCF REF vs Ensembl FASTA | **5,181 / 5,181 records match, 0 mismatches** across all four contigs |
| Contig-name resolution across FASTA/BAM/VCF | FASTA `20`/`19`/`4`/`1` vs VCF+BAM `chr20`/`chr19`/`chr4`/`chr1`; resolved automatically by `PileupCountsProvider._normalise_contig` |
| Region window alignment | 1,000,000 bp / `SEQ_LEN` 64 = 15,625 windows exactly, for every cell |
| End-to-end smoke extraction, 64 kb per contig, new `--contig` path | **4/4 wrote a valid cache**; 0 non-finite LLRs; depths as intended; LLR ranges non-degenerate (binomial −32.8…+113.9, PB −29.1…+204.1) |

The smoke caches were written to a scratch directory outside the experiment
namespace and are not benchmark artefacts.

One structural observation worth recording in advance, because it is a property
of the *data* and not of any result: `chr1_segdup` has **8.8%** zero-depth loci
in the smoke slice, versus 0.0–0.2% for the other three. That is the expected
signature of segmental duplication — reads that cannot be placed uniquely are
simply absent — and it confirms the region is the mapping regime it was selected
to be.

### 8.3 Code change required to reach a second chromosome, and its blast radius

`pileup_counts.load_counts` and `read_level_pileup.load_reads` hardcoded
`contig="chr21"`. Both now take `contig: str = "chr21"` and pass it through;
`extract_bench_v12.py` gained `--contig` (same default) and records the contig in
the cache. `PileupCountsProvider` already supported arbitrary contigs, so no
provider logic changed. The default preserves every previous invocation
byte-for-byte, and the full suite passes unchanged (§18). **No router file, no
threshold and no caller formula was touched.**

### 8.4 Extraction command (exactly what runs)

`run_v14_extract.sh`, three batches of four concurrent jobs (native, then 30x,
then 15x), each cell run exactly once:

```
python3 extract_bench_v12.py \
  --fasta data/reference/<contig>_full.fa \
  --bam   data/giab_hg002_v14/<name>_<tag>.bam \
  --vcf   data/giab_hg002_v14/<name>.vcf.gz \
  --bed   data/giab_hg002_v14/<name>_highconf.bed \
  --contig <contig> --region <start> <stop> \
  --features none \
  --out   cache/bench_v14/<name>_<tag>.npz
```

`--features none` because this benchmark has no neural arm. No `--legacy-pb`:
the PB depth-clipping A/B is settled in devlog 13, and every cell here uses the
repaired PB and only the repaired PB.

### 8.5 Expected artefacts and completion criteria

12 files `cache/bench_v14/{chr20_neutral,chr19_gcrich,chr4_atrich,chr1_segdup}_{full,30x,15x}.npz`.

Extraction is declared complete only when, for all 12: the file exists; its
`labels` array has exactly 1,000,000 / 64 × 64 = 15,625 × 64 = 1,000,000 entries
minus the windows dropped as uncallable; `labels.size % 64 == 0`; the stored
`contig` and `region` match the pre-registered values; and no LLR is non-finite.
These are asserted mechanically by `verify_v14_caches` before any metric is
computed, and the audit is written to
`results/bench_v14/independence_audit.json`.

### 8.6 How much of the data the benchmark actually uses

The scoring frame is `label ∈ {0, SNP}` — indel and no-call loci are excluded, as
in every devlog since 9. Loci outside the high-confidence BED carry no truth and
are excluded upstream by the provider. So the benchmark uses fewer loci than are
extracted, and the exact scored count per cell is reported in §11 rather than
assumed here.

### 8.7 Stopping criteria (restated, unchanged from §3.5)

Twelve cells, extracted once, scored once. No cell is added, dropped,
re-extracted or re-scored after a result is read, except to fix an outright
crash, which would be disclosed. No threshold is changed on the basis of any
number produced today. The experiment ends when §3.5's checklist is satisfied,
whatever the verdict.

## 9. Extraction — what was produced

All 12 cells extracted once each, 09:53–10:25 (32 min wall-clock, four
concurrent jobs on 8 cores). No cell was re-run.

| Region | loci extracted | SNP | native mean depth | 30x | 15x |
|---|---:|---:|---:|---:|---:|
| `chr20_neutral` | 947,904 | 765 | 70.2x | 30.04x | 14.99x |
| `chr19_gcrich` | 886,848 | 1,083 | 66.2x | 29.75x | 14.92x |
| `chr4_atrich` | 966,016 | 1,111 | 73.7x | 29.99x | 15.07x |
| `chr1_segdup` | 656,512 | 1,128 | 47.6x | 30.16x | 15.10x |

Locus counts are identical across the three depths of a region, as they must be:
window tiling depends on the FASTA and the high-confidence BED, not on the reads.
`chr1_segdup` yields 31% fewer loci per Mb than `chr20_neutral` because its
high-confidence coverage is 72.2% rather than 96.1% — that is the segdup regime
showing up in callable territory, not a defect.

**Total scored: 10,368,078 loci, 12,261 SNP** across 4 chromosomes, 4 Mb of
distinct genomic sequence, 3 depth regimes.

### 9.1 Independence audit (`results/bench_v14/independence_audit.json`)

All 12 cells pass, mechanically: stored contig equals the pre-registered contig;
stored coordinates equal the pre-registered span; locus count is window-aligned
and equals the count obtained by independently rebuilding the window tiling;
no non-finite LLR; every reconstructed coordinate inside the region; and **no
cell on chr21 or any other previously used contig**. Because the four contigs
are disjoint from chr21, coordinate overlap with any prior training, tuning or
test region is impossible rather than merely improbable.

Coordinate reconstruction was separately validated against the truth VCF
(`test_bench_v14.py::test_reconstructed_coordinates_land_on_real_vcf_variants`):
every SNP-labelled locus sits exactly on a VCF SNV coordinate and every tiled
VCF SNV is labelled — zero misplacements in either direction, on all four
contigs. This matters because the structure strata and the audit both rest on
these coordinates, and the caches do not store them.

## 10. Main results — frozen cascade vs PB-only

Full table in `results/bench_v14/summary.md`; artefact in
`results/bench_v14/robustness_results.json`.

| Cell | depth | loci | SNP | F1 PB-only | F1 cascade | ΔF1 | ΔF1 95% CI | routed→PB | disagree | verdict |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---|
| `chr20_neutral` | full | 947,626 | 765 | 0.983290 | 0.982028 | −0.001262 | [−0.003849, +0.001228] | 0.0053% | 4 | **DEGRADED** |
| `chr20_neutral` | 30x | 947,626 | 765 | 0.987718 | 0.987718 | +0.000000 | [0, 0] | 0.0325% | 0 | PRESERVED |
| `chr20_neutral` | 15x | 947,626 | 765 | 0.984375 | 0.984375 | +0.000000 | [0, 0] | 0.1743% | 0 | PRESERVED |
| `chr19_gcrich` | full | 886,344 | 1,083 | 0.972997 | 0.972122 | −0.000875 | [−0.002668, +0.000856] | 0.0255% | 4 | PRESERVED |
| `chr19_gcrich` | 30x | 886,344 | 1,083 | 0.983097 | 0.983097 | +0.000000 | [0, 0] | 0.0599% | 0 | PRESERVED |
| `chr19_gcrich` | 15x | 886,344 | 1,083 | 0.974860 | 0.974407 | −0.000454 | [−0.001415, 0] | 0.2184% | 1 | PRESERVED |
| `chr4_atrich` | full | 965,687 | 1,111 | 0.996413 | 0.995966 | −0.000447 | [−0.001393, 0] | 0.0030% | 1 | PRESERVED |
| `chr4_atrich` | 30x | 965,687 | 1,111 | 0.998203 | 0.998203 | +0.000000 | [0, 0] | 0.0133% | 0 | PRESERVED |
| `chr4_atrich` | 15x | 965,687 | 1,111 | 0.987799 | 0.987799 | +0.000000 | [0, 0] | 0.0936% | 0 | PRESERVED |
| `chr1_segdup` | full | 656,369 | 1,128 | 0.882949 | 0.888167 | +0.005218 | [+0.001418, +0.008972] | 0.0507% | 19 | **IMPROVED** |
| `chr1_segdup` | 30x | 656,369 | 1,128 | 0.875817 | 0.878873 | +0.003056 | [+0.000211, +0.006140] | 0.0597% | 12 | **IMPROVED** |
| `chr1_segdup` | 15x | 656,369 | 1,128 | 0.858105 | 0.858105 | +0.000000 | [0, 0] | 0.1233% | 0 | PRESERVED |

Every cell is **powered** by the §3.4 rule: all have ≥ 765 SNP (threshold 100),
and the widest CI half-width is 0.00254 (threshold 0.01). So no cell escapes
into UNDERPOWERED, and 12/12 cells carry a real verdict.

### 10.1 Pooled

| Pool | loci | SNP | F1 PB-only | F1 cascade | ΔF1 | ΔF1 95% CI | routed→PB | disagree | verdict |
|---|---:|---:|---:|---:|---:|---|---:|---:|---|
| **all** | 10,368,078 | 12,261 | 0.955726 | 0.956306 | +0.000581 | [+0.000075, +0.001089] | 0.0704% | 41 | IMPROVED |
| `chr20_neutral` | 2,842,878 | 2,295 | 0.985126 | 0.984702 | −0.000425 | [−0.001284, +0.000417] | 0.0707% | 4 | PRESERVED |
| `chr19_gcrich` | 2,659,032 | 3,249 | 0.976978 | 0.976532 | −0.000447 | [−0.001186, +0.000150] | 0.1013% | 5 | PRESERVED |
| `chr4_atrich` | 2,897,061 | 3,333 | 0.994152 | 0.994003 | −0.000149 | [−0.000457, 0] | 0.0366% | 1 | PRESERVED |
| `chr1_segdup` | 1,969,107 | 3,384 | 0.872477 | 0.875236 | +0.002759 | [+0.001182, +0.004402] | 0.0779% | 31 | IMPROVED |
| depth full | 3,456,026 | 4,087 | 0.957447 | 0.958374 | +0.000927 | [−0.000304, +0.002146] | 0.0185% | 28 | PRESERVED |
| depth 30x | 3,456,026 | 4,087 | 0.959773 | 0.960702 | +0.000929 | [+0.000180, +0.001764] | 0.0393% | 12 | IMPROVED |
| depth 15x | 3,456,026 | 4,087 | 0.949850 | 0.949730 | −0.000119 | [−0.000364, 0] | 0.1534% | 1 | PRESERVED |
| ordinary structure only (no segdup) | 8,398,971 | 8,877 | 0.985504 | 0.985173 | −0.000331 | [−0.000672, 0] | 0.0686% | 10 | PRESERVED |

**Observation, not interpretation:** no pooled result is DEGRADED. One cell of
twelve is. The pooled `all` result is positive and its CI excludes zero, driven
entirely by `chr1_segdup`.

### 10.2 The one DEGRADED cell, examined

`chr20_neutral_full` fails the pre-registered rule, and it is worth being exact
about *which clause* fails. Its bootstrap CI is [−0.003849, +0.001228] — it
**contains zero**, so the difference is not statistically significant. What fires
is the second clause of the DEGRADED definition: |ΔF1| ≥ 0.001.

The physical event behind ΔF1 = −0.001262 is: PB-only makes 26 false positives,
the cascade makes 28, recall is 1.000 for both (765/765 TP). **Two false
positives.** At 765 SNP, one FP moves F1 by ≈ 0.00065, so the pre-registered
0.001 equivalence bound sits at roughly 1.5 false calls — below the granularity
at which this cell can resolve anything.

That is a limitation of the rule as pre-registered for a region this size, and it
is recorded as such in §17 and §19. It does **not** license changing the verdict:
the rule was fixed in advance, it fires, and §21 applies it as written.

## 11. All three arms

Excerpt; full table in `results/bench_v14/summary.md`. The binomial-only arm is
included to price the cheap caller on its own.

| Cell | arm | F1 | precision | recall | TP | FP | FN |
|---|---|---:|---:|---:|---:|---:|---:|
| `chr20_neutral` full | binomial | 0.980141 | 0.961055 | 1.000000 | 765 | 31 | 0 |
| `chr20_neutral` full | PB | 0.983290 | 0.967130 | 1.000000 | 765 | 26 | 0 |
| `chr20_neutral` full | cascade | 0.982028 | 0.964691 | 1.000000 | 765 | 28 | 0 |
| `chr20_neutral` 30x | binomial | 0.973835 | 0.951372 | 0.997386 | 763 | 39 | 2 |
| `chr20_neutral` 30x | PB | 0.987718 | 0.976982 | 0.998693 | 764 | 18 | 1 |
| `chr20_neutral` 30x | cascade | 0.987718 | 0.976982 | 0.998693 | 764 | 18 | 1 |

Accuracy and balanced accuracy are reported per arm per cell in the artefact.
They are uninformative at this class balance (12,261 SNP in 10.37 M loci, so
accuracy is ≥ 0.9999 for every arm everywhere) and no claim rests on them —
which is exactly why the pre-registration made F1 the endpoint and required
TP/FP/FN beside it.

## 12. Compute and wall-clock

Throughput measured in a quiet single process on the v14 BAMs after every
extraction had exited:

| Regime | binomial | PB | ratio |
|---|---:|---:|---:|
| ~70x | 2,579,314 loci/s | 7,760 loci/s | 332x |
| 30x | 3,118,544 loci/s | 7,462 loci/s | 418x |
| 15x | 2,651,321 loci/s | 9,688 loci/s | 274x |

| Cell | loci | PB loci | PB compute fraction | PB-only s | cascade s | speedup |
|---|---:|---:|---:|---:|---:|---:|
| `chr20_neutral` full | 947,626 | 50 | 0.0053% | 122.1 | 0.4 | **327x** |
| `chr20_neutral` 30x | 947,626 | 308 | 0.0325% | 127.0 | 0.3 | 368x |
| `chr20_neutral` 15x | 947,626 | 1,652 | 0.1743% | 97.8 | 0.5 | 185x |
| `chr19_gcrich` full | 886,344 | 226 | 0.0255% | 114.2 | 0.4 | 306x |
| `chr19_gcrich` 30x | 886,344 | 531 | 0.0599% | 118.8 | 0.4 | 334x |
| `chr19_gcrich` 15x | 886,344 | 1,936 | 0.2184% | 91.5 | 0.5 | 171x |
| `chr4_atrich` full | 965,687 | 29 | 0.0030% | 124.4 | 0.4 | 329x |
| `chr4_atrich` 30x | 965,687 | 128 | 0.0133% | 129.4 | 0.3 | 396x |
| `chr4_atrich` 15x | 965,687 | 904 | 0.0936% | 99.7 | 0.5 | 218x |
| `chr1_segdup` full | 656,369 | 333 | 0.0507% | 84.6 | 0.3 | 284x |
| `chr1_segdup` 30x | 656,369 | 392 | 0.0597% | 88.0 | 0.3 | 334x |
| `chr1_segdup` 15x | 656,369 | 809 | 0.1233% | 67.7 | 0.3 | 205x |

**PB compute fraction: 0.0030%–0.2184%. Speedup: 171x–396x.** chr21's range was
0.0013%–0.17% and 171x–369x. The compute claim transfers to four new
chromosomes essentially unchanged — the routed fraction on the *hardest* region
(`chr1_segdup`, 0.05–0.12%) is no worse than on the easiest chr21 regions.

The routed fraction rises as depth falls (0.003% at 70x → 0.09–0.22% at 15x),
which is the expected direction: shallower pileups leave the binomial less
certain more often. It also rises with GC content and repeat density
(`chr19_gcrich` is highest at every depth), as H2 predicted.

Only caller compute is measured. BAM I/O and pileup construction are not
included in the speedup and are paid by both arms; extraction timings in the
caches were taken with four concurrent jobs and are upper bounds, not clean
measurements.

## 13. Accuracy-vs-compute (diagnostic only — the cutoff stays frozen)

`chr1_segdup_full`, the hard region:

| cutoff | routed→PB | F1 | ΔF1 vs PB | disagreements | projected speedup |
|---:|---:|---:|---:|---:|---:|
| 2.500 | 0.0102% | 0.887959 | +0.005010 | 23 | 322x |
| 4.000 | 0.0261% | 0.887651 | +0.004702 | 20 | 306x |
| 5.000 | 0.0413% | 0.888167 | +0.005218 | 19 | 292x |
| **5.412 ← frozen** | **0.0507%** | **0.888167** | **+0.005218** | **19** | **284x** |
| 6.000 | 0.0516% | 0.888167 | +0.005218 | 19 | 284x |
| 7.000 | 8.6716% | 0.888683 | +0.005733 | 16 | 11x |
| 10.000 | 12.5987% | 0.886627 | +0.003677 | 11 | 8x |
| 20.000 | 21.9849% | 0.884686 | +0.001737 | 4 | 5x |

The frozen cutoff sits on a **plateau** (5.0–6.0 give bit-identical F1) that ends
in a cliff at 7.0, where the routed fraction jumps from 0.05% to 8.7% and the
speedup collapses from 284x to 11x for +0.0005 F1. The cliff is structural, not
accidental: `|binomial_llr − 7.0| ≤ 7.0` first admits `binomial_llr ≥ 0`, and an
enormous mass of loci with no non-reference read sit at LLR ≈ 0. The frozen
point sits just below that cliff on every region tested.

This answers **H4**: the knee does not move between chr21 and four new
chromosomes, including the GC-rich and segdup regimes. The cutoff behaves like a
property of the statistic, not a chr21 artefact. **No cutoff was changed.**

One negative diagnostic worth stating plainly: on `chr20_neutral_full`, the
DEGRADED cell, **no cutoff between 4.5 and 8.0 removes the −0.001262**. The two
extra false positives come from loci the binomial is *confidently wrong* about
(router margins 6.6, 13.7 and 20.0), so widening the band does not reach them
until it has already surrendered the speedup. Widening the router is not a fix
for this failure mode — which is a reason to report it, not to attempt it.

## 14. Depth-stratified results (pooled over all 12 cells)

| depth bin | loci | SNP | F1 binomial | F1 PB | F1 cascade | ΔF1 | routed→PB | disagree |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1–9x | 580,741 | 1,021 | 0.830285 | 0.877713 | 0.878307 | +0.000594 | 0.8766% | 1 |
| 10–29x | 4,405,610 | 5,186 | 0.964860 | 0.978463 | 0.978271 | −0.000191 | 0.0448% | 2 |
| 30–47x | 1,983,724 | 2,218 | 0.977313 | 0.980209 | 0.980846 | +0.000637 | 0.0077% | 7 |
| 48–59x | 630,754 | 805 | 0.970357 | 0.965144 | 0.966265 | +0.001121 | 0.0070% | 10 |
| 60–79x | 1,748,118 | 1,957 | 0.989623 | 0.988122 | 0.989873 | +0.001751 | 0.0015% | 13 |
| ≥80x | 831,167 | 691 | 0.978030 | 0.971147 | 0.973888 | +0.002741 | 0.0013% | 8 |

Two observations:

1. **Below 48x, PB is the better caller** (0.978 vs 0.965 binomial at 10–29x) and
   the cascade tracks PB. **Above 48x the ordering inverts**: the fixed-ε
   binomial *outscores* the repaired PB (0.9896 vs 0.9881 at 60–79x; 0.978 vs
   0.971 at ≥80x). This independently reproduces devlog 13 §23's finding — PB is
   no longer the accuracy ceiling at high depth — now on four new chromosomes
   rather than on chr21 alone.
2. The cascade is at or above both arms in every bin above 30x. It is not
   "approximating PB" there; it is mixing two callers whose errors are
   anti-correlated, and the mixture beats both.

The disagreement count rises with depth (1–2 events below 30x, 8–13 above 48x)
even though the routed fraction *falls* by two orders of magnitude — the
disagreements are concentrated in the loci the router does **not** send to PB.

## 15. Structure-stratified results (pooled over all 12 cells)

Strata are external GIAB v3.1 annotations, fixed before the experiment.

| stratum | loci | SNP | F1 PB | F1 cascade | ΔF1 | routed→PB | disagree |
|---|---:|---:|---:|---:|---:|---:|---:|
| in `alldifficult` | 3,919,953 | 5,781 | 0.912150 | 0.913236 | +0.001086 | 0.0934% | 40 |
| not in `alldifficult` | 6,448,125 | 6,480 | 0.993676 | 0.993752 | +0.000077 | 0.0564% | 1 |
| in `lowmap_segdup` | 2,119,047 | 3,675 | 0.881653 | 0.884217 | +0.002564 | 0.0817% | 31 |
| not in `lowmap_segdup` | 8,249,031 | 8,586 | 0.985596 | 0.985254 | −0.000342 | 0.0675% | 10 |
| in `tandemrepeats` | 574,524 | 744 | 0.905329 | 0.901373 | **−0.003956** | 0.2560% | 9 |
| not in `tandemrepeats` | 9,793,554 | 11,517 | 0.959273 | 0.960194 | +0.000922 | 0.0595% | 32 |

**40 of 41 disagreements fall inside `alldifficult`.** In easy sequence the
cascade reproduces PB almost exactly (1 disagreement in 6.45 M loci).

The tandem-repeat/homopolymer stratum is the one place where the cascade is
**worse** than PB by more than the equivalence bound (ΔF1 −0.0040 on 744 SNP).
It is also the stratum with the highest routed fraction (0.256%). This is
reported, not fixed.

## 16. Controls

Matched-coverage random routing (5 seeds) and reverse-confidence routing, both
spending exactly the frozen router's PB budget:

| Cell | PB budget | frozen router ΔF1 | random routing ΔF1 | reverse-confidence ΔF1 |
|---|---:|---:|---:|---:|
| `chr20_neutral` full | 50 | −0.001262 | −0.003150 ± 0.000000 | −0.003150 |
| `chr20_neutral` 30x | 308 | +0.000000 | −0.013883 ± 0.000000 | −0.013883 |
| `chr20_neutral` 15x | 1,652 | +0.000000 | −0.022220 ± 0.000000 | −0.022220 |
| `chr19_gcrich` full | 226 | −0.000875 | −0.009077 ± 0.000000 | −0.009077 |
| `chr19_gcrich` 30x | 531 | +0.000000 | −0.020267 ± 0.000000 | −0.020267 |
| `chr19_gcrich` 15x | 1,936 | −0.000454 | −0.036018 ± 0.000172 | −0.036103 |
| `chr4_atrich` full | 29 | −0.000447 | **+0.000447** ± 0.000000 | +0.000447 |
| `chr4_atrich` 30x | 128 | +0.000000 | −0.005359 ± 0.000000 | −0.005359 |
| `chr4_atrich` 15x | 904 | +0.000000 | −0.015752 ± 0.000000 | −0.015752 |
| `chr1_segdup` full | 333 | +0.005218 | +0.002938 ± 0.000000 | +0.002938 |
| `chr1_segdup` 30x | 392 | +0.003056 | +0.002370 ± 0.000000 | +0.002370 |
| `chr1_segdup` 15x | 809 | +0.000000 | +0.000285 ± 0.000000 | +0.000285 |

**Both controls collapse onto the binomial-only arm**, exactly (their ΔF1 equals
F1(binomial) − F1(PB) in every cell, and the random arm has zero variance across
5 seeds). The reason is mechanical and worth stating rather than hiding: with a
budget of 29–1,936 loci out of ~1 M, a random or maximally-confident subset
essentially never contains a locus where the two callers disagree, so routing it
to PB changes no call.

That is precisely the control result the design was after. The frozen router
recovers **all** of PB's advantage — up to the disagreements analysed in §17 —
with 0.003–0.22% of PB's compute, whereas the same budget spent randomly, or
spent on the loci the binomial is *most* sure about, recovers essentially none
of it (e.g. `chr19_gcrich` 15x: frozen −0.00045 vs random −0.036, an 80-fold
difference). The router is exploiting structured disagreement, not merely
spending a PB budget.

Two honest caveats: the controls are *degenerate* rather than merely worse, so
they establish "the routed set is where the callers disagree" and nothing
finer; and on `chr4_atrich` full the random control (+0.000447) beats the frozen
router (−0.000447), because there PB is the worse caller and touching fewer loci
with it is better.

## 17. Disagreement and failure-boundary analysis

41 disagreements in 10,368,078 loci — rate **3.95×10⁻⁶**. chr21's rate over
24.58 M loci was 6.1×10⁻⁷, so the rate is ~6.5x higher on the new data. It
remains far below the 10⁻⁴ rate that §3.6 pre-registered as falsifying.

### 17.1 What the disagreements do

| effect | count | meaning |
|---|---:|---|
| `router_avoids_fp` | 26 | PB called a false positive; the cascade did not |
| `router_fp` | 8 | the cascade called a false positive; PB did not |
| `router_fn` | **5** | **the cascade lost a true SNP that PB found** |
| `router_recovers_fn` | 2 | PB missed a true SNP; the cascade found it |

Net 28 favourable, 13 unfavourable — which is why the pooled ΔF1 is positive.

Distribution across cells:

| cell | disagreements | breakdown |
|---|---:|---|
| `chr20_neutral_full` | 4 | 3 `router_fp`, 1 `router_avoids_fp` |
| `chr19_gcrich_full` | 4 | 3 `router_fp`, 1 `router_avoids_fp` |
| `chr19_gcrich_15x` | 1 | 1 `router_fp` |
| `chr4_atrich_full` | 1 | 1 `router_fp` |
| `chr1_segdup_full` | 19 | 14 `router_avoids_fp`, 3 `router_fn`, 2 `router_recovers_fn` |
| `chr1_segdup_30x` | 12 | 10 `router_avoids_fp`, 2 `router_fn` |
| all other 6 cells | 0 | — |

**All 5 lost true SNP are in `chr1_segdup`**, at two depths of that one region.
On chr21, devlog 13 found zero lost true SNP in 24.58 M loci. This is a genuinely
new failure mode, and it is confined to one structural regime.

### 17.2 The failure boundary, by axis

| axis | bin | loci | disagreements | rate |
|---|---|---:|---:|---:|
| VAF | 0.15–0.25 | 5,677 | 8 | **1.41×10⁻³** |
| VAF | 0.00–0.15 | 566,773 | 31 | 5.47×10⁻⁵ |
| VAF | ≥0.75 | 6,274 | 0 | 0 |
| quality | Q0–20 | 5,400 | 4 | **7.41×10⁻⁴** |
| quality | Q20–25 | 26,126 | 5 | 1.91×10⁻⁴ |
| quality | Q30–35 | 2,145,274 | 0 | 0 |
| quality | Q35–60 | 7,770,003 | 30 | 3.86×10⁻⁶ |
| depth | 48–59x | 630,754 | 10 | 1.59×10⁻⁵ |
| depth | 10–29x | 4,405,610 | 2 | 4.54×10⁻⁷ |
| truth | SNP | 12,261 | 7 | **5.71×10⁻⁴** |
| truth | non-SNP | 10,355,817 | 34 | 3.28×10⁻⁶ |

Disagreement is ~350x more likely at VAF 0.15–0.25 than in the high-VAF bins,
~190x more likely below Q20 than above Q35, ~35x more likely above 48x than at
10–29x, and ~175x more likely at a true SNP than at a reference locus.

### 17.3 The mechanism — one axis explains everything

Every disagreement record was dumped with its evidence
(`results/bench_v14/disagreements.json`). The three populations separate cleanly:

**`router_fp` (8) — the binomial is confidently wrong, positively.**
Signature: `binomial_llr` 12.5 to 27.0 (a confident variant call), `pb_llr` −5.5
to +9.3 (PB says no), VAF 0.15–0.21, and **mean base quality Q15–Q24**. Example:
`chr20_neutral_full`, depth 72, VAF 0.164, Q22.0, binomial +13.61, PB −5.46,
router margin 6.61 → not routed.
*Cause:* the binomial's fixed ε = 0.01 is **optimistic** when the reads are
Q15–Q24 (true error rate 0.004–0.03). It over-calls; PB, which reads the actual
qualities, does not. The router does not intervene because the binomial is
confident — just wrong.

**`router_avoids_fp` (26) and `router_fn` (5) — the binomial is confidently
wrong, negatively, and that is usually right.**
Signature for both: VAF 0.07–0.13, **high base quality Q27–Q39**, `binomial_llr`
−21 to +1.3 (confidently no variant), `pb_llr` +11 to +24 (PB calls it).
*Cause:* this is the PB high-depth false-positive mode identified in devlog 13
§18.1 — a handful of high-quality non-reference reads at low VAF is overwhelming
evidence to a per-read Poisson-binomial with no term for mapping error. In
segmental duplication those reads come from a paralogous copy. The binomial's
pessimism suppresses them. **26 times the binomial is right and the cascade
dodges a PB false positive; 5 times the low-VAF variant is real and the cascade
loses it.**

So the failure boundary is a **single axis**, not four:

> The frozen cascade departs from PB exactly at loci where the fixed-ε binomial
> is **confidently wrong** — and the router, which measures *uncertainty* rather
> than *error*, cannot see them by construction.

Low VAF (≲0.25) is where "confidently wrong" lives, because that is where the
fixed-ε model and the per-read quality model disagree most. Base quality
determines the *direction*: below ~Q25 the binomial over-calls (`router_fp`),
above ~Q30 it under-calls (`router_avoids_fp` / `router_fn`). Depth amplifies
both, which is why disagreements concentrate above 48x. Segdup regions
concentrate the high-quality/low-VAF corner, which is why all 5 lost variants
are there.

This is why widening the cutoff does not help (§13): these loci have large router
margins *because the binomial is confident*. A router built on binomial
confidence is structurally blind to binomial error. Recording that is the main
scientific result of this experiment.

**No change was made to the router in response to any of this.**

## 18. Statistical analysis

* **Estimator.** Paired bootstrap on ΔF1, 10,000 resamples, module-default seed,
  identical code path to devlogs 11–13 (`evaluate_quality_error.paired_bootstrap`).
  Both arms are scored on the same resample, so the interval is for the paired
  difference and absorbs the fact that the arms agree on 99.9996% of loci.
* **Statistically significant** (CI excludes 0): `chr1_segdup` full
  (+0.005218 [+0.001418, +0.008972]) and 30x (+0.003056 [+0.000211, +0.006140]);
  the `chr1_segdup` pool; the `all` pool (+0.000581 [+0.000075, +0.001089]); the
  30x depth pool. All are *positive* — the cascade beating the arm it
  approximates, for the mechanistic reason in §17.3.
* **Statistically indistinguishable from zero** (CI contains 0): the remaining
  10 cells and every other pool, including `chr20_neutral_full`.
* **Practically negligible.** Seven cells have ΔF1 exactly 0.000000 — the two
  arms make literally identical calls on every one of ~950,000 loci. The
  remaining differences outside `chr1_segdup` are 1–4 calls per million loci.
* **Underpowered.** No cell is underpowered by the pre-registered rule
  (all ≥ 765 SNP; widest CI half-width 0.00254 vs the 0.01 threshold).
  However, power is adequate only for the *endpoint*: with 765–1,128 SNP per
  cell, one false call moves F1 by ~0.0005–0.0007, so the 0.001 equivalence bound
  is roughly 1.5 calls wide. Distinguishing "equivalent" from "worse by two
  calls" is beyond this design, and §10.2's DEGRADED verdict sits exactly there.
* **Explicitly underpowered sub-analyses.** Every stratified cell in §17.2 rests
  on 1–10 events. The *directions* are consistent across axes and across
  regions, but no individual bin is powered, and the mechanism in §17.3 is
  supported by the evidence dumps, not by the bin counts.
* **Multiplicity.** Twelve cells × one endpoint were tested with no correction.
  Under the pre-registered rule one DEGRADED cell drives the verdict; with 12
  tests at 95% coverage, seeing one cell outside its interval is unsurprising.
  This is stated because it is true, **not** to argue the verdict away — the rule
  was fixed in advance and is applied as written in §21.
* **Controls** are degenerate rather than merely worse (§16), so they support a
  qualitative claim ("the routed set is where the callers disagree") and not a
  quantitative one.

## 19. Limitations

1. **One sample, one platform.** HG002, Illumina 2×250, novoalign, GRCh38, held
   constant deliberately so that chromosome/structure was the only variable.
   Nothing here supports a claim about other samples, library preparations,
   read lengths, aligners or sequencing technologies.
2. **Four 1 Mb windows.** 4 Mb of distinct sequence, 12,261 SNP. Wide in
   structure, narrow in quantity next to chr21's 24.6 M loci.
3. **One window per regime.** `chr19_gcrich` is one GC-rich window, not the
   GC-rich regime. Region-level effects cannot be separated from
   regime-level effects with n = 1 window each.
4. **Depth realisations are downsamples of one library**, not independent
   sequencing runs. The three depths of a region are three regimes of one
   region.
5. **The 0.001 equivalence bound is ~1.5 false calls wide** at these SNP counts
   (§18). The rule was inherited unchanged from devlog 13 for consistency, but
   it is coarser than the effects being measured.
6. **High sequencing-error data is still a proxy.** `chr1_segdup` is a genuine
   *mapping*-error regime, which is a real improvement on devlog 13's
   quality-stratum proxy, but no BAM with an elevated *substitution* error rate
   was tested.
7. **SNP only.** Indels, MNPs and structural variants are outside the label set.
   In segdup regions this matters more than elsewhere.
8. **PB above 48x still sees at most `MAX_READS = 48` reads** (devlog 13 §18.4).
   Its ≥48x numbers are for a correct estimator on a subsample.
9. **Only caller compute is timed.** BAM I/O and pileup construction are excluded
   from the speedup and are paid by both arms.
10. **`chr1_segdup` truth is itself hardest there.** GIAB high-confidence calls in
    segmental duplication are the calls most likely to be wrong, so the F1 ≈ 0.88
    level in that region should not be read as an absolute accuracy.
11. **No correction for multiple comparisons** (§18).

## 20. Reproducibility

* **Pre-registration:** §1–§7 written before any caller ran on these loci; §8
  written after the interrupted-run audit but before the relaunch and before any
  metric existed. Two selection-rule revisions are disclosed in §5.2, both made
  while only annotation statistics existed. No section of §1–§7 was edited after
  results were read.
* **Data:** `fetch_v14_data.sh` (public GIAB URLs, recorded in §5.5),
  `prepare_v14_depths.sh` (`samtools view -s 42.<fraction>`, fractions derived
  from measured native depth).
* **Regions:** `select_v14_regions.py` → `results/bench_v14/region_selection.json`.
  Deterministic function of public annotation.
* **Extraction:** `run_v14_extract.sh`, 12 cells, once each, command in §8.4.
* **Scoring:** `bench_v14_crosschrom.py` → `results/bench_v14/robustness_results.json`,
  `disagreements.json`, `independence_audit.json`. Tables rendered by
  `summarize_v14.py` → `results/bench_v14/summary.md` (pure formatting; it
  computes nothing the benchmark did not).
* **Tests:** `test_bench_v14.py` (29 tests, one integration test requiring the
  caches). Full relevant suite: **241 passed**, before and after the code change.
* **Code change:** `contig` threaded through `pileup_counts.load_counts`,
  `read_level_pileup.load_reads` and `extract_bench_v12.py --contig`, all
  defaulting to `chr21` so every prior invocation is reproduced unchanged
  (§8.3). **`cascade.py`, `cheap_router.py`, `hybrid_router.py`,
  `robustness_benchmark.py`, `bench_v12_stage1.py` and `binomial_baseline.py`
  are byte-identical to their state at session start** (`git diff` empty).
* **Frozen constants**, read from `cascade.py` and asserted by the test suite:
  binomial threshold 7.0, PB threshold 10.5, router cutoff 5.411872376933351.
* **Wall-clock:** downloads ~12 min; region selection ~1 min; extraction 32 min
  (4 concurrent); benchmark 5 min 56 s.
* **Git:** no commit, no push, no branch operation, no history rewrite. HEAD
  remains `720050e` on `main`. All new work is uncommitted.
* **Timings caveat:** extraction timings inside the caches were produced with
  four concurrent jobs and are upper bounds; every speedup claim uses the
  quiet-process throughput probes of §12.
* **No historical artefact was modified.** `results/robustness_v11/*`,
  `results/bench_v12/*`, `results/bench_v13/*` and `History/3–13_DEVLOG.md` are
  untouched.

## 21. Final verdict

### 21.1 Applying the pre-registered rule as written

§3.4: *"A single powered DEGRADED cell makes the verdict NEGATIVE or MIXED."*
§3.4 mapping: *"NEGATIVE — at least one powered cell DEGRADED."*

`chr20_neutral_full` is powered (765 SNP, CI half-width 0.00254) and DEGRADED
(|ΔF1| = 0.001262 ≥ 0.001). Therefore:

> ## Pre-registered verdict: **NEGATIVE**

The STRONG NEGATIVE clause does **not** fire: no pooled result is DEGRADED, and
the cascade loses true SNP in exactly **one** region (`chr1_segdup`), not "more
than one region".

This verdict is recorded as the rule produces it. It is not softened, and the
rule is not being renegotiated after the fact. What it rests on is stated
exactly in §10.2: two false positives out of 765 SNP, in a cell whose confidence
interval contains zero.

### 21.2 On the requested four-level generalisation scale

> ## **PARTIALLY GENERALIZES**

Justification, and its relation to §21.1: the pre-registered NEGATIVE stands and
governs. On the coarser generalisation question it maps to *partial* because the
evidence splits cleanly in two:

* **What transfers.** The compute claim transfers intact (0.0030%–0.2184% of
  loci routed, 171x–396x speedup, on four new chromosomes including the hardest
  one). Accuracy equivalence holds in 11 of 12 cells, in all four region pools,
  in all three depth pools, and in the no-segdup pool. Seven cells reproduce
  PB's calls *exactly*. The frozen cutoff's knee does not move.
* **What does not.** One powered cell exceeds the pre-registered equivalence
  bound. A failure mode absent from 24.6 M chr21 loci appears here: 5 true SNP
  lost, all in segmental duplication, all at VAF ≈ 0.11–0.13 with high base
  quality. The tandem-repeat stratum is worse than PB by −0.0040 F1. The
  disagreement rate is ~6.5x chr21's.

### 21.3 The six pre-registered questions, answered

1. **Does frozen router → PB preserve PB-level accuracy on genuinely unseen
   data?** *Mostly, with one documented exception.* 11/12 cells and 9/9 pools
   are PRESERVED or IMPROVED; 1 cell is DEGRADED on a two-false-positive effect
   whose CI contains zero.
2. **How much PB compute does it eliminate?** *99.78%–99.997%.* PB runs on
   0.0030%–0.2184% of loci; projected speedup 171x–396x. Indistinguishable from
   chr21.
3. **Does the result generalise to genuinely different genomic data?** *Across
   chromosome identity and GC/isochore structure, yes* — chr20, chr19 (GC 0.591)
   and chr4 (GC 0.363) all PRESERVED at every depth. *Across mapping difficulty,
   it changes character rather than failing*: in segdup the cascade **beats** PB
   (+0.0028 pooled, CI excludes 0) while simultaneously losing 5 true variants
   PB found.
4. **What regimes cause failures?** A single axis: **loci where the fixed-ε
   binomial is confidently wrong** — VAF ≲ 0.25, with base quality determining
   direction (< Q25 over-calls, > Q30 under-calls), amplified by depth > 48x and
   concentrated in segmental duplication. 40 of 41 disagreements fall inside
   GIAB `alldifficult`. The router cannot see these loci by construction: it
   measures uncertainty, and these loci are confident.
5. **Is the frozen cutoff robust?** *Yes, on this data.* It sits on a flat
   plateau (5.0–6.0 bit-identical) ending in a structural cliff at 7.0 where
   speedup collapses 284x → 11x. The knee is in the same place on all four new
   chromosomes. It was not re-tuned, and §13 shows re-tuning would not fix the
   one DEGRADED cell anyway.
6. **Bonus finding, unprompted.** Above 48x the fixed-ε binomial *outscores* the
   repaired PB (§14) on four new chromosomes — an independent replication of
   devlog 13 §23's finding that PB is no longer the accuracy ceiling at high
   depth.

### 21.4 One paragraph, conservatively

On four chromosomes the project had never touched, spanning GC 0.363–0.591 and
one region where 100% of high-confidence bases lie in segmental duplication, the
frozen binomial → router → PB cascade reproduced PB-only accuracy in 11 of 12
region/depth cells and in every pooled result, disagreeing with PB on 41 of
10,368,078 loci (3.95×10⁻⁶), while running PB on 0.0030%–0.2184% of loci for a
projected 171x–396x speedup — so the *compute* claim generalises off chr21
essentially unchanged, and the *accuracy* claim generalises with one exception
that the pre-registered rule scores as a failure. That exception
(`chr20_neutral_full`, ΔF1 −0.001262, CI containing zero, two false positives out
of 765 SNP) makes the pre-registered verdict **NEGATIVE**, and it is reported as
such rather than reasoned away. The scientifically important result is not the
verdict but the failure boundary: the cascade departs from PB exactly where the
fixed-ε binomial is *confidently wrong* — low VAF, extreme base quality, high
depth, segmental duplication — and a router that scores binomial *uncertainty*
is structurally incapable of routing those loci, which is why no cutoff between
4.5 and 8.0 repairs the one failing cell. In segmental duplication this cuts both
ways in the same region: the cascade avoided 24 PB false positives and lost 5
true variants, netting a statistically significant *improvement* over the caller
it is supposed to approximate. Nothing in this document licenses a claim about
another sample, another platform, indels, or any chromosome beyond the four
1 Mb windows tested, and no threshold, cutoff or router file was changed at any
point.

## 22. Recommended next experiment

**Route on disagreement risk, not on binomial uncertainty — and test it against
this frozen cascade as the baseline.**

§17.3 shows the current router's blind spot is structural: it routes loci where
the binomial is *unsure*, but the errors live where the binomial is *confidently
wrong*. Those loci have a common, cheap signature already present in the count
matrix — VAF in roughly 0.05–0.25 combined with mean base quality far from Q25 —
and none of it requires the read tensor or a neural model.

The concrete experiment: add a second, disjoint cheap predicate
(`0.05 ≤ VAF ≤ 0.25` and (`meanBQ < 25` or `meanBQ > 32`)) as an **additional**
route-to-PB trigger, leaving `|binomial_LLR − 7.0| ≤ 5.412` exactly as it is, and
measure on this same 12-cell benchmark. It is a single-variable change against a
frozen, now cross-chromosome-validated baseline and an existing harness, so it is
a clean A/B. The pre-registered question would be whether it recovers the 8
`router_fp` and 5 `router_fn` events, and what it costs in routed fraction —
§12 suggests the budget is available, since the VAF 0.15–0.25 band holds only
5,677 of 10.37 M loci.

Two secondary follow-ups, in priority order:

1. **Cross-sample transfer** (HG003/HG004, same platform and the same four
   windows) — the last untested axis of the "does it generalise" question, and
   cheap now that the extraction path is contig-general.
2. **PB's low-VAF false-positive mode**, still the largest accuracy defect in the
   system (devlog 13 §25, independently reproduced here at ≥48x and in segdup).
   The 26 `router_avoids_fp` events are a ready-made test set for a mapping-error
   term in the PB error floor.

What should **not** happen next: widening or re-tuning the router cutoff (§13
shows it does not address the failure mode and costs an order of magnitude of
speedup), or revisiting Mamba (settled negative in devlogs 10 and 12).
