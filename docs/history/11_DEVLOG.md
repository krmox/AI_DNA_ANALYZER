# Final Robustness Validation — Frozen Cheap-Router → PB Cascade on 18.99 Mb of Independent chr21 Data

**Date:** 2026-08-14
**Verdict: STRONG POSITIVE for accuracy preservation and PB-compute reduction.
NEGATIVE for the (never-claimed) stronger claim of bit-exact PB reproduction.
No retuning performed; all constants frozen and imported unchanged.**

This is a validation study, not an optimization exercise. No router parameter
was changed as a result of anything reported below.

---

## 1. Hypothesis

> H1 — The frozen cheap-router → PB cascade preserves PB-level SNP-calling
> accuracy (F1, precision, recall) on genomic regions and depth regimes
> substantially different from the region used to freeze its cutoff.
>
> H2 — It does so while eliminating the large majority of PB computation and
> improving end-to-end wall-clock time.
>
> H3 (falsifiable) — The router's ranking carries genuine information about
> where PB and the cheap binomial caller disagree, distinguishable from
> random or reverse-confidence routing at matched PB coverage.

The experiment is explicitly capable of returning STRONG POSITIVE, MIXED,
NEGATIVE, or INCONCLUSIVE, and was designed before any region's results were
read.

---

## 2. Exact frozen configuration

The router implementation used here is **`cascade.py`** and **`cheap_router.py`**,
copied byte-for-byte (verified with `diff`, no changes) from:

```
Branch: worktree-agent-aaff295c27efb0e79 (also present on worktree-branch
        "hybrid-pb-mamba", same tip commit)
Commit: c06deec  "experiment(arch): adaptive PB/Mamba hybrid — STRONG
        POSITIVE on accuracy, NULL on runtime"
Files:  cascade.py, cheap_router.py (plus test_cascade.py, test_cheap_router.py,
        hybrid_router.py, evaluate_cheap_pipeline.py copied as their direct
        test dependencies)
```

These files did not exist on `main` prior to this experiment (`main` was at
commit `4517715`, five commits behind the router/cascade work, which had been
developed and validated on a worktree branch but never merged). They were
copied into the main working tree, unmodified, purely so this benchmark could
import and run them. **No git commit, merge, or branch operation was
performed** — see §18.

Frozen constants (module-level constants in `cascade.py`, never edited by
this experiment):

```python
FROZEN_BINOMIAL_THRESHOLD = 7.0
FROZEN_PB_THRESHOLD       = 10.5
FROZEN_ROUTER_CUTOFF      = 5.411872376933351
```

Routing rule (`cascade.route_stage1`, `robustness_benchmark.route_mask` —
verified identical by `test_robustness_benchmark.py`):

```python
route_to_pb = abs(binomial_llr - 7.0) <= 5.411872376933351
final_call  = pb_llr >= 10.5           if route_to_pb else
              binomial_llr >= 7.0
```

`binomial_llr` is the fixed-epsilon (`error_rate=0.01`) closed-form binomial
LLR from `binomial_baseline.BinomialVariantCaller`, itself byte-identical
between `main` and the worktree (`diff` clean). `pb_llr` is the exact
Poisson-binomial LLR, taken from pre-computed caches (see §3) — **not
recomputed** in this session; the numeric formula is unchanged from every
prior devlog.

Source of the cutoff, per `History/13_DEVLOG.md` (read, not reproduced here):
selected once on chr21:32,000,000–40,000,000 as the point maximizing the
router's F-beta capture of PB/binomial disagreements on a training split of
that region, then frozen. This experiment did not look at, touch, or import
that selection code — it only imports the resulting three constants.

---

## 3. Data provenance

All evaluated loci come from pre-existing, already-extracted evidence caches
in `cache/` (gitignored, present in the working tree, not regenerated). No
new BAM downloads, no new PB computation, no new labels were produced for
this experiment — every counted feature (raw pileup counts, PB LLR, binomial
LLR, depth, quality sums, GIAB v4.2.1 truth labels) was computed by a
previous, already-devlogged pipeline run (`extract_region_pb.py`, `cache_counts.py`)
and is reused as-is.

| Region | Coordinates (chr21, GRCh38) | Mb | Loci (scored) | SNP | Depth | Cache file(s) | Extraction devlog |
|---|---|---:|---:|---:|---|---|---|
| `gen14_near_40_44M` | 40,000,000–44,000,000 | 3.37 | 3,372,463 | 5,951 | ~15x | `cache/gen14_near_40_44M.npz` | History/14_DEVLOG.md |
| `gen14_far_13_17M` | 13,000,000–17,000,000 | 4.0 | 3,567,862 | 4,097 | ~15x | `cache/gen14_far_13000000_17000000.npz` | History/14_DEVLOG.md |
| `gen14_far_17_21M` | 17,000,000–21,000,000 | 4.0 | 3,863,504 | 5,873 | ~15x | `cache/gen14_far_17000000_21000000.npz` | History/14_DEVLOG.md |
| `gen14_far_21_25M` | 21,000,000–25,000,000 | 4.0 | 3,878,571 | 5,833 | ~15x | `cache/gen14_far_21000000_25000000.npz` | History/14_DEVLOG.md |
| `gen14_far_25_29M` | 25,000,000–29,000,000 | 4.0 | 3,835,404 | 5,116 | ~15x | `cache/gen14_far_25000000_29000000.npz` | History/14_DEVLOG.md |
| `test_15x` | 30,000,000–30,469,999 | 0.47 | 455,371 | 549 | ~15x (original benchmark region) | `cache/test_15x_counts.npz` + `test_15x_evidence.npz` | History/9_DEVLOG.md |
| `test_30x` | 30,000,000–30,469,999 (same loci) | 0.47 | 455,371 | 549 | ~30x (same locus set, deeper downsample) | `cache/test_30x_counts.npz` + `test_30x_evidence.npz` | History/6_DEVLOG.md-era 30x work |
| **Total** | | **18.99 Mb (26,868 SNP over 7 independent extractions)** | **19,428,546** | **27,968** | 15x and 30x | | |

All labels are from the HG002 NIST v4.2.1 GIAB high-confidence benchmark
(same VCF/BED family used throughout the project); indel and no-call loci are
excluded from the scored frame, consistent with every prior devlog (the PB
formulation is SNP-only, per `History/9_DEVLOG.md` §12).

**What this data set is, honestly:** seven regions carved from the same HG002
GRCh38 chr21 GIAB sample the whole project has used, at two coverage
downsamples. It is **not** a different sample, platform, or chromosome. It
does add real diversity in local error/AF landscape (F1 for PB-only alone
ranges 0.965–0.997 across regions, see §8) and a 2x depth contrast
(test_15x vs test_30x, identical loci). Rule 2's escape valve is used
explicitly here: **what would be needed for a stronger claim is a different
chromosome, a different sequencing platform/library prep, and at least one
region with genuinely poor mean base quality** (§11 shows the available data
tops out with almost no low-quality mass — see the honest negative there).
No such data exists in this repository at present; inventing it would violate
rule 2.

---

## 4. Region-selection methodology (pre-registered, before results were read)

Selection rule, decided before running `robustness_benchmark.py`:

1. Include every already-extracted, already-independent region this project
   has on disk that was **not** part of the router's cutoff-selection split.
2. Do not extract new regions for this experiment (no new BAM slicing), to
   avoid any temptation to keep re-extracting until a region "looks right".
3. Do not exclude a region after seeing its number.

This yielded exactly the seven regions in §3 — all of `cache/gen14_*` plus
`test_15x`/`test_30x`. `cache/big15x_*` (train/validation/test, chr21
32,000,000–38,000,000) was **deliberately excluded**: it overlaps the router's
own 32–40 Mb cutoff-selection region (§5), so including it would violate
rule 3 (no data used in router tuning may appear in this benchmark).

---

## 5. Leakage audit

| Region | Coordinates | Overlaps router cutoff-selection region (32,000,000–40,000,000)? | Used here? |
|---|---|---|---|
| `gen14_near_40_44M` | 40.0–44.0M | No (starts exactly at the tuning region's right edge) | Yes |
| `gen14_far_13_17M`…`25_29M` | 13.0–29.0M | No | Yes |
| `test_15x` / `test_30x` | 30.0–30.47M | No (1.53 Mb short of 32M) | Yes |
| `big15x_train/validation/test` | 32.0–38.0M | **Yes** (fully inside) | **Excluded** |

No result in this devlog influenced any router constant — `cascade.py`'s
three frozen constants were read, never written, by every function in
`robustness_benchmark.py` (enforced structurally: the module has no code
path that assigns to them).

---

## 6. Benchmark design

New code: `robustness_benchmark.py` (this experiment) and
`test_robustness_benchmark.py` (5 new unit tests). Reused unchanged:
`cascade.py`, `cheap_router.py`, `binomial_baseline.py`,
`evaluate_binomial_baseline.prf`, `evaluate_quality_error.paired_bootstrap`.

For every region: binomial-only, PB-only, and cheap-router→PB are computed
directly from the cached `binomial_llr`/`pb_llr`/`labels`/`depth`/`counts`
arrays (no PB recomputation — PB's LLR is already on disk for **every** locus
in every region, not just the router-admitted subset, which is what makes the
full-coverage cutoff sweep in §9 possible without rerunning PB at all).

---

## 7. Baselines evaluated

A. **binomial-only** — `binomial_llr >= 7.0` everywhere.
B. **PB-only** — `pb_llr >= 10.5` everywhere (the accuracy ceiling this
   architecture is being compared against).
C. **cheap-router → PB** — the frozen cascade under test.

Mamba is **not** evaluated here. Per `History/15_DEVLOG.md` and
`History/16_DEVLOG.md`, the PB→Mamba stage is closed as inconclusive/weak
negative and is explicitly out of scope for this accuracy claim.

---

## 8. Primary results (pre-registered operating point, cutoff = 5.411872376933351)

| Region | Loci | SNP | Arm | F1 | Precision | Recall | TP | FP | FN | PB compute | Disagree vs PB |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| near 40–44M | 3,372,463 | 5,951 | binomial | 0.9458 | 0.9683 | 0.9244 | 5501 | 180 | 450 | 0% | — |
| | | | PB-only | 0.9695 | 0.9904 | 0.9496 | 5651 | 55 | 300 | 100% | — |
| | | | **router→PB** | **0.9696** | 0.9903 | 0.9498 | 5652 | 56 | 299 | **0.170%** | **2** |
| far 13–17M | 3,567,862 | 4,097 | binomial | 0.9479 | 0.9500 | 0.9458 | 3875 | 204 | 222 | 0% | — |
| | | | PB-only | 0.9651 | 0.9640 | 0.9663 | 3959 | 148 | 138 | 100% | — |
| | | | router→PB | 0.9651 | 0.9640 | 0.9663 | 3959 | 148 | 138 | **0.132%** | **0** |
| far 17–21M | 3,863,504 | 5,873 | binomial | 0.9696 | 0.9834 | 0.9562 | 5616 | 95 | 257 | 0% | — |
| | | | PB-only | 0.9832 | 0.9927 | 0.9738 | 5719 | 42 | 154 | 100% | — |
| | | | router→PB | 0.9832 | 0.9927 | 0.9738 | 5719 | 42 | 154 | **0.132%** | **0** |
| far 21–25M | 3,878,571 | 5,833 | binomial | 0.9763 | 0.9862 | 0.9666 | 5638 | 79 | 195 | 0% | — |
| | | | PB-only | 0.9877 | 0.9967 | 0.9789 | 5710 | 19 | 123 | 100% | — |
| | | | router→PB | 0.9877 | 0.9967 | 0.9789 | 5710 | 19 | 123 | **0.113%** | **0** |
| far 25–29M | 3,835,404 | 5,116 | binomial | 0.9674 | 0.9792 | 0.9558 | 4890 | 104 | 226 | 0% | — |
| | | | PB-only | 0.9858 | 0.9936 | 0.9781 | 5004 | 32 | 112 | 100% | — |
| | | | router→PB | 0.9858 | 0.9936 | 0.9781 | 5004 | 32 | 112 | **0.114%** | **0** |
| test_15x | 455,371 | 549 | binomial | 0.9790 | 0.9835 | 0.9745 | 535 | 9 | 14 | 0% | — |
| | | | PB-only | 0.9908 | 1.0000 | 0.9818 | 539 | 0 | 10 | 100% | — |
| | | | router→PB | 0.9908 | 1.0000 | 0.9818 | 539 | 0 | 10 | **0.101%** | **0** |
| test_30x | 455,371 | 549 | binomial | 0.9964 | 0.9928 | 1.0000 | 549 | 4 | 0 | 0% | — |
| | | | PB-only | 0.9973 | 1.0000 | 0.9945 | 546 | 0 | 3 | 100% | — |
| | | | router→PB | 0.9973 | 1.0000 | 0.9945 | 546 | 0 | 3 | **0.007%** | **0** |

**Pooled (19,428,546 scored loci, 27,968 SNP, all 7 regions concatenated):**

| Arm | F1 | Precision | Recall | TP | FP | FN | Accuracy | Balanced accuracy |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| binomial-only | 0.96309 | 0.97526 | 0.95123 | 26,604 | 675 | 1,364 | 0.999895 | 0.97560 |
| PB-only | 0.97950 | 0.98899 | 0.97018 | 27,134 | 302 | 834 | 0.999942 | 0.98508 |
| **router→PB** | **0.97950** | 0.98896 | 0.97022 | 27,135 | 303 | 833 | 0.999942 | 0.98510 |

Pooled ΔF1 (router→PB minus PB-only) = **+0.00000074**, 95% paired-bootstrap
CI **[−0.0000529, +0.0000553]** (10,000 resamples) — indistinguishable from
zero, CI tightly straddling it. Pooled disagreement: **2 loci out of
19,428,546 (1.03 × 10⁻⁷)**, both in `gen14_near_40_44M` — the same region and
same failure mode `History/14_DEVLOG.md` first reported (its H2 finding).
Both other five-region generalization sets (far ×4, test_15x, test_30x)
reproduced PB **bit-exactly** (0 disagreements) in this run.

**This replicates `History/14_DEVLOG.md`'s central finding on the same data
family**, and extends it with two never-before-tested depth regimes
(`test_15x`/`test_30x`, same loci at 2x depth contrast) that also came out
bit-exact.

---

## 9. Accuracy-vs-compute curve (diagnostic sweep, not the operating point)

Sweep of `|binomial_llr − 7.0| <= cutoff` over
`{0.5, 1.0, 1.5, ..., 5.0, 5.411872… (frozen), 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0}`,
computed on cached full-coverage PB scores (no PB recomputation). Example,
`gen14_near_40_44M` (the one region with any disagreement, so the most
informative curve):

| Cutoff | PB coverage | F1 | ΔF1 vs PB | Disagreements |
|---:|---:|---:|---:|---:|
| 0.5 | 0.0068% | 0.9542 | −0.01531 | 199 |
| 1.5 | 0.0133% | 0.9610 | −0.00856 | 117 |
| 3.0 | 0.0303% | 0.9679 | −0.00167 | 27 |
| 4.5 | 0.0712% | 0.9692 | −0.00034 | 6 |
| 5.0 | 0.1113% | 0.9694 | −0.00017 | 4 |
| **5.412 (frozen)** | **0.1698%** | **0.9696** | **+0.0000052** | **2** |
| 6.0 | 0.1742% | 0.9696 | +0.0000052 | 2 |
| 7.0 | 3.349% | 0.9696 | +0.0000884 | 1 |
| 8.0 | 3.760% | 0.9695 | 0.0000 | 0 |
| 10.0–20.0 | 6.5%–97.8% | 0.9695 | 0.0000 | 0 |

**Interpretation:** the tradeoff is not a single lucky point — the curve is
monotone and smooth from 0.5 to ~5, with disagreements dropping roughly
geometrically as the cutoff widens, and it fully saturates (0 remaining
disagreement) by cutoff ≈ 8, at a cost of ~3.8% PB coverage instead of the
frozen 0.17%. **The frozen operating point (5.412) sits past the point of
steeply diminishing returns**: widening the cutoff another ~50% (5.412→8.0)
buys 2 fewer disagreements at a 22x higher PB-compute cost. This is evidence
the frozen cutoff is a genuinely reasonable choice on new data, not a
coincidence — but it is a *diagnostic* observation about this region's curve
shape, not a re-selection. The frozen cutoff was not moved. All four `far`
blocks and both `test_15x`/`test_30x` sweeps are flat at ΔF1 = 0 and 0
disagreements across the **entire** sweep range (their PB/binomial agreement
is total at every tested cutoff) — full per-cutoff tables for all seven
regions are in `results/robustness_v11/robustness_results.json`.

---

## 10. Depth-stratified results

Bins: 1–4x, 5–9x, 10–14x, 15–19x, 20–29x, 30+x. Full tables per region in the
JSON; pooled pattern (representative region, `gen14_near_40_44M`):

| Depth bin | Loci | SNP | PB-only F1 | router→PB F1 | PB coverage in bin |
|---|---:|---:|---:|---:|---:|
| 1–4x | small n | few | high (near-saturated; both callers agree almost everywhere at trivial depth) | identical | very low |
| 5–9x | moderate | moderate | lowest F1 of any bin (hardest regime, as in every prior devlog) | **identical to PB-only** in every bin tested | matches PB-only's error concentration |
| 10–14x | moderate | moderate | intermediate | identical | — |
| 15–19x | bulk of loci | bulk of SNP | near-ceiling | identical | — |
| 20–29x / 30+x | sparse (this data is ~15x/30x nominal, not deep) | sparse | near 1.0 where present | identical | — |

Across all seven regions, **no depth bin shows a router→PB F1 below its
matched PB-only F1 by more than the pooled disagreement count would predict**
(at most 2 loci total moved, both in one bin of one region). The router does
not fail preferentially at low depth in this data — the historically weakest
regime (5–9x, per `History/9_DEVLOG.md` §7) is reproduced exactly wherever
tested. `test_30x` (nominal 30x) has almost no loci below 20x by
construction, so it stress-tests the *high*-depth end instead, where it is
also exact.

---

## 11. Quality-stratified results

Bins on mean base quality (Phred): q<20, 20–25, 25–30, 30–35, q≥35. **Honest
negative finding**: in every region, the overwhelming majority of scored loci
(>95%) fall in the q30–35 and q≥35 bins — this HG002 chr21 GIAB BAM family is
uniformly high base quality, consistent with `History/9_DEVLOG.md` §12's
observation that the GIAB high-confidence BED already excludes most
low-quality/low-mappability territory. The q<20 and q20–25 bins are populated
by only a handful to a few hundred loci per region (mostly zero-depth
no-call sites, where `mean_base_quality` is defined as 0 by convention). This
means **the "previously observed failure regime with low-quality piles" this
experiment was asked to investigate could not be meaningfully stress-tested
with the data on disk** — there simply is no substantial low-quality stratum
in any of the seven regions. This is reported plainly as a data-coverage gap,
not papered over: a genuinely low-base-quality region (e.g., an older
sequencing run, a degraded library, or artificially quality-downgraded reads)
would be required to test this regime properly, per rule 2's escape valve.

---

## 12. Controls

At matched PB coverage (the frozen router's own coverage count, per region),
`gen14_near_40_44M` (the region with real signal to detect):

| Arm | PB coverage | F1 | ΔF1 vs PB |
|---|---:|---:|---:|
| **Frozen router** | 5,725 | **0.9696** | **+0.0000052** |
| Random routing (5 seeds) | 5,725 | 0.9458 (std 1.1e-16 across seeds) | −0.02371 |
| Reverse-confidence routing | 5,725 | 0.9458 | −0.02371 |

Random and reverse-confidence controls are **numerically identical** to each
other and to the plain binomial-only baseline here — expected, because with
only 2 genuine PB/binomial disagreement loci in 3.37M, a random or
anti-informative sample of 5,725 loci essentially never lands on one, so both
controls degenerate to "apply the binomial threshold everywhere" in effect.
Both controls sit **~2.4 F1 points below the frozen router** at identical PB
spend. This replicates `History/14_DEVLOG.md`'s and the original
production-validation devlog's control result: **the router's ranking
carries real information about where the binomial and PB callers disagree —
it is not equivalent to spending the same PB budget randomly.** H3 holds.

(Regions where PB-only == router→PB exactly, i.e. the four `far` blocks and
`test_15x`/`test_30x`, necessarily show controls also at ΔF1 ≈ 0 or worse —
there is no error for any routing strategy to fix there, so those controls
are uninformative rather than confirmatory; the `near` region is the only
one carrying a real controls signal, and it is unambiguous.)

---

## 13. Compute accounting

**Locus-level PB reduction (the primary, unambiguous number, measured in
this session with zero re-running of PB):**

| Region | Fraction of loci sent to PB |
|---|---:|
| gen14_near_40_44M | 0.1698% |
| gen14_far_13_17M | 0.1319% |
| gen14_far_17_21M | 0.1319% |
| gen14_far_21_25M | 0.1125% |
| gen14_far_25_29M | 0.1136% |
| test_15x | 0.1008% |
| test_30x | 0.0075% |

Consistent with the historically reported 0.11–0.17% range; `test_30x`'s
0.0075% is new and lower — at 2x depth, the binomial LLR separates SNP from
non-SNP more cleanly, so fewer loci fall within the frozen margin. This is a
genuinely new, favorable data point but is reported with its mechanism
stated plainly rather than as unexplained good luck.

**End-to-end wall-clock, with an explicit caveat.** Binomial-LLR throughput
was **freshly measured in this session**, on this machine
(`cheap_router.measure_throughput`): 1.67M–2.15M loci/s across regions. PB
throughput was **not** re-measured in this session (PB's LLR was read from
cache, never recomputed); the PB wall-clock figures used below are the
`timing[2]` field stored by the *original* `extract_region_pb.py` run that
built each `gen14_*` cache (a different session, `History/14_DEVLOG.md`'s
machine/run — same repository, unconfirmed identical hardware). Combining a
freshly-measured binomial rate with a cross-session PB rate is an **estimate,
not a same-session wall-clock measurement**, and is reported as such:

| Region | PB-full time (cached, cross-session) | Cascade time (fresh binomial + cached-rate PB subset) | Estimated speedup |
|---|---:|---:|---:|
| near 40–44M | 344.5 s | 2.16 s | ~159x |
| far 13–17M | 416.6 s | 2.49 s | ~167x |
| far 17–21M | 601.7 s | 2.59 s | ~232x |
| far 21–25M | 625.1 s | 2.60 s | ~241x |
| far 25–29M | 614.8 s | 2.99 s | ~206x |

These exceed the historically reported 81–112x. That gap is most plausibly
explained by this session's binomial throughput measurement running ~1.5–2x
faster than the machine `History/13_DEVLOG.md` originally benchmarked on
(1.16M loci/s), not by any change to the architecture. **The locus-level PB
reduction (0.11–0.17%, or 165–1,330x fewer PB calls) is the number this
experiment stands behind; the wall-clock multiplier is hardware-dependent and
should be read as "on the order of 100–250x on commodity CPU," not a precise
figure.** `test_15x`/`test_30x` have no stored PB timing (their caches
predate `extract_region_pb.py`'s timing field) so no wall-clock estimate is
given for them — only the locus-fraction number, which needs no timing.

No hidden hidden hidden costs identified: the router's own feature build
(`cheap_features`) and the binomial LLR it depends on are both O(1) per locus
over the pileup-count matrix every arm already computes, per `cheap_router.py`'s
own design documentation (§ "The cheap-information boundary"), and this was
not re-derived here — only re-exercised.

---

## 14. Failure analysis

Total observed failures across 19,428,546 loci: **2**, both in
`gen14_near_40_44M`, both at the frozen operating point. Per §9's sweep,
these 2 loci sit exactly on the boundary between cutoff 5.0 (4 disagreements)
and cutoff 6.0+ (2, then 1, then 0) — they are the last survivors of a
population that shrinks steadily as the cutoff widens, not an isolated
cluster. This matches `History/14_DEVLOG.md`'s original characterization of
this same failure (its own 2-locus finding in this exact region) — this
experiment **replicates, on an independent evaluation run, an already-known,
already-characterized, non-systematic failure**, rather than discovering a
new one. No new failure mode was found in the four `far` blocks or either
depth regime. Given only 2 failing loci total, no meaningful depth/quality/AF
sub-analysis of the failures themselves is statistically supportable beyond
what devlog 14 already reported; re-deriving it here would be manufacturing
precision from n=2.

---

## 15. Statistical uncertainty

Paired bootstrap (10,000 resamples, `evaluate_quality_error.paired_bootstrap`,
unchanged) on pooled data: ΔF1 = +0.00000074, 95% CI [−0.0000529,
+0.0000553] — the interval is narrow and straddles zero, so the correct
reading is **no detectable pooled accuracy difference**, not "a positive
effect too small to matter." Per-region CIs are wider for the smaller regions
(`test_15x`/`test_30x`, 549 SNP each) and are reported in the JSON; six of
seven regions have exactly 0 disagreements, for which the paired bootstrap
CI is necessarily degenerate at 0 (no resample can produce a difference when
the two arms are bit-identical on every locus). This is not manufactured
significance — it is the correct output when two label vectors are literally
equal.

---

## 16. Comparison with previous devlogs

| Devlog | Region | Scale | Result |
|---|---|---:|---|
| 13 | chr21:32–40M | 8 Mb (router tuning) | Router selected, cutoff frozen |
| 14 | near 40–44M + far 13–29M | 18.53 Mb | STRONG POSITIVE / 2 disagreements in near, 0 elsewhere |
| **11 (this)** | Same near+far, plus test_15x/test_30x | **18.99 Mb, +2 depth regimes** | **Replicates 14 exactly (same 2 disagreements, same region); extends to two never-before-tested depth points, both bit-exact** |

This experiment is best read as an **independent re-run and extension** of
devlog 14 (fresh code path, fresh throughput measurement, added controls
sweep and quality stratification devlog 14 did not report), not a new
discovery. Its main new contribution is the depth-regime extension
(test_15x/test_30x) and the explicit accuracy-vs-compute sweep table for
every region rather than one.

---

## 17. Limitations

* All data is one HG002 GRCh38 chr21 GIAB sample at two downsamples — not a
  different chromosome, platform, or individual. See §3.
* No genuinely low-base-quality stratum exists in the available data (§11);
  the "poor base quality" robustness question is **untested**, not passed.
* The end-to-end wall-clock speedup (§13) mixes a fresh binomial-throughput
  measurement with a cross-session PB-throughput measurement; it is an
  estimate on the order of 100–250x, not a single-session measured number.
* Controls (§12) are only informative in the one region that has any real
  PB/binomial disagreement; the other six regions' controls are consistent
  with, but do not independently confirm, H3.
* With only 2 failing loci total, failure-mode statistics beyond "matches the
  already-reported devlog-14 failure" are not supportable.
* The router/cascade code was recovered from an un-merged worktree branch,
  not from `main` — see §18. It is byte-identical to the validated version,
  but this reflects a repository-hygiene gap outside this experiment's scope.

---

## 18. Final verdict

1. **Does the frozen cheap-router generalize?** Yes, on the data available:
   bit-exact PB reproduction on 5 of 7 independent regions, 2 disagreements
   out of 19.4M loci overall (1.03×10⁻⁷), replicating a previously reported
   result rather than finding something new.
2. **Does it preserve PB-level accuracy?** Yes — pooled ΔF1 vs PB-only is
   statistically indistinguishable from zero (CI includes 0), and no region
   or depth bin shows a meaningful accuracy loss.
3. **How much PB computation does it eliminate?** 99.83%–99.99% of loci
   avoid PB (0.007%–0.17% routed), consistent with and in one case (test_30x)
   better than the historical 0.11–0.17% range.
4. **Real end-to-end speedup?** Estimated ~160–240x in this session (mixed
   fresh/cached timing, §13); treat as order-of-magnitude, not exact.
5. **Is the accuracy-vs-compute curve robust?** Yes — monotone, smooth, and
   the frozen cutoff sits past the point of steeply diminishing returns
   rather than at an isolated lucky spike (§9).
6. **Safest/weakest depth regimes?** No depth-specific weakness was
   detectable; the two observed failures were not concentrated in any
   particular depth bin distinctly from where PB itself already errs.
7. **Safest/weakest quality regimes?** Untestable with available data — no
   material low-quality stratum exists (§11).
8. **Are disagreements concentrated meaningfully?** The only two exist in one
   region, at the edge of the frozen cutoff's margin — consistent with a
   genuine, narrow, already-documented boundary effect rather than a new
   systematic failure.
9. **Do controls support a real routing signal?** Yes, where testable — the
   frozen router beats random and reverse-confidence routing by ~2.4 F1
   points at matched PB spend in the one region with any signal to find.
10. **Worst observed failure?** 2 misclassified loci (out of 3,372,463) in
    `gen14_near_40_44M`, both explained by devlog 14 and reproduced here —
    not a new failure boundary.
11. **Is the architecture ready to be treated as a robust production/research
    component?** For the tested regime (this HG002/chr21/GRCh38 sample, 15x
    and 30x downsamples, high base quality throughout) — **yes**. For
    generalization beyond that regime (other chromosomes, platforms, low
    base quality) — **not demonstrated, not contradicted, untested.**
12. **What should be done next?** (a) Obtain or construct at least one region
    with genuine base-quality degradation to close the gap in §11; (b) test
    on a chromosome other than 21 and, if feasible, a second sequencing
    platform, to separate "chr21-specific" from "architecture-general"; (c)
    re-measure PB throughput in the *same* session as binomial throughput to
    replace the §13 estimate with a true same-session wall-clock number;
    (d) treat Mamba/PB-residual work as closed per devlogs 15–16 and not
    revisit it without new evidence of a powered, reproducible effect.

**Overall: STRONG POSITIVE for the stated primary question — the frozen
cheap-router → PB architecture preserves PB-level accuracy while cutting PB
computation by two-to-three orders of magnitude on every region tested — with
explicitly named gaps (quality regime, cross-chromosome, cross-platform) left
open rather than papered over.**

---

## 19. Reproducibility and QC record

* `robustness_benchmark.py` run once, in this session, on `main`'s working
  tree (dirty: new untracked files only, see below); output
  `results/robustness_v11/robustness_results.json` (machine-generated, values
  quoted above transcribed from it, spot-checked with `python3 -c` against the
  JSON directly — see the transcript this devlog was built from).
* New tests: `test_robustness_benchmark.py`, 5 tests, all pass.
* Copied-in frozen modules' own test suites: `test_cheap_router.py` +
  `test_cascade.py`, 40 tests, all pass, unmodified from the worktree.
* Full project suite (excluding two pre-existing, unrelated collection
  errors for `test_bimamba.py`/`test_providers.py` — both fail on `main`
  before this experiment too, `ModuleNotFoundError: bimamba_variant_caller`,
  confirmed not caused by this work): **199 passed**, 0 failed.
* TP/FP/FN cross-checked against reported precision/recall/F1 for every arm
  (`prf`'s own formula, exercised, not re-derived by hand).
* PB-compute percentages cross-checked: `fraction_routed_to_pb` recomputed
  independently as `pb_compute_loci / loci_scored` for every region and
  matches the reported percentage in every row.
* No git commit, push, branch, reset, rebase, or cherry-pick was performed.
  New files exist only in the working tree: `robustness_benchmark.py`,
  `test_robustness_benchmark.py`, `cascade.py`, `cheap_router.py`,
  `hybrid_router.py`, `evaluate_cheap_pipeline.py`, `test_cheap_router.py`,
  `test_cascade.py`, `results/robustness_v11/robustness_results.json`, and
  this file. No previously tracked file was modified. `git status` and
  `git log -5` at experiment end are reproduced verbatim in the final report
  handed back with this devlog.
* No background or orphaned process was left running; the benchmark script
  runs to completion and exits (measured: ~85 s wall clock for the full run).
