# bench_v12 — Three-Stage Cascade Benchmark on a Newly Extracted chr21 Region

**Date:** 2026-08-14
**Status of this file when sections 1–6 were written: PROTOCOL ONLY.**
Sections 1–6 below were written **before any Stage-1/2/3 metric was computed or
read**. The only thing that had been run at that point was the *extraction* of
the new region's evidence caches (`extract_bench_v12.py`), which produces
features/LLRs/labels but no accuracy metric, and a 62,592-locus smoke test of
that extractor. Results sections (7+) were appended afterwards and no protocol
sentence above them was edited after results were seen.

---

## 1. Research question

Three questions, per the commissioning brief:

1. **Stage 1** — Does the frozen cheap-router → PB component remain accurate and
   computationally advantageous on genomic data that is *new* relative to
   `History/11_DEVLOG.md` (which scored 19,428,546 loci over 7 chr21 blocks)?
2. **Stage 2** — Does PB → Mamba provide measurable benefit when evaluated with
   sufficient statistical power and an appropriately-defined candidate pool?
3. **Stage 3** — Does the complete three-stage cascade provide a meaningful
   accuracy/compute tradeoff versus PB-only, Mamba-only, and simpler baselines?

This is a validation benchmark. No threshold, cutoff, or architecture parameter
is refit anywhere in it.

## 2. Hypotheses (falsifiable, stated before results)

* **H1** — On chr21:31.0–32.0 Mb (never previously scored) the frozen router →
  PB arm reproduces PB-only F1 to within a paired-bootstrap CI that includes
  zero, at ≤1% PB locus coverage.
* **H2** — H1 holds at all three depth regimes extracted (≈15x, ≈30x, full
  ≈70–80x) and in the newly-defined *hard* strata (low depth, low VAF,
  PB-uncertain, lower base quality) cut from the 7 already-cached regions.
* **H3** — PB → Mamba, evaluated on a PB-uncertain candidate pool large enough
  to have power and on a region genuinely held out from the Mamba
  checkpoints' training split, produces ΔF1 ≥ 0 versus PB-only.
  *Prior expectation from `History/16_DEVLOG.md` (worktree line): this will
  fail or be null. The experiment is designed to be able to say so.*
* **H4** — The complete cascade lies on the accuracy/compute Pareto frontier,
  i.e. no arm gives strictly more accuracy at strictly less compute.

## 3. Dataset provenance

| Item | Value |
|---|---|
| Sample | HG002 (NA24385), GIAB |
| Reference | GRCh38, `data/reference/chr21_full.fa` |
| Truth | NIST GIAB HG002 v4.2.1 benchmark VCF + high-confidence BED, region-sliced: `data/giab_hg002_real_train/hg002_chr21_31_33M.vcf.gz`, `..._highconf.bed` |
| New region | chr21:31,000,000–32,000,000 (1.0 Mb, window-aligned to SEQ_LEN=64) |
| BAM @15x | `data/giab_hg002_real_15x/hg002_chr21_31_33M_15x.bam` |
| BAM @30x | `data/giab_hg002_real_30x/hg002_chr21_31_33M_30x.bam` |
| BAM @full | `data/giab_hg002_real_train/hg002_chr21_31_33M.bam` |
| Extraction | `extract_bench_v12.py` (new; delegates all statistics to existing modules) |
| Output caches | `cache/bench_v12/new31_32M_{15x,30x,full}.npz` |
| Reused caches | the 7 regions of `results/robustness_v11` (`cache/gen14_*.npz`, `cache/test_{15,30}x_*.npz`), **read-only, not regenerated** |
| Mamba checkpoints | `.claude/worktrees/agent-aaff295c27efb0e79/checkpoints_big15x/best_big15x_seed{A,B,C}.pt` and `best_big15x_shuffled.pt` (read-only; trained on chr21:32–38 Mb), plus `checkpoints_pb_residual/best_v1_repr_seed{A,B}.pt`, `best_v2_gate_seedA.pt` (trained on chr21:31–33 Mb → **leaked w.r.t. this region**, reported separately and excluded from every verdict) |

**Honest scope statement, carried forward from devlog 11 and not weakened:**
everything here is still one HG002 GRCh38 chr21 sample. No second chromosome,
individual, or sequencing platform exists in this repository, and none was
downloaded (explicit instruction). "Independent" below always means
*independent genomic coordinates and/or an independent depth realisation*,
never an independent sample.

## 4. Dataset selection rationale

The commissioning brief asks for data that is new and harder relative to
devlog 11. Search of every BAM on disk (read-span measured with pysam, not
assumed) gave:

| BAM | chr21 span covered |
|---|---|
| `giab_chr21/chr21_slice.bam` | 0–5.0 Mb (acrocentric/unmappable) |
| `giab_chr21_heldout_30M/chr21_heldout.bam` | 0–0.47 Mb (mock coordinates) |
| `giab_hg002_real_{30M,15x,30x}/…30M…` | 30.00–30.47 Mb → **already scored** (test_15x/test_30x) |
| `giab_hg002_real_{train,15x,30x}/…31_33M…` | 31.00–33.00 Mb |
| `giab_hg002_chr21_12Mb/…32_44M_15x.bam` | 32.00–44.00 Mb |

Leakage constraint: the router cutoff was selected once on **chr21:32–40 Mb**
(devlog 13), so 32–40 Mb is forbidden. `robustness_v11` already scored
13–29 Mb, 30.00–30.47 Mb and 40–44 Mb. Intersecting "on disk", "not router
tuning territory", "not already scored" leaves exactly **chr21:31.0–32.0 Mb**.
That is the region used. It is 1.0 Mb — small — so the design compensates for
power in two pre-registered ways:

1. **Three depth realisations of the same locus set** (15x / 30x / full). Full
   depth (~70–80x) is a regime devlog 11 never tested at all (it topped out at
   30x), and each realisation is an independent read sample.
2. **New hard strata cuts** on the 7 already-cached regions of devlog 11 — the
   brief's explicitly sanctioned fallback. These are *new cuts of previously
   pooled data*, not new regions, and are labelled as such everywhere. They buy
   large-n statements about exactly the regimes devlog 11 only reported in
   aggregate.

Note on 31–33 Mb: it is the historical *training* region for the older
`checkpoints_pb_residual` models and for `cache/train_15x_*`. It is **not**
part of the router's cutoff-selection region and contains no router-fitted
parameter, so it is valid for Stage 1. For Stage 2 it is only valid for the
`checkpoints_big15x` models (trained on 32–38 Mb); the older checkpoints are
treated as leaked (§3).

## 5. Frozen configuration

```python
FROZEN_BINOMIAL_THRESHOLD = 7.0
FROZEN_PB_THRESHOLD       = 10.5
FROZEN_ROUTER_CUTOFF      = 5.411872376933351
route_to_pb = abs(binomial_llr - 7.0) <= 5.411872376933351
final_call  = pb_llr >= 10.5 if route_to_pb else binomial_llr >= 7.0
```

`cascade.py`, `cheap_router.py`, `binomial_baseline.py`, `model_pb_residual.py`,
`residual_prior.py`, `locus_evidence.py`, `test_cascade.py`,
`test_cheap_router.py` were verified byte-identical (`diff -q`) to the
worktree source of the `robustness_v11` provenance record.
`md5sum`: `cascade.py 7b55630dcd7d7a7531f38ce5a3a6a3f4`,
`cheap_router.py e23b0a0c75b77c709fcdf5f738ba93ef`. **No frozen file was
edited by this experiment.**

Stage-2 Mamba decision rule is also frozen, not fitted: `mamba_score >=
FROZEN_PB_THRESHOLD` (the residual head is zero-initialised onto the PB LLR
axis by construction, `model_pb_residual.py`), exactly as devlog 16 used it.

## 6. Experimental protocol (pre-registered)

### 6.1 Stage 1
Arms: **binomial-only**, **PB-only**, **frozen router → PB**. Scoring frame:
`label ∈ {0 (Normal), 1 (SNP)}`, indels/no-calls excluded (project-wide
convention). Datasets: the three new depth realisations, plus the four
pre-registered hard strata applied to each of the 7 devlog-11 regions and
their pool:

* `hard_low_depth` — `depth <= 9`
* `hard_pb_uncertain` — `|pb_llr - 10.5| <= 13.5`
* `hard_low_vaf` — `0 < vaf <= 0.35` (`residual_metrics.vaf_from_counts`)
* `hard_low_quality` — `mean_base_quality < 30` at `depth > 0`

Metrics per arm: F1, precision, recall, TP, FP, FN, TN, accuracy, balanced
accuracy; plus PB locus fraction, disagreement count/rate vs PB, ΔF1 vs PB
with 10,000-resample paired bootstrap (`evaluate_quality_error.paired_bootstrap`).
Depth bins `(1-4, 5-9, 10-14, 15-19, 20-29, 30+)` and quality bins
`(<20, 20-25, 25-30, 30-35, >=35)`. Diagnostic cutoff sweep over
`{0.5 … 20.0}` including the frozen point — **diagnostic only, the production
cutoff is not re-selected from it.**

Compute: same-session measurement of *both* binomial throughput and PB
throughput on the new region (devlog 11 could only mix a fresh binomial rate
with a cross-session PB rate; closing that gap is a pre-registered goal here).

### 6.2 Stage 2
Candidate pool, defined on caller output only, never on labels:
`|pb_llr - FROZEN_PB_THRESHOLD| <= 13.5` (identical to devlog 16's
`CANDIDATE_MARGIN`; adopted unchanged rather than re-searched). Because
`extract_bench_v12.py` retained 44-channel features for **every window
containing a pool locus**, the pool is *feature-unconstrained* here — the
limitation devlog 16 had to declare on its gen14 arm.

Stage-2 routing sweep (pre-registered, identical to devlog 16, not extended):
`t ∈ {0.5, 1, 2, 3, 5, 8, 13.5}` with `mamba_routed = |pb_llr - 10.5| < t`.

**Two pre-declared reporting points, both fixed before results:**
* `t = 0.5` — the cheap, production-like point (what devlog 15/16's selection
  rule converges to).
* `t = 13.5` — the whole pool, the only point in the sweep with the power the
  brief demands. Declared here as the **primary test of H3**, because devlog
  16 demonstrated that a "cheapest non-degrading threshold" selection rule
  mechanically re-derives a ~30-locus operating point and can therefore never
  answer the scientific question. No selection rule is run in this study; the
  whole sweep is reported.

**Decision-threshold addendum (added while writing the Stage-2 script, still
before any Stage-2 metric was computed):** loading a checkpoint revealed it
carries its own `best_threshold` (7.0 for `big15x_seedA`), selected on *its
own* validation split during training and therefore not leaking this region.
Devlog 16 used `FROZEN_PB_THRESHOLD = 10.5` for the Mamba decision. Rather
than pick one after seeing results, **both are declared now**: primary =
10.5 (unchanged from devlog 16); secondary diagnostic = each checkpoint's own
stored validation threshold, reported for every checkpoint in
`own_threshold_diagnostic`. Neither is refit here.

Checkpoints: **all four** `big15x` checkpoints (seedA/B/C + shuffled control),
none discarded regardless of result, plus the three leaked
`checkpoints_pb_residual` models reported in a separate, explicitly-marked
block. Per checkpoint: decisions changed, PB errors corrected, new errors
introduced, ΔTP/ΔFP/ΔFN, ΔF1 + 95% CI, precision, recall, balanced accuracy,
depth stratification.

Controls at matched Mamba coverage: random routing (5 seeds), reverse-
confidence routing (route the *most* PB-confident loci), and the
shuffled-prior checkpoint.

### 6.3 Stage 3
Arms on the new region's full frame: binomial-only; PB-only; router→PB;
PB→Mamba; router→PB→Mamba (`cascade.run_cascade`, sweeping the same frozen
`t` grid); Mamba-on-pool/PB-elsewhere as the "Mamba-only" stand-in (the
44-channel features exist only for pool windows, so a literal whole-frame
Mamba-only arm is **not computable** and is reported as such rather than
approximated silently). Per arm: full accuracy block, fraction of loci
reaching each stage, windows touched, window amplification factor, PB compute,
Mamba compute, wall-clock, speedup vs PB-only, ΔF1 vs PB + CI, disagreement
rate vs PB. Accuracy-vs-compute curve over the sweep.

### 6.4 Stopping criteria (pre-registered)
Stop after: Stage 1 on the three new depth realisations + the hard strata;
Stage 2 over the full pre-registered sweep on all four non-leaked checkpoints
plus controls; Stage 3 over the same sweep; compute accounting; test suite.
**No arm is re-run with a modified threshold after its result is seen.** An
inconclusive or negative result terminates the study — it does not trigger
another cycle.

---

## 7. Stage 1 results

### 7.1 The new region, three depth realisations (chr21:31.0–32.0 Mb, 925,418 scored loci, 924 SNP each)

| Depth arm | Arm | F1 | Precision | Recall | TP | FP | FN | Balanced acc. | PB locus fraction | Disagreements vs PB |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **15x** | binomial-only | 0.96164 | 0.9601 | 0.9632 | 890 | 37 | 34 | 0.98158 | 0% | — |
| | PB-only | 0.98206 | 0.9869 | 0.9773 | 903 | 12 | 21 | 0.98863 | 100% | — |
| | **router→PB** | **0.98206** | 0.9869 | 0.9773 | 903 | 12 | 21 | 0.98863 | **0.1885%** | **0** |
| **30x** | binomial-only | 0.97861 | 0.9672 | 0.9903 | 915 | 31 | 9 | 0.99511 | 0% | — |
| | PB-only | 0.99192 | 0.9871 | 0.9968 | 921 | 12 | 3 | 0.99837 | 100% | — |
| | **router→PB** | **0.99192** | 0.9871 | 0.9968 | 921 | 12 | 3 | 0.99837 | **0.0335%** | **0** |
| **full (~69x)** | binomial-only | 0.98560 | 0.9716 | 1.0000 | 924 | 27 | 0 | 0.99999 | 0% | — |
| | PB-only | **0.68407** | 0.9361 | **0.5390** | 498 | 34 | **426** | 0.76946 | 100% | — |
| | **router→PB** | **0.98876** | 0.9778 | 1.0000 | 924 | 21 | 0 | 0.99999 | **0.0070%** | **439** |

**H1 and H2 hold at 15x and 30x: the router→PB arm is bit-exact with PB-only
(0 disagreements out of 925,418 loci at each depth; paired-bootstrap ΔF1 CI
degenerate at [0, 0] because the two label vectors are literally equal), at
0.19% and 0.033% PB locus coverage.** This is a clean generalization of
devlog 11's result onto coordinates no previous devlog had scored.

**At full depth the comparison inverts, and not because of the router.**
PB-only collapses to F1 0.684 (426 false negatives out of 924 SNP) while the
cascade, which consults PB at only 0.007% of loci, stays at 0.989. The 439
"disagreements vs PB" here are 439 loci where the router was *right* and PB
was wrong. The mechanism is diagnosed in §15 — it is a genuine defect in the
PB implementation outside the frozen router, not a router failure, and it was
**not fixed** (per the brief's rules), only characterised.

Pooled over the three depth realisations (2,776,254 locus-evaluations,
2,772 SNP) the numbers are dominated by the full-depth failure and are
reported for completeness only: binomial-only F1 0.97534, PB-only 0.90140,
router→PB **0.98760**, ΔF1 vs PB **+0.0862** [95% CI +0.0785, +0.0942].
**Pooling across a regime where the reference arm is broken is not a
meaningful accuracy statement about the router**; the per-depth rows above are
the result.

### 7.2 New hard strata on the seven devlog-11 regions (new cuts of already-scored data)

Pre-registered strata, pooled across all seven regions (19.4 M loci):

| Stratum | Loci | SNP | binomial F1 | PB-only F1 | router→PB F1 | Disagreements | PB locus fraction | ΔF1 95% CI |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| depth ≤ 9x | 3,006,471 | 5,004 | 0.88415 | 0.92893 | 0.92905 | 1 | 0.659% | [0.0, +0.00035] |
| \|PB−10.5\| ≤ 13.5 | 801,727 | 2,863 | 0.62398 | 0.81006 | 0.81013 | 2 | 3.032% | [−0.00046, +0.00069] |
| VAF ≤ 0.35 | 387,246 | 2,881 | 0.70774 | 0.86610 | 0.86594 | 1 | 5.995% | [−0.00049, 0.0] |
| mean BQ < 30 | 290,310 | 570 | 0.62958 | 0.81463 | 0.81494 | 2 | 2.975% | [−0.00204, +0.00300] |

These are the regimes where the cheap binomial caller is weakest (F1 0.62–0.88
versus PB's 0.81–0.93) — i.e. exactly where a router that mis-triages would be
punished. **The router remains statistically indistinguishable from PB-only in
every stratum** (1–2 disagreeing loci per stratum, all CIs straddling or
touching zero), while paying 0.66%–6.0% PB coverage inside the hard strata
(higher than the 0.11–0.19% whole-region figure, which is the router
behaving as designed: it spends more PB where evidence is genuinely marginal).

### 7.3 Accuracy-vs-compute curve on the new region (diagnostic sweep)

| Cutoff | 15x PB coverage | 15x F1 | 15x disagreements | 30x PB coverage | 30x F1 | 30x disagreements |
|---:|---:|---:|---:|---:|---:|---:|
| 0.5 | 0.0040% | 0.96966 | 29 | 0.0020% | 0.98443 | 16 |
| 2.0 | 0.0143% | 0.97710 | 9 | 0.0066% | 0.98977 | 6 |
| 3.5 | 0.0372% | 0.98099 | 2 | 0.0118% | 0.99137 | 3 |
| 4.0 | 0.0412% | **0.98206** | **0** | 0.0156% | 0.99191 | 2 |
| 5.0 | 0.1104% | 0.98206 | 0 | 0.0251% | **0.99192** | **0** |
| **5.412 (frozen)** | **0.1885%** | **0.98206** | **0** | **0.0335%** | **0.99192** | **0** |
| 8.0 | 0.7089% | 0.98206 | 0 | 0.1040% | 0.99192 | 0 |
| 20.0 | 98.19% | 0.98206 | 0 | 31.45% | 0.99192 | 0 |

The curve is monotone and saturates (0 disagreements) at cutoff ≈ 4.0 at 15x
and ≈ 5.0 at 30x. **The frozen cutoff of 5.412 sits just past saturation at
both depths** — it was selected on a different region (32–40 Mb) and lands in
the correct place here, with margin, without being re-selected. The frozen
cutoff was not moved.

### 7.4 Depth and quality stratification on the new region

15x: F1 rises monotonically with depth (1–4x 0.872, 5–9x 0.964, 10–14x 0.985,
15–19x 0.992, 20–29x 0.995) and **router→PB equals PB-only in every bin.**
30x: same pattern, all bins identical. Full depth: the 30x+ bin (918,369 loci,
918 SNP) is where PB collapses to 0.683 while the cascade holds 0.991.

Base-quality bins reproduce devlog 11 §11's honest negative: >97% of loci sit
in q30–35 and q≥35, the q<20 bin holds 13–545 loci with no SNP at all. The
low-quality regime **still cannot be tested** with data in this repository.
The `mean BQ < 30` hard stratum in §7.2 (290,310 loci, 570 SNP pooled across
the seven old regions) is the closest available proxy and shows no router
degradation.

### 7.5 Controls

At matched PB coverage on the new region (15x: 1,744 loci; 30x: 310):

| Arm | 15x F1 | 30x F1 |
|---|---:|---:|
| Frozen router | **0.98206** | **0.99192** |
| Random routing (5 seeds, mean) | 0.96186 | 0.97861 |
| Reverse-confidence routing | 0.96164 | 0.97861 |

Both controls collapse to the binomial-only baseline, ~2.0 F1 points (15x) and
~1.3 points (30x) below the frozen router at identical PB spend. H3 replicates
on new coordinates.

## 8. Stage 2 results

### 8.1 Candidate pool

Pool rule `|pb_llr − 10.5| ≤ 13.5`, feature-unconstrained (this benchmark's
extractor keeps the 44-channel features for every window containing a pool
locus, removing devlog 16's gen14 limitation).

| Dataset | Region loci | Pool loci | Pool SNP |
|---|---:|---:|---:|
| chr21:31–32 Mb @15x | 925,952 | 37,873 | **69** |
| chr21:40–44 Mb @15x (powered replication) | see §8.4 | see §8.4 | see §8.4 |

**The 31–32 Mb pool is 37,873 loci but only 69 SNP.** Locus count is not the
binding constraint on F1 power — SNP count is — and 69 positives is far below
the "hundreds to 1000+" the brief asks for. This was visible from the pool
size alone, so a second, larger region (40–44 Mb, the largest block on disk
that is outside the router's 32–40 Mb tuning region and outside the Mamba
checkpoints' 32–38 Mb training split) was extracted with the same frozen
protocol purely for power. That decision was made on a sample-size number, not
on any effect direction, and both results are reported.

### 8.2 Sweep on chr21:31–32 Mb (PB-only pool F1 = 0.76336, TP 50 / FP 12 / FN 19)

| Seed | t=0.5 (n=10) | t=2 (n=33) | t=3 (n=69) | t=8 (n=1,446) | t=13.5 (n=37,873) |
|---|---:|---:|---:|---:|---:|
| seedA | +0.0023 | +0.0163 | +0.0125 | +0.0019 | −0.0090 [−0.076, +0.055] |
| seedB | +0.0023 | +0.0466 | +0.0533 | +0.0533 | **+0.0533 [+0.012, +0.098]** |
| seedC | +0.0023 | +0.0163 | +0.0086 | −0.0023 | −0.0134 [−0.085, +0.054] |
| **shuffled control** | **+0.0366** | **+0.0959** | +0.0108 | −0.5725 | −0.6322 [−0.706, −0.547] |

(values are ΔF1 vs PB-only within the pool)

**The control is the headline here.** The shuffled-prior checkpoint — a model
whose prior was deliberately destroyed, included precisely as a negative
control — produces the **largest positive ΔF1 of any arm** at t ≤ 2
(+0.037 to +0.096, CI excluding zero) before collapsing catastrophically at
larger coverage. Any "improvement" observed at small t therefore cannot be
attributed to learned signal: at 10–33 routed loci out of 69 pool SNP, a model
merely needs to shift borderline calls in the direction that happens to help
this particular label vector. **Small-t positive deltas in this pool are noise
that the control reproduces or exceeds.** Devlog 16 could not see this because
it did not run the shuffled control at the informative coverages.

At the well-powered end of the sweep (t = 13.5, whole pool), seedA and seedC
are negative-but-not-significant, seedB is significantly positive (+0.053).
With 69 SNP, this is **inconclusive/underpowered**, and the three seeds
disagree in sign.

### 8.3 Secondary diagnostic: each checkpoint's own validation threshold

| Checkpoint | own threshold | pool F1 | ΔF1 vs PB | 95% CI |
|---|---:|---:|---:|---|
| big15x_seedA | 7.0 | 0.80000 | +0.0366 | [−0.0011, +0.0782] |
| big15x_seedB | 9.5 | 0.82540 | +0.0620 | [+0.0245, +0.1073] |
| big15x_seedC | 7.5 | 0.79070 | +0.0273 | [0.0000, +0.0620] |
| big15x_shuffled_control | 9.5 | 0.11723 | −0.6461 | [−0.7196, −0.5618] |
| LEAKED_pb_residual_v1_seedA | 10.5 | 0.79032 | +0.0270 | [−0.0081, +0.0626] |
| LEAKED_pb_residual_v1_seedB | 12.5 | 0.77966 | +0.0163 | [−0.0380, +0.0686] |
| LEAKED_pb_residual_v2_gate_seedA | 12.5 | 0.77966 | +0.0163 | [−0.0380, +0.0686] |

Under each checkpoint's own threshold the three real seeds are positive
(+0.027 to +0.062) and the shuffled control is catastrophic — a different
picture from the frozen-10.5 rule, which suggests **most of what the frozen
rule measures at t=13.5 is threshold mismatch, not model quality**. This is
reported as the pre-declared secondary diagnostic (§6.2 addendum); it is not
promoted to the primary result, and no threshold was selected from it. The
three `LEAKED_*` rows are trained on this very region and carry no evidential
weight; they are shown only so that no seed is silently dropped.

### 8.4 Powered replication on chr21:40–44 Mb @15x (the primary Stage-2 evidence)

Extracted fresh with `extract_bench_v12.py --features scores` (3,372,544 loci,
5,943 SNP). Held out from **every** checkpoint used: `big15x_*` were trained on
32–38 Mb, `pb_residual_*` on the 31–33 Mb training split — so on this region
all seven checkpoints are legitimately out-of-training (the `LEAKED_` prefix in
the JSON keys refers to the 31–32 Mb dataset and does **not** apply here; the
names were fixed before this region was added and were deliberately not
renamed after the fact).

**Candidate pool: 232,673 loci, 854 SNP** — the "hundreds to 1000+" power
target of the brief is met for the first time in this project's history.
PB-only pool F1 = 0.77524 (TP 576 / FP 56 / FN 278).

ΔF1 vs PB-only *within the pool*, frozen decision rule (Mamba ≥ 10.5):

| t | routed | seedA | seedB | seedC | shuffled control |
|---:|---:|---|---|---|---|
| 0.5 | 53 | −0.0056 [−0.0126,+0.0008] | −0.0005 [−0.0060,+0.0046] | −0.0056 [−0.0126,+0.0008] | **+0.0070 [+0.0001,+0.0143]** |
| 1.0 | 89 | **−0.0180 [−0.0279,−0.0084]** | +0.0008 [−0.0056,+0.0072] | **−0.0180 [−0.0279,−0.0084]** | **+0.0108 [+0.0028,+0.0195]** |
| 2.0 | 171 | **−0.0295 [−0.0424,−0.0173]** | +0.0012 [−0.0066,+0.0090] | **−0.0268 [−0.0392,−0.0150]** | **+0.0162 [+0.0048,+0.0277]** |
| 3.0 | 308 | **−0.0486 [−0.0648,−0.0334]** | +0.0005 [−0.0083,+0.0092] | **−0.0433 [−0.0586,−0.0286]** | +0.0029 [−0.0106,+0.0171] |
| 5.0 | 930 | **−0.0559 [−0.0739,−0.0389]** | +0.0002 [−0.0090,+0.0093] | **−0.0550 [−0.0719,−0.0387]** | **−0.0693** |
| 8.0 | 4,899 | **−0.0601 [−0.0786,−0.0426]** | +0.0002 [−0.0090,+0.0093] | **−0.0668 [−0.0849,−0.0494]** | **−0.3353** |
| 13.5 | 232,673 | **−0.0644 [−0.0833,−0.0466]** | −0.0007 [−0.0099,+0.0085] | **−0.0735 [−0.0925,−0.0556]** | **−0.5219** |

Paired change breakdown at t = 13.5: seedA fixes 45 PB errors and introduces
99 (net −54); seedC fixes 35, introduces 103 (net −68); seedB fixes 28,
introduces 20 (net +8, not significant). **The direction is consistent with
the F1 deltas and with `History/16_DEVLOG.md`'s ceiling diagnostic, now
replicated on a region that study only evaluated in feature-constrained form,
with a decision-relevant sample size.**

Under each checkpoint's own validation threshold (secondary diagnostic):
seedA +0.0063 [−0.0018,+0.0146] (ns), seedB +0.0110 [+0.0023,+0.0199] (small
significant positive), seedC −0.0012 (ns), pb_residual_v1_seedA −0.0131
(significant negative), pb_residual_v1_seedB / v2_gate_seedA −0.0313
(significant negative), shuffled −0.5524. **Across six effectively
independent checkpoints, one shows a +0.011 F1 gain, two show nothing, and
three show significant losses.** That is not a benefit; it is a coin flip with
a downward bias.

**Depth stratification (seedA, whole pool):** Mamba is worse than PB in every
populated bin except 15–19x (1–4x 0.356 vs 0.519; 5–9x 0.727 vs 0.797;
10–14x 0.815 vs 0.880; 15–19x 0.839 vs 0.831; 20–29x 0.923 vs 0.963). The
damage is concentrated at low depth, which is exactly where a neural refiner
would have to earn its place.

## 9. Stage 3 results — the complete cascade

Whole-frame arms on chr21:40–44 Mb (3,370,737 scored loci, 5,943 SNP):

| Arm | F1 | Precision | Recall | TP | FP | FN | Loci reaching PB | Loci reaching Mamba | ΔF1 vs PB (95% CI) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| binomial-only | 0.94541 | 0.9681 | 0.9238 | 5490 | 181 | 453 | 0 | 0 | −0.0238 |
| **PB-only** | 0.96925 | 0.9900 | 0.9494 | 5642 | 57 | 301 | 100% | 0 | — (reference) |
| **router→PB (2-stage)** | **0.96925** | 0.9898 | 0.9495 | 5643 | 58 | 300 | **0.169%** | 0 | +0.000004, 1 disagreeing locus |
| 3-stage, seedB, t=0.5 | 0.96939 | — | — | — | — | — | 0.169% | 51 | +0.00015 [−0.00055,+0.00084] |
| 3-stage, seedB, t=13.5 | 0.96942 | — | — | — | — | — | 0.169% | 5,595 | +0.00017 [−0.00097,+0.00128] |
| 3-stage, seedA, t=13.5 | 0.96311 | — | — | — | — | — | 0.169% | 5,595 | **−0.00614 [−0.00814,−0.00427]** |
| 3-stage, seedC, t=13.5 | 0.96241 | — | — | — | — | — | 0.169% | 5,595 | **−0.00684 [−0.00884,−0.00490]** |
| Mamba-on-pool/PB-elsewhere (seedA) | 0.96399 | — | — | — | — | — | 100% | 6.9% of frame | −0.00526 |
| Mamba-on-pool/PB-elsewhere (seedB) | 0.96974 | — | — | — | — | — | 100% | 6.9% of frame | +0.00049 |

A literal **Mamba-only whole-frame arm is not computable** with the artefacts
on disk (44-channel features exist only for windows containing a PB-uncertain
locus; materialising them for 3.37 M loci would be a fresh full-region
extraction with no scientific question attached). This is stated rather than
approximated silently; the `Mamba-on-pool/PB-elsewhere` rows are the closest
computable stand-in and are labelled as such.

**Verdict for Stage 3: adding the Mamba stage moves the complete cascade
nowhere useful.** The best case (seedB) is a ΔF1 of +0.0002 with a CI
straddling zero; two of the three checkpoints are significantly *negative*;
and the third-stage costs 31.7–61.0x window amplification for the privilege.
The two-stage `router→PB` cascade is the Pareto-optimal arm in this
experiment: it matches PB-only accuracy (1 disagreeing locus in 3.37 M) at
0.169% of PB's locus cost, and nothing that adds Mamba improves on it.

## 10. Accuracy-vs-compute curves

**Stage-1 curve (new region, §7.3):** monotone, saturating at cutoff ≈4.0
(15x) / ≈5.0 (30x), frozen cutoff 5.412 sits just past saturation with margin
at both depths, costing 0.19% / 0.03% of loci.

**Stage-3 curve (40–44 Mb, Mamba coverage vs accuracy):**

| Mamba loci | Windows touched | Amplification | seedA ΔF1 | seedB ΔF1 | seedC ΔF1 |
|---:|---:|---:|---:|---:|---:|
| 0 (router→PB) | 0 | — | +0.000004 | +0.000004 | +0.000004 |
| 51 | 48 | 60.2x | −0.00031 | +0.00015 | −0.00031 |
| 297 | 272 | 58.6x | −0.00468 | +0.00027 | −0.00413 |
| 2,970 | 1,809 | 39.0x | −0.00578 | +0.00026 | −0.00630 |
| 5,595 | 2,769 | 31.7x | −0.00614 | +0.00017 | −0.00684 |

Every unit of Mamba compute buys either nothing (seedB) or a loss
(seedA/seedC). **The frontier is a single point at zero Mamba compute.**

## 11. Controls

| Control | Stage | Result |
|---|---|---|
| Random routing (5 seeds), matched PB coverage | 1 | Collapses to binomial-only (0.9619 vs router 0.9821 at 15x); ~2 F1 points below the router at identical spend |
| Reverse-confidence routing, matched PB coverage | 1 | Same collapse (0.9616) |
| Random routing, matched Mamba coverage | 2 | Degenerate — at 53/232,673 the random pick lands on no decision-relevant locus, so ΔF1 is exactly 0 (same degeneracy devlogs 15/16 report) |
| Reverse-confidence routing to Mamba | 2 | Degenerate for the same reason |
| **Shuffled-prior checkpoint** | 2, 3 | **The informative control, and it fires.** Significantly *positive* at t ≤ 2 (+0.007 to +0.016 on 40–44 Mb; +0.037 to +0.096 on 31–32 Mb), then catastrophically negative at full coverage (−0.52). Any small-coverage "improvement" from a real checkpoint is therefore not distinguishable from a model with a destroyed prior. |

The Stage-1 controls confirm the router carries real information. The Stage-2
controls confirm the opposite for the neural stage: at the coverages where
Mamba looks good, a deliberately broken model looks better.

## 12. Statistical analysis

* All comparisons paired: every arm scored on identical loci, differences
  taken locus-wise, 10,000-resample paired bootstrap
  (`evaluate_quality_error.paired_bootstrap`, unchanged, seed fixed).
* **Statistically significant positive:** none of the Mamba arms under the
  frozen decision rule. One (seedB, own-threshold diagnostic, +0.011 on
  40–44 Mb) is significant but is not replicated by any other checkpoint and
  is not the pre-registered primary rule.
* **Statistically indistinguishable from PB:** `router→PB` on every dataset at
  15x/30x (CI [0,0] where the arms are bit-identical; [−0.0005,+0.0007] on the
  largest hard stratum); Mamba seedB at every sweep point.
* **Statistically significant negative:** Mamba seedA and seedC at t ≥ 1.0 on
  the powered 40–44 Mb pool (up to −0.073), and the three `pb_residual`
  checkpoints under their own thresholds (−0.013 to −0.031).
* **Underpowered / inconclusive:** the whole 31–32 Mb Stage-2 arm (69 pool
  SNP); the seedA-vs-seedC "replication" (their scores are numerically
  identical, so they are one replicate, reproducing `History/16_DEVLOG.md` §10
  — verified again here: identical F1 to 5 decimals at every sweep point).
* Sample sizes are stated next to every number rather than hidden behind
  aggregates over millions of easy loci: the decision-relevant n for Stage 2
  is **854 pool SNP / 232,673 pool loci** (40–44 Mb) and **69 / 37,873**
  (31–32 Mb), not the 3.37 M-locus frame.

## 13. Depth and quality stratification

Covered in §7.4 (Stage 1, new region), §7.2 (hard strata over 19.4 M loci) and
§8.4 (Stage 2 by depth). Summary: the router is depth- and
quality-insensitive within the tested envelope (identical to PB in every
populated bin at 15x/30x); Mamba's damage is concentrated at 1–14x. The
low-base-quality regime remains untestable with the data in this repository —
the `mean BQ < 30` proxy stratum (290,310 loci) is as close as the GIAB
high-confidence BED allows, and shows no router degradation.

## 14. Compute accounting

**Same-session throughput** (both rates measured in one process, on a
124,032-locus slice of the region being scored — this closes
`History/11_DEVLOG.md` §13's acknowledged gap of mixing a fresh binomial rate
with a cross-session PB rate):

| Dataset | Binomial loci/s | PB loci/s | PB-only projected | Cascade projected | Speedup |
|---|---:|---:|---:|---:|---:|
| 31–32 Mb @15x | 2,130,328 | 8,099 | 114.3 s | 0.65 s | **175.9x** |
| 31–32 Mb @30x | 1,900,291 | 6,320 | 146.4 s | 0.54 s | **273.2x** |
| 31–32 Mb @full | 2,111,127 | 5,753 | 160.9 s | 0.45 s | **357.7x** |

Cross-check: the projected PB-only time (114.3 s at 15x) is within 30% of the
independently measured PB time in the extraction run that built the same cache
(87.7 s), and the discrepancy is in the conservative direction (the throughput
sample ran while another extraction competed for CPU). The full-depth speedup
figure is arithmetically correct but **should not be quoted as a benefit**,
since PB is broken in that regime (§15).

**Locus-level PB reduction (the number this benchmark stands behind):**
0.0070%–0.1885% of loci routed to PB on the new region (530x–14,000x fewer PB
evaluations), 0.169% on 40–44 Mb, 0.66%–6.0% inside the deliberately hard
strata.

**Neural compute is accounted in windows, not loci**, as the brief requires:

| Stage-2 coverage | Mamba loci | Windows touched | Window-amplified loci | Amplification |
|---:|---:|---:|---:|---:|
| t = 0.5 | 51 | 48 | 3,072 | 60.2x |
| t = 3.0 | 297 | 272 | 17,408 | 58.6x |
| t = 8.0 | 2,970 | 1,809 | 115,776 | 39.0x |
| t = 13.5 | 5,595 | 2,769 | 177,216 | 31.7x |

At sparse routing, amplification saturates near `SEQ_LEN = 64` (each routed
locus sits in its own window); it falls to 31.7x only when routing is dense
enough for windows to overlap. **A "0.0015% of loci reach Mamba" claim would
be a 60x understatement of the real cost.** Measured Mamba inference
throughput on this machine (CPU, 441,088 loci scored per checkpoint):
0.33–0.48 s, i.e. ~1.0–1.3 M window-loci/s — cheap per locus, but it is the
amplified locus count that must be multiplied.

Preprocessing is the honest elephant: producing the read tensor PB needs took
117 s / Mb at 15x and 590 s / Mb at full depth, dwarfing both callers. The
cascade only needs read-level evidence at routed loci, so the true end-to-end
saving is larger than the caller-only arithmetic above; that larger claim is
**not** quantified here because the extractor materialises reads chunk-wise for
all loci and would need restructuring to measure it honestly.

## 15. Failure analysis — a real defect found in the PB implementation

**The most consequential finding of this benchmark is not about the router.**

At full depth (~69x mean) PB-only's recall collapses to 0.539: 426 of 924 true
SNP are missed, and every one of them has `pb_llr` **exactly 0.0**. The
distribution of that failure by depth is unambiguous:

| Depth bin | Loci | Fraction with `pb_llr == 0` | SNP in bin | Fraction of SNP with `pb_llr == 0` |
|---|---:|---:|---:|---:|
| 10–30x | 7,048 | 0.0000 | 6 | 0.000 |
| 30–48x | 52,440 | 0.0000 | 73 | 0.000 |
| **48–60x** | 152,232 | 0.0007 | 197 | **0.487** |
| **60–80x** | 511,186 | 0.0005 | 497 | **0.497** |
| **80x+** | 203,014 | 0.0005 | 151 | **0.550** |

The boundary is exactly 48, and 48 is `read_level_pileup.MAX_READS`. The
mechanism, traced through the code:

1. `read_level_pileup` retains at most `MAX_READS = 48` reads per locus, so the
   quality-evidence tensor has 48 columns regardless of true depth.
2. `candidate_alt` computes the ALT count `k` from the **uncapped** pileup
   count matrix, so at 69x a homozygous SNP yields `k ≈ 54`.
3. `poisson_binomial_llr` clips `k` to the tensor width
   (`block_k = np.clip(k, 0, block_p.shape[1])` → 48), i.e. "all 48 observed
   reads support ALT".
4. Under H0 that has log-probability −inf, so the LLR is `+inf`, and the final
   `np.nan_to_num(llr, posinf=0.0)` converts the **strongest possible evidence
   for a variant into exactly 0.0** — a score below the 10.5 threshold, hence a
   confident false negative.
5. Empirically: of the 426 missed SNP, 76% have `k ≥ 48`, and none of the
   PB-scored SNP with `pb_llr ≠ 0` do.

**This was not fixed.** It sits in `quality_error_model.poisson_binomial_llr`
/ `read_level_pileup`, outside the frozen router, and the brief's rules say to
document rather than repair a defect that is itself a benchmark finding.
Fixing it (either by not capping reads, by scaling `k` to the retained sample,
or by mapping `+inf` to a large positive rather than 0) is a concrete,
testable next task, and every result in devlogs 3–11 that used ≤30x data is
unaffected because the failure requires depth > 48.

**Consequence for the architecture, stated plainly:** the frozen cascade was
*accidentally* protected. It routes only 0.007% of full-depth loci to PB
(because at high depth the binomial LLR is far from its threshold), so PB's
defect reaches almost nothing, and the cascade scores 0.98876 where its own
"accuracy ceiling" scores 0.68407. That is luck, not design — the router has
no depth awareness — and it is reported as luck.

Other failures found: **1 disagreeing locus** between `router→PB` and PB-only
in 3.37 M loci on the freshly re-extracted 40–44 Mb (devlog 11 reported 2 on
its own extraction of the same coordinates; the SNP counts differ slightly,
5,943 here vs 5,951 there, because this is an independent extraction from a
different BAM slice, so this is a qualitative replication, not a bit-exact
one). No router failure was found at 15x or 30x on the new region, or in any
hard stratum beyond the 1–2 loci already reported.

## 16. Unexpected findings

1. **PB breaks above 48x depth** (§15) — the benchmark set out to test the
   router and found a defect in the component the router is supposed to be
   approximating.
2. **The shuffled-prior control beats the real checkpoints at small Mamba
   coverage** (§11) — an artefact-detector that fires, and that retroactively
   explains why several previous devlogs saw small, unreplicable positive
   deltas at tiny routing volumes.
3. **`best_big15x_seedA.pt` and `best_big15x_seedC.pt` again produce
   numerically identical scores** on two datasets neither had seen, confirming
   `History/16_DEVLOG.md` §10's anomaly on fresh data. Three "seeds" are two.
4. **The frozen router's PB coverage falls with depth** (0.19% → 0.033% →
   0.007%), extending devlog 11's `test_30x` observation to a third depth: the
   cascade gets *cheaper* as data gets better, which is the opposite of most
   cost models.
5. `positions` written by `extract_region_pb.py` (and inherited by this
   benchmark's extractor) are *cache-order indices offset by the chunk start*,
   not true genomic coordinates, because BED filtering drops loci. They are
   monotone and therefore safe for the index-mapping this benchmark does, but
   they must not be read as coordinates. Documented, not changed.

## 17. Reproducibility information

```bash
# 1. extract the new region (three depth realisations)
python3 extract_bench_v12.py --fasta data/reference/chr21_full.fa \
  --bam data/giab_hg002_real_15x/hg002_chr21_31_33M_15x.bam \
  --vcf data/giab_hg002_real_train/hg002_chr21_31_33M.vcf.gz \
  --bed data/giab_hg002_real_train/hg002_chr21_31_33M_highconf.bed \
  --region 31000000 32000000 --out cache/bench_v12/new31_32M_15x.npz
#   ... same with the 30x / full BAMs and --features none
# 2. extract the powered Stage-2 region (Mamba scored in-process)
python3 extract_bench_v12.py --fasta data/reference/chr21_full.fa \
  --bam data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_15x.bam \
  --vcf data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz \
  --bed data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed \
  --region 40000000 44000000 --features scores --out cache/bench_v12/pool40_44M_15x.npz
# 3. evaluate
python3 bench_v12_stage1.py
python3 bench_v12_stage23.py
python3 bench_v12_stage23.py --region-npz cache/bench_v12/pool40_44M_15x.npz \
  --out results/bench_v12/stage23_pool40_44M_results.json --score-cache results/bench_v12/unused_40_44M.npz
python3 -m pytest test_bench_v12.py test_cascade.py test_cheap_router.py -q
```

| Artefact | Path |
|---|---|
| Stage-1 results | `results/bench_v12/stage1_results.json` |
| Stage-2/3, 31–32 Mb | `results/bench_v12/stage23_results.json` |
| Stage-2/3, 40–44 Mb (powered) | `results/bench_v12/stage23_pool40_44M_results.json` |
| Cached Mamba scores, 31–32 Mb | `results/bench_v12/mamba_scores_new31_32M_15x.npz` |
| Evidence caches | `cache/bench_v12/*.npz` (gitignored) |
| New code | `extract_bench_v12.py`, `bench_v12_stage1.py`, `bench_v12_stage23.py`, `test_bench_v12.py` |

Software state: Python 3.14, numpy/torch/pysam as installed in the working
environment; CPU inference (no CUDA device present). Bootstrap seeds are the
module defaults (`BOOTSTRAP_SEED` in `evaluate_quality_error`, 20260811 in
`residual_metrics`), unchanged. Every frozen constant is imported, never
redefined. **No historical artefact was modified**: `results/robustness_v11/*`
and `History/3..11_DEVLOG.md` are byte-unchanged, and this benchmark writes
only into `results/bench_v12/`, `cache/bench_v12/` and this file.

**Test status:** `test_bench_v12.py` 9 new tests pass. Frozen-module suites
`test_cascade.py` + `test_cheap_router.py` (40 tests) pass unmodified. Full
project suite, excluding the two pre-existing collection errors on
`test_bimamba.py` / `test_providers.py` (`ModuleNotFoundError:
bimamba_variant_caller`, present on `main` before this work and unrelated to
it): **208 passed, 0 failed** (199 before this experiment + 9 new).

## 18. Limitations

* Still one sample (HG002), one reference (GRCh38), one chromosome (21). No
  cross-sample, cross-platform or cross-chromosome claim is made or implied.
* The genuinely-unscored, non-tuning territory available on disk was **1.0 Mb**
  (924 SNP). Statistical strength for Stage 1 comes from three depth
  realisations of it plus 19.4 M loci of new hard-stratum cuts, not from a
  large new region.
* The 40–44 Mb Stage-2 region was previously scored by devlog 11 at Stage 1; it
  is new only as a *Stage-2, feature-unconstrained* dataset. It was added after
  the 31–32 Mb pool turned out to hold 69 SNP — a decision driven by sample
  size, not by any observed effect direction, and disclosed here rather than
  presented as part of the original plan.
* `seedA` and `seedC` are one effective replicate (§16.3), so Stage 2 rests on
  two independent `big15x` checkpoints plus three `pb_residual` ones, not six.
* No low-base-quality stratum exists in the data (devlog 11 §11's gap is
  unclosed).
* Wall-clock projections combine measured throughputs with locus counts; the
  throughput sample ran under CPU contention from a concurrent extraction, so
  the speedups are conservative estimates rather than isolated-machine numbers.
* End-to-end preprocessing savings (§14) are argued but not measured.
* The full-depth arm's PB numbers describe a broken component, so any pooled
  statistic that mixes it with 15x/30x is uninterpretable and is flagged as
  such wherever it appears.

## 19. Final verdict

**Stage 1 — POSITIVE, and generalises.** On coordinates no previous devlog had
scored (chr21:31–32 Mb) the frozen cheap-router → PB cascade is **bit-exact
with PB-only at 15x and 30x** (0 disagreements in 925,418 loci at each depth)
while sending 0.19% and 0.033% of loci to PB, and it beats random and
reverse-confidence routing by ~1.3–2.0 F1 points at matched PB spend. Across
19.4 M loci of newly-cut hard strata (low depth, low VAF, PB-uncertain, low
base quality — the regimes where the cheap caller is weakest, F1 0.62–0.88) it
remains statistically indistinguishable from PB-only, with 1–2 disagreeing loci
per stratum. The frozen cutoff, selected once on a different region, lands just
past the saturation point of the accuracy-vs-compute curve on this new data.
Same-session throughput measurement gives **176x–357x** projected speedup,
replacing devlog 11's mixed-session estimate.

**Stage 2 — NEGATIVE (not merely "no evidence").** With the power target met
for the first time (232,673 pool loci, **854 pool SNP**, on a region held out
from every checkpoint), PB → Mamba is significantly *harmful* for two of the
three `big15x` checkpoints (ΔF1 −0.018 to −0.073, CIs strictly below zero at
every t ≥ 1.0, introducing 2–3x more errors than it fixes) and null for the
third. Under the checkpoints' own thresholds, six effectively-independent
models give one small significant gain (+0.011), two nulls and three
significant losses. The shuffled-prior control **outperforms** the real models
at low routing volume, which disqualifies small-coverage positive deltas as
evidence. This is evidence of no benefit, not absence of evidence — while
remaining a statement about *these* checkpoints and *this* residual
architecture, not about neural refinement in general.

**Stage 3 — the third stage does not earn its place.** The complete cascade's
best case is ΔF1 +0.0002 (CI straddling zero) and its typical case is
−0.006 significant, at 31.7–61.0x window amplification. The Pareto frontier of
this system is a single point: **binomial → frozen router → PB**, at PB-level
accuracy for 0.17% of PB's cost.

**Unplanned but most important: PB itself is broken above 48x depth** (§15),
silently converting the strongest variant evidence into a zero LLR and losing
~50% of SNP at high depth. Every previous result in this project is unaffected
(all used ≤30x), but any move to deeper data must fix this first.

## 20. Recommendation for the next research phase

1. **Adopt the two-stage cascade (binomial → frozen router → PB) as the
   production component** for ≤30x short-read data on this sample/reference.
   It is now validated on 19.4 M previously-scored loci (devlog 11) plus
   1.85 M loci of genuinely new coordinates at two depths plus 3.37 M loci of
   independent re-extraction, with controls, at 100–350x compute reduction.
2. **Close the PB depth defect before any high-coverage work.** Concrete fix
   candidates in §15; whichever is chosen needs a unit test asserting that a
   locus with `k = n = MAX_READS` produces a large positive LLR, not 0.
3. **Stop investing in the PB-residual Mamba stage.** Three independent studies
   (devlogs 10, 16 and this one) now agree it ranges from null to significantly
   harmful, and this one had the sample size to say so. Any revival should
   start from a different hypothesis about *what information PB lacks* (§15
   suggests one: PB is blind to depth beyond 48 reads), not from retraining the
   same residual head.
4. **The binding constraint on every remaining question is data, not method.**
   The highest-value next acquisition is a second chromosome and a
   low-base-quality / degraded-library region; without them, "does this
   generalise beyond chr21 HG002" cannot be answered by any amount of further
   analysis of what is on disk.
5. Treat the biological/bioinformatics phase as **unblocked for ≤30x
   variant-calling on this data family**, and blocked for high-coverage work
   until item 2 is done.
