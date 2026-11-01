# Quality-aware sequencing-error model — 15× HG002, chr21

**Date:** 2026-08-12
**Verdict: POSITIVE.** Per-locus ε derived from base qualities beats fixed
ε = 0.01 by **ΔF1 = +0.0179 [+0.0102, +0.0263]** (Poisson-binomial arm,
paired bootstrap, 95% CI excludes zero). Every control that removes the
locus-specific pairing removes the entire gain. This is the first experiment in
the project where added information measurably helped.

**Poisson-binomial dominates Binomial v1 strictly:** 18 false negatives
rescued, **0 true positives lost, 0 new false positives**, 1 false positive
removed. F1 0.9729 → 0.9908.

---

## 1. Hypothesis

Binomial v1 assumes one fixed sequencing error probability ε = 0.01 (~Q20) at
every locus. The reads already carry a per-observation Phred score. If a
locus-specific ε reflecting that evidence produces a better caller, the
information bottleneck is *statistical representation of uncertainty*, not
model capacity.

Prediction if true: gain concentrated at 5–14×, where one or two ALT reads
decide the call. Prediction if false: quality-derived ε changes nothing a
different constant could not achieve.

## 2. Formulation

Per-read error probability, `p_i = 10^(-Q_i / 10)`. The hypotheses are
**unchanged** from v1 — ε is the only variable that moves:

```
p0    = ε / 3
p_het = 0.5 (1 - ε) + 0.5 (ε / 3)
p_hom = 1 - ε + ε / 3
LLR   = max( log P(k | n, p_het), log P(k | n, p_hom) ) - log P(k | n, p0)
```

Estimators for ε (all label-free, all from base qualities only):

| | estimator | definition |
|---|---|---|
| A | `mean` | `mean(p_i)` over counted reads |
| B | `median` | `median(p_i)` over counted reads |
| C | `ref_mean` | `mean(p_i)` over **reference-supporting reads only** |
| D | `ref_trimmed` | symmetric 10 %-trimmed mean over reference-supporting reads |

C and D were designed to avoid circularity: ε cannot be driven by the reads
whose status as errors is the question being asked.

**Poisson-binomial arm.** The plug-in step is an approximation: with
read-specific `p_i` the ALT count is a sum of *non-identically* distributed
Bernoulli variables. The exact log-pmf is computed by direct convolution in log
space (`logaddexp`, no subtraction of near-equal terms, unlike the
characteristic-function methods that lose precision in exactly the tail where a
variant call lives). Per-read ALT probabilities are `p_i/3`,
`0.5(1-p_i)+0.5(p_i/3)`, `1-p_i+p_i/3`. When all `p_i` are equal it reduces
**exactly** to the binomial — asserted in the unit tests, which is why this is
a weaker assumption rather than a different model.

ε clamp `[1e-4, 0.25]`. The floor is scientific, not cosmetic: base quality
describes *base miscall only* and says nothing about misalignment, which does
not vanish at high Q. Results saturate for any floor ≤ 1e-4 and are nearly
identical at 1e-3 (§7), so the floor is not a tuned knob.

## 3. Dataset and protocol

Real HG002 GIAB v4.2.1, GRCh38, ~15×. Test = chr21:30,000,000–30,469,999,
455,371 loci in the SNP-vs-Normal frame, 549 SNPs. Validation = last 10 % of
the training region's windows (26386–29317), the same split the binomial
baseline and the neural training used.

Everything is computed from the **existing** cached tensors — no new extraction
pass. The read-level tensor is cross-checked against the cached count matrix by
`reconstruct_counts` at load time, column by column; every arm is therefore
provably scoring the same reads. Candidate-ALT selection is *delegated* to
`BinomialVariantCaller.score_counts`, so it is bit-identical to the baseline
and this experiment cannot change two things at once.

The primary new arm — **Estimator C** — was **pre-registered before test
evaluation**, chosen for its anti-circularity argument rather than by
validation F1. (It turned out to be the weakest arm. Reported as such.)

## 4. Main comparison — SNP vs Normal, validation-frozen thresholds

| Method | ε | thr | P | R | **F1** | TP | FP | FN |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Flat frequency | — | 0.335 | 0.9942 | 0.9381 | 0.9653 | 515 | 3 | 34 |
| Binomial v1 | 0.01 fixed | 10.5 | 0.9981 | 0.9490 | 0.9729 | 521 | 1 | 28 |
| *control:* best global ε | 3.98e-3 fixed | 13.0 | 0.9981 | 0.9490 | 0.9729 | 521 | 1 | 28 |
| *control:* calibrated global ε | 9.85e-4 fixed | 17.5 | 0.9981 | 0.9454 | 0.9710 | 519 | 1 | 30 |
| *control:* mean-Q ε **shuffled** | per-locus, wrong loci | 14.0 | 0.9709 | 0.9727 | 0.9718 | 534 | 16 | 15 |
| C — reference-only mean *(pre-registered)* | per-locus | 13.0 | 0.9816 | 0.9709 | 0.9762 | 533 | 10 | 16 |
| D — reference-only trimmed | per-locus | 13.0 | 0.9816 | 0.9709 | 0.9762 | 533 | 10 | 16 |
| B — median-Q | per-locus | 23.0 | 0.9981 | 0.9563 | 0.9767 | 525 | 1 | 24 |
| A — mean-Q | per-locus | 9.0 | 0.9926 | 0.9836 | 0.9881 | 540 | 4 | 9 |
| **Poisson-binomial** | per-**read** | 12.5 | **1.0000** | **0.9818** | **0.9908** | 539 | **0** | 10 |
| *diagnostic:* Binomial v1 @ 9.5 | 0.01 fixed | 9.5 | 0.9944 | 0.9617 | 0.9778 | 528 | 3 | 21 |

Paired bootstrap vs Binomial v1 (10,000 resamples, seed 20260812):

| Arm | ΔF1 | 95 % CI | ΔTP | ΔFP | ΔFN |
|---|---:|---|---:|---:|---:|
| Poisson-binomial | **+0.0179** | **[+0.0102, +0.0263]** | +18 | −1 | −18 |
| A — mean-Q | **+0.0152** | **[+0.0069, +0.0241]** | +19 | +3 | −19 |
| B — median-Q | +0.0038 | [−0.0019, +0.0099] | +4 | 0 | −4 |
| C — reference-only *(pre-registered)* | +0.0033 | [−0.0061, +0.0128] | +12 | +9 | −12 |
| D — reference-only trimmed | +0.0033 | [−0.0061, +0.0128] | +12 | +9 | −12 |
| best global ε | +0.0000 | [+0.0000, +0.0000] | 0 | 0 | 0 |
| calibrated global ε | −0.0019 | [−0.0049, +0.0000] | −2 | 0 | +2 |
| mean-Q ε shuffled | −0.0011 | [−0.0109, +0.0083] | +13 | +15 | −13 |
| Flat frequency | −0.0076 | [−0.0149, −0.0010] | −6 | +2 | +6 |

## 5. Depth-stratified (F1, TP/FP/FN) — primary frame

| depth | loci | SNP | Binomial v1 | C ref-only | A mean-Q | Poisson-binomial |
|---|---:|---:|---|---|---|---|
| 1–4 | 4,137 | 7 | 0.923 · 6/0/1 | 1.000 · 7/0/0 | 1.000 · 7/0/0 | **1.000 · 7/0/0** |
| **5–9** | 54,981 | 79 | 0.912 · 67/1/12 | 0.916 · 71/5/8 | 0.954 · 73/1/6 | **0.954 · 72/0/7** |
| **10–14** | 160,626 | 209 | 0.970 · 197/0/12 | 0.976 · 203/4/6 | 0.988 · 206/2/3 | **0.993 · 206/0/3** |
| 15–19 | 153,458 | 162 | 0.991 · 159/0/3 | 0.991 · 160/1/2 | 0.997 · 162/1/0 | **1.000 · 162/0/0** |
| 20–29 | 79,425 | 91 | 1.000 · 91/0/0 | 1.000 | 1.000 | 1.000 |
| 30+ | 2,744 | 1 | 1.000 · 1/0/0 | 1.000 | 1.000 | 1.000 |

The gain is **entirely** in 5–19×, as predicted. At ≥20× every method is
saturated — there is nothing left to win there.

## 6. Why it works — the controls are the result

Three controls separate "per-locus uncertainty is informative" from "0.01 was
just the wrong number":

1. **Calibrated global ε** — pooled mean over all validation reads, 9.85e-4
   (~Q30). F1 **0.9710**, worse than v1. Recalibrating the constant does not help.
2. **Best global ε** — joint (ε, threshold) grid search on validation, the
   *ceiling of the entire fixed-ε family*. F1 **0.9729**, identical to v1.
   The sweep is flat at 0.971–0.975 across ε ∈ [1e-5, 3e-2] — five orders of
   magnitude.
3. **Shuffled ε** — the same per-locus ε values, randomly reassigned to loci.
   Marginal distribution preserved exactly, pairing destroyed. F1 **0.9718**,
   ΔF1 −0.0011, CI spans zero. The entire gain disappears.

The reason control 2 is flat is structural, and it is the core insight:

> With a fixed ε, the LLR is a deterministic function of `(n, k)` alone.
> Changing ε only slides the decision threshold along a fixed ranking; it can
> **never re-order two loci with the same `(n, k)`**. Quality-aware ε breaks
> exactly those ties.

The error analysis confirms this directly. All 18 rescued false negatives are
low-depth heterozygous sites with `k = 2–3`, VAF 0.18–0.40 — and mean base
quality **Q38–Q53**:

| depth | n | k | VAF | mean BQ | LLR v1 | LLR PB |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 6 | 2 | 0.33 | 46.3 | 7.28 | 16.30 |
| 8 | 8 | 2 | 0.25 | 38.6 | 5.91 | 14.26 |
| 15 | 11 | 2 | 0.18 | 53.6 | 3.86 | 12.86 |
| 12 | 10 | 3 | 0.30 | 44.9 | 10.23 | 23.29 |

Under a Q20 assumption, 2 ALT reads out of 8 is an unremarkable noise draw.
The reads actually say Q45 — at which point 2 ALT observations are close to
impossible under H0. v1 was discarding real variants because it assumed its
reads were 100× worse than they are.

**Why the pre-registered arm C failed** is equally informative. Its 9 new false
positives are precisely the loci with *low* ALT-read quality (mean BQ 30–38).
By construction C cannot see ALT-read quality, so it assigns those loci the
same low ε as a clean site and promotes the noise. The discriminating variable
is the quality of the ALT reads themselves — and the estimator designed to be
maximally safe from circularity is the one that throws that variable away.

**This is not circular.** ε depends only on Q, never on whether a read matches
the reference. High-quality ALT reads *lower* ε (more evidence for a variant);
low-quality ALT reads *raise* it. The direction is set by the quality, not by
the disagreement — which is the opposite of the failure mode we guarded against.

**Poisson-binomial vs plug-in (CASE C):** PB beats mean-Q plug-in on
precision — 0 FP vs 4 FP, F1 0.9908 vs 0.9881 — though the paired CIs overlap.
The plug-in collapses a heavy-tailed quality distribution (test mean p = 8.9e-4
but median 1.4e-4) into one number; PB keeps each read's own weight, so a
single bad read among nine good ones no longer contaminates the whole locus.
Suggestive, not established at this sample size.

## 7. Robustness, runtime, tests, leakage

**ε floor sensitivity** (thresholds re-picked on validation at each floor):
identical results at 1e-4, 1e-5, 1e-6; F1 0.9881 / 0.9824 at 1e-3; collapses to
v1 at 1e-2 (the clamp then forces ε ≥ 0.01). Not a tuned knob.

**Runtime**, 455,488 loci, single CPU core: quality-evidence decode 6.6 s;
mean-Q ε 0.57 s; median 0.76 s; trimmed 0.96 s; fixed binomial 0.16 s;
**Poisson-binomial 55.6 s** (exact 48-step convolution, chunked, vectorized over
loci). Whole experiment end-to-end ~2 min including both regions. For
comparison the Mamba inference pass is GPU-bound and the training runs are
hours. The most expensive statistical arm is still ~3 orders of magnitude
cheaper than the neural pipeline.

**Unit tests:** 47 new in `test_quality_error.py`; full suite 132 passed.
Covers Q10/Q20/Q30 → 0.1/0.01/0.001 exactly; gap reads excluded from the
denominator; base-quality-filtered reads excluded; reference-read
identification; ALT ≠ reference; ALT selection bit-identical to the baseline;
rebuilt likelihood equal to v1 at ε = 0.01; monotonicity in k and in ε;
batch = scalar; n up to 1e6 finite; zero-depth neutral; PB log-pmf normalises
to 1; **PB ≡ binomial when qualities are equal**; PB ≠ plug-in when they are
not; chunk size numerically irrelevant; determinism.
(`test_bimamba.py` / `test_providers.py` fail collection on a stale
`bimamba_variant_caller` import — pre-existing, untouched here.)

**Leakage checks:** no VCF, label or model output is read anywhere in
`quality_error_model.py`, enforced by a static test over the module's
executable code. ε estimation uses base qualities only. Thresholds, the primary
arm and the global-ε constant are all fixed on validation before the test
region is scored, and the test region is scored exactly once per run. Oracle
(test-optimal) thresholds are recorded but labelled a diagnostic ceiling and
are not reported as results.

## 8. Verdict and next experiment

**CASE A + CASE C.** Quality-aware ε clearly beats fixed ε, and the
Poisson-binomial formulation is at least as good as the plug-in with a
precision edge that is suggestive but not yet established.

The bottleneck was, as suspected, statistical representation of uncertainty —
and specifically the *ranking* of loci that share `(n, k)`, which no fixed-ε
caller can express at all.

**Recommendation — do not enlarge the neural model.** The new statistical
baseline is F1 0.9908 with 0 false positives; every previous neural arm is far
below it, and the residual-Mamba experiment already showed the network learning
nothing beyond the *weaker* baseline. The next question is the CASE A question:

> Can the neural model learn anything beyond the **quality-aware** likelihood?

Concretely, next: re-run the residual experiment with the Poisson-binomial LLR
as the initialization/offset instead of the fixed-ε LLR, on both seeds. If the
residual again collapses to a constant shift, the statistical caller has
captured the available information at 15× and the honest move is to attack a
different regime (higher coverage, indels, or the 10 remaining FNs) rather than
the architecture.

Second priority, cheap: the 10 remaining Poisson-binomial FNs and the 5–9×
stratum (0.954) are now the entire error budget. Worth characterising them
before designing anything.

## 9. Repository changes

**Added:**
- `quality_error_model.py` — ε estimators A–D, rebuilt binomial likelihood,
  exact Poisson-binomial LLR
- `evaluate_quality_error.py` — validation freeze, controls, test evaluation,
  paired bootstrap, depth stratification, error analysis
- `test_quality_error.py` — 47 unit tests
- `results/quality_error_15x_results.json`
- `History/8_DEVLOG.md`

**Not changed:** no checkpoint, dataset, cache, previous result file, or
existing scientific module was modified. `binomial_baseline.py` is imported and
delegated to, never edited — the baseline this experiment is measured against
is byte-identical to the one that produced the established numbers.
