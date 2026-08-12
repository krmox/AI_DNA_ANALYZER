# A regime with a real error budget — PB-residual Mamba at 15× over 3.81 Mb

## WEAK POSITIVE

The neural residual carries information the Poisson-binomial caller does not
model. That information is real, its mechanism is identified, and it is not a
threshold move. But its effect size is not reliably above the noise floor: two
of three seeds beat the statistical caller with confidence intervals excluding
zero, and the third does not.

```
paired bootstrap vs Poisson-binomial, 3.81 Mb test, 5,704 SNPs
  seed A (20260811)  dF1 = +0.0012  [+0.0001, +0.0024]   dFP -12
  seed B (424242)    dF1 = +0.0019  [+0.0009, +0.0030]   dFP -19
  seed C (13371337)  dF1 = +0.0003  [-0.0007, +0.0015]   dFP  -4
```

Seed C was run as an unconditional replication check *after* seeds A and B had
been scored on test, and it is reported for that reason. With only A and B the
result would have been called STRONG POSITIVE. It would have been wrong.

What survives all three seeds: the direction is never negative, every seed
removes false positives, every seed's residual is locus-specific rather than
constant, every seed correlates with `alt_strand_bias`, and every seed beats
its own constant-shift null. What does not survive: the claim that the
improvement is reliably detectable on a 271-error budget.

---

## 1. Why this regime, and what was ruled out first

The 15× and 30× SNP benchmarks were exhausted: PB scored 539/0/10 and 549/2/0,
leaving error budgets of 10 and 2 loci. The instruction was to find a regime
where the statistical caller has a substantial, measurable error budget.

### Indels — audited, and rejected before any code was written

Priority 1 was indels, and the audit killed it. A *trivial* fixed-epsilon
binomial on the raw indel counts, at oracle thresholds:

| regime | class | F1 | TP | FP | FN |
|---|---|---:|---:|---:|---:|
| 15× | DELETION | **1.0000** | 77 | 0 | 0 |
| 15× | INSERTION | 0.9873 | 39 | 0 | 1 |
| 30× | DELETION | **1.0000** | 77 | 0 | 0 |
| 30× | INSERTION | **1.0000** | 40 | 0 | 0 |

One error at 15×, zero at 30×. True indels sit at VAF 0.70–0.86 while normal
loci carry an indel read 0.005–0.05% of the time; the classes barely overlap.
The indel infrastructure (indel candidate generation, a gap-channel prior, an
indel-centric feature builder) was therefore **not built**.

Recorded for the future, since the audit was done properly: truth labelling
already handles indels correctly — insertions on the anchor base, deletions on
the actually-deleted bases with the anchor skipped — and the read tensor already
carries `FLAG_GAP`, `FLAG_INSERTION` and `FLAG_DELETION_NEXT`. What is missing
is that `candidate_alt` is SNP-only, the 44-channel representation is
SNP-centric (3 of 44 channels touch indels), and **the Poisson-binomial prior is
mathematically inapplicable to indels**: its per-read error probability is the
base-quality Phred, and a gap has no base and therefore no Phred. An indel
experiment would need a different statistical baseline, not a reused one.

### Priorities 2–3 — structurally absent from the old region

* **Low-mappability: none.** Minimum mean MAPQ across 455,488 loci is **66.6**,
  median 70, zero loci below 50.
* **Homopolymers: real noise, still separable.** Deletion background
  contamination rises 40× from 0.00051 (homopolymer ≤4) to 0.02226 (≥10), and
  21 of 77 deletions sit in ≥10-mers — but peak background VAF there is 0.105
  against true-variant VAF 0.86.

The common cause is the GIAB high-confidence BED: it **excludes hard regions by
construction**, so the difficulty is filtered out before any model sees the
data. This is the single most important structural fact about this benchmark.

### The regime that did have headroom

Measured PB error density: **32.5 errors/Mb at 15×**, 6.9 at 30×, ~0 for indels.
The binding constraint was never coverage or variant class — it was evaluation
size. A 0.47 Mb test set yields ~10 errors no matter what.

## 2. Dataset

New region acquired, identical protocol to all previous data.

| | |
|---|---|
| Region | chr21:32,000,000–40,000,000 (8 Mb of a 12 Mb download) |
| Source BAM | GIAB HG002 NIST Illumina 2×250, novoalign, GRCh38 |
| Truth | NIST v4.2.1 benchmark VCF — 20,715 records in 32–44 Mb |
| Confidence | v4.2.1 `benchmark_noinconsistent.bed`, 93.06% of the span |
| REF concordance | **20,715 / 20,715 (100%)** against `chr21_full.fa` |
| Downsample | `--subsample 0.203 --subsample-seed 20260811` (the recorded 15× recipe) |
| Measured depth | **13.79×**; zero loci above 48 reads, so `MAX_READS=48` unchanged |
| Materialized | 7,611,136 loci, 10,658 SNP |

Depth distribution (test blocks): 1–4× 39,923 loci; 5–9× 540,188; 10–14×
1,435,554; 15–19× 1,180,977; 20–29× 519,522; 30–39× 10,287; 40–49× 17; 50+ 0.

### Split — designed for validation power

1 Mb blocks assigned by index mod 4 → train / validation / test / test.

| role | Mb | loci | SNP | blocks |
|---|---:|---:|---:|---|
| train | 1.89 | 1,894,656 | 2,544 | 0, 4 |
| validation | 1.91 | 1,907,840 | **2,410** | 1, 5 |
| test | 3.81 | 3,808,640 | **5,704** | 2, 3, 6, 7 |

Two deliberate properties. **Train stays at 1.89 Mb**, the size every previous
experiment used, so evaluation power is the only moved variable and no result
can be attributed to more training data. **Roles interleave** rather than taking
contiguous halves, because difficulty drifts along a chromosome and contiguous
blocks would confound role with region.

Assignment is per 64-locus *window*, so no window straddles two roles. The only
contact between roles is at block boundaries, where a ≤250 bp read can overlap
the neighbour: at most 250 bp per 1,000,000 bp block (0.025%), no shared locus,
and every locus labelled from the VCF rather than from its neighbours.

Validation now holds 2,410 SNPs against 251 previously, and test 5,704 against
549.

## 3. Statistical baselines — the headroom check

Thresholds selected on the validation blocks, frozen, then applied to test.

| arm | thr | P | R | F1 | 95% CI | TP | FP | FN |
|---|---:|---:|---:|---:|---|---:|---:|---:|
| frequency | 0.255 | 0.9201 | 0.9386 | 0.9293 | [0.9243, 0.9341] | 5354 | 465 | 350 |
| binomial v1 | 7.0 | 0.9747 | 0.9385 | 0.9562 | [0.9523, 0.9600] | 5353 | 139 | 351 |
| **Poisson-binomial** | **10.5** | **0.9899** | **0.9623** | **0.9759** | **[0.9730, 0.9787]** | **5489** | **56** | **215** |

**271 test errors**, against 10 in the old 15× set and 2 at 30× — a 27× and 135×
larger budget. Validation carries 103 errors. The caller hierarchy is also
restored: PB beats binomial v1 by +0.0197 F1, where at 30× the two were
indistinguishable. The regime is not saturated, so the experiment proceeded.

## 4. The information gap

Error analysis of PB on test:

**False negatives (215).** Depth median 7. Alternate-read histogram: 0 reads →
78 loci, 1 → 91, 2 → 31, 3 → 15, ≥4 → 0. **169 of 215 (79%) carry ≤1 alternate
read and 78 carry none at all.** For those the variant is absent from the reads
and no architecture reading this pileup can recover them. Contestable: the 46
loci with 2–3 alternate reads.

**False positives (56).** Depth median 11, VAF median 0.40, median 4 alternate
reads, 12 with ≥6. Loci with real alternate support that are not variants — the
artefact class.

Complementary signal, measured at **fixed alternate-read count** (which removes
PB's own variable, so any remaining separation is something PB cannot express),
on **validation only**:

| feature | k=2 | k=3 | k=4 | k=5 |
|---|---:|---:|---:|---:|
| **alt_strand_bias** (inverted) | **0.674** | **0.660** | **0.862** | **0.739** |
| alt_mean_phred | 0.939 | 0.984 | 1.000 | 1.000 |
| alt_mean_mapq | 0.502 | 0.500 | 0.500 | 0.500 |
| alt_near_end_fraction | 0.505 | 0.511 | 0.517 | 0.515 |
| alt_mean_read_position | 0.495 | 0.532 | 0.507 | 0.453 |

`alt_mean_phred` scores highly but is **not complementary**: PB already consumes
per-read base quality, and at fixed k its LLR is a monotone function of exactly
those numbers. MAPQ and read position carry nothing (AUC ≈ 0.50) — unsurprising,
since MAPQ is ~70 everywhere in this data.

**`alt_strand_bias` is the one genuinely complementary signal.** It argues
*against* a variant, so it can only remove false positives — and for the first
time there were 56 to remove rather than 2. Predicted ceiling for a perfect
artefact detector: F1 0.9759 → 0.9808, ΔF1 +0.0049, against a CI half-width of
0.0028. Detectable in principle.

## 5. Model, gate and training

`final_logit = PB_logit + residual(44-channel evidence)`, residual
zero-initialized. Architecture, representation, prior, optimizer and loss all
unchanged from the previous experiments — the only new element is the split.

### Initialization gate — PASSED, all seeds, over all 1,907,840 validation loci

| property | measured |
|---|---|
| max abs deviation | 1.068e-05 (float32 ULP here: 1.526e-05) |
| Pearson | 0.999999999999997 |
| Spearman | 0.999999974101353 |
| **max rank-swap gap** | **9.295e-07** vs ULP 1.526e-05 (16× margin) |
| call disagreements | **0** |
| confusion vs PB | **2,338 / 31 / 72 — identical** |

One tolerance was corrected **before any test scoring**: `SPEARMAN_TOLERANCE`
1e-8 → 1e-6. A fixed absolute bound on that statistic is not well posed, since
its deviation grows with locus count purely through float32 tie resolution
(3.4e-08 at 190k loci, 2.0e-10 at 4k) — a property of the sample size, not the
model. The binding criterion is the exact, sample-size-independent rank-swap-gap
check, which passes with a 16× margin.

### Configuration

epochs 8 · batch 16 windows · seq_len 64 · lr 3e-4 · AdamW · focal γ=3.0 ·
class weights α=0.75 cap 20 · mutation-presence sampler · residual λ=0 ·
selection = max validation average precision · no early stopping.

Seeds 20260811, 424242, and 13371337 (post-hoc replication check).

### Validation selection was a real operation this time

Unlike 30×, where F1, AP and AUC were all exactly 1.000000 for all nine
candidates, here AP spans 0.9828–0.9929 across epochs. Selection discriminates.

| arm | selected epoch | val F1 | val AP | TP | FP | FN |
|---|---:|---:|---:|---:|---:|---:|
| PB (epoch 0) | — | 0.9784 | 0.9916 | 2338 | 31 | 72 |
| seed A | 3 | 0.9807 | 0.9929 | 2339 | 21 | 71 |
| seed B | 5 | 0.9813 | 0.9929 | 2341 | 20 | 69 |
| seed C | 2 | 0.9801 | 0.9925 | 2342 | 27 | 68 |

Selection stayed on the pre-registered average-precision rule. That mattered:
seed A's epoch 8 had the best validation F1 (0.9832) but not the best AP, and
switching to F1 after seeing the numbers would have been choosing the rule to
suit the result.

## 6. Frozen test results

| arm | thr | P | R | F1 | 95% CI | TP | FP | FN |
|---|---:|---:|---:|---:|---|---:|---:|---:|
| Poisson-binomial | 10.5 | 0.9899 | 0.9623 | 0.9759 | [0.9730, 0.9787] | 5489 | 56 | 215 |
| **PB-residual seed A** | 7.0 | 0.9921 | 0.9627 | **0.9771** | [0.9744, 0.9799] | 5491 | 44 | 213 |
| **PB-residual seed B** | 9.5 | 0.9933 | 0.9628 | **0.9778** | [0.9750, 0.9806] | 5492 | 37 | 212 |
| **PB-residual seed C** | 7.5 | 0.9906 | 0.9623 | **0.9763** | [0.9733, 0.9791] | 5489 | 52 | 215 |
| PB + constant shift | 10.5 | 0.9899 | 0.9623 | 0.9759 | [0.9730, 0.9787] | 5489 | 56 | 215 |
| *control:* shuffled prior | 9.5 | 0.5518 | 0.9862 | 0.7077 | [0.6997, 0.7160] | 5625 | 4568 | 79 |

### Paired bootstrap vs Poisson-binomial (10,000 resamples)

| arm | ΔF1 | 95% CI | ΔTP | ΔFP | ΔFN |
|---|---:|---|---:|---:|---:|
| seed A | +0.0012 | **[+0.0001, +0.0024]** | +2 | −12 | −2 |
| seed B | +0.0019 | **[+0.0009, +0.0030]** | +3 | −19 | −3 |
| seed C | +0.0003 | **[−0.0007, +0.0015]** | 0 | −4 | 0 |
| binomial v1 | −0.0197 | [−0.0226, −0.0170] | −136 | +83 | +136 |
| frequency | −0.0466 | [−0.0509, −0.0425] | −135 | +409 | +135 |
| shuffled prior | −0.2682 | [−0.2765, −0.2600] | +136 | +4512 | −136 |
| constant shift | +0.0000 | [+0.0000, +0.0000] | 0 | 0 | 0 |

**Two of three seeds exclude zero; seed C does not.** This is the finding that
sets the verdict.

## 7. Error overlap

| arm | PB FPs removed | new FPs | PB FNs recovered | FNs lost | disagreements |
|---|---:|---:|---:|---:|---:|
| seed A | 15 | 3 | 14 | 12 | 44 |
| seed B | **19** | **0** | 12 | 9 | 40 |
| seed C | 8 | 4 | 13 | 13 | 38 |

Precision rises in every seed (0.9899 → 0.9921 / 0.9933 / 0.9906) while recall
never falls. **Recall was not bought with false positives** — the failure mode
the protocol warns about does not occur. Seed B is the cleanest case: 19 of PB's
56 false positives removed, none created.

The false-negative ledger is close to a wash in every seed (14/12, 12/9, 13/13),
which is expected: 79% of PB's false negatives carry ≤1 alternate read and are
unreachable in principle.

## 8. Residual diagnostics

| | seed A | seed B | seed C |
|---|---:|---:|---:|
| mean | −0.747 | +0.833 | +0.527 |
| std | 1.210 | 1.198 | 0.916 |
| mean on SNP | −11.131 | −4.596 | — |
| **dispersion / location** | **1.62** | **1.44** | **1.74** |
| standalone ROC-AUC | 0.0195 | 0.0576 | — |
| corr. with **alt_strand_bias** | **+0.335** | **+0.174** | **+0.344** |
| corr. with PB LLR | −0.525 | −0.336 | −0.510 |
| corr. with VAF | −0.537 | −0.311 | — |
| corr. with alt_mean_phred | −0.461 | −0.331 | — |
| corr. with depth | +0.213 | +0.199 | — |

Against the four hypotheses the protocol requires testing:

1. **Constant shift — rejected.** Dispersion/location is 1.44–1.74; a constant
   would be ≈0. The 30× experiment's residual, by contrast, was flat.
2. **Threshold movement — rejected.** See §9.
3. **Shrinkage of the prior — partially present.** Correlation −0.34 to −0.53
   with the PB LLR means the residual does pull confident scores toward the
   boundary. This is stated rather than glossed: the residual is *not* purely
   complementary information.
4. **Genuinely complementary evidence — present.** Correlation +0.17 to +0.34
   with `alt_strand_bias` in every seed, the one feature identified as
   complementary on validation *before* training, and the mechanism (false
   positive removal) matches the prediction exactly.

## 9. Constant-shift null

| adjustment | thr | F1 | TP | FP | FN |
|---|---:|---:|---:|---:|---:|
| PB, unadjusted | 10.500 | 0.9759 | 5489 | 56 | 215 |
| PB + best free constant (validation-fitted) | 10.500 | 0.9759 | 5489 | 56 | 215 |
| PB at seed A's own mean residual (−0.747) | 11.247 | 0.9755 | 5474 | 45 | 230 |
| PB at seed B's own mean residual (+0.833) | 9.667 | 0.9760 | 5503 | 70 | 201 |
| PB at seed C's own mean residual (+0.527) | 9.973 | 0.9757 | — | — | — |
| **neural, seeds A / B / C** | — | **0.9771 / 0.9778 / 0.9763** | | | |

**Every neural arm beats its own constant-shift equivalent**, and the best free
additive constant fitted on validation is exactly 0.0 — PB's threshold was
already validation-optimal. The improvement is therefore not a moved decision
boundary. This is the sharpest distinction from the 30× result, where the
residual *was* reproducible only as a shift.

## 10. Depth-stratified

| depth | loci | SNP | PB F1 | PB TP/FP/FN | seed B F1 | seed B TP/FP/FN |
|---|---:|---:|---:|---|---:|---|
| 1–4 | 39,923 | 71 | 0.8033 | 49/2/22 | 0.8099 | 49/1/22 |
| 5–9 | 540,188 | 946 | 0.9470 | 866/17/80 | 0.9493 | 862/**8**/84 |
| 10–14 | 1,435,554 | 2,284 | 0.9840 | 2239/28/45 | 0.9870 | 2246/**21**/38 |
| 15–19 | 1,180,977 | 1,683 | 0.9952 | 1676/9/7 | 0.9958 | 1676/7/7 |
| 20–29 | 519,522 | 647 | 1.0000 | 647/0/0 | 1.0000 | 647/0/0 |
| 30–39 | 10,287 | 12 | 1.0000 | 12/0/0 | 1.0000 | 12/0/0 |

The gain is concentrated at **5–14×**, where the false positives live (FP 17→8
and 28→21). Above 20× both callers are perfect, which is the 30× result
reappearing as a stratum rather than as a whole benchmark.

## 11. Controls

**Shuffled prior** (PB LLR permuted across loci, marginal preserved): test F1
collapses to 0.7077 with 4,568 false positives, ΔF1 −0.2682
[−0.2765, −0.2600]. Validation F1 never exceeds 0.70 across 8 epochs against
0.98 with the true prior. The model depends on the genuine locus-specific
statistical evidence, not on merely receiving another numerical channel. The
residual head meanwhile learns large corrections (means −14 to −28 in the
earlier 30× run under the same control), so it has capacity and gradient to
spare.

## 12. Runtime

7,611,136 loci. CPU is a 12-thread desktop part; GPU an RTX 3070.

| stage | seconds | loci/s | device |
|---|---:|---:|---|
| pileup count extraction | 595.7 | 12,777 | CPU |
| read-level extraction (48 slots) | 1,162.7 | 6,546 | CPU |
| 44-channel feature build | 103.3 | 73,680 | CPU |
| **Poisson-binomial scoring** | **890.5** | **8,547** | CPU |
| binomial v1 scoring | 6.77 | 1,124,276 | CPU |
| **Mamba inference** | **6.29** | **1,209,780** | GPU |

The practical conclusion is unchanged from 30× and is worth restating because it
inverts the usual expectation: **the statistical prior is the expensive
component and the neural network is almost free.** The Mamba costs 6.3 s on top
of an 890 s prior — a 0.7% surcharge — and runs 142× faster per locus than the
PB convolution it corrects. If the residual's improvement is worth having, its
compute cost is not the obstacle.

## 13. Incidental engineering result

`GiabAlignmentProvider._is_confident` scanned every BED interval for every
window — on this region ~125,000 windows × 2,934 intervals ≈ 367 M Python
comparisons per extraction pass, paid twice. Since `load_bed_regions` already
returns sorted, merged, non-overlapping intervals, a `bisect_right` finds the
only candidate interval in O(log n).

**One extraction stage went from ~35 minutes to 51 seconds.** Output verified
identical, not assumed: the 30× test region was re-extracted with the new code
and its `counts` and `labels` arrays compared byte-for-byte against the
committed cache, plus 10 equivalence tests against a brute-force reference
covering exact fills, both off-by-one edges, inter-interval gaps, before/after
all intervals, the empty list, randomised interval sets, and 4,000 random
windows against the real GIAB BED.

This cost had been paid by every extraction in the project's history.

## 14. Comparison across all three regimes

| | 15× (0.47 Mb) | 30× (0.47 Mb) | **15× (3.81 Mb)** |
|---|---:|---:|---:|
| test SNPs | 549 | 549 | **5,704** |
| validation SNPs | 251 | 251 | **2,410** |
| PB F1 | 0.9908 | 0.9982 | 0.9759 |
| PB TP/FP/FN | 539/0/10 | 549/2/0 | 5489/56/215 |
| PB error budget | 10 | 2 | **271** |
| validation errors | 2 | **0** | **103** |
| validation AP range | — | 1.000000 exactly | 0.9828–0.9929 |
| neural ΔF1 | +0.0000 | +0.0000 | **+0.0012 / +0.0019 / +0.0003** |
| call disagreements | 0 | 0 | **38–44** |
| residual dispersion/location | — | ~0 (constant) | **1.44–1.74** |
| verdict | NULL (underpowered) | SATURATED | **WEAK POSITIVE** |

The progression is coherent. At 0.47 Mb the neural model made *zero* different
calls from PB at either coverage — there was nothing to learn and no way to
measure it. At 3.81 Mb it makes 38–44 different calls, most of them correct, and
its residual is locus-specific rather than constant. **The variable that
unlocked a measurable effect was evaluation size, not coverage, not variant
class, and not model capacity.**

## 15. Verdict

**WEAK POSITIVE.**

Against the stated criteria:

* *Reproducibly beats the baseline* — in direction yes, all three seeds ΔF1 ≥ 0
  and all three remove false positives; in magnitude, two of three.
* *Survives paired bootstrap* — seeds A and B yes, seed C no.
* *Replicates across seeds* — **weakly**. This is the binding limitation.
* *Cannot be explained by threshold movement* — **yes, decisively.** Every seed
  beats its own constant-shift null, and the best free constant is 0.0.
* *Residual diagnostics show genuine locus-specific information* — **yes.**
  Dispersion/location 1.44–1.74 and correlation +0.17 to +0.34 with
  `alt_strand_bias` in every seed.

That is the WEAK POSITIVE definition almost word for word: a small improvement
appears, seed replication is weak, and the residual diagnostics show genuine
locus-specific information.

### What this establishes

There **is** a regime where the learned representation provides information the
Poisson-binomial caller does not model. It is the artefact-detection regime:
loci with real alternate-read support that are not variants, distinguished by
strand asymmetry that a per-read independence model structurally cannot express.
The effect is worth 4–19 false positives out of 56, at 5–14× depth, for 0.7% of
the pipeline's compute.

### What it does not establish

That the effect is large enough to matter operationally, or that it is stable
enough to deploy. ΔF1 of +0.0003 to +0.0019 on a 271-error budget is a handful
of loci, and one seed in three cannot distinguish it from noise.

### Where to go next, in order

1. **More evaluation, not more model.** The effect is at the edge of
   detectability with 5,704 test SNPs. The obvious next step is 30–50 Mb —
   pooled chromosomes or a whole high-confidence genome — which would settle
   the magnitude question outright. Nothing about the architecture needs to
   change to run it, and the `_is_confident` fix makes it affordable.
2. **More seeds before any claim is strengthened.** Three seeds gave verdicts of
   "significant, significant, not significant". Five to ten would establish
   whether the mean effect is +0.0011 with seed variance, or whether A and B
   were fortunate.
3. **Target the mechanism directly.** The residual is currently a mixture of
   genuine strand-bias detection and prior shrinkage (correlation −0.34 to −0.53
   with the LLR). A residual constrained to be a function of artefact features
   only — with the count/quality channels ablated — would test whether the
   complementary component alone carries the gain, and would be a cleaner
   scientific object than the current mixture.
4. **Do not pursue indels, low-mappability, or higher coverage on GIAB
   high-confidence data.** All three were measured and are saturated or absent,
   for the same structural reason: the high-confidence BED removes difficulty
   by construction.
