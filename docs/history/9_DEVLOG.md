# Neural residual above the Poisson-binomial caller — 15× HG002, chr21

**Date:** 2026-08-12
**Verdict: NULL.** All three architecture variants produce **byte-identical test
calls to the Poisson-binomial caller**: TP 539, FP 0, FN 10, F1 0.9908,
ΔF1 = 0.0000 with 95 % CI [0.0000, 0.0000]. Two of the three selected the
*untrained* checkpoint on validation, i.e. training never beat its own starting
point.

The result is not "the architecture is too small". Error analysis identifies a
structural reason, stated in §7: **93.7 % of the remaining error budget is
recall, half of it at loci carrying ≤1 alternate read, while the only
complementary signal the network can add is a precision signal — and precision
is already 1.000.**

---

## 1. Bottleneck identified before building anything

The 14-channel `raw_pileup` representation was condemned as a *representation*
by two measurements, not by intuition:

1. **It cannot see the variable the previous experiment proved decisive.**
   `History/8_DEVLOG.md` established that ALT-read base quality is what
   separates a real low-depth heterozygote from noise: estimators that saw it
   gained ΔF1 +0.015…+0.018, the reference-only estimator that could not was
   null and added 9 FPs. `raw_pileup` supplies exactly one quality channel,
   `mean_base_quality`, **pooled over reference and alternate reads together**
   (`raw_pileup.py:162`). At k=2 of n=10 the decisive signal arrives diluted
   5:1 and confounded with the background.
2. **Integer `n` and `k` are not recoverable from it.** The fractions divide by
   `called + gap_reads` while `log_depth` counts reads admitted before the
   base-quality filter. Measured: the two disagree at **427,737 of 455,488**
   test loci (93.9 %).

That is a sufficient explanation for the earlier residual model collapsing to a
constant shift — it had no locus-specific information the prior lacked, so a
constant was the correct thing to learn.

## 2. The three variants

| | Variant | Changed vs the previous residual experiment | Hypothesis |
|---|---|---|---|
| **1** | Representation | 44-channel evidence + PB prior; backbone, loss, split, optimizer, sampler all unchanged | The bottleneck is representation |
| **2** | Gated residual | Variant 1 + `final = prior + gate(x)·residual(x)` | Lets the model say *where* the prior is reliable instead of shifting it everywhere |
| **3** | Read-level set encoder | *not run* — gated on Variant 1 showing life | — |

Variant 1 first, because it is the minimal change: the only moved variable.

**The 44 channels** (`locus_evidence.py`), each group with an explicit hypothesis:
`count` (12) — integer n, k, ref/other counts, depth, gap, VAF;
`alt_quality` (6) and `ref_quality` (7) — mean/min/max Phred, mean error, Q30
fraction **on each side of the ref/alt split**, plus the alt−ref contrast;
`strand_position` (14) — strand balance, read position, end distance, MAPQ,
again split by ref/alt; `reference` (5) — the base one-hot.

`strand_position` is the group that could in principle beat PB: the
Poisson-binomial model treats reads as independent Bernoulli draws indexed only
by quality, so it *structurally cannot* express "these 3 high-quality ALT reads
are all on one strand".

**Capacity is controlled.** Only the input projection widens: +3,900 parameters
(415,964 vs 412,064), asserted by a unit test. A gain could not be capacity.

**Initialization is exact.** The residual head is zeroed, so the untrained model
reproduces the PB decision score to 5.5e-6 (float32 rounding). Epoch 0 *is* the
PB baseline, measured on this experiment's own split with its own threshold
grid.

## 3. A documented protocol change: selection metric

The validation split holds **251 SNPs**, and the PB prior — where the model
*starts* — already scores 249 TP / 0 FP / 2 FN on it. Validation F1 therefore
has a total dynamic range of **two loci**. Selecting checkpoints or
architectures on it would be selecting on noise.

Selection was keyed on validation **average precision** instead, which ranks all
188k validation loci and is threshold-free. F1 is still logged every epoch. The
change follows from the split's size and used no test label.

## 4. Pre-flight, before any GPU time

A logistic residual probe (`score = PB_llr + w·features`, ~45 parameters, fitted
on the fit split, evaluated on validation) moved validation from 249 TP/0 FP/2 FN
to 250/0/1 — **one locus**. Reported at the time as too thin to be evidence.
It was correct.

## 5. Training result — validation

| Run | Selected epoch | Val AP | vs PB AP 0.99931 | Val F1 |
|---|---|---|---|---|
| V1 representation, seed 20260811 | 3 | 0.99939 | **+0.00008** | 0.9960 (= PB) |
| V1 representation, seed 424242 | **0 (untrained)** | 0.99931 | 0.00000 | 0.9960 (= PB) |
| V2 gated, seed 20260811 | **0 (untrained)** | 0.99931 | 0.00000 | 0.9960 (= PB) |

Validation F1 was 0.9960 — exactly the PB value — at **every epoch of every
run** bar one (0.9940). Per-epoch AP oscillated *below* the untrained baseline
in 21 of 24 trained epochs. The seed-A epoch-3 blip of +0.00008 did not
replicate on seed B.

The residual did become locus-specific, unlike the previous experiment: its
SNP-gap standard deviation reached 1.4–2.8 against a mean of 0.6–1.7 (the old
model was std ≈ 0). The gate also varied across loci (0.83 ± 0.07). **Neither
converted into accuracy.**

## 6. Test result — evaluated once, thresholds frozen

| Method | thr | P | R | F1 | TP | FP | FN |
|---|---:|---:|---:|---:|---:|---:|---:|
| Binomial v1 | 10.5 | 0.9981 | 0.9490 | 0.9729 | 521 | 1 | 28 |
| Mamba 14ch *(previous)* | — | 0.9693 | 0.9217 | 0.9449 | 506 | 16 | 43 |
| Residual binomial Mamba *(previous)* | 10.5 | 0.9981 | 0.9490 | 0.9729 | 521 | 1 | 28 |
| **Poisson-binomial** | 12.5 | 1.0000 | 0.9818 | 0.9908 | 539 | 0 | 10 |
| V1 representation, seed A | 10.5 | 1.0000 | 0.9818 | **0.9908** | 539 | 0 | 10 |
| V1 representation, seed B | 12.5 | 1.0000 | 0.9818 | **0.9908** | 539 | 0 | 10 |
| V2 gated, seed A | 12.5 | 1.0000 | 0.9818 | **0.9908** | 539 | 0 | 10 |

Paired bootstrap vs PB: **ΔF1 = 0.0000, CI [0.0000, 0.0000], ΔTP = ΔFP = ΔFN = 0**
for all three. Depth-stratified: identical to PB in every stratum including
5–9× and 10–14×.

Seed A is the interesting case. Its residual is genuinely non-zero on test
(mean +1.32, std 1.42, 1st–99th percentile −1.80…+3.98) — it moves scores, but
**never across the decision boundary at a single locus out of 455,371**.

*Oracle diagnostic only, not a result:* at test-optimal thresholds seed A scores
0.9927 (541/0/8) against PB's 0.9918 (541/1/8) — one false positive, at an
operating point chosen with test labels. One locus is not evidence.

## 7. Why — the error budget is the answer

Measured on the fit split (1,852 SNPs, no test data involved):

```
PB @ 12.5:  TP=1793  FP=4  FN=59     precision 0.99777   recall 0.96814
error budget:  false negatives 59/63 = 93.7%
               false positives  4/63 =  6.3%
```

PB's false negatives by alternate-read count:

| k | loci | share | median n |
|---:|---:|---:|---:|
| 0 | 8 | 13.6 % | 5 |
| 1 | 22 | 37.3 % | 6 |
| 2 | 27 | 45.8 % | 10 |
| 3 | 2 | 3.4 % | 11 |

**Half of the remaining misses (30/59) carry ≤1 alternate read; 8 carry none at
all.** For those the variant is absent from the reads: no architecture reading
this pileup can recover them. That is a sampling limit at 15×, not a modelling
limit.

Which features carry information *complementary* to PB? Testing at **fixed k**
(which removes PB's own variable, so any remaining separation is something PB
cannot express), on contested loci:

| feature | AUC at k=2 | AUC at k=3 |
|---|---:|---:|
| `alt_mean_phred` | 0.986 | 0.980 |
| `alt_max_phred` | 0.974 | 0.981 |
| **`alt_strand_bias`** | **0.276 (0.724 inverted)** | **0.251 (0.749 inverted)** |
| `alt_near_end_fraction` | 0.506 | 0.512 |
| `alt_mean_mapq` | 0.509 | — |

The quality features score highly but are *not* complementary — PB already
consumes per-read quality, and at fixed k its LLR is a monotone function of
exactly those numbers. Read position and mapping quality carry **nothing**
(AUC ≈ 0.51).

**`alt_strand_bias` is the one genuinely complementary signal**, worth
AUC ≈ 0.72–0.75 at fixed k, and PB structurally cannot represent it. But strand
bias, near-end clustering and low MAPQ are all *artefact* signals: they argue
**against** a variant, so they can only remove false positives.

> Maximum conceivable gain from a perfect artefact detector: **4 loci**.
> Recall-side errors it cannot touch: **59 loci**.

The one thing the neural model can add is a precision tool, and PB's precision
is already 1.0000 on test and 0.9978 on 1.69 M fit loci. There is no headroom
for it to recover. The null result is what this error budget predicts.

## 8. Honest limitations

* **The experiment is underpowered.** Validation resolution is 2 loci; test
  headroom above PB is 10 FN of which ~5 are unreachable. A real improvement of
  1–2 loci could not be distinguished from noise here. "No effect detected" is
  the correct reading — not "no effect exists".
* Variant 3 (read-level set encoder) was not run: it was gated on Variant 1
  showing life, and it did not.
* Only 8 epochs, one learning rate, one loss. The loss was *not* re-tuned,
  deliberately — changing representation and objective together would have made
  the outcome unattributable. The systematic positive residual mean (+2 to +4 in
  early epochs) is consistent with the focal objective preferring a global shift,
  and remains an untested confound.
* The 44-channel representation reorders counts into ref/alt slots, which
  `raw_pileup` deliberately refused to do. This makes the model's task easier by
  construction and is a genuine change of philosophy, documented in
  `locus_evidence.py`.

## 9. Recommendation — stop optimizing at 15×

**Do not enlarge the network, and do not run Variant 3 on this data.** Neither
addresses the binding constraint, which is that the reads do not contain the
missing variants.

Priorities, in order:

1. **Change the regime, not the model.** The complementary signal that exists —
   strand/artefact structure — only pays where there are false positives to
   remove. Move to a setting with FP pressure: higher coverage (30–50×, where
   recall saturates and precision becomes the binding constraint), low-mappability
   regions, or the indel classes, which the whole project has so far excluded
   from the primary metric.
2. **Fix the evaluation before the model.** Any future architecture claim needs
   a validation split with enough positives to resolve it — a larger region, or
   pooled chromosomes. At 251 validation SNPs, architecture selection is not
   measurable, and that is now the first blocker.
3. **If the neural direction continues, target recall at k=2**, the only
   contestable stratum (27 of 59 fit-split misses). That is a question about
   the *prior* over variant sites — local sequence context, homopolymers,
   known-variant density — not about read-level evidence, which PB has already
   exhausted.

## 10. Repository changes

**Added:** `locus_evidence.py`, `model_pb_residual.py`, `train_pb_residual.py`,
`cache_locus_evidence.py`, `sanity_pb_residual.py`, `evaluate_pb_residual.py`,
`test_locus_evidence.py` (22 tests), `checkpoints_pb_residual/`,
`cache/{train,test}_15x_evidence.npz`,
`results/pb_residual_{sanity,15x_results,v1_seedA_train,v1_seedB_train,v2_gate_seedA_train}.json`,
this devlog.

**Not changed:** no previous checkpoint, dataset, cache, result file or
scientific module was modified. `RawPileupMamba` is imported as the backbone and
`BinomialVariantCaller` is delegated to for ALT selection — neither is edited,
so every baseline in the comparison table is byte-identical to the one that
produced its published number. Full suite: 154 passed.
