# 30× coverage regime — PB-residual Mamba vs the Poisson-binomial caller

## SATURATED REGIME

At approximately 30× the Poisson-binomial caller reaches **recall 1.0000** with
**two false positives in 455,371 loci**. Doubling coverage converted every one
of the 15× false negatives into a true positive — the sampling limit the 15×
error budget predicted — and in doing so removed the error budget that any
comparison would need in order to resolve a difference.

The frozen PB-residual model is, at both seeds, **exactly the Poisson-binomial
caller**: validation could not distinguish any checkpoint from any other, so
the honest selection rule returned epoch 0, whose residual is identically zero
by the initialization gate. ΔF1 = +0.0000, 95% CI [+0.0000, +0.0000], **zero
call disagreements at every one of 455,371 loci**.

This is not a failure of the experiment. It is the benchmark reporting that it
has no headroom left: the maximum improvement any model could demonstrate over
a threshold-tuned PB at 30× is **one locus**.

The question asked — *"can RawPileupMamba learn stable locus-specific
information that improves variant calling beyond the Poisson-binomial caller at
30×?"* — cannot be answered by this benchmark, because the benchmark cannot
resolve the answer. It is not that the answer is "no"; it is that 30× on this
region destroys the measurement rather than enabling it.

---

## 1. Dataset

Identical to the 15× experiment in every respect except sequencing depth.

| | |
|---|---|
| Sample | HG002 NA24385, NIST Illumina 2×250bp, novoalign, GRCh38 |
| Truth | GIAB NIST v4.2.1 benchmark VCF |
| Confidence regions | v4.2.1 `benchmark_noinconsistent.bed` |
| Reference | `data/reference/chr21_full.fa` (Ensembl 110) |
| Test region | chr21:30,000,000–30,469,999 → 455,488 loci |
| Train region | chr21:31,000,000–33,000,000 → 1,876,352 loci |
| Split | last 10% of training-region window batches = 187,776 loci |

Construction, the single changed variable:

```bash
samtools view -b --subsample 0.406 --subsample-seed 20260811   # 30x
samtools view -b --subsample 0.203 --subsample-seed 20260811   # 15x, unchanged
```

The fraction is exactly doubled and the seed held. `samtools` keeps reads by
hash threshold, so the 15× read set is a strict subset of the 30× set —
verified directly: all 2,996 15× reads in a test window are present among the
5,980 30× reads, zero exceptions. The two regimes are the same underlying data.

**Truth labels are byte-identical between regimes** (`np.array_equal` on both
regions): 455,488 test loci / 549 SNP / 40 insertion / 77 deletion;
1,876,352 train loci / 2,103 SNP / 234 insertion / 761 deletion. No change to
the truth set, region, candidate-ALT logic, split methodology, evaluation
definition, model, objective, optimizer or loss.

### The one forced change, and why keeping the old value would have been worse

The read tensor is capped at `MAX_READS`. At 15× the cap **never bound**: max
depth 38 against a cap of 48. At 30× it binds on **4,234 test loci (0.93%)**,
max depth 62 (79 in train).

This is not cosmetic. The Poisson-binomial prior takes its per-read error
probabilities from the read tensor but its alternate count `k` from the
full-depth count matrix. At a truncated locus the ratio is inflated — a
heterozygous SNP at depth 62 would present as VAF ≈ 0.65 rather than 0.5 —
which manufactures exactly the high-depth false positives this experiment
exists to measure.

The codebase's own integrity assertion settles it. `load_region` rebuilds the
count matrix from the read tensor and demands exact equality. Measured:

| cap | loci violating the read/count integrity gate |
|---|---|
| 48 | **4,234** |
| 96 | **0** |

`MAX_READS` was therefore raised 48 → 96 for the 30× dataset. The invariant
that actually held at 15× is *"the read tensor contains every read at the
locus"*; keeping 48 would have been the silent change of feature semantics,
not raising it. No scientific code changed — every module is shape-generic in
the read axis, so this is one command-line argument. `k > 48` never occurs, so
the cap-48 defect would have been a bias, not a crash.

The cap-48 arm remains exactly recoverable for free: the extractor keeps the 48
smallest read hashes in sorted order, so `reads[:, :48]` reproduces the cap-48
tensor bitwise (verified).

## 2. Measured depth distribution

Measured from the pipeline's own count matrix, not from the BAM filename or the
subsample parameter.

| | 15× | 30× |
|---|---:|---:|
| mean (test) | 15.03 | **29.97** |
| median (test) | 15 | **30** |
| min / max (test) | 1 / 38 | 4 / 62 |
| mean (train) | ~14.4 | 28.80 |

| depth | 15× loci | 30× loci |
|---|---:|---:|
| 1–4 | 4,139 | 4 |
| 5–9 | 54,997 | 1,172 |
| 10–14 | 160,672 | 6,504 |
| 15–19 | 153,494 | 28,446 |
| 20–29 | 79,442 | 183,411 |
| 30–39 | 2,744 | 189,811 |
| 40–49 | 0 | 43,050 |
| 50+ | 0 | 3,090 |

Confirmed ≈30×.

## 3. Statistical baselines at 30×

Thresholds selected on the 30× validation split; test labels not consulted.

| arm | thr | P | R | F1 | TP | FP | FN |
|---|---:|---:|---:|---:|---:|---:|---:|
| frequency (VAF) | 0.255 | 0.9982 | 0.9982 | 0.9982 | 548 | 1 | 1 |
| binomial v1 | 11.5 | 0.9964 | 1.0000 | 0.9982 | 549 | 2 | 0 |
| binomial mean-Q | 12.5 | 0.9964 | 1.0000 | 0.9982 | 549 | 2 | 0 |
| **Poisson-binomial** | **14.0** | **0.9964** | **1.0000** | **0.9982** | **549** | **2** | **0** |

**Recall is exactly 1.0000.** Every one of the 549 SNPs is recovered.

Two consequences, both fatal to the comparison this experiment set out to make:

**(a) The error budget is two false positives.** Their LLRs are 16.19 and
24.25; the lowest-scoring true SNP is 16.83. The first sits below every true
SNP and is removable by threshold alone — an oracle threshold of 16.5 yields
549/1/0, F1 0.99909. The second sits *inside* the true-SNP LLR range and no
threshold on the PB score can remove it. So the ceiling for any model over a
threshold-tuned PB is **one locus in 455,371**, and reaching it requires
information PB structurally cannot express. A one-locus difference cannot
survive a paired bootstrap.

**(b) PB's 15× advantage has evaporated.** At 15× PB beat binomial v1 by
+0.0179 F1 (18 fewer FN). At 30× all three callers score identically. The
quality-aware modelling that was the previous experiment's positive result is
no longer distinguishable once depth is sufficient.

## 4. Neural baseline — 14-channel RawPileupMamba at 30×

Established architecture and training methodology, not enlarged, not retuned.
Validation-frozen threshold 0.80 (epoch 1).

| regime | P | R | F1 | TP | FP | FN |
|---|---:|---:|---:|---:|---:|---:|
| 30× | 0.9874 | 0.9982 | **0.9928** | 548 | 7 | 1 |
| 30× PB, same test set | 0.9964 | 1.0000 | **0.9982** | 549 | 2 | 0 |

The 14-channel model **scales but does not catch up**: it is worse than the
statistical caller at 30× (7 FP vs 2, 1 FN vs 0). Its validation is *not*
saturated (F1 0.944–0.988, AP 0.69–0.999 across epochs), which is itself
informative — the saturation is a property of the PB prior's separation of this
data, not of the dataset alone. Runtime: 8 epochs, ~13 min including extraction.

## 5. Primary experiment — PB-residual Mamba

`final_logit = PB_logit + residual(44-channel evidence)`, residual
zero-initialized.

### Initialization sanity gate — PASSED, both seeds

Run over the whole validation split, both experiment seeds:

| property | measured | required |
|---|---|---|
| max abs deviation | 1.511e-05 | ≤ float32 scale (1 ULP = 3.052e-05) |
| Pearson | 0.999999999999998 | 1.0 |
| Spearman | 0.999999999798215 | 1.0 |
| call disagreements | **0** | 0 |
| confusion vs PB | **251 / 0 / 0, identical** | identical |

Exact 1.0 Spearman is not a float32-achievable criterion: values closer than
one ULP may have their ranks swapped by rounding. The exact form of the
requirement was checked instead — the largest true gap between any two values
whose ranks moved is **1.893e-06**, against a float32 ULP of **3.052e-05**.
Every reordering is between values float32 cannot distinguish, and no decision
changes. This tolerance was widened once, before any test number was read, and
the reason is recorded here rather than quietly applied.

### Training configuration

Inherited unchanged from the 15× experiment; nothing retuned for coverage.

| | |
|---|---|
| epochs | 8 |
| batch size | 16 windows |
| sequence length | 64 loci |
| learning rate | 3e-4 |
| optimizer | AdamW, default weight decay |
| loss | focal, γ = 3.0 |
| class weighting | calibrated, α = 0.75, cap 20.0 |
| sampler | mutation-presence weighted, with replacement |
| split | sequential 90/10 by window batch |
| seeds | 20260811, 424242 |
| residual λ | 0.0 |
| checkpoint selection | max validation average precision |
| early stopping | none; all 8 epochs run, best retained |

### Validation selection — the metric is also saturated

Section 9 anticipated that F1 might lack dynamic range and named average
precision as the fallback. At 30× **the fallback is degenerate too**:

```
epoch  thr    F1        AP        AUC       TP  FP  FN   residual mean ± std
  0   14.0  1.000000  1.000000  1.000000  251   0   0   +0.000 ± 0.000  (== PB)
  1   12.0  1.000000  1.000000  1.000000  251   0   0   +0.955 ± 0.558
  2   13.5  1.000000  1.000000  1.000000  251   0   0   +1.386 ± 0.462
  4   17.0  1.000000  1.000000  1.000000  251   0   0   +2.811 ± 0.829
  5   10.5  1.000000  1.000000  1.000000  251   0   0   −1.764 ± 0.935
  8    9.5  1.000000  1.000000  1.000000  251   0   0   −1.593 ± 1.309
```

F1, AP and ROC-AUC are all exactly 1.000000 for all nine candidates, at both
seeds, because PB separates the validation split perfectly (251/0/0). **No
threshold-free selection metric remains.** The selection rule therefore
returned epoch 0 — the untrained model — at both seeds.

No third metric was sought. Choosing a selection rule after observing that the
pre-registered ones tie would be selecting a rule to manufacture a result, and
is the specific failure mode this protocol forbids.

The residual also does not converge. Seed A's mean across epochs runs
+0.96, +1.39, +0.84, +2.81, −1.76, +2.70, −0.22, −1.59; seed B's runs
+1.62, +4.01, +1.12, −0.52, −1.82, +0.73, +0.05, +4.04. Both flip sign
repeatedly, dispersion is the same order as location, and the two trajectories
share no structure.

## 6. Frozen test results

Evaluated once, after freezing. Validation thresholds throughout.

| arm | thr | P | R | F1 | TP | FP | FN |
|---|---:|---:|---:|---:|---:|---:|---:|
| frequency | 0.255 | 0.9982 | 0.9982 | 0.9982 | 548 | 1 | 1 |
| binomial v1 | 11.5 | 0.9964 | 1.0000 | 0.9982 | 549 | 2 | 0 |
| **Poisson-binomial** | 14.0 | 0.9964 | 1.0000 | 0.9982 | 549 | 2 | 0 |
| **PB-residual seed A** | 14.0 | 0.9964 | 1.0000 | **0.9982** | 549 | 2 | 0 |
| **PB-residual seed B** | 14.0 | 0.9964 | 1.0000 | **0.9982** | 549 | 2 | 0 |
| PB + constant shift | 19.5 | 0.9982 | 0.9982 | 0.9982 | 548 | 1 | 1 |
| *diagnostic:* seed A last epoch | 9.5 | 0.9964 | 1.0000 | 0.9982 | 549 | 2 | 0 |
| *diagnostic:* seed B last epoch | 15.5 | 0.9964 | 1.0000 | 0.9982 | 549 | 2 | 0 |
| *control:* shuffled prior | 9.0 | 0.6362 | 1.0000 | 0.7776 | 549 | 314 | 0 |

95% CI on PB's F1 and on every tied arm: the paired bootstrap below is the
correct interval, since the arms are not independent.

### Paired bootstrap vs Poisson-binomial, 10,000 resamples

| arm | ΔF1 | 95% CI | ΔTP | ΔFP | ΔFN |
|---|---:|---|---:|---:|---:|
| **PB-residual seed A** | **+0.0000** | **[+0.0000, +0.0000]** | 0 | 0 | 0 |
| **PB-residual seed B** | **+0.0000** | **[+0.0000, +0.0000]** | 0 | 0 | 0 |
| seed A last epoch | +0.0000 | [+0.0000, +0.0000] | 0 | 0 | 0 |
| seed B last epoch | +0.0000 | [+0.0000, +0.0000] | 0 | 0 | 0 |
| binomial v1 | +0.0000 | [+0.0000, +0.0000] | 0 | 0 | 0 |
| frequency | −0.0000 | [−0.0037, +0.0036] | −1 | −1 | +1 |
| shuffled-prior control | −0.2206 | [−0.2450, −0.1971] | 0 | +312 | 0 |

The interval is degenerate for the primary arms because they make *identical
calls*, not because the bootstrap failed.

## 7. Depth-stratified results

| depth | loci | SNP | PB F1 | PB TP/FP/FN | seed A F1 | seed A TP/FP/FN |
|---|---:|---:|---:|---|---:|---|
| 1–4 | 4 | 0 | — | 0/0/0 | — | 0/0/0 |
| 5–9 | 1,172 | 2 | 1.0000 | 2/0/0 | 1.0000 | 2/0/0 |
| 10–14 | 6,501 | 12 | 1.0000 | 12/0/0 | 1.0000 | 12/0/0 |
| 15–19 | 28,432 | 54 | 1.0000 | 54/0/0 | 1.0000 | 54/0/0 |
| 20–29 | 183,358 | 218 | 0.9954 | 218/2/0 | 0.9954 | 218/2/0 |
| 30–39 | 189,775 | 208 | 1.0000 | 208/0/0 | 1.0000 | 208/0/0 |
| 40–49 | 43,039 | 54 | 1.0000 | 54/0/0 | 1.0000 | 54/0/0 |
| 50+ | 3,090 | 1 | 1.0000 | 1/0/0 | 1.0000 | 1/0/0 |

Identical in every stratum. Note the prediction being tested — that neural
gains would appear specifically at higher coverage — is refuted in a way that
was not anticipated: **PB is perfect in every stratum above 20×**, and its only
two errors fall in the 20–29× band. There is no high-coverage regime where
precision becomes the binding constraint, because at these depths precision
does not degrade at all.

## 8. Error overlap

| arm | PB FN | arm FN | shared FN | PB FN recovered | arm FP not made by PB | disagreements |
|---|---:|---:|---:|---:|---:|---:|
| PB-residual seed A | 0 | 0 | 0 | 0 | **0** | **0** |
| PB-residual seed B | 0 | 0 | 0 | 0 | **0** | **0** |
| seed A last epoch | 0 | 0 | 0 | 0 | 0 | **0** |
| seed B last epoch | 0 | 0 | 0 | 0 | 0 | **0** |
| shuffled control | 0 | 0 | 0 | 0 | 312 | 312 |

PB has **zero** false negatives, so the question "does Mamba recover PB's false
negatives?" has no content at 30× — there are none to recover. The failure mode
the protocol warns about (buying recall with false positives) cannot arise
either, since recall is already 1.0.

The striking entry is the *diagnostic* row: even the fully trained last-epoch
checkpoints, evaluated at their own validation thresholds, produce **identical
calls at all 455,371 loci**. The trained residual moves scores substantially
but never across the decision boundary at a single locus — the same phenomenon
the 15× experiment recorded, reproduced here at double the coverage.

## 9. Residual diagnostics

| | frozen seed A | seed A last epoch | seed B last epoch |
|---|---:|---:|---:|
| mean (all) | −4.9e-10 | −1.5869 | +4.0628 |
| std (all) | 0.0000 | 1.1097 | 0.8854 |
| mean on SNP | +0.0000 | −23.1283 | −9.3692 |
| mean on Normal | −0.0000 | −1.5608 | +4.0790 |
| standalone ROC-AUC | 0.5012 | 0.0187 | 0.0140 |
| corr. with depth | −0.0009 | +0.0100 | +0.0018 |
| corr. with VAF | −0.0008 | **−0.8148** | **−0.6686** |
| corr. with PB LLR | +0.0031 | **−0.6997** | **−0.5728** |
| corr. with alt mean phred | +0.0014 | −0.4472 | −0.4411 |
| corr. with alt strand bias | −0.0056 | +0.6309 | +0.4963 |
| corr. with alt mean MAPQ | +0.0004 | −0.5034 | −0.4993 |
| dispersion / location | 971 | 0.70 | 0.22 |

The frozen model's residual is exactly zero — category **(A) constant**, at the
value zero, as the initialization gate guarantees.

The trained residuals are the interesting object, and they are **category (D)
but in the wrong direction**. They correlate −0.70 and −0.57 with the PB LLR
and −0.81 and −0.67 with VAF: the residual is largely *subtracting a scaled
copy of the prior*. Its standalone ROC-AUC is 0.019 and 0.014 — far below 0.5,
meaning it ranks true SNPs *lower* than normal loci (mean −23.1 on SNPs against
−1.6 on normal loci for seed A). This is shrinkage toward the decision
boundary, which is what a focal objective produces when the prior is already
confident and correct, and it is not additional evidence. Its apparent
"information" is a negative function of the LLR the model was handed, not
something read from the reads.

## 10. Constant-shift null

| adjustment | thr | F1 | TP | FP | FN |
|---|---:|---:|---:|---:|---:|
| PB, unadjusted | 14.000 | 0.9982 | 549 | 2 | 0 |
| PB + best validation-fitted constant | 19.500 | 0.9982 | 548 | 1 | 1 |
| PB + seed A last epoch's own mean residual (−1.587) | 15.587 | 0.9982 | 549 | 2 | 0 |
| PB + seed B last epoch's own mean residual (+4.063) | 9.937 | 0.9973 | 549 | 3 | 0 |
| neural, both seeds | 14.000 | 0.9982 | 549 | 2 | 0 |

**constant-shift PB ≥ neural model.** Per the pre-stated rule, the neural model
has not demonstrated additional information. The free-parameter shift trades
the removable false positive for a false negative and lands on the same F1;
the neural arms land on the same F1 by making literally the same calls.

## 11. Shuffled-prior control

The most informative arm in the experiment.

```
        thr    F1        AP       TP   FP   FN   residual mean ± std
ep 0   68.5  0.004505  0.00141    1  192  250    0.00 ± 0.00
ep 3    6.5  0.667568  0.33620  247  242    4   −21.62 ± 3.07
ep 7    9.0  0.665742  0.33792  240  230   11   −27.29 ± 2.89
test    9.0  0.777620      —    549  314    0        —
```

Permuting the PB LLR across loci collapses validation F1 from 1.0000 to ~0.67
and test F1 to 0.7776 (ΔF1 −0.2206, CI [−0.2450, −0.1971]). The model depends
on the genuine locus-specific statistical evidence, not merely on receiving an
additional numerical feature.

The load-bearing observation is the **residual magnitude**: against a corrupted
prior the head learns means of −14 to −28 logits with std 2.7–5.4, sustained
across epochs. The residual head therefore has ample capacity, gradient and
dynamic range. **That is what makes the primary null substantive rather than an
artefact of a dead head**: with the true prior the identical machinery settles
at exactly zero, not because it cannot learn a correction but because there is
no correction left to make.

## 12. Runtime

Test region, 455,488 loci. CPU is a 12-thread desktop part; GPU is an
RTX 3070.

| stage | seconds | loci/s | device |
|---|---:|---:|---|
| pileup count extraction | 57.6 | 7,908 | CPU |
| read-level extraction (96 slots) | 147.2 | 3,094 | CPU |
| 44-channel feature build | 11.2 | 40,668 | CPU |
| **Poisson-binomial scoring** | **461.9** | **986** | CPU |
| binomial v1 scoring | 0.219 | 2,076,518 | CPU |
| **Mamba inference** | **0.499** | **912,693** | GPU |

The practical finding inverts the usual expectation: **the statistical prior is
the expensive component and the neural network is nearly free.** The exact
Poisson-binomial convolution costs 462 s; adding the entire Mamba on top of the
same evidence costs 11.7 s (features + inference), a 4% surcharge on a
610 s pipeline. The neural model is ~925× faster per locus than the prior it
is built on.

The genuinely cheap option is binomial v1 at 2.08 M loci/s — 2,100× faster than
PB — which at 30× **scores identically to PB** (549/2/0). For this depth and
region, the correct engineering choice is the simplest caller, not the most
sophisticated one.

## 13. Comparison with the 15× regime

Reported side by side, not pooled.

| | 15× | 30× |
|---|---:|---:|
| mean depth | 15.03 | 29.97 |
| PB precision | 1.0000 | 0.9964 |
| PB recall | 0.9818 | **1.0000** |
| PB F1 | 0.9908 | 0.9982 |
| PB TP / FP / FN | 539 / 0 / 10 | **549 / 2 / 0** |
| binomial v1 F1 | 0.9729 | 0.9982 |
| PB advantage over binomial v1 | +0.0179 | **+0.0000** |
| 14-channel Mamba F1 | — | 0.9928 |
| PB-residual F1 | 0.9908 (tie) | 0.9982 (tie) |
| ΔF1 vs PB | +0.0000 | +0.0000 |
| call disagreements vs PB | 0 | 0 |
| validation SNPs | 251 | 251 |
| validation PB errors | 2 FN | **0** |
| validation F1 range | 2 loci | **0 loci** |

Depth did exactly what the 15× analysis predicted for *recall*: the 10 residual
false negatives, of which roughly half carried ≤1 alternate read, were a
sampling limit and disappeared when the reads arrived. That prediction is
confirmed.

The accompanying prediction — that higher coverage would create false-positive
pressure and hence a regime where the artefact signals (strand bias, read
position) that PB cannot express would pay — is **refuted**. Precision did not
become the binding constraint. PB gained two false positives, one of them
threshold-removable, and remains perfect in every depth stratum above 20×.
Moving from 15× to 30× reduced the total error budget from 10 loci to 2, and
reduced the *contestable* budget to 1.

Both regimes return the same verdict for the same underlying reason, reached
from opposite directions: at 15× there was no precision headroom for the
model's only complementary signal to recover; at 30× there is no headroom of
any kind.

## 14. Test integrity

* Frozen before evaluation: checkpoints, model configuration, thresholds,
  feature transformation, statistical prior, seed-selection rule.
* The test set was evaluated **once** per frozen configuration. No checkpoint,
  threshold, feature or architecture was modified after any test number was
  read.
* Test labels were not used in feature construction, training, architecture
  selection, hyperparameter selection, threshold selection or early stopping.
  Every threshold traces to the 30× validation split; the neural selection
  metric was pre-registered in `train_pb_residual.py` during the 15× work.
* The diagnostic arms (last-epoch checkpoints) and the oracle thresholds are
  labelled as diagnostics throughout and are not substituted for the frozen
  result anywhere.
* One tolerance was widened during the run — the initialization gate's Spearman
  criterion — before any test number was read, for a stated numerical reason,
  and replaced with a stricter exact check (§5).

## 15. Unit tests and file integrity

**179 tests pass**, including 26 new ones (`test_analysis_30x.py` 19,
`test_30x_init_gate.py` 7). No pre-existing test was modified.

All **38** 15× artifacts — caches, checkpoints, result JSONs and BAMs —
verified byte-identical by `md5sum -c` after the experiment. No previous
dataset, checkpoint, cache or result file was written to. New artifacts are
isolated in `checkpoints_30x_pb_residual/`, `results/30x_pb_residual/`,
`data/giab_hg002_real_30x/`, `data/giab_hg002_30x_readlevel/` and
`cache/*_30x_*.npz`.

## 16. Verdict

**SATURATED REGIME.**

The Poisson-binomial caller reaches recall 1.0000 with 2 false positives in
455,371 loci at 30×. One of those two is removable by moving the threshold;
the other is not removable by any threshold on the PB score. The contestable
error budget is therefore **one locus**, which no paired bootstrap can resolve.

Within that budget:

* the PB-residual Mamba is, at both seeds, **exactly the Poisson-binomial
  caller** — validation could not separate any checkpoint from any other on F1,
  average precision *or* ROC-AUC, so selection returned the zero-residual epoch;
* the trained checkpoints, evaluated as diagnostics, make **identical calls at
  all 455,371 loci** despite carrying large residuals;
* the constant-shift null matches the neural model, so no additional
  information is demonstrated;
* the shuffled-prior control shows the machinery works and would learn a large
  correction if one were available.

This is a negative result about the *benchmark*, not about the architecture.
The experiment was designed to test whether depth creates headroom for neural
modelling; it establishes that on this region depth removes headroom instead.

### What this rules out, and what it does not

It does **not** show that read-level structure is uninformative for variant
calling. It shows that on 470 kb of high-confidence chr21 at 30×, a
well-calibrated statistical caller is already at the ceiling, so the question
cannot be posed here. The measurement instrument, not the hypothesis, is what
failed.

### Where the question is still answerable

The binding constraint is now unambiguous and it is the **evaluation set**, not
the model. 549 test SNPs and 251 validation SNPs with zero validation errors
cannot resolve anything. Before any further architecture work:

1. **Move to a regime with a real error budget.** The candidates PB actually
   struggles with are outside this experiment's frame: low-mappability and
   segmental-duplication regions, homopolymers and tandem repeats, and the
   indel classes — which the primary metric has excluded from the start, and
   where the 30× test set still holds 40 insertions and 77 deletions untouched.
2. **Enlarge the evaluation before enlarging the model.** Pooled chromosomes or
   a whole-genome high-confidence set, sized so validation carries enough
   errors to separate checkpoints at all. At present architecture selection is
   not a measurable operation, and no amount of GPU time changes that.
3. **Stop testing on saturated benchmarks.** Both 15× and 30× now return the
   same tie for different reasons. A third coverage point on this region would
   answer nothing.

The honest summary of the two experiments together: the Poisson-binomial caller
with per-read quality is an extremely strong baseline for germline SNP calling
on high-confidence Illumina data, and this project has not yet constructed a
setting where it can be beaten.
