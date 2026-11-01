# Neural Residual Above the Poisson-Binomial Caller — Powered 15× HG002 Benchmark

**Date:** 2026-08-12
**Verdict: WEAK POSITIVE**

This experiment answers a more difficult question than the original 15× benchmark:

> **Can a neural residual model extract additional information from the pileup that is not represented by the Poisson-binomial statistical caller?**

The answer is **possibly yes, but the evidence is not yet strong enough to claim a robust improvement**.

The key change was not a larger neural network.

It was a larger evaluation region.

---

## 1. Why the previous benchmark was insufficient

The original experiment used:

* HG002 NIST v4.2.1
* GRCh38
* chr21:30,000,000–30,469,999
* approximately 15× coverage
* 455,371 callable test loci
* 549 SNPs

The Poisson-binomial caller already achieved:

```text
TP = 539
FP = 0
FN = 10
F1 = 0.9908
```

This left almost no measurable error budget.

The validation split contained only 251 SNPs, making architecture selection extremely low-power.

Therefore a neural model could not demonstrate a meaningful improvement even if a small amount of complementary information existed.

The correct response was not to enlarge the network.

It was to enlarge the benchmark.

---

## 2. Powered benchmark

The new split increased the evaluation region while preserving the training scale.

| Role       |  Region |  SNPs | Change       |
| ---------- | ------: | ----: | ------------ |
| Train      | 1.89 Mb | 2,544 | Same scale   |
| Validation | 1.91 Mb | 2,410 | 9.6× larger  |
| Test       | 3.81 Mb | 5,704 | 10.4× larger |

The expected Poisson-binomial error budget increased from approximately 10 test errors to approximately 271.

This makes model selection and paired comparison measurable for the first time.

The test set was kept completely isolated from training and model selection.

---

## 3. Statistical baseline

At 15× coverage on the powered benchmark:

| Method                     | Precision | Recall |     F1 |    TP | FP |  FN |
| -------------------------- | --------: | -----: | -----: | ----: | -: | --: |
| Poisson-binomial           |    0.9899 | 0.9623 | 0.9759 | 5,489 | 56 | 215 |
| PB-residual Mamba — seed A |    0.9921 | 0.9627 | 0.9771 | 5,491 | 44 | 213 |
| PB-residual Mamba — seed B |    0.9933 | 0.9628 | 0.9778 | 5,492 | 37 | 212 |
| PB-residual Mamba — seed C |         — |      — |      — |     — |  — |   — |

The important observation is that the neural model did not merely trade recall for precision.

Seed B:

```text
PB:
56 FP / 215 FN

Mamba:
37 FP / 212 FN
```

It removed 19 false positives while also recovering 3 false negatives.

---

## 4. Paired bootstrap

The primary comparison is against the Poisson-binomial caller.

### Seed A

```text
ΔF1 = +0.0012
95% CI = [+0.0001, +0.0024]

ΔTP = +2
ΔFP = -12
ΔFN = -2
```

### Seed B

```text
ΔF1 = +0.0019
95% CI = [+0.0009, +0.0030]

ΔTP = +3
ΔFP = -19
ΔFN = -3
```

Both confidence intervals exclude zero.

However, this is not sufficient to call the result a strong positive because replication across seeds is still limited.

---

## 5. The third seed

Seed C was run as an unconditional replication check.

Final result:

```text
ΔF1 = +0.0003
95% CI = [-0.0007, +0.0015]

ΔFP = -4
```

The confidence interval crosses zero.

Therefore:

```text
2 / 3 seeds → positive
1 / 3 seeds → statistically inconclusive
```

The correct verdict is:

# WEAK POSITIVE

Calling this a strong positive would overstate the evidence.

---

## 6. Why the improvement is scientifically interesting

The improvement was predicted before training from the information gap between the two models.

The Poisson-binomial caller models per-read base quality, but it does not model every property of the read evidence.

The most important complementary feature discovered was:

```text
ALT strand bias
```

At fixed alternate-read count, strand bias showed meaningful discriminative power.

This is important because it represents information that the Poisson-binomial model cannot directly encode.

The neural residual therefore has access to a genuinely different signal.

---

## 7. Error-overlap analysis

The important question is not simply:

> Did Mamba find more variants?

The important question is:

> Can Mamba remove Poisson-binomial errors without creating comparable new errors?

Seed B provides the strongest example.

```text
Poisson-binomial:
56 FP

Mamba:
37 FP
```

The model removed:

```text
19 PB false positives
```

while creating:

```text
0 new false positives
```

and additionally recovered:

```text
3 PB false negatives
```

This is fundamentally different from simply moving the classification threshold.

---

## 8. Constant-shift null

A major concern from previous experiments was that the neural residual might simply learn:

```text
PB score + constant
```

That would not demonstrate new information.

The constant-shift control rejected this explanation.

For each seed, the best free constant fitted on validation was:

```text
shift = 0
```

The neural models still achieved higher test F1 than the corresponding constant-shift PB caller.

Therefore the observed effect cannot be explained solely by moving the decision threshold.

---

## 9. Residual diagnostics

The learned residual is genuinely locus-specific.

Residual dispersion:

```text
Seed A: ~1.62
Seed B: ~1.44
Seed C: ~1.74
```

This is fundamentally different from the earlier 15× experiment where the residual behaved approximately like a constant shift.

The residual also correlates with ALT strand bias:

```text
r ≈ +0.17 ... +0.34
```

This supports the hypothesis that the model is learning information related to read-level artefacts.

However, the residual is not purely complementary.

It also correlates negatively with the Poisson-binomial LLR:

```text
r ≈ -0.34 ... -0.53
```

Therefore part of the learned correction appears to shrink or recalibrate the statistical prior.

The model is learning a mixture of:

1. complementary evidence;
2. statistical-score correction.

---

## 10. Depth-stratified behaviour

The improvement is concentrated where the benchmark actually contains errors.

At approximately:

```text
5–9×
10–14×
```

the model reduces false positives.

Above approximately:

```text
20×
```

both methods approach saturation.

This is important because it explains why simply increasing coverage did not solve the original benchmark.

At higher depth, the statistical caller becomes too accurate for the neural model to demonstrate meaningful additional information.

---

## 11. Controls

### Shuffled-prior control

The Poisson-binomial prior was shuffled between loci while preserving its distribution.

Result:

```text
F1 ≈ 0.7077
FP = 4,568
```

The model therefore depends heavily on the actual locus-specific statistical evidence.

The PB prior is not just an arbitrary numerical feature.

### Constant-shift control

Best validation-fitted constant:

```text
shift = 0
```

The neural improvement therefore cannot be reduced to a simple threshold displacement.

---

## 12. What was ruled out

Several seemingly promising directions were investigated and found unsuitable for this benchmark.

### Indels

The current Poisson-binomial formulation is not directly applicable to indels because a gap does not have a base-quality Phred score.

The indel classes were therefore not used as the primary neural-vs-PB comparison.

Interestingly, simple statistical baselines were already extremely strong:

```text
deletions: F1 = 1.0000
insertions: F1 = 0.9873
```

There was insufficient error budget.

### Low-mappability regions

The GIAB high-confidence BED largely excludes these difficult regions.

Observed mean MAPQ remained extremely high:

```text
minimum mean MAPQ ≈ 66.6
```

Therefore this dataset cannot meaningfully test the low-mappability hypothesis.

### Homopolymers

Homopolymers generate substantially more background noise, but true variants remained strongly separated:

```text
true variants: VAF ≈ 0.86
background:    VAF ≈ 0.105
```

They did not provide the desired competitive regime.

---

## 13. The real bottleneck

The major discovery of this experiment is methodological.

The original problem was not:

```text
Mamba is too small
```

It was:

```text
the benchmark was saturated
```

At 15× on the original region:

```text
PB:
0 FP
10 FN
```

There was simply not enough measurable headroom.

After increasing the evaluation region:

```text
5,704 SNPs
≈271 PB errors
```

the neural correction became detectable.

The model did not need to become dramatically larger.

The benchmark needed enough statistical power to distinguish models.

---

## 14. Runtime

The computational profile is also interesting.

For the powered 15× benchmark:

```text
Poisson-binomial prior:
~890 s CPU

Mamba inference:
~6.3 s GPU
```

The neural model is approximately:

```text
~142× faster per locus
```

than computing the Poisson-binomial prior on CPU.

This does NOT mean the complete neural system is 142× faster than the complete caller, because the neural model currently depends on the PB prior.

The correct interpretation is:

> Once the statistical prior is available, the neural correction is extremely cheap relative to computing the prior itself.

This distinction matters for future architecture design.

---

## 15. Scientific conclusion

The experiment provides evidence that:

> **A neural residual above the Poisson-binomial caller can learn additional locus-specific information that improves variant calling in a sufficiently powered 15× benchmark.**

However:

> **The magnitude and reproducibility of the improvement are not yet established strongly enough for a strong claim.**

Current evidence:

```text
Seed A: +0.0012 F1
Seed B: +0.0019 F1
Seed C: +0.0003 F1
```

Therefore the correct project status is:

# WEAK POSITIVE

This is the first experiment in which the neural model has demonstrated measurable additional information above the strongest statistical baseline.

---

## 16. What this means for the project

The project has passed an important milestone.

The question is no longer:

> "Can Mamba beat a simple binomial caller?"

It already can under some regimes.

The harder question is now:

> "Can the neural component provide a robust advantage large enough to matter against a strong statistical caller?"

That requires more evaluation rather than immediately increasing model capacity.

Recommended next steps:

1. Increase the evaluation region further, approximately 30–50 Mb.
2. Run 5–10 seeds.
3. Preserve the same statistical and leakage controls.
4. Measure whether the ΔF1 magnitude stabilizes.
5. Only after the effect is reproducible, investigate architectural improvements.
6. Test a hybrid architecture where the statistical component handles high-confidence local evidence while Mamba handles longer-range/contextual information.

Do not interpret the current +0.1–0.2 percentage-point improvement as a final benchmark victory.

The result is promising evidence, not a finished competitive claim.

---

## 17. Repository state

The powered benchmark was committed and pushed.

```text
Commit: 6c31876
Branch: worktree-30x-pb-residual
```

Main report:

```text
History/11_DEVLOG.md
```

Results:

```text
results/big15x/
```

Engineering verification:

```text
202 tests passed
```

All previous 15× and 30× artifacts were verified byte-identical.

No stale experiment or monitoring processes remained after cleanup.

---

# Final status

```text
15× original benchmark       → NULL
30× benchmark                → SATURATED
15× powered benchmark       → WEAK POSITIVE

Current best evidence:
PB-residual Mamba
      ↓
additional locus-specific information
      ↓
small F1 improvement
      ↓
replication not yet conclusive
```

The next bottleneck is no longer neural capacity.

It is **statistical power and reproducibility**.
