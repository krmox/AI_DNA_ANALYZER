# Residual Mamba over a binomial prior — 15× HG002, chr21

**Date:** 2026-08-11
**Verdict: NULL.** The residual model, selected the way it would actually be
deployed, is *identical* to the binomial v1 caller. The one arm that shows a
gain does not replicate across seeds and is reproducible by adding a scalar to
the LLR.

---

## 1. Question

Can the Mamba architecture extract information from the raw pileup that
**improves upon** the binomial LLR caller, rather than merely relearning or
distorting the same statistic?

The two previous experiments could not answer this. `mamba_14ch` (F1 0.9449)
and `mamba_llr_15ch` (F1 0.9390) both lost to binomial v1 (F1 0.9729), and
because each learned its decision function from scratch, a loss was
uninformative: it could mean "no extra information exists" or "the model
failed to rediscover a statistic it was handed as one channel among fifteen."
Initializing *at* the binomial caller removes that ambiguity. Any change from
the starting point is, by construction, the model's own contribution.

## 2. Method

### Architecture

```
final_logit = binomial_prior_logit + residual(x)
```

The binomial LLR is computed from raw integer pileup counts by the existing,
unmodified `BinomialVariantCaller` (v1: fixed ε = 0.01, `max(H_het, H_hom)`
over `H0`), then placed in a 4-class prior logit vector
(`residual_prior.py`):

| class | prior logit |
|---|---|
| Normal | 0 |
| **SNP** | **LLR** |
| Insertion | 0 |
| Deletion | 0 |

Indel channels are set equal to Normal, not to a large negative number: the
binomial caller has no opinion about indels, and pinning them would have
smuggled in a non-binomial prior. Consequently
`prior[SNP] − prior[Normal] = LLR` exactly, so the model's SNP decision score
*is* the binomial LLR and the baseline's threshold transfers unchanged.

The residual branch is the unmodified `RawPileupMamba` backbone over the same
14-channel raw-pileup representation (same width 128, depth 6, dropout 0.15;
412,064 parameters). Its final classifier weight **and** bias are zeroed, so
`residual(x) ≡ 0` at initialization for any input and any backbone state. The
layer is not dead — its gradient is `δ · final_norm(h)`, non-zero as soon as a
prediction is wrong — so training escapes zero on the first step.

### Scoring in logit space (documented protocol deviation)

`softmax(prior)[SNP]` saturates to exactly 1.0 in float32 above LLR ≈ 17,
which would tie together most true SNPs, and the binomial's operating point
(LLR 10.5 ⇒ p = 0.99992) lies beyond the old 0.50–0.99 probability grid
entirely. This is a representational impossibility, not a preference. All
scores and thresholds therefore use the SNP−Normal **logit gap**, swept over
`arange(-20, 200, 0.5)` — the *same grid* `evaluate_binomial_baseline` used to
select binomial v1's own 10.5. The selection objective is SNP F1, matching how
the binomial threshold was chosen. The residual model and the baseline are now
selected by an identical procedure on an identical axis. Nothing else changed:
same train/val/test regions, seed, batch size, sequence length, optimizer, LR,
epochs, class weighting, focal γ = 3.0, and mutation-presence sampler.

### Initialization gate (run before any training)

`sanity_residual_init.py` on the 455,488 held-out test loci, two backbone
seeds:

| quantity | result |
|---|---|
| Pearson (model logit, LLR) | **1.000000000000** |
| Spearman | **1.000000000000** |
| max abs logit difference | **7.16 × 10⁻⁶** (float32 storage only) |
| ROC-AUC | 0.9999805 — identical to binomial |
| PR-AUC | 0.9938519 — identical to binomial |
| F1 @ frozen threshold 10.5 | 0.972923 — identical |
| TP / FP / FN | 521 / 1 / 28 — identical |
| call-level disagreements | **0 of 455,488** |

The gate passed. The untrained model is not "close to" the binomial baseline;
it is the binomial baseline.

### Training and λ sweep

λ ∈ {0, 1e-4, 1e-3, 1e-2} on `λ · mean(residual²)`, 8 epochs each, selection
on validation only. Two arms are reported, both validation-only:

* **`residual_selected`** — best over all candidates *including epoch 0*, the
  untrained model. Epoch 0 is a legitimate candidate: it is what early
  stopping returns when no amount of training improves on the prior.
* **`residual_best_trained`** — best over epochs 1–8 only, needed because a
  residual pinned at exactly zero has nothing to diagnose.

**No trained epoch at any λ beat the untrained model on validation.**
Validation F1: epoch 0 = 0.9960; the best trained epoch anywhere in the sweep
= 0.9940 (λ = 1e-2, epoch 1). Selection returned epoch 0, λ = 0.

## 3. Results — SNP vs Normal, 15× test region (455,371 loci, 549 SNPs)

Every threshold frozen from validation; the test set was never used for
selection.

| method | P | R | F1 | 95% CI | ROC-AUC | PR-AUC | TP | FP | FN |
|---|---|---|---|---|---|---|---|---|---|
| frequency | 0.9942 | 0.9381 | 0.9653 | [0.9534, 0.9761] | 0.99998 | 0.9920 | 515 | 3 | 34 |
| **binomial v1** | 0.9981 | 0.9490 | **0.9729** | [0.9624, 0.9824] | 0.99998 | 0.9939 | 521 | 1 | 28 |
| mamba 14ch | 0.9693 | 0.9217 | 0.9449 | [0.9308, 0.9578] | 0.99997 | 0.9855 | 506 | 16 | 43 |
| mamba+LLR 15ch | 0.9257 | 0.9526 | 0.9390 | [0.9235, 0.9530] | 0.99997 | 0.9805 | 523 | 42 | 26 |
| **residual (selected)** | 0.9981 | 0.9490 | **0.9729** | [0.9624, 0.9824] | 0.99998 | 0.9939 | 521 | 1 | 28 |
| residual (best trained) | 0.9981 | 0.9563 | 0.9767 | [0.9671, 0.9851] | 0.99997 | 0.9936 | 525 | 1 | 24 |
| *shuffled-prior control* | 0.4021 | 0.8871 | 0.5534 | [0.5242, 0.5800] | 0.99606 | 0.2805 | 487 | 724 | 62 |

Paired bootstrap ΔF1 against binomial v1 (2000 resamples, identical draws):

| method | ΔF1 | 95% CI | P(better) |
|---|---|---|---|
| residual (selected) | **+0.0000** | [+0.0000, +0.0000] | — |
| residual (best trained) | +0.0038 | [+0.0009, +0.0081] | 0.978 |
| residual (best trained, **seed 424242**) | **+0.0000** | [−0.0029, +0.0029] | 0.341 |
| mamba 14ch | −0.0280 | [−0.0409, −0.0153] | 0.000 |
| mamba+LLR 15ch | −0.0340 | [−0.0495, −0.0197] | 0.000 |
| shuffled-prior control | −0.4195 | [−0.4506, −0.3919] | 0.000 |

### Depth-stratified (F1, primary frame)

| depth | loci | SNP | binomial v1 | residual (best trained) | mamba 14ch | mamba+LLR |
|---|---|---|---|---|---|---|
| 1–4 | 4,137 | 7 | 0.9231 | 0.9231 | 0.8000 | 0.7692 |
| 5–9 | 54,981 | 79 | 0.9116 | 0.9116 | 0.8333 | 0.7922 |
| 10–14 | 160,626 | 209 | 0.9704 | 0.9780 | 0.9580 | 0.9486 |
| 15–19 | 153,458 | 162 | 0.9907 | 0.9938 | 0.9649 | 0.9730 |
| 20–29 | 79,425 | 91 | 1.0000 | 1.0000 | 0.9889 | 0.9891 |
| 30+ | 2,744 | 1 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |

The entire best-trained difference lives in the 10–19× bins (+3 and +1 TP).
Below 10×, where the binomial actually struggles (F1 0.91), the residual
changes nothing at all.

## 4. Diagnostics — what did the residual actually learn?

This is the part that matters, and it is unambiguous.

**The residual carries essentially no locus-specific information.**

| diagnostic | best-trained residual |
|---|---|
| residual on all loci | mean **+0.273**, sd 0.589 |
| residual on SNP loci | mean **+0.244**, sd 0.548 |
| residual on Normal loci | mean **+0.273**, sd 0.589 |
| **residual alone as a SNP score, ROC-AUC** | **0.530** |
| residual alone, PR-AUC | 0.0013 (base rate 0.0012) |
| corr(residual, LLR) | +0.016 |
| corr(residual, depth) | −0.007 |
| corr(residual, VAF) | ≈ 0 |

The residual's distribution on true SNPs is *indistinguishable* from its
distribution on Normal loci — if anything it is very slightly **lower** on
SNPs. Mean residual by depth is flat (+0.29 → +0.24 across 1× to 30+×); by
VAF flat (+0.20 → +0.40); by LLR flat (+0.25 → +0.20). Its standalone ROC-AUC
of 0.530 against a base rate of 0.12% is noise. What the model learned is a
**near-constant positive bias of ≈ +0.27 logits**, which is not evidence — it
is a threshold slide.

**Error-level accounting (best-trained arm, vs binomial v1):**

| | count |
|---|---|
| binomial FNs recovered | **4** of 28 |
| binomial FPs corrected | 0 of 1 |
| new FPs introduced | **0** |
| binomial TPs overturned | **0** |
| total call changes | 4 |

Four recovered FNs, no cost. That looks good in isolation. It does not
survive two controls:

**Control A — the constant-shift null.** Replacing the whole network with its
own mean residual (+0.272), i.e. evaluating the plain binomial at the
equivalent threshold 9.728, recovers 2 of those 4 loci (F1 0.9748 vs the
model's 0.9767). And simply scanning the binomial's own single knob does
better than the neural model outright: **plain binomial @ 9.5 gives F1 0.9778**
(TP 528, FP 3), and at its test-optimal 6.0 gives **F1 0.9790**. The residual
model's oracle is 0.9808 — a 0.002 edge over a baseline with no network at
all, obtained by 412k parameters.

**Control B — seed replicate.** Re-running λ = 1e-2 with seed 424242 gives a
best-trained arm with **ΔF1 = +0.0000, CI [−0.0029, +0.0029]**, recovering 1
FN while overturning 1 TP. Its mean residual is **−0.607** — the *opposite
sign* from the first seed's +0.272. The direction of the "learned correction"
is not even stable across initializations. The +0.0038 does not replicate.

**Ablation — shuffled prior.** Training with the LLR permuted across loci
(marginal preserved, association destroyed) collapses to F1 0.5534 (P 0.402,
724 FP). This confirms the prior carries essentially all the signal and that
the residual branch on its own cannot substitute for it. Note the honest
limitation: because the shuffle destroys the prior rather than the residual,
this control demonstrates that the *binomial statistic* is doing the work — it
is not, by itself, evidence about the residual's quality. The residual-alone
AUC of 0.530 is the direct evidence there, and it is negative.

## 5. Conclusion

**Classification: NULL.**

Under the selection rule that would actually be deployed — validation-only,
with the untrained model as a legitimate candidate — the residual Mamba is
*bit-for-bit the binomial v1 caller*: 521 TP, 1 FP, 28 FN, ΔF1 = 0.0000. Given
a model initialized exactly at the binomial caller, 1.69 M training loci, four
regularization strengths and 8 epochs each, gradient descent found **no**
modification of the raw-pileup evidence that validation preferred over leaving
the statistic alone.

The forced-training arm classifies as *weak positive* in isolation
(+0.0038 F1, CI [+0.0009, +0.0081], 4 FNs recovered at zero FP cost), and this
is reported honestly — but it cannot carry the headline. It fails all three
follow-ups: validation did not select it, it does not replicate across seeds
(the residual's sign flips), and its effect is reproduced by adding a scalar
to the LLR — a manoeuvre the baseline can perform on its own, and which at
threshold 9.5 beats the neural model outright.

This is a **negative answer to the scientific question, not a failure of the
experiment**. The experiment is now clean: the previous ambiguity — "did the
Mamba lose because there is no extra information, or because it never found
the statistic?" — is resolved. Handed the statistic for free and asked only to
improve on it, the model had nothing to add. At 15× on GIAB HG002 chr21, the
14-channel raw pileup appears to contain no SNP-relevant signal beyond what
the depth-aware binomial likelihood ratio already extracts.

The most valuable practical finding is negative in a second way: the binomial
caller's threshold is under-tuned. Moving it from 10.5 to 9.5 gains more
(F1 0.9729 → 0.9778) than the entire neural residual did, at a cost of 2 FP.
That is where the next hour of work belongs, not in a larger network.

### What would change the verdict

The residual's flatness against depth, VAF and LLR suggests the 14 channels
genuinely lack per-locus discriminative signal beyond `(n, k, ε)`. Candidates
worth testing before concluding the pileup is exhausted: strand-bias and
position-in-read features (absent from the current 14 channels, and the
classic source of the FP/FN the binomial model cannot see), and the low-depth
1–9× regime where the binomial is weakest (F1 0.91) but only 86 test SNPs
exist — too few to resolve a difference here, so a larger low-depth test
region is needed first.

## 6. Reproduction

```bash
# 1. gate (must pass before training)
python sanity_residual_init.py --test-counts cache/test_15x_counts.npz \
    --binom-v1-threshold 10.5

# 2. lambda sweep
for L in 0 1e-4 1e-3 1e-2; do
  python train_residual_binomial.py --train-counts cache/train_15x_counts.npz \
    --residual-lambda $L --epochs 8 --tag residual_lambda$L \
    --out results/residual_lambda${L}_train.json
done

# 3. shuffled-prior control + seed replicate
python train_residual_binomial.py --train-counts cache/train_15x_counts.npz \
  --residual-lambda 1e-2 --epochs 8 --shuffle-prior --tag residual_shuffled \
  --out results/residual_shuffled_train.json
python train_residual_binomial.py --train-counts cache/train_15x_counts.npz \
  --residual-lambda 1e-2 --epochs 8 --seed 424242 \
  --tag residual_lambda1e-2_seed424242 \
  --out results/residual_lambda1e-2_seed424242_train.json

# 4. validation-only selection, then test evaluation, then verdict
python select_residual_lambda.py --training-json results/residual_lambda*_train.json
python evaluate_residual_experiment.py --test-counts cache/test_15x_counts.npz \
  --freq-threshold 0.335 --binom-v1-threshold 10.5 \
  --old-checkpoint checkpoints_rawpileup_15x/best_rawpileup_15x.pt \
  --llr-checkpoint checkpoints_rawpileup_llr/best_llr_real.pt \
  --residual-selected-checkpoint checkpoints_residual_binomial/best_residual_lambda0.pt \
  --residual-best-trained-checkpoint checkpoints_residual_binomial/epoch1_residual_lambda1e-2.pt \
  --shuffled-checkpoint checkpoints_residual_binomial/best_residual_shuffled.pt
python summarize_residual_experiment.py
```

Artefacts: `results/residual_experiment_summary.json` (machine-readable
verdict), `results/residual_experiment_results.json`,
`results/residual_init_sanity.json`, `results/residual_lambda_selection.json`,
`results/residual_seed424242_test.json`. New checkpoints live in
`checkpoints_residual_binomial/`; no existing checkpoint, dataset, BAM, VCF,
BED or prior experiment script was modified.

Tests: `test_residual_binomial.py` (17 tests — exact zero-initialization,
gradient escape, prior fidelity, monotonicity, and byte-equality of the LLR
against `BinomialVariantCaller`). Full suite: 56 passed.
`test_bimamba.py` and `test_providers.py` fail to collect on `main` for an
unrelated pre-existing reason (they import a `bimamba_variant_caller` package
that does not exist in this layout); untouched by this work.
