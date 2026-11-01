# Read-level pileup representation — 15× HG002, chr21

**Date:** 2026-08-11
**Verdict: NEGATIVE.** The read-level representation did not provide measurable
information beyond the statistical baseline under this dataset and regime. It
was worse — and, uniquely among the experiments so far, the extra information
was *actively harmful*: every control that destroys read-level detail
**improved** the model.

---

## 1. Hypothesis

An aggregate caller sees only `(n, k, ε)`. It cannot distinguish "10 alternate
reads of 30, all low-quality, clustered at read ends, all on one strand" from
"10 alternate reads of 30, high-quality, well-mapped, strand-balanced". If that
distinction carries signal, a model that observes individual reads should beat
the binomial caller, especially at 5–14× where the binomial is weakest.

Prediction if true: higher recall at similar precision, concentrated at low
depth. Prediction if false: no gain, and the controls that destroy read-level
structure cost nothing.

## 2. Representation

`[loci, 48 reads, 14 per-read features]`, uint8-packed, plus the locus
reference one-hot.

Per-read features: base one-hot (A/C/G/T/gap), base quality, mapping quality,
strand, normalized read position, normalized distance from nearest read end,
near-end flag (<5 bp), insertion-next, deletion-next, and **low-quality-base**
— a read whose base the aggregate caller discards entirely.

**R = 48 chosen from the data, not by convention.** Maximum observed depth is
43 (train) / 38 (test) under the project's MAPQ and flag filters, so
**truncation rate = 0.000000** and the deterministic sampling path never
activates. No observation is discarded anywhere in either region. (The
sampling path is implemented and unit-tested so the representation stays
correct on deeper data.)

**Read order** is a deterministic BLAKE2b hash of the read name — stable across
processes (Python's salted `hash` would not be), and statistically independent
of every feature, so BAM traversal order cannot become a hidden feature.

`mismatch_density` (NM/read length) is extracted and cached but **excluded from
the primary feature set**: it is alignment-derived rather than truth-derived,
but it partially encodes reference-discordance, which the protocol restricts.

## 3. Architecture

```
[B, L, 48, 14] reads
  → shared per-read MLP encoder (14→64→128)
  → masked attention pooling + masked mean pooling   (permutation-invariant)
  + explicit depth term  log1p(n_valid)/5
  + locus reference one-hot
  → [B, L, 128] locus embedding
  → the UNMODIFIED 6-layer MambaBlock stack (width 128, dropout 0.15)
  → per-locus 4-class logits
```

Flattening `[R, F]` was rejected: it would make the model order-dependent —
the exact thing being controlled for — and scale parameters with `R`.

## 4. Parameter count

| model | parameters |
|---|---|
| aggregate `RawPileupMamba` (14ch) | 412,064 |
| **`ReadLevelMamba`** | **420,961 (+2.2%)** |

The read encoder is shared across all 48 reads and all 64 loci, so widening the
input from 14 numbers to 48×14 costs ~9k parameters. The backbone is identical.

## 5. Dataset

Unchanged from every prior 15× experiment. Train/val: chr21 31–33 Mb (1,876,352
loci, 2,103 SNP), sequential 90/10 split by window. Test: chr21 30.0–30.47 Mb,
455,488 loci, 549 SNP / 454,822 Normal / 40 ins / 77 del, median depth 15×.
Same BAMs, same GIAB truth VCFs, same high-confidence BEDs.

## 6. Leakage controls

* No VCF is opened during feature construction; labels come from a separate
  inherited path and are used only for evaluation.
* No feature encodes "differs from reference". The extractor never receives the
  reference base while building read rows (unit-tested: identical reads yield
  identical features regardless of the reference).
* The candidate alternate allele, used **only in the post-hoc error analysis**,
  is the most-supported non-reference base from observed counts — the same rule
  the binomial caller uses.
* The locus reference one-hot *is* supplied to the model (see §7 — bug found).
  It is reference-genome information, not truth, and every baseline in the
  comparison already has it (channels 9–13 of the aggregate representation).
* Feature names are asserted free of truth-derived tokens by unit test.

## 7. Sanity checks — and one real bug

All 13 gate checks passed on **both** regions before training
(`results/readlevel/sanity_{test,train}.json`):

| check | result |
|---|---|
| **aggregate counts reconstructed exactly** | **0 mismatched loci of 455,488; max abs difference 0.0 in all 10 columns** |
| labels/region/BAM identical to aggregate cache | pass |
| features finite, bounded, padded rows all-zero | pass |
| extraction deterministic; matches cache byte for byte | pass |
| read order not sorted by any feature | 0 of 5,000 deep loci monotone |
| truncation rate / zero-depth loci | 0.000 / 0 |
| strand balance (forward fraction) | 0.4998 |
| fraction of admitted reads the binomial **discards** | **19.5%** |

**Bug found and fixed mid-experiment (documented, not hidden).** The first
training run reached only SNP AUC 0.74 / F1 0.21. Cause: the read-level model
had no access to the locus reference base, so "every read says G" was
uninterpretable — it cannot be a variant call without knowing what the
reference says. Every baseline (binomial, 14-channel Mamba) already receives
this. I added the reference one-hot at the locus level (+768 parameters), which
restored AUC to 0.999, and re-ran. All reported results use the fixed model.
The uncorrected run is *not* reported as a result; it was a defect in my
design, not a finding.

## 8. Training configuration

Identical to the previous 15× runs, imported from `train_raw_pileup` rather
than restated: seq_len 64, batch 16, AdamW lr 3e-4, focal γ = 3.0, calibrated
class weights (α = 0.75, cap 20), mutation-presence sampler, 8 epochs, seed
20260811, sequential 90/10 coordinate-order split. Checkpoint and threshold
selected on validation `mutation_macro_f1` via the same
`sweep_thresholds` grid the aggregate model used — deliberately the *aggregate
model's* protocol, so "read-level vs aggregate" is a one-variable comparison.
45 s/epoch on an RTX 3070.

## 9. Main results — SNP vs Normal, 455,371 loci, 549 SNPs

All thresholds frozen from validation. Test never used for selection.

| method | P | R | F1 | 95% CI | ROC-AUC | PR-AUC | TP | FP | FN |
|---|---|---|---|---|---|---|---|---|---|
| frequency | 0.9942 | 0.9381 | 0.9653 | [0.9534, 0.9761] | 0.99998 | 0.9920 | 515 | 3 | 34 |
| **binomial v1** | 0.9981 | 0.9490 | **0.9729** | [0.9624, 0.9824] | 0.99998 | 0.9939 | 521 | **1** | 28 |
| binomial v2 | 1.0000 | 0.9454 | 0.9719 | [0.9612, 0.9816] | 0.99998 | 0.9940 | 519 | 0 | 30 |
| mamba 14ch | 0.9693 | 0.9217 | 0.9449 | [0.9308, 0.9578] | 0.99997 | 0.9855 | 506 | 16 | 43 |
| **read-level** | 0.8848 | 0.9089 | **0.8967** | [0.8780, 0.9148] | 0.99809 | 0.9433 | 499 | **65** | 50 |
| read-level (seed 424242) | 0.8938 | 0.9199 | 0.9066 | [0.8874, 0.9242] | 0.99649 | 0.9508 | 505 | 60 | 44 |

Paired bootstrap ΔF1 (2000 identical resamples):

| comparison | ΔF1 | 95% CI | P(better) |
|---|---|---|---|
| read-level − binomial v1 | **−0.0762** | [−0.0949, −0.0578] | 0.000 |
| read-level − binomial v1 (seed 2) | **−0.0663** | [−0.0847, −0.0487] | 0.000 |
| read-level − mamba 14ch | **−0.0482** | [−0.0673, −0.0288] | 0.000 |

The read-level model is significantly worse than the binomial baseline **and**
significantly worse than the aggregate neural model it was meant to improve on.
Both CIs exclude zero and the result replicates on a second seed.

This is not threshold miscalibration. At **oracle** (test-optimal, diagnostic
only) thresholds the ordering is unchanged: binomial v1 0.9790, mamba 14ch
0.9565, read-level 0.8999.

## 10. Depth-stratified results

F1, with TP/FP/FN:

| depth | SNP | binomial v1 | mamba 14ch | read-level | read-level decoupled (C6) |
|---|---|---|---|---|---|
| 1–4× | 7 | **0.923** 6/0/1 | 0.800 6/2/1 | 0.600 6/7/1 | 0.737 7/5/0 |
| **5–9×** | 79 | **0.912** 67/1/12 | 0.833 65/12/14 | **0.667 63/47/16** | 0.824 68/18/11 |
| **10–14×** | 209 | **0.970** 197/0/12 | 0.958 194/2/15 | 0.939 191/7/18 | 0.948 192/4/17 |
| 15–19× | 162 | **0.991** 159/0/3 | 0.965 151/0/11 | 0.952 150/3/12 | 0.968 153/1/9 |
| 20–29× | 91 | **1.000** 91/0/0 | 0.989 89/0/2 | 0.978 88/1/3 | 0.994 90/0/1 |
| 30+× | 1 | 1.000 | 1.000 | 1.000 | 1.000 |

**The hypothesis predicted the read-level model would help most at 5–14×. It
did the opposite**: at 5–9× it produced 47 false positives against the
binomial's 1, the single worst stratum in the experiment.

## 11. Ablations and controls

All trained with the identical pipeline; only the input information changes.
Tensor shape and parameter count are held fixed (ablations *zero* a feature
group rather than removing it), so these measure information, not capacity.

| arm | what is destroyed | F1 | Δ vs full read-level |
|---|---|---|---|
| **full read-level** | — | **0.8967** | — |
| C2 permute reads | read order only | 0.8967 | **+0.0000 (bit-identical)** |
| C1 aggregate-only | all per-read variation | 0.9117 | **+0.0150** |
| C4 no strand | strand | 0.9316 | **+0.0349** |
| C5 no read-position | position, end-distance, near-end | 0.9361 | **+0.0394** |
| C3 no quality | base + mapping quality, low-qual flag | 0.9378 | **+0.0411** |
| **C6 decouple base↔attributes** | which read carries which base | **0.9393** | **+0.0427** |

Three things follow, and they are unusually clean:

**CONTROL 2 is exactly satisfied.** Permuting whole read rows reproduces the
primary result bit for bit (identical TP/FP/FN/AUC/AP). The aggregation is
provably permutation-invariant (unit-tested) and the deployed pipeline confirms
it numerically. Read order is not a hidden feature.

**CONTROL 6 is the decisive one, and it is inverted.** The decoupling permutes
each read's attribute bundle across reads *stratified by counted status*, which
preserves **all ten aggregate count columns byte-for-byte** — allele counts,
depth, gaps, insertions, `quality_sum`, `mapping_sum`. Everything the binomial
caller and the 14-channel model can see is untouched; only "which read carried
which base under what conditions" is destroyed. That is the exact quantity the
hypothesis is about. Destroying it **improves F1 by 0.043**, and brings the
model to statistical parity with the aggregate Mamba (Δ = −0.0056, CI
[−0.0213, +0.0105]).

**Every attribute group is individually harmful.** Removing quality, strand, or
position each improves the model by 0.035–0.041. There is no subset of the
read-level attributes whose inclusion helps.

## 12. Error analysis

Read-level model vs binomial v1 (frozen thresholds):

| | count |
|---|---|
| binomial FNs recovered by the model | **11** of 28 |
| binomial FPs corrected | 1 of 1 |
| **new FPs introduced** | **65** |
| **binomial TPs overturned** | **33** |
| shared FNs | 17 |

Against the declared success criterion ("530 TP / 2 FP would be interesting;
530 TP / 30 FP is not"), this is 499 TP / 65 FP — a clear failure.

Profiling the alternate-supporting reads (candidate ALT from counts only):

| locus set | n | alt BQ median | ref BQ median | alt fwd-strand | alt reads |
|---|---|---|---|---|---|
| **model-only recovered** | 11 | **40** | 39 | 0.50 | 2 |
| **model-only false positives** | 65 | **16** | 44 | 0.00 | 1 |
| shared FN | 17 | 40 | 40 | 0.33 | 2 |
| TP (both methods) | 488 | 40 | 45 | 0.50 | 10 |

This is the mechanism, and it cuts both ways. The 11 genuine recoveries are
exactly what the hypothesis predicts: low depth (median 9×), VAF 0.30, LLR 7.96
— just under the binomial's 10.5 cutoff — with **high-quality** (BQ 40),
strand-balanced alternate reads. The model really did use read-level evidence
to rescue them.

But the 65 false positives have alternate-allele base quality **median 16**
against reference-read quality 44, typically **one** marginal-quality alternate
read at depth 7, and zero strand balance. The model learned to treat
marginal-quality, single-read, strand-biased evidence as *support* rather than
as noise — the opposite of the correct inference, and precisely the failure
mode the extra information was supposed to prevent. Net: +11 / −65 / −33.

## 13. Runtime

Extraction 84 s (test region) / 294 s (train region), cached once and consumed
byte-identically by every arm. Storage 721 MB uint8 (train), 29 MB compressed.
Training 45 s/epoch. Inference 1.67 s for 455,488 loci = 273k loci/s, versus
0.51 s for the aggregate model and 0.21 s for the binomial caller — i.e. the
read-level model costs ~3× the aggregate model and ~8× the binomial, for worse
accuracy.

## 14. Statistical uncertainty

All F1 values carry percentile bootstrap 95% CIs over 2000 resamples;
method-to-method comparisons use a **paired** bootstrap on identical resamples,
because the methods share most of their errors and independent intervals badly
overstate the uncertainty in a difference. The primary negative result
(Δ = −0.0762, CI [−0.0949, −0.0578]) excludes zero by a wide margin and
replicates on a second seed (−0.0663, CI [−0.0847, −0.0487]). With 549 test
SNPs, differences below ~0.01 F1 are not resolvable here; the observed
difference is 8× that.

## 15. Scientific interpretation

The evidence does **not** support the hypothesis. Three distinct readings can
be separated cleanly:

1. **Is read-level information useful here? No — it is harmful.** CONTROL 6
   holds every aggregate statistic byte-identical and destroys only base↔
   attribute association; performance *improves* by 0.043. If the association
   carried exploitable signal, destroying it could not help. Every individual
   attribute ablation points the same way.

2. **Is this an architecture failure or an information failure?** Partly
   architecture, but that does not rescue the hypothesis. With read-level
   detail destroyed (CONTROL 6), the model reaches parity with the aggregate
   14-channel Mamba (Δ = −0.0056, CI spanning zero) — so the encoder/pooling
   design is sound and reproduces the aggregate model when given aggregate
   information. Adding genuine read-level detail on top of that *degrades* it.
   The representation is not being under-exploited; it is being over-fitted.

3. **Does the neural approach beat the statistical baseline? Still no.** Every
   neural variant in this project now sits below binomial v1: 0.9449 (14ch),
   0.9390 (LLR), 0.9729 (residual, by construction), 0.8967 (read-level). The
   binomial caller remains undefeated at 0.9729, and a simple threshold move to
   9.5 takes it to 0.9778.

The mechanism is visible in the error analysis: at 15× coverage, most loci
offer 1–2 alternate reads, so "the alternate reads are low-quality and
strand-biased" is a statement about a sample of size one or two. There is
essentially no per-locus read-level *distribution* to characterize at this
depth. The model fits noise in those one-or-two-read patterns, and the fit does
not transfer — which is exactly why removing the information helps.

## 16. Limitations

* **Mapping quality is nearly constant in this BAM** (mean 69.997, min 26, max
  70). The "alternate reads have unusual mapping quality" hypothesis is
  effectively **untestable** on this data and is not refuted by this experiment.
* **549 test SNPs**, of which only 86 lie below 10× — the regime of interest.
  A larger low-depth truth set is needed to resolve small effects there.
* Base qualities in this BAM reach 115 with a mean of 37, which is not standard
  Illumina encoding; quality-dependent conclusions may not transfer to other
  data.
* Only one architecture and one training budget (8 epochs, 45 s each) were
  tried. The read-level model's validation trajectory was noticeably noisier
  than the aggregate model's (SNP F1 swinging 0.73–0.92 across epochs), which is
  consistent with over-fitting but also means a longer or regularized schedule
  was not explored. This was deliberate — the protocol forbids changing
  multiple variables — but it bounds the claim.
* Read-level indel evidence was captured but indels (117 test loci) were not
  the primary metric.
* The conclusion is specific to **15× Illumina HG002 chr21**. It says nothing
  about somatic calling at low VAF, long reads, or higher depth, where
  per-locus read distributions are actually populated.

## 17. Does the evidence support the hypothesis?

**No.** Stated plainly, as required:

> The read-level representation did not provide measurable information beyond
> the statistical baseline under this dataset and regime.

Stronger than a null: the read-level association is *actively harmful* at this
depth, and every control that destroys it improves the model. The one piece of
the hypothesis that survived is narrow and real — 11 low-depth variants were
recovered on genuinely high-quality alternate evidence the binomial rejected —
but it is swamped by 65 false positives from the same mechanism running in
reverse.

## 18. Exact next experiment

Not a bigger network, and not more read-level features. The result says the
per-locus read sample at 15× is too small to characterize. Two candidates, in
priority order:

**(a) Fix the binomial's error model instead of replacing it.** The single
most valuable number in this report is that binomial v1's threshold is
untuned: 10.5 → 9.5 gives F1 0.9778, beating every neural model ever trained
here, and its oracle is 0.9790. A per-locus ε estimated from the *observed*
base-quality distribution of the alternate-supporting reads (not the flat 0.01,
and not v2's mean-quality version, which underperformed) is a 20-line change
that directly encodes the "are the alternate reads trustworthy" intuition this
experiment tried to learn — as a statistic, where 1–2 reads is enough, rather
than as a learned pattern, where it is not.

**(b) If read-level is to be revisited, change the regime, not the model.**
Run the identical pipeline at 30× and 74× on the same loci. The prediction is
explicit and falsifiable: if the failure is "too few alternate reads to
characterize a distribution", the read-level model's deficit against the
binomial should shrink monotonically with depth, and CONTROL 6's *inversion*
should disappear. If the deficit persists at 74×, read-level structure is
simply not informative for germline SNP calling in this data and the line
should be closed.

## 19. Reproduction

```bash
# PHASE 1-3: extract, gate, verify exact count reconstruction
python cache_read_level.py --fasta data/reference/chr21_full.fa \
  --bam data/giab_hg002_real_15x/hg002_chr21_30M_15x.bam \
  --vcf data/giab_hg002_real_30M/hg002_chr21_30M.vcf.gz \
  --bed data/giab_hg002_real_30M/hg002_chr21_30M_highconf.bed \
  --region 29999999 30469999 --out data/giab_hg002_15x_readlevel/test_15x_reads.npz
python sanity_read_level.py --reads data/giab_hg002_15x_readlevel/test_15x_reads.npz \
  --counts cache/test_15x_counts.npz --out results/readlevel/sanity_test.json

# PHASE 4: primary model (+ seed replicate)
python train_read_level.py --train-reads data/giab_hg002_15x_readlevel/train_15x_reads.npz \
  --train-counts cache/train_15x_counts.npz --epochs 8 --tag readlevel_primary \
  --out results/readlevel/readlevel_primary_train.json

# PHASE 6: controls and ablations
for C in aggregate_only decouple permute_reads; do
  python train_read_level.py ... --control $C --tag readlevel_$C; done
for G in quality strand position; do
  python train_read_level.py ... --drop-groups $G --tag readlevel_drop_$G; done

# PHASE 5: evaluation + verdict
python evaluate_read_level.py --test-reads ... --test-counts cache/test_15x_counts.npz \
  --freq-threshold 0.335 --binom-v1-threshold 10.5 --binom-v2-threshold 17.5 \
  --aggregate-checkpoint checkpoints_rawpileup_15x/best_rawpileup_15x.pt \
  --readlevel-checkpoint checkpoints_readlevel/best_readlevel_primary.pt \
  --extra-checkpoint <name>=<path> ...
python summarize_read_level_experiment.py --training-json results/readlevel/*_train.json
```

Artefacts: `results/readlevel/readlevel_experiment_summary.json` (verdict),
`readlevel_experiment_results.json`, `sanity_{test,train}.json`, and one
`*_train.json` per arm. Data in `data/giab_hg002_15x_readlevel/`, checkpoints in
`checkpoints_readlevel/`. No existing dataset, checkpoint, or prior experiment
script was modified.

Tests: `test_read_level.py`, 29 tests — exact count reconstruction, hash-order
determinism, mask correctness, permutation invariance, ablation shape
invariance, control-preservation properties, torch/numpy feature-builder bit
equality, and absence of truth-derived features. Full suite: **85 passed**
(`test_bimamba.py` / `test_providers.py` fail to collect for a pre-existing
unrelated reason — they import a `bimamba_variant_caller` package absent from
this layout).
