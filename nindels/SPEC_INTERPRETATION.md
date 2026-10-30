# Indel Pattern Engine — spec interpretation (v2, after author review)

v1 was reviewed by the author on 2026-10-30. Their 10 corrections are folded in below and
marked **[C1]…[C10]**. **ASSUMPTION** = a default I chose; it may be wrong and is not
scientifically established. **SPEC ISSUE** = a place where the original math does not do what
its text claims. Every arbitrary constant is listed in §12 together with where it came from.

## 0. Scope: an independent caller, not a re-scorer  [C1, C9]

`nindels-call` is a complete indel caller with two stages:

1. **Candidate discovery from BAM** (C++, htslib). Every I/D operation in every read CIGAR is
   collected, left-normalized against the reference and grouped by allele (pos, REF, ALT).
2. **Pattern Engine scoring** of every candidate (spec §2–§12) → PASS / LowProb / UNCERTAIN VCF.

bcftools (mpileup | call) is an **external comparison caller only**. Optional shadow-mode run:
the same Pattern Engine scores bcftools' indel set (`--candidates-vcf`). It is always labelled
as such, and its recall is never reported as the module's recall. The candidate stage's own
**recall ceiling** (truth indels present among candidates) is reported separately, so losses
from discovery and losses from scoring are not mixed up.

Known structural limits of CIGAR-based discovery, stated up front:
- Indels the aligner did not represent as I/D are missed. This includes most > 50 bp events and
  soft-clipped breakpoints. Long-indel recall is reported, and it is expected to be poor.
- Complex or clustered events are scored one allele at a time (no multi-variant haplotype
  assembly).

## 1. Pattern P_i = (t, L, m, c, h, r)  — spec §2

| Field | Implementation |
|---|---|
| t | INS / DEL after left-normalization |
| L | \|len(ALT) − len(REF)\| |
| m | inserted/deleted bases |
| c | reference, 20 bp left of the event + 20 bp right of it (deleted bases excluded) — **ASSUMPTION** |
| h | length and base of the longest homopolymer run touching the event in the reference |
| r | smallest period p ∈ 1..6 that tiles m (0 if none); copies of that unit in the reference at the site; Shannon entropy of c |

## 2. Evidence E_i — every feature reaches the model  [C4]

Per candidate, using every primary, non-duplicate, QC-passed read that overlaps the event ±3 bp
(**ASSUMPTION**: no MAPQ floor at counting, so that MAPQ itself can be a feature):

| Group | Features (all enter X_i, §8) |
|---|---|
| counts | n_alt, n_ref, n_other (reads carrying a *different* indel here), f = n_alt/(n_alt+n_ref), L_depth = ln(1+n_alt+n_ref) |
| strand | alt forward fraction; strand-bias = \|log-odds ratio\| of (alt fwd, alt rev, ref fwd, ref rev) with +0.5 pseudo-counts |
| read position | mean of min(dist to read start, dist to read end)/read length over alt reads, and its difference from the same value over ref reads |
| soft clipping | fraction of alt reads with soft clip ≥ 10 bp; same fraction over all overlapping reads |
| MAPQ | mean MAPQ of alt reads; fraction of overlapping reads with MAPQ = 0; fraction with MAPQ < 20 |
| base quality | mean BQ of the ±5 flanking bases (plus inserted bases) in alt reads |
| alignment consistency | mismatch rate of alt reads in a ±20 bp window excluding the event; n_other/(depth); after PairHMM, the fraction of reads whose CIGAR-based allele agrees with their likelihood-best haplotype |

f undefined (denominator 0) → `missing`, never 0. Features that cannot be computed get an
explicit missingness indicator (§8).

## 3. Distances and similarity — spec §5 (as written)

d_type ∈ {0,1}; d_len = min(1, \|ln(1+L_i)−ln(1+L_j)\|/s_L); d_ctx = Lev(c_i,c_j)/max(len), with
empty-empty = 0. Additional distances (**ASSUMPTION**, because the spec names h and r but gives no
metric for them):
- d_motif = Lev(m_i,m_j)/max(len)
- d_hp = min(1, \|ln(1+h_i)−ln(1+h_j)\|/s_h)
- d_rep = ½·[period_i ≠ period_j] + ½·min(1, \|ln(1+copies_i)−ln(1+copies_j)\|/s_r)

D = Σw·d/Σw, S = 1 − D. The s_* values are set to the 95th percentile of the raw log-difference
over random training pairs. w = 1 for all features in the baseline. Variant `V-w-grid` picks w on
**validation** log-loss of P_context.

## 4. Memory classes — spec §6, §14  [C3]

- Index bucket = (t, length bin {1,2,3–5,6–10,11–20,21–50,>50}, hp bin {≤3,4–7,8–12,>12}, period 0–6).
- Inside a bucket: deterministic leader clustering in genomic order. A pattern joins the nearest
  leader if D < θ = 0.15, otherwise it becomes a new leader. Medoid = the member with minimal summed
  D over a reservoir of ≤ 64 members chosen by deterministic hash.
- n_pos and n_neg come **only from GIAB truth labels on the training split**. The module's own calls
  are never fed back into memory.
- The memory is a **prior over contexts**. It is one feature group in a model whose read-evidence
  features stay in place. Ablation is required: every reported model is run with and without the
  memory features.

## 5. Retrieval and aggregation — spec §7, §7.1  [C5]

- Retrieval = own bucket plus adjacent length and hp bins (same t), then the top K = 32 classes by S.
  N_pos, N_neg and S_max are computed from the K retrieved classes. P_context = (N_pos+α)/(N_pos+N_neg+α+β).
  α and β are estimated by beta-binomial method of moments over training classes.
- **N_eff_kish** (spec formula): W²/ΣS²·n. This is Kish's effective sample size. It measures how
  *evenly the similarity weight is spread across classes*. If one class takes nearly all the weight
  it shrinks toward the size of that class; if the weight is even it equals the total count. It does
  **not** discount repeated observations: 1000 near-identical records with S ≈ 1 give N_eff = 1000.
- **N_eff_loci** (alternative): the same Kish formula, but n_c is replaced by the number of *distinct
  training loci* in the class (events within 50 bp on the same chromosome are merged into one locus).
  The same locus seen at several depths, or in several reads, counts once. This measures how many
  independent genomic sites stand behind the prior.
- Both are features. They are compared as variants, never used together.

## 6. Local read-level verification — diploid genotype likelihoods  [C2]

- Haplotypes: H_ref = reference window [event − 2·readlen, event + 2·readlen]; H_alt = the same
  window with the candidate allele applied.
- P(R|H): a 3-state PairHMM (match/ins/del) with BQ-derived emissions and gap-open Q45 / gap-extend
  Q10 (GATK defaults, **ASSUMPTION**), computed in log space. The read is free to start anywhere
  on H (global in the read, local in the haplotype).
- Per-read MAPQ cap: log P' = log[(1−ε)·P(R|H) + ε·P_floor], with ε = 10^(−MAPQ/10). A
  mis-mapped read then cannot dominate the result.
- **Diploid GL** for g ∈ {0/0, 0/1, 1/1}:
  `GL(g) = Σ_r log[ (2−g)/2 · P(R_r|H_ref) + g/2 · P(R_r|H_alt) ]`.
  Under 0/1, each read is allowed to come from either haplotype.
- Features passed to the model: `QUAL_dip = log-sum-exp(GL01, GL11) − GL00` and the genotype
  posterior. GT = argmax over g of GL + log prior. Priors 0/0 : 0/1 : 1/1 = 1 : 1e-4 : 0.5e-4
  (arbitrary; the decision threshold is learned on top of them anyway).
- **Original Δℓ** (spec §8, kept as an intentionally limited baseline variant `V-dl-orig`):
  Δℓ = Σ log P(R|H_alt) − Σ log P(R|H_ref) = GL(1/1) − GL(0/0) without the MAPQ mixture. It is a
  hom-alt vs hom-ref test.
- Read cap: at most 250 reads per candidate (deterministic hash of the read name), because the
  BAMs are 300× before downsampling. Arbitrary, documented.
- **Count-based diploid GL** (always computed, no HMM): a binomial on n_alt/n_ref with an indel
  error rate e = 0.01. This gives every candidate a GT, and the stage-1 P_base has a genotype-aware
  input even when the HMM is skipped.

## 7. Nucleotide-level inspection — spec §11

Computed during reassessment, for each position k in the event ±10 bp:
Q_k = ¼(B_k + M_k + A_k + S_k), each term in [0,1]:
- B = mean BQ/40
- M = mean MAPQ/60
- A = fraction of reads whose local alignment at k agrees with their best haplotype
- S = 1 − \|support at k from HMM-best haplotypes − f\|

min_k Q_k and the position of the minimum are passed as features (with a missingness indicator)
[C4], and are also written to INFO for diagnosis. Q_k never classifies a base on its own (spec §11).

## 8. Final model — spec §12

- X_i = evidence features (§2) + memory features (§5) + Δℓ/QUAL_dip + GT posterior + C_r + min Q_k.
- Missing values: each group that can be absent (HMM features, Q_k, f) has an indicator m ∈ {0,1}.
  The value is set to the training mean whenever m = 1. The indicator carries the information.
- Logistic regression with standardized features, IRLS, L2 λ = 1e-2 (arbitrary; the variant
  `V-lambda` picks it on validation), fitted in C++.
- τ_call is the value that maximizes indel F1 on **validation**. Calibration (reliability bins,
  Brier score, ECE) is measured on validation and reported on test without refitting.
- FILTER: PASS if P ≥ τ_call; LowProb otherwise; UNCERTAIN if f is missing, depth < 4, or ≥ 2
  features fall outside the training 0.5–99.5 % range (out of distribution). Only PASS is scored.
- GT: from the HMM diploid GL when reassessed, otherwise from the count GL.

## 9. Aggressiveness and complexity — spec §9, §10

- C_r = ¼(H + R + L + G): H = min(1, hp/15), R = min(1, repeat_span/50), L = 1 − entropy(c)/2,
  G = fraction of overlapping reads with MAPQ < 20 or soft clip ≥ 10.
- P_base = a stage-1 logistic with X_i minus the HMM/Q_k features, calibrated on validation. Then
  U = 1 − \|2·P_base − 1\|.
- D_i = \|P_context − P_evid\|, where P_evid is a stage-1 logistic on the evidence features only.
- a_U = a_C = a_D = ⅓. τ_strict and τ_loose = the 90th and 30th percentiles of R on training.
- A ∈ {0, .25, .5, .75, 1} is swept on validation, reporting the fraction reassessed, runtime and F1.
  A **never** enters X_i. The final model is trained with every candidate reassessed (A = 1) and
  applied with missingness indicators when A < 1. **ASSUMPTION**; the alternative of training one
  model per A is listed as a variant.

## 10. SPEC ISSUES (carried from v1 and resolved by the corrections)

1. Δℓ is a hom test → replaced by diploid GL [C2]; the original is kept as `V-dl-orig`.
2. Summing over all classes → top-K retrieval.
3. Kish N_eff ≠ independence → both versions kept and explained [C5].
4. Memory is a context prior → ablated in every result table [C3].
5. E_i was unused → all features enter X_i [C4].
6. U needs a calibrated baseline → stage-1 model.
7. (new) A model trained only at A = 1 sees no missing HMM features in training. The missingness
   indicators would then be unidentifiable. Fix: during training, a deterministic 50 % of
   candidates have their HMM features masked (missingness augmentation). **ASSUMPTION.**

## 11. Data and evaluation protocol  [C7, C10]

| Split | Data (GRCh38, NIST v4.2.1 truth + confident BED) | Used for |
|---|---|---|
| train | HG002 chr20, NHGRI HiSeq 2×148 NovoAlign 300× downsampled to 15/30/50× | memory, s_*, α/β, every coefficient, τ_strict/τ_loose |
| validation | HG002 chr21, same platform and depths | τ_call, λ, w grid, variant selection, A choice, calibration |
| test (touched once, after the config is frozen) | HG001 chr19, same platform, 15/30/50× | final tables only |

- Platform: same platform for train and test, because of what is available. "Performance by
  platform" is therefore **unknown**. A second platform would be needed to measure it.
- HG002 and HG001 are unrelated individuals (AJ son vs CEPH). A different chromosome gives a second
  independent axis.
- A frozen config (`config/frozen.json` + memory file + model file, all with sha256) is committed
  **before** any test run. The test runs once.
- Scoring: **rtg vcfeval** (the haplotype-aware engine inside hap.py; bundled JRE, no Docker). Run in
  genotype mode and in allele mode (`--squash-ploidy`). Stratification uses the annotated
  TP/FP/FN VCFs plus a Python script:
  - length {1, 2–5, 6–15, 16–50, > 50}
  - genotype {het, hom-alt}
  - depth {15, 30, 50×}
  - repeat context: GIAB v3.3 stratification BEDs for homopolymer ≥ 7 and tandem repeats, plus "other"
  - clustered: another truth or call variant within 10 bp
- Raw VCFs, logs, `/usr/bin/time -v` output and commands are kept under `/mnt/hdd/nindels/runs/<run_id>/`.

## 12. Register of arbitrary parameters

| Parameter | Value | Origin |
|---|---|---|
| context flank | 20 bp | arbitrary |
| span padding for counting | 3 bp | arbitrary |
| min n_alt to keep a candidate | 2 | arbitrary; sets the recall ceiling (measured) |
| min alt fraction to keep a candidate | 0.03 | arbitrary; recall ceiling (measured) |
| w_k | 1 | arbitrary (spec says "learn"); `V-w-grid` |
| s_L, s_h, s_r | 95th percentile on training | data-derived (training) |
| θ (cluster radius) | 0.15 | arbitrary |
| K (retrieved classes) | 32 | arbitrary |
| reservoir per class | 64 | arbitrary |
| α, β | method of moments | data-derived (training) |
| PairHMM gap open / extend | Q45 / Q10 | GATK default, not re-fitted |
| reads per HMM | 250 | arbitrary (cost cap) |
| count-GL error rate e | 0.01 | arbitrary |
| genotype priors | 1 : 1e-4 : 5e-5 | arbitrary |
| v_1..v_4, a_U, a_C, a_D | equal | arbitrary (spec gives no values) |
| Q_k weights | equal | arbitrary |
| τ_strict / τ_loose | 90th / 30th percentile of R | data-derived (training), percentiles arbitrary |
| L2 λ | 1e-2 | arbitrary; `V-lambda` on validation |
| τ_call | max validation F1 | data-derived (validation) |
| OOD rule | ≥ 2 features outside 0.5–99.5 % | arbitrary |
