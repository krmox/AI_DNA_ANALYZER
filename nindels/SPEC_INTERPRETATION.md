# Indel Pattern Engine — spec interpretation (draft for confirmation)

Status: **not yet confirmed by the author**. Nothing is implemented. Items marked
**ASSUMPTION** are my choices where the spec is silent or ambiguous; items marked
**SPEC ISSUE** are places where I think the math does not do what the text claims.
Faithful versions are implemented as the baseline; my alternatives are benchmarked only
as named variants.

## 0. What the spec is (and is not)

The spec describes a **re-scorer of indel candidates**, not a caller. It assumes some
upstream process emits candidates plus a baseline score (the "AI_DNA_ANALYZER baseline",
§15 shadow mode). But AI_DNA_ANALYZER is **SNP-only** (README: "scope SNP-only"), so no
baseline indel caller exists to shadow.

**ASSUMPTION A0:** `nindels-call` = (1) a permissive C++ candidate generator that collects
every indel seen in read CIGARs, left-normalizes it and keeps it if n_alt ≥ 2, then
(2) the Pattern Engine scores every candidate. This is the "nindels" caller compared
against bcftools. In shadow mode the engine is also run over bcftools' own indel
candidates (variant `V-shadow-bcf`). In that mode bcftools' QUAL is *not* treated as a
calibrated probability (rule 7). Instead it is calibrated with a fitted logistic before
being used as P_base.

## 1. Pattern P_i = (t, L, m, c, h, r)  — spec §2

| Field | Implementation |
|---|---|
| t | INS / DEL after left-normalization against the reference |
| L | \|len(ALT) − len(REF)\| |
| m | inserted/deleted bases (after normalization) |
| c | **ASSUMPTION:** reference ±20 bp around the normalized position (41 bp, plus the deleted bases for DEL) |
| h | length of the homopolymer run that contains/abuts the event (in reference), plus its base |
| r | shortest tandem-repeat unit (1–6 bp) that tiles m; repeat unit length and copy count in the reference at the site; local entropy of c |

## 2. Evidence E_i — spec §2.1, §3, §4

n_alt and n_ref are counted per read (not per base) at the normalized locus. A read is
*informative* if it spans the event ±3 bp with MAPQ ≥ 20 (**ASSUMPTION**). n_alt = reads
whose CIGAR carries an equivalent (normalized) indel. n_ref = spanning reads without it.
Reads with a *different* indel at the locus count toward neither. f = n_alt/(n_alt+n_ref).
If the denominator is 0, f is set to NaN and explicitly flagged as missing.
L_depth = ln(1+n_alt+n_ref). q (mean BQ around the event, mean MAPQ), b (alt reads on the
forward strand vs reverse, Fisher-style SB), and a (mean distance to read end, soft-clip
fraction, MQ0 fraction) are recorded as diagnostics. Per the spec, they are **not** in
X_i (§12 lists only seven features), but see V-evidence below.

## 3. Distances and similarity — spec §5

- d_type ∈ {0,1}; d_length = min(1, \|ln(1+L_i) − ln(1+L_j)\| / s_L).
- d_seq(c) = Levenshtein(c_i, c_j) / max(len). Empty-empty = 0.
- **ASSUMPTION:** d_motif = normalized Levenshtein on m. d_hp = min(1, \|ln(1+h_i)−ln(1+h_j)\|/s_h).
  d_rep = mean of [unit-length mismatch, min(1, \|ln(1+copies_i)−ln(1+copies_j)\|/s_r)].
  The spec defines only type, length and context distances explicitly.
- D_ij = Σ w_k d_k / Σ w_k, S = 1 − D.
- **ASSUMPTION (weights):** the spec says w_k, s_L, s_h, s_r must be "learned" but gives no
  objective. Baseline: all w_k = 1, and s_* = the 95th percentile of the corresponding raw
  log-difference on training pairs. This is documented as arbitrary per rule 8. Variant
  `V-w-grid`: w chosen by grid search maximizing validation log-loss of P_context.

## 4. Memory classes — spec §6, §14

- **ASSUMPTION (clustering is unspecified):** bucket key = (t, length bin {1,2,3–5,6–10,11–20,21–50,>50},
  hp bin {0–3,4–7,8–12,>12}, repeat-unit length 0–6). Within a bucket, greedy leader
  clustering in a deterministic input order (sorted by locus): a pattern joins the nearest
  class if D < θ = 0.15, otherwise it starts a new class. The medoid is recomputed at the
  end over a bounded reservoir (≤ 64 members per class, deterministic hash sampling).
- Each class stores n_pos and n_neg from **GIAB truth labels on training data only**. These are
  independent labels, never the engine's own predictions.
- Labels: a candidate is positive if it matches a truth indel after normalization inside the
  GIAB confident BED. It is negative if it lies inside the BED with no truth match. Candidates
  outside the BED are not used for memory and are not scored.

## 5. Retrieval and aggregation — spec §7, §7.1

- **ASSUMPTION:** retrieval = the candidate's own bucket plus neighbouring length/hp bins of the
  same type, then the top K = 32 classes by S. The spec says "retrieved classes" without defining
  the set. This matters: summing S·n over *all* classes means a huge dissimilar class with
  S = 0.2 can dominate (see SPEC ISSUE 2).
- N_pos, N_neg, P_context = (N_pos+α)/(N_pos+N_neg+α+β), with α and β fitted on training data
  (beta-binomial method of moments over classes, **ASSUMPTION**). S_max = max S.
- N_eff = W² / Σ S²(n_pos+n_neg), exactly as written.

## 6. Local verification Δℓ — spec §8, §11

- H_ref/H_alt = reference window ±(read length) around the event, without or with the indel.
- P(R|H): a standard 3-state PairHMM (match/insert/delete), using base qualities for emissions
  and fixed gap-open Q45 / gap-extend Q10 (GATK defaults, **ASSUMPTION**). It is computed in
  log space, with the read realigned locally against each haplotype. MAPQ caps the per-read
  contribution (the per-read log-likelihood is mixed with a uniform-error floor ε = 10^(−MAPQ/10)).
- Δℓ = Σ_r log P(R_r|H_alt) − Σ_r log P(R_r|H_ref), **exactly as specified**.
- Q_k (§11): computed per position in the event window during reassessment, using the
  normalized BQ, MAPQ, fraction of reads whose local alignment agrees with the best haplotype,
  and support fraction. The weights are equal (arbitrary, documented). The spec's X_i does not
  include Q_k, so min_k Q_k is emitted as a diagnostic INFO field and does not change the call.

## 7. Aggressiveness and complexity — spec §9, §10

- C_r = Σ v·(H, R, L, G) with H = min(1, hp_len/15), R = min(1, repeat_span/50), L = 1 − Shannon
  entropy(c)/2, G = fraction of informative reads with MAPQ < 20 or soft clip > 10 bp.
  v = 0.25 each (arbitrary, documented).
- U_i needs a calibrated P_base, but there is no baseline caller. **ASSUMPTION:** P_base is a stage-1
  logistic fitted on training data using (f, L_depth, P_context, N_eff, S_max, C_r), i.e. X_i
  without Δℓ. It is calibrated on validation data before U is computed.
- D_i (disagreement) is undefined in the spec. **ASSUMPTION:** D_i = \|P_context − P_evidence\|,
  where P_evidence is a stage-1 logistic on (f, L_depth) only.
- a_U = a_C = a_D = 1/3. τ_strict and τ_loose are set to the 90th and 30th percentiles of R on
  training data. A is swept over {0, 0.25, 0.5, 0.75, 1}, and the cost and F1 are reported per A.
- A never enters X_i (rule 11). C_r enters X_i only as a feature with a learned coefficient
  (rule 12).

## 8. Final model — spec §12

- X_i = (f, L_depth, P_context, N_eff, S_max, Δℓ, C_r). Missing Δℓ (not reassessed) is handled with
  an explicit missingness indicator and the value set to 0 *only together with* the indicator. This
  is the standard missing-indicator method, so missing is not silently zero. Missing f gives
  UNCERTAIN.
- Logistic regression is fitted in C++ (IRLS, L2 = 1e-3) on the training split. τ_call is chosen on
  the validation split to maximize indel F1. Calibration is checked on test with a reliability
  table and Brier score.
- Output: PASS when P ≥ τ_call. FILTER=UNCERTAIN when f is missing, depth < 4, or the candidate is
  outside the training support (any feature outside the training 0.5–99.5 percentile range). Other
  candidates get FILTER=LowProb. Only PASS records are scored.
- **ASSUMPTION (genotype, not in spec):** GT = 1/1 if f ≥ 0.8, else 0/1. The spec has no genotype
  model, so F1 is reported both at the allele level (match ignoring GT) and at the genotype level.
- SNPs: nindels does not call SNPs, so SNP F1 = N/A for nindels.

## 9. SPEC ISSUES found before coding

1. **Δℓ (§8) is a hom-alt vs hom-ref test.** Under ℓ(H) = Σ log P(R|H), a heterozygous indel has
   about half its reads strongly favouring H_ref, so Δℓ is near 0 or negative for true hets, which
   are most germline indels. The faithful version will be implemented and measured. Variant
   `V-dl-diploid` uses Δℓ = Σ log(½P(R|ref)+½P(R|alt)) − Σ log P(R|ref) (plus a hom-alt term).
2. **Unbounded retrieval + S = 1−D weights.** S is bounded below by 0 but is rarely near 0, so
   with all classes summed, a large unrelated class contributes heavily. Restricting to top-K is
   necessary in practice. I treat that as part of "indexed retrieval" (§14), not a math change.
3. **N_eff (§7.1) is Kish's ESS and does not do what the text says.** If 1000 records are all
   near-identical (S ≈ 1), N_eff = 1000. Kish ESS only discounts *unequal weights*, not
   redundancy. The text's goal ("a thousand highly similar records should not count as a
   thousand") is not achieved by this formula. Implemented as written. Variant `V-neff-loci`
   counts distinct training loci per class with ±50 bp deduplication.
4. **The pattern memory carries no locus-specific evidence.** P_context is a prior per pattern
   class, essentially an estimate of the error rate for homopolymer/repeat contexts. Expect it to
   help precision in homopolymers and do little elsewhere. This is the central hypothesis (§17), and
   the benchmark should answer it with an ablation: the X_i model with vs without the three
   memory features.
5. **The evidence vector E_i is computed but mostly unused.** Strand bias, read position and soft
   clipping are classic indel-FP signals, yet only f and L_depth reach X_i. Variant `V-evidence`
   adds SB and mean read-end distance to X_i.
6. **U requires a calibrated baseline that does not exist** (see §7 above). This adds a stage-1 model.

## 10. Data splits (no leakage)

| Split | Sample / region | Use |
|---|---|---|
| train | HG002 chr20 (GIAB Illumina 2×250, GRCh38) | memory, s_*, α/β, stage-1 and final logistic coefficients, τ_strict/τ_loose |
| validation | HG002 chr21 | τ_call, calibration check, variant selection |
| test (touched once) | HG001 chr19, NIST v4.2.1 truth (BAM already on disk, ~10 GB) | final table, depth strata via downsampling to 15×/30×/50× |

The HG001 BAM's depth and platform still need to be checked. Difficult-region strata come from
GIAB v3.3 stratifications (homopolymer, tandem repeat, low-mappability BEDs). Long indels > 50 bp
are reported separately; short-read CIGAR-based candidates will mostly miss them, and the report
will say so.

## 11. Evaluation tool

**vcfdist** or **rtg vcfeval** will be chosen once the tools are available. hap.py needs Docker,
which is not installed. bcftools stats on its own does not do haplotype-aware matching, which is
essential for indels in repeats.
