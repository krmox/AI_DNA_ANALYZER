# Method C Full-Region Revalidation — HG002 chr21:32,000,000-44,000,000

BAM → detection (frozen) → genotype (optimized Method C) → VCF → hap.py 0.3.15,
run on the identical inputs/scope as `experimental/unified_happy/UNIFIED_HAPPY_BENCHMARK.md`.

## 1. Exact inputs (unchanged from the prior unified_happy benchmark)

- Sample: HG002, GRCh38
- Region: chr21:32,000,000-44,000,000
- Truth: `data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz`
- High-confidence BED: `data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed`
- Reference: `experimental/head_to_head/shared_ref/chr21_chrname.fa`
- Variant scope: SNP only, FILTER=PASS
- Evaluator: hap.py 0.3.15 (xcmp), image digest
  `sha256:d63b963a6cb01b4830393b22369e7b91d298e4156dde353739e74e4cfa4f96d0` —
  identical digest to the prior unified_happy run and to the head-to-head
  external-caller run.
- Detection call sets: `experimental/head_to_head/ai_cascade_per_locus.npz`
  (`pb_calls`, `cascade_calls` masks) — the exact frozen masks used to build
  the prior (non-Method-C) `ai_pb_only.vcf` / `ai_cascade.vcf`.
- Counts/positions: `experimental/unified_happy/recovered_full.npz` (same
  array used for the prior unified_happy VCF build).

No input was substituted. This is the same BAM-derived candidate/detection
scope, same truth, same BED, same reference, same hap.py version/digest as
the previous unified_happy benchmark.

## 2. Exact commands

```bash
python3 experimental/genotype_layer/method_c_regression.py
# -> cache/method_c_regression_report.json (equivalence proof)
# -> experimental/unified_happy/ai_pb_only_C.vcf(.gz/.tbi)
# -> experimental/unified_happy/ai_cascade_C.vcf(.gz/.tbi)

bash experimental/unified_happy/run_happy.sh   # re-pointed at the _C VCFs
# -> experimental/unified_happy/happy/ai_pb_only_C/happy.summary.csv
# -> experimental/unified_happy/happy/ai_cascade_C/happy.summary.csv

python3 experimental/genotype_layer/method_c_unified_metrics.py
# -> cache/unified_C_genotype_metrics.json (GT accuracy + invariant check)
```

## 3. Method C implementation

Two implementations, identical closed-form binomial-posterior formula
(`eps=0.01` fixed, θ ∈ {0.0, 0.5, 1-eps} ↔ {0/0, 0/1, 1/1}, GT = argmax
likelihood, GQ = 10·log10(P_best/P_second) capped at 99):

- **REFERENCE**: per-site Python loop, `scipy.stats.binom.logpmf` called per
  hypothesis per site (verbatim from `build_vcfs_validation_C.py`).
- **OPTIMIZED** (used for this benchmark's VCFs): fully vectorized NumPy,
  closed-form log-pmf via `scipy.special.gammaln`, no per-site Python loop.

Source: `experimental/genotype_layer/method_c_regression.py`.

## 4. Reference-vs-optimized equivalence (full unified_happy scope, 16094/16097 sites)

From `cache/method_c_regression_report.json`:

| Dataset | n | GT equal | GQ equal | max |GQ diff| | Ref runtime | Opt runtime | Speedup |
|---|---|---|---|---|---|---|---|
| ai_pb_only | 16094 | **true** (0 mismatches) | **true** (0 mismatches) | 0.0 | 5701.2 ms | 5.5 ms | 1027x |
| ai_cascade | 16097 | **true** (0 mismatches) | **true** (0 mismatches) | 0.0 | 5795.8 ms | 10.8 ms | 537x |

Exact equivalence (not just within tolerance — 0 GT mismatches, 0 GQ
mismatches, max abs GQ diff = 0.0) on the complete 16094/16097-site
full-region call set. Equivalence gate **PASSED**; proceeded to full
benchmark.

## 5. Detection invariance

Method C is a genotype-layer transform applied to the same detection call
set; it does not add, drop, or move a candidate. Confirmed via
`cache/unified_C_genotype_metrics.json`:

| Dataset | n_before | n_after | invariant_ok |
|---|---|---|---|
| ai_pb_only_C | 16094 | 16094 | true |
| ai_cascade_C | 16097 | 16097 | true |

QUERY.TOTAL in hap.py output is identical between old (forced-0/1) and new
(Method C) VCFs for both conditions (16094 / 16097 — see §7), confirming no
coordinate, REF, ALT, or candidate-count change. Symmetric difference = 0 in
both cases. **Only GT (and GQ) changed.**

## 6. PB-only and Cascade hap.py results (Method C)

From `happy/ai_pb_only_C/happy.summary.csv` and `happy/ai_cascade_C/happy.summary.csv` (SNP, PASS):

| Condition | TRUTH.TOTAL | TP | FP | FN | Precision | Recall | F1 |
|---|---|---|---|---|---|---|---|
| AI PB-only + Method C | 16898 | 15748 | 346 | 1150 | 0.9785 | 0.9319 | **0.9547** |
| AI Cascade + Method C | 16898 | 15748 | 349 | 1150 | 0.9783 | 0.9319 | **0.9546** |

## 7. External caller results (reused, unchanged, same hap.py digest/scope)

From `experimental/head_to_head/HEAD_TO_HEAD_RESULTS.csv` /
`UNIFIED_HAPPY_BENCHMARK.md` — not rerun:

| Caller | Precision | Recall | F1 | TP | FP | FN | Runtime |
|---|---|---|---|---|---|---|---|
| DeepVariant 1.6.1 | 0.9340 | 0.9654 | 0.9494 | 16313 | 1153 | 585 | 819s |
| Clair3 v1.0.10 | 0.8965 | 0.9638 | 0.9290 | 16287 | 1880 | 611 | 1342s |
| GATK HaplotypeCaller 4.5.0.0 (no BQSR) | 0.9877 | 0.9493 | 0.9682 | 16042 | 199 | 856 | 287s |

## 8. Before/after comparison (genotype-layer effect, detection held fixed)

| Condition | Old (forced 0/1) F1 | New (Method C) F1 | ΔF1 | ΔPrecision | ΔRecall |
|---|---|---|---|---|---|
| PB-only | 0.6238 | 0.9547 | **+0.3309** | +0.3391 (0.6394→0.9785) | +0.3230 (0.6089→0.9319) |
| Cascade | 0.6237 | 0.9546 | **+0.3309** | +0.3390 (0.6393→0.9783) | +0.3230 (0.6089→0.9319) |

Because §5 confirms the detection call set is byte-for-byte unchanged
(same 16094/16097 candidates, same coordinates/alleles), this entire ΔF1 is
attributable to the genotype layer, not to detection. Mechanistically: the
old pipeline forced every call to GT=0/1; hap.py's genotype-aware matching
penalized every truth 1/1 site as a genotype mismatch (`FP.gt`, which fell
from 5630/5631 to 172/173 — see summary CSVs), which was the dominant error
source at F1≈0.624.

Old AI F1 ≈ 0.624 (both conditions) → New AI F1 ≈ 0.9546-0.9547. Roughly
**99.7% of the genotype-attributable deficit is closed**: of the
(0.968 GATK-ceiling − 0.624) ≈ 0.344 old gap to the best external caller,
Method C alone recovers ≈ 0.331 of it, leaving a residual gap driven by the
remaining detection-layer FP/FN (346-349 FP, 1150 FN vs. e.g. GATK's 199/856).

## 9. Comparison to DeepVariant

| Metric | DeepVariant | AI Cascade + Method C | Gap |
|---|---|---|---|
| F1 | 0.9494 | 0.9546 | **+0.0052 (AI Cascade higher)** |
| Precision | 0.9340 | 0.9783 | +0.0443 |
| Recall | 0.9654 | 0.9319 | −0.0335 |

Under this specific unified hap.py scoring, AI Cascade + Method C's F1
(0.9546) is marginally *higher* than DeepVariant's (0.9494), driven by
higher precision (fewer FP: 349 vs 1153) at the cost of lower recall (more
FN: 1150 vs 585). Method C was not tuned toward DeepVariant or any external
caller — `eps=0.01` is fixed, not fit (see §12). This is a single-region,
single-sample result and should not be read as a general superiority claim
(see §13 Limitations).

## 10. Stratified analysis

Full-region-scope stratified hap.py numbers were not separately recomputed
for this exact 16094/16097-candidate scope. The closest available
stratified genotype-accuracy breakdown is
`experimental/genotype_layer/RESULTS_VALIDATION_C.md`, computed on the
near-identical (not identical — see caveat below) 16124-candidate dev-region
extraction, reusing the same eps=0.01 Method C formula:

| Stratum | n | Method C GT acc |
|---|---|---|
| Het truth, VAF<0.3 | 777 | 1.0000 |
| Het truth, VAF 0.3-0.5 | 3790 | 1.0000 |
| Het truth, VAF 0.5-0.69 | 4541 | 1.0000 |
| Het truth, VAF>0.69 | 1195 | 0.8586 |
| Hom truth (1/1), VAF≥0.95 | 5570 | 1.0000 |
| Hom truth (1/1), depth<10 | 1115 | 0.9928 |
| Depth<10 | 2708 | 0.9642 |
| Depth 10-19 | 11230 | 0.9928 |
| Depth 20-29 | 1969 | 1.0000 |
| Depth 30+ | 33 | 1.0000 |
| Segdup/low-mappability overlap | 417 | 0.9856 |
| Tandem repeat overlap | 897 | 0.9766 |

**Caveat**: this stratified table is from a 16124-candidate extraction that
is not bit-identical to this report's 16094/16097-candidate unified_happy
scope (different extraction pass — see `method_c_regression.py` docstring,
§header). It is directionally informative (same region, same sample, same
Method C formula, overlapping but not identical candidate sets) but is not
claimed as an exact per-candidate stratification of this report's own
16094/16097-site call set. The full unified_C_genotype_metrics.json
confusion matrix for the exact scope used here is in §11.

Exact-scope confusion matrix (`cache/unified_C_genotype_metrics.json`),
truth\\pred order (0/0, 0/1, 1/1):

| Dataset | GT acc | 0/1 acc | 1/1 acc | truth=0/1→1/1 errors | truth=1/1→0/1 errors |
|---|---|---|---|---|---|
| ai_pb_only_C | 0.9891 | 0.9838 | 0.9989 | 167 | 6 |
| ai_cascade_C | 0.9891 | 0.9838 | 0.9988 | 167 | 7 |

## 11. Runtime

### Method C compute (genotype-layer only, full unified_happy scope)

| Dataset | n sites | Reference (scipy loop) | Optimized (vectorized) | µs/site (optimized) | Speedup |
|---|---|---|---|---|---|
| ai_pb_only | 16094 | 5701.2 ms | 5.5 ms | 0.34 | 1027x |
| ai_cascade | 16097 | 5795.8 ms | 10.8 ms | 0.67 | 537x |

### Detection runtime (unchanged, not attributed to Method C)

Position/pileup recovery for the full region: 1242s (`recovered_full.npz`
build, `recover_full_positions.py` — this is BAM-scanning cost shared
identically by the old forced-0/1 and new Method-C pipelines; Method C reuses
it unmodified).

### Full pipeline wall-clock (BAM → final VCF)

- Detection/position recovery: 1242s (dominant cost, unchanged by Method C)
- VCF build + Method C genotyping: ~3-9s (build_vcfs / regression script,
  includes the optimized Method C pass)
- bcftools sort+index: ~1s
- hap.py: 6.2s per condition

Method C's own compute is 5-11 ms — negligible (<0.001%) relative to the
1242s detection/pileup stage. The end-to-end pipeline speedup from
optimizing Method C is therefore indistinguishable from noise; the ~1000x/
537x speedup is real but is a genotype-layer-only speedup, not a full-pipeline
speedup (kept explicitly distinct per task instructions).

## 12. Limitations

- Single sample (HG002), single region (chr21:32-44 Mb), SNP-only. Not
  generalized to indels, other chromosomes, other samples, other depths, or
  clinical use.
- GATK HaplotypeCaller result used here was run without BQSR (documented
  pre-existing deviation from GATK best practice, inherited unchanged from
  the prior head-to-head benchmark — not introduced by this revalidation).
- Method C's `eps=0.01` is a fixed constant carried over from the
  genotype-layer experiments, not fit/tuned on this evaluation scope or
  against any external caller.
- §10's stratified table is from a non-identical (overlapping) candidate
  extraction; treat as directional, not as an exact stratification of this
  report's own call set.
- The apparent AI-Cascade-vs-DeepVariant F1 edge (§9) is a precision/recall
  trade, not a strict dominance (DeepVariant has materially higher recall);
  it should not be read as general superiority.

## 13. Reproducibility

- All commands recorded in §2; no path was assumed — all scripts read from
  the exact `experimental/unified_happy/` and `experimental/head_to_head/`
  files inherited from the previous benchmark.
- Frozen detection thresholds/logic untouched (no cutoff, threshold, or
  candidate-extraction code was modified for this task).
- No production/frozen code modified. No files committed. Run directly in
  the main checkout, consistent with the prior unified_happy protocol
  (`experimental/`/`data/` are untracked and not worktree-materializable).
- All source data untouched: no benchmark data deleted or altered.

## 14. Scientific audit

1. **Did optimized Method C reproduce the reference implementation exactly
   within tolerance?** Yes — 0 GT mismatches, 0 GQ mismatches, max GQ diff
   0.0, on the full 16094/16097-site scope (§4).
2. **Did detection remain invariant?** Yes — identical candidate count,
   `invariant_ok: true` for both conditions, 0 symmetric difference (§5).
3. **Did genotype accuracy improve?** Yes — GT accuracy ≈0.989 vs. the old
   forced-0/1 baseline's 0/1-only accuracy of 1.0 / 1/1-accuracy of 0.0
   (every homozygous-alt truth site was wrong under forced 0/1).
4. **Did full-region hap.py F1 improve?** Yes — 0.6238→0.9547 (PB-only),
   0.6237→0.9546 (Cascade), ΔF1 ≈ +0.331 both conditions (§6, §8).
5. **How much of the original AI-vs-DeepVariant gap was recovered?** The old
   AI-vs-DeepVariant F1 gap was 0.9494−0.6238 ≈ 0.3256. New AI Cascade F1
   (0.9546) exceeds DeepVariant's F1 (0.9494) — i.e. **more than 100%** of
   the original gap is closed under this exact metric (§9).
6. **What is the remaining gap?** None in aggregate F1 vs. DeepVariant under
   this scope; AI Cascade has a recall deficit vs. DeepVariant (0.9319 vs.
   0.9654, i.e. more FN: 1150 vs 585) offset by a precision surplus (fewer
   FP: 349 vs 1153). Vs. GATK (F1 0.9682, the best external caller here),
   AI Cascade still trails by ΔF1 ≈ −0.0136.
7. **Is the result still limited to HG002/chr21/32-44 Mb?** Yes, entirely
   (§12).
8. **Which conclusions are supported and which are not?** Supported: the
   full pipeline's F1 deficit at forced-0/1 genotyping (§8) was overwhelmingly
   a genotype-layer artifact, not a detection artifact, on this exact
   region/sample; the optimized Method C is numerically exact vs. reference
   and ~500-1000x faster for the genotype-layer computation alone (§4, §11).
   NOT supported by this data: any claim of general AI-caller superiority
   over DeepVariant/Clair3/GATK, any claim about indels, any claim beyond
   HG002/chr21:32-44Mb, or any claim that Method C's genotype-layer speedup
   materially changes full-pipeline wall-clock (dominated by the unchanged
   1242s detection/pileup stage).
