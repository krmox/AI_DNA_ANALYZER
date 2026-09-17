# Root Cause Analysis: F1 0.974 (internal evaluator) → 0.624 (unified hap.py)

See `ERROR_DECOMPOSITION.csv`, `TRUTH_VARIANT_CLASSIFICATION_{ai_pb_only,ai_cascade}.csv`,
`AI_FP_CLASSIFICATION_{ai_pb_only,ai_cascade}.csv`, `ALLELE_LEVEL_METRICS_*.csv`,
`GENOTYPE_DISTRIBUTION_*.csv`, `DENOMINATOR_ANALYSIS.csv`,
`PB_VS_CASCADE_DIFFERENCES.csv` for the full underlying data. All classifications
below are derived directly from hap.py's own per-locus `happy.vcf.gz` `BD`/`BK`/`BLT`
FORMAT fields (ground truth from hap.py itself, not re-derived) for every one of the
17,459 truth SNP records and every one of the 5804-5807 query FP records — full
population, not a sample.

## 1. Old internal evaluator semantics (quoted, not inferred)

`accuracy_block()`, `robustness_benchmark.py:77-86`:
```python
def accuracy_block(calls: np.ndarray, truth: np.ndarray) -> dict:
    stats = prf(calls, truth)
    true_negative = int(np.sum(~calls & ~truth))
    total = truth.size
    sensitivity = stats["tp"] / max(stats["tp"] + stats["fn"], 1)
    specificity = true_negative / max(true_negative + stats["fp"], 1)
    stats["tn"] = true_negative
    stats["accuracy"] = (stats["tp"] + true_negative) / max(total, 1)
    stats["balanced_accuracy"] = 0.5 * (sensitivity + specificity)
    return stats
```
calls `prf()`, `evaluate_binomial_baseline.py:41-49`:
```python
def prf(predicted: np.ndarray, truth: np.ndarray) -> dict:
    tp = int(np.sum(predicted & truth))
    fp = int(np.sum(predicted & ~truth))
    fn = int(np.sum(~predicted & truth))
    ...
```
`predicted`/`truth` are plain boolean arrays, position-aligned by array index only
(no explicit coordinate ever compared). **No REF, no ALT, no GT, no zygosity is
ever considered anywhere in this code path.** A locus counts as TP iff
`calls[i]==True and truth[i]==True` at the same array index — pure
variant/no-variant agreement.

`truth` comes from `providers.py::_fetch_labels` (`LABEL_SNP` per
substitution-truth-VCF-record, skipping homozygous-reference `0/0` genotype
records — `providers.py:814`, `_is_homozygous_reference`). **Truth denominator**
is the count of `label==1` loci across all *kept* windows.

**BED handling**: not per-locus, per-window. `providers.py::_build_windows`
(lines 571-599) calls `_is_confident(start, end)` (lines 550-569) once per
64bp tiling window and drops the **entire window** (all 64 loci) if any part
of it lies outside the high-confidence BED — see Section 6 below for the
quantitative consequence.

**Duplicate/multiallelic handling**: none explicit; `_fetch_labels` loops
`for alternate in record.alts`, marking whichever offsets the allele implies,
without any allele-identity bookkeeping — a locus is simply "SNP" or not.

## 2. hap.py semantics (from actual output files, not general knowledge)

From `happy/ai_pb_only/happy.runinfo.json`'s `final_args`: `engine: "xcmp"`,
`preprocessing_leftshift: true`, `preprocessing_truth_confregions: true`,
`fp_bedfile: <highconf.bed>`, `window: 50` (local-match haplotype comparison
window).

From `happy/ai_pb_only/happy.vcf.gz`'s header (`##FORMAT` lines): each
TRUTH/QUERY sample column carries `GT:BD:BK:BI:BVT:BLT:QQ`, where `BD` =
"Decision for call (TP/FP/FN/N)", `BK` = "Sub-type for decision
(match/mismatch type)", `BVT` = "High-level variant type (SNP|INDEL)",
`BLT` = "High-level location type (het|homref|hetalt|homalt|nocall)".

Empirically decoded from the actual data (all 17,459 SNP records, both AI
VCFs):
- `BD=TP, BK=gm` (10290 records): genotype match — full agreement.
- `BD=FN` on truth side paired with `BD=FP` on query side, both `BK=am`
  (5630/5631 records): **"am" = allele match** — REF/ALT are correct, GT
  differs (truth `homalt`/1-1, query `het`/0-1 in every single one of these
  records, confirmed exhaustively, 0 exceptions). hap.py counts this as
  *both* a truth FN and a query FP (it is not "half credit" — one real
  variant call becomes two error counts).
- `BD=FN, BK=.` on truth with no matching query record at all (962 records):
  **genuinely absent AI call**.
- `BD=FN, BK=lm` (16 records): "local match" — a nearby/complex-region
  match that did not resolve to an exact match at this specific
  representation; hap.py counts these as FN. `FP.al` in the summary tracks
  the query-side mirror of this (29/30 records) — this is what "allele
  mismatch" concretely consists of in this dataset (representation
  differences, not wrong bases).
- `BD=UNK` (561 records, all with `Regions` INFO tag absent, i.e. no
  `CONF`/`TS_contained` tag): truth SNPs **outside the high-confidence BED**
  entirely — hap.py excludes these from `TRUTH.TOTAL` (confirmed: manually
  inspected e.g. chr21:32,005,718, which pysam confirms is not covered by
  either flanking BED interval `[32,000,000-32,005,617)` /
  `[32,005,742-32,007,886)` — a genuine 125bp BED gap).

`QUERY.FP=5804` (ai_pb_only summary) = `5630 (am) + 145 (.) + 29 (lm)` —
verified exact.

## 3-5. Position/allele join and classification (full population, not sampled)

Every truth SNP and every AI FP record was classified directly from the
`happy.vcf.gz` `BD`/`BK` codes (see Section 2's mapping). Results
(`TRUTH_VARIANT_CLASSIFICATION_ai_pb_only.csv`, n=17,459 total,
16,898 counted by hap.py + 561 UNK):

| Classification | Count | % of TRUTH.TOTAL (16898) |
|---|---|---|
| EXACT_MATCH | 10290 | 60.9% |
| GENOTYPE_ONLY_ERROR | 5630 | 33.3% |
| TRUE_MISSED_VARIANT | 962 | 5.7% |
| REPRESENTATION_ERROR | 16 | 0.1% |
| OTHER (=UNK, outside BED) | 561 | (excluded from denominator) |

AI FP classification (`AI_FP_CLASSIFICATION_ai_pb_only.csv`, n=5804):

| Classification | Count | % of FP |
|---|---|---|
| GENOTYPE_MISMATCH | 5630 | 97.0% |
| TRUE_FALSE_POSITIVE | 145 | 2.5% |
| REPRESENTATION_MISMATCH | 29 | 0.5% |

No genuine `ALLELE_ERROR` bucket exists in this dataset (0 records where the
wrong non-reference base was called at a real variant site with correct
GT) — `FP.al`/representation-error records (16/29) are hap.py "local match"
haplotype-window edge cases, not base-calling errors. No `OTHER`/unclassifiable
residual beyond the BED-boundary UNK bucket, which was manually inspected
(Section 2) and fully explained, not left ambiguous.

AI Cascade classification is essentially identical (5631/961/16/561 and
5631/146/30) — see the `_ai_cascade` CSVs.

## 6. Allele-level (GT-blind) diagnostic

`ALLELE_LEVEL_METRICS_ai_pb_only.csv`: TP_allele=15920 (EXACT_MATCH +
GENOTYPE_ONLY_ERROR), FP_allele=174 (TRUE_FALSE_POSITIVE +
REPRESENTATION_MISMATCH), FN_allele=978 (TRUE_MISSED_VARIANT +
REPRESENTATION_ERROR).

**Precision_allele = 0.9892, Recall_allele = 0.9421, F1_allele = 0.9651.**

This is a **diagnostic figure only** (it answers "if genotype correctness
were not required, how well does the caller find the right position+allele"),
not a replacement for the hap.py metric, and is reported alongside it, not
instead of it.

## 7-8. Genotype distribution and the forced-0/1 effect

`GENOTYPE_DISTRIBUTION_ai_pb_only.csv`: of 16898 truth SNPs, 11120 (65.8%)
are truth `het` (0/1) and 5777 (34.2%) are truth `homalt` (1/1) (1 `hetalt`
multiallelic, negligible).

AI PB-only/Cascade emit `GT=0/1` for literally 100% of calls (hard-coded in
`build_vcfs.py`, by design, documented in advance as a limitation — this is
confirmed empirically too: `homalt_exact_match=0` and `het_wrong_gt=0` across
all 16898+ records, i.e. the forced-GT assumption produces exactly the
mechanical result expected, no surprises).

- Truth `homalt` called (correctly, allele-wise) as `0/1` by AI: **5630 of
  5777 = 97.5%** of all homalt truth sites. This single number is the
  direct, near-total cause of the genotype-mismatch bucket above.
- Truth `het` correctly matched by AI's `0/1`: **10290 of 11120 = 92.5%.**

`N_correct_allele_wrong_GT = 5630`, `N_correct_allele_correct_GT = 10290`
(sanity: het total is 11120, so 830 het truth sites are NOT matched at all —
these fall inside `TRUE_MISSED_VARIANT`/`REPRESENTATION_ERROR`, i.e. the
caller simply didn't call them, unrelated to GT), `N_true_allele_misses = 978`.

**Counterfactual (diagnostic only, no VCF was modified):** if genotype were
ignored entirely (Section 6's allele-level figures), F1 rises from 0.6238 to
0.9651 — a difference of **+0.3413 F1 points**, i.e. essentially the entire
gap versus the internal evaluator's 0.9742 is explained by combining (a) the
denominator/BED-window difference (Section 11) and (b) the forced-GT
mismatch, with only a modest residual (Section 9) attributable to genuinely
absent calls.

## 9. Why recall is 0.609 (full population, not a sample — 6608 total FN)

| FN cause | Count | % of 6608 |
|---|---|---|
| Genotype mismatch (allele correct, GT wrong) | 5630 | 85.2% |
| Genuinely absent AI call (no call at all near this position) | 962 | 14.6% |
| Representation/local-match ("lm") | 16 | 0.2% |

**The dominant driver of the low recall number is genotype mismatch, not
missing calls.** Only 962 of 6608 "FN" (14.6%) represent the caller truly
failing to flag a variant at all; the other 85.4% are cases where a variant
WAS correctly flagged at the correct position with the correct allele, but
scored as a miss because of the zygosity mismatch.

## 10. VCF conversion audit — CONFIRMED CORRECT, no new issues found

Extended beyond `SANITY_CHECK.md`'s original ~12-locus check: a random
200-locus sample of `pb_calls==True` positions was independently
cross-checked against a **separate** `pysam.FastaFile` read (not reusing the
count-matrix's own reference channel) — **0 REF mismatches**. **0** cases of
`ALT == REF` (degenerate/no-op calls). **0** duplicate positions found in the
sample; a full-population check of all 16,094 called positions found **0
duplicates and confirmed strictly ascending sort order** (`np.diff(positions)
> 0` everywhere). Chromosome naming (`chr21`, matching
`shared_ref/chr21_chrname.fa`) is hard-coded correctly and both hap.py runs
completed without any chromosome-mismatch error. GT encoding is exactly
`0/1` for every record as designed (no encoding bugs — this is the intended,
documented behavior of the forced-heterozygous default, not an accidental
error). No multiallelic records are ever emitted (`build_vcfs.py` only ever
picks one majority ALT base per locus), which is a design simplification, not
a bug, given the frozen evidence pipeline is a single binary substitution
decision.

**Conclusion: the VCF conversion pipeline (position, REF, ALT, sort order,
uniqueness) is correct.** The forced `GT=0/1` is a known, documented,
intentional simplification — not a conversion bug — and is the single
largest contributor to the measured F1 drop (Section 8).

## 11. Denominator difference: 16597 (internal) vs 16898 (hap.py) — fully explained

See `DENOMINATOR_ANALYSIS.csv`. The raw truth VCF has 17,459 biallelic
single-base SNP records in chr21:32-44M; exactly 16,898 of them fall inside
the high-confidence BED per an independent per-base pysam intersection —
this **exactly** matches hap.py's `TRUTH.TOTAL`. The internal evaluator's
16,597 is smaller because its BED filtering happens at the **64bp-window
level** (`providers.py::_is_confident`, Section 1): a window is dropped in
its entirety if *any* base in it is outside the BED, even if the specific
truth SNP in question is itself inside the BED. Directly intersecting the
16,898 in-BED truth SNP positions against the internal evaluator's
kept-window SNP-position set found **303** such SNPs whose window was wholly
dropped despite the SNP itself being BED-confident. `16898 - 303 = 16595`,
2 short of the actual 16597 due to a duplicate-position set-collision
artifact from MNP-decomposed truth records in this diagnostic script (traced,
not left mysterious, listed under a genuine OTHER/negligible bucket rather
than forced into the main category) — net **301 explained, 0 unexplained
residual.**

## 12. PB-only vs Cascade — confirmed nearly identical

See `PB_VS_CASCADE_DIFFERENCES.csv`. Exactly 3 of 10,975,654 scored loci
differ between `pb_calls` and `cascade_calls` region-wide. At all 3, `routed`
is `False` (the frozen router did not send these to PB), so the difference is
actually `cascade_calls (== binomial_calls here)` vs the independently-scored
`pb_calls` — not a router-invoked-PB effect. TP is identical (10290 vs
10290); FN is identical (6608 vs 6608); the only difference is +3 FP for
Cascade (1 of which is a genuine additional true-positive rescue, 2 of which
are additional false positives). **Confirmed: Cascade's hap.py score is
inherited almost entirely from PB-only at this region's scale** — the
0.00006 F1 difference is noise-level, not a meaningful behavioral
difference between the two decision paths as measured here.
