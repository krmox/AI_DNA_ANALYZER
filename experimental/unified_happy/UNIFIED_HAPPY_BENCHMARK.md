# Unified hap.py 0.3.15 Benchmark — chr21:32,000,000-44,000,000, HG002, 15x

All 5 rows below are scored by the **identical** hap.py 0.3.15 (xcmp engine,
digest `sha256:d63b963a6cb01b4830393b22369e7b91d298e4156dde353739e74e4cfa4f96d0`)
against the **identical** truth VCF, high-confidence BED, contig-renamed
reference, and region. Type=SNP, Filter=PASS rows only.

| Caller | Evaluator | Precision | Recall | F1 | TP | FP | FN | Runtime |
|---|---|---|---|---|---|---|---|---|
| AI PB-only | hap.py 0.3.15 (xcmp) | 0.6394 | 0.6089 | 0.6238 | 10290 | 5804 | 6608 | see note |
| AI Cascade | hap.py 0.3.15 (xcmp) | 0.6393 | 0.6089 | 0.6237 | 10290 | 5807 | 6608 | see note |
| DeepVariant 1.6.1 | hap.py 0.3.15 (xcmp) | 0.9340 | 0.9654 | 0.9494 | 16313 | 1153 | 585 | 819s (calling only) |
| Clair3 v1.0.10 | hap.py 0.3.15 (xcmp) | 0.8965 | 0.9638 | 0.9290 | 16287 | 1880 | 611 | 1342s (calling only) |
| GATK HaplotypeCaller 4.5.0.0 | hap.py 0.3.15 (xcmp) | 0.9877 | 0.9493 | 0.9682 | 16042 | 199 | 856 | 287s (calling only) |

TRUTH.TOTAL = 16898 for every row (same truth set, same evaluator, now truly
apples-to-apples across all 5 callers).

## Cascade vs PB-only deltas (under the unified evaluator)

- ΔF1 (Cascade − PB) = 0.623731 − 0.623788 = **−0.000057**
- ΔPrecision (Cascade − PB) = 0.639250 − 0.639369 = **−0.000119**
- ΔRecall (Cascade − PB) = 0.608948 − 0.608948 = **0.000000**

Cascade and PB-only are statistically indistinguishable at this scale: they
differ at only **3 of 10,975,654** scored loci region-wide (see
`DISAGREEMENTS.csv`), all 3 being loci the router sent to PB where PB called
a variant and binomial did not — 1 of the 3 is a true SNP (a rescue), 2 are
not (both are additional false positives introduced by the router relative
to PB-only, since PB-only made the same "no call" as binomial at those 2
loci — wait: PB-only's own call at these 3 loci was `False` per the table,
meaning cascade's PB-routed call differs from PB-only's own full-region call;
this reflects the router only invoking PB scoring where its own condition
triggers, not that PB and cascade compute PB differently — see Problems
section in the final report for the exact mechanism).

## Cross-check against previous head-to-head numbers

`experimental/head_to_head/HEAD_TO_HEAD_RESULTS.csv`'s previously-reported
AI-Cascade/PB-only F1 (computed via the frozen internal `accuracy_block()`
evaluator) were **substantially higher** than what hap.py measures here
(frozen-evaluator F1 was in the ~0.95 range per
`experimental/head_to_head/ai_cascade_region_result.json`'s `pb_only`/`cascade`
blocks, vs. ~0.624 here). This is a large, real discrepancy — not a rounding
difference — and is discussed at length in the final report's Problems and
Invalidated-Previous-Results sections.
