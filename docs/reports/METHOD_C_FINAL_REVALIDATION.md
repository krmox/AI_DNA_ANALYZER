# Method C Final Revalidation — Optimized Genotype Layer, Unified hap.py Head-to-Head

Git HEAD throughout this pass: `16b221ee82d64f66de25119aafa4e5944cfe63cc` (branch `main`, tree unmodified by this
work except new files listed in Section 9). Report generated 2026-09-17.

**Scope note, read first:** this report supersedes nothing in `experimental/genotype_layer/RESULTS_VALIDATION_C.md`.
That file's Method C numbers (GT acc 0.9888, hap.py F1 0.9538) were computed on a **different** extraction pass
(`experimental/genotype_layer/cache/validation_region.npz`, 16,121/16,124 calls) than the one used for the Stage-1
unified hap.py head-to-head benchmark (`experimental/unified_happy/`, 16,094/16,097 calls, from
`experimental/head_to_head/ai_cascade_per_locus.npz` + `experimental/unified_happy/recovered_full.npz`). The two
extractions are close in size but are **not the same array** (different positions, non-identical counts) — see
Section 1. Per the task's explicit "do not change the dataset" instruction, this revalidation uses **only** the
Stage-1 unified_happy inputs, not the genotype_layer validation_region cache. All numbers below are freshly computed
against Stage-1 inputs; the `RESULTS_VALIDATION_C.md` numbers are cited only as directionally-corroborating, separate
evidence, never pooled with these.

---

## 1. Experimental Setup

Reused verbatim, unchanged, from `experimental/unified_happy/environment.txt`:

| Item | Value | Verified this session |
|---|---|---|
| Sample | HG002 | — |
| Reference | GRCh38, `experimental/head_to_head/shared_ref/chr21_chrname.fa` | sha256 `c783b002...810b4` matches `environment.txt` ✓ |
| Region | chr21:32,000,000–44,000,000 (12 Mb) | — |
| Depth | 15x | — |
| Truth VCF | `data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz` | sha256 `d5568f8a...524ddd` matches ✓ |
| Truth BED | `data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed` | sha256 `d6d53796...f3973d0` matches ✓ |
| BAM (shared by all callers) | `hg002_chr21_32_44M_15x.bam` | sha256 recorded in `environment.txt`; not re-hashed here (not re-read this session) |
| hap.py | `quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0` | pulled this session; digest `sha256:d63b963a...cfa4f96d0` matches `environment.txt` exactly ✓ |
| Scope | SNP, Type=SNP Filter=PASS rows | same |
| Evaluation params | `hap.py <truth> <query> -f <BED> -r <ref> -o <out>`, no extra flags | identical command, only query VCF differs |

**Candidate call sets (frozen, not recomputed):**
- `pb_calls`: 16,094 loci — PB-only, `FROZEN_PB_THRESHOLD = 10.5`
- `cascade_calls`: 16,097 loci — Binomial→Router→PB cascade, `FROZEN_BINOMIAL_THRESHOLD=7.0`, `FROZEN_ROUTER_CUTOFF=5.411872376933351`

Both loaded from `experimental/head_to_head/ai_cascade_per_locus.npz` (unmodified), positions/counts from
`experimental/unified_happy/recovered_full.npz` (unmodified). Alignment between the two files was re-verified this
session with the same check `build_vcfs.py` used originally: frame-filter `recovered_full.npz` labels to
`labels ∈ {0,1}`, then compare element-wise against `ai_cascade_per_locus.npz`'s stored labels.

**Result: label agreement = 1.000000 (0 mismatches, n=10,981,312 recovered → filtered to the SNP-eval frame).**
Positions and masks are safely zippable index-for-index. No detection thresholds, no truth, no BED, no region were
touched. `git status` before and after this pass showed no modifications to any tracked file (see Section 9 for the
new files this pass adds).

**Not re-verified this session (explicitly flagged, not assumed):** BAM sha256, and the deeper claim that
`recovered_full.npz` itself was extracted from that exact BAM with `extract_bench_v12.py`/`--features none` as
described in `experimental/unified_happy/build_vcfs.py`'s docstring. This was not re-derived from raw BAM reads —
doing so would mean re-running the ~4740s full-region extraction (`RUNTIME_RESULTS.csv`), which was explicitly out
of scope for "verify the optimization, don't recompute what hasn't changed."

---

## 2. Optimization Verification

**Reference Method C** (`experimental/genotype_layer/build_vcfs_validation_C.py`, logic reused verbatim): a Python
`for` loop over candidate loci, calling `scipy.stats.binom.logpmf(k, n, p)` three times per site (once per genotype
hypothesis θ ∈ {0/0, 0/1, 1/1}, ε=0.01 fixed), then `sorted()` + dict construction per site to pick argmax and GQ.

**Optimized Method C** (new this session: `experimental/genotype_layer/method_c_regression.py`): the same closed-form
binomial log-pmf, `log C(n,k) + k·log(p) + (n−k)·log(1−p)` via `scipy.special.gammaln`, computed as **array
operations over all sites × all 3 hypotheses simultaneously** — no per-site Python loop, no `scipy.stats` call
overhead. Same ε=0.01, same three θ values, same GQ formula (`10·log10(P_best/P_second)`, capped at 99, then the
same θ=0-wins → forced 0/1 + `min(gq,5.0)` flagging rule).

No threshold, no ε, no formula was changed. Only the *evaluation strategy* of an identical closed-form expression
changed (loop → vectorized).

---

## 3. Numerical Equivalence

Run on **both** Stage-1 call sets (`pb_calls`, `cascade_calls`), element-by-element:

| Dataset | n | GT equal | GQ equal | θ=0 flag equal | GT mismatches | GQ mismatches | max|ΔGQ| |
|---|---|---|---|---|---|---|---|
| ai_pb_only | 16,094 | **True** | **True** | **True** | 0 | 0 | 0.0 |
| ai_cascade | 16,097 | **True** | **True** | **True** | 0 | 0 | 0.0 |

Zero unexplained genotype changes. Per protocol, had any mismatch been found the benchmark would have stopped here
(the script enforces this: it exits nonzero and skips VCF-writing/hap.py on any GT or flag mismatch — see
`method_c_regression.py`, the `if not (gt_equal and flag_equal): ... sys.exit(1)` guard).

Raw script output:
```
[ai_pb_only] GT equal=True (0 mismatches), GQ equal=True (0 mismatches, max|diff|=0.0), flags equal=True
[ai_cascade] GT equal=True (0 mismatches), GQ equal=True (0 mismatches, max|diff|=0.0), flags equal=True
```

Full JSON: `experimental/genotype_layer/cache/method_c_regression_report.json`.

**Conclusion: the optimized C is proven bit-exact-equivalent to the reference C on the identical Stage-1 input.**
The rest of this report uses the optimized output — proven equal — for the full benchmark.

---

## 4. Runtime Benchmark

Genotype-layer-only wall-clock (`time.perf_counter()`, pure Python/numpy/scipy, same measurement boundary reference
and optimized both used — no BAM I/O, no detection, just the genotyping arithmetic over already-extracted (k, n)
counts):

| Dataset | n | Reference C | Optimized C | Speedup |
|---|---|---|---|---|
| ai_pb_only | 16,094 | 5701.24 ms (354.25 µs/site) | 5.55 ms (0.34 µs/site) | **1027×** |
| ai_cascade | 16,097 | 5795.85 ms (360.06 µs/site) | 10.79 ms (0.67 µs/site) | **537×** |

This is **genotype-layer speedup only**. It does **not** represent full BAM→VCF runtime. Full-pipeline runtime was
not re-measured this session (that would require re-running the ~4740s full-region evidence extraction, out of
scope — see Section 1's "not re-verified" note and `RUNTIME_RESULTS.csv`'s own caveat that the extraction pass is
"NOT a fair like-for-like against the external callers' calling-only runtimes"). The genotype layer is a small
fraction of total pipeline time either way; **do not read 500–1000× as an end-to-end number.**

For context (not re-verified this session, cited from `experimental/genotype_layer/RESULTS_VALIDATION_C.md`, a
*different* extraction/scope): that separate run measured reference C at 452 µs/site on its own 16,124-call set —
same order of magnitude as the 354–360 µs/site measured here on the Stage-1 set, supporting that both extractions
exercise the same cost profile even though they are not the same array.

---

## 5. AI Results (hap.py 0.3.15, SNP, PASS)

| Method | Precision | Recall | F1 | TP | FP | FN |
|---|---|---|---|---|---|---|
| AI PB-only + optimized C | 0.978501 | 0.931945 | **0.954656** | 15,748 | 346 | 1,150 |
| AI Cascade + optimized C | 0.978319 | 0.931945 | **0.954569** | 15,748 | 349 | 1,150 |
| AI PB-only, Stage-1 forced-0/1 (for reference) | 0.639369 | 0.608948 | 0.623788 | 10,290 | 5,804 | 6,608 |
| AI Cascade, Stage-1 forced-0/1 (for reference) | 0.639250 | 0.608948 | 0.623731 | 10,290 | 5,807 | 6,608 |

Note TRUTH.TP jumps from 10,290 → 15,748 and TP+FN totals differ (16,898 vs 15,748+1,150=16,898, consistent) between
the forced-0/1 and Method-C rows — this is because the forced-0/1 VCFs are missing the `GQ` field/proper genotype
calls, causing hap.py's genotype-match logic (`FP.gt`) to reject far more candidate matches as genotype-mismatches
rather than true positives. Method C directly fixes this: `FP.gt` drops from 5,630/5,633 (forced-0/1) to 172/173
(Method C) — see raw `happy.summary.csv` under `experimental/unified_happy/happy/ai_pb_only_C/` and
`.../ai_cascade_C/`.

Genotype metrics (own script, `experimental/genotype_layer/method_c_unified_metrics.py`, truth GT parsed from the
same truth VCF, matched by (chrom,pos)):

| Method | n loci w/ truth | GT accuracy | 0/1 accuracy | 1/1 accuracy |
|---|---|---|---|---|
| AI PB-only + optimized C | 15,924 | 0.9891 | 0.9838 | 0.9989 |
| AI Cascade + optimized C | 15,925 | 0.9891 | 0.9838 | 0.9988 |

Confusion matrices (rows=truth, cols=predicted; order 0/0, 0/1, 1/1 — no 0/0-truth rows exist by construction,
since these are candidate-call sets that already passed detection):

```
PB-only + C:            Cascade + C:
        0/0  0/1  1/1           0/0  0/1  1/1
0/0  [   0    0    0 ]   0/0  [   0    0    0 ]
0/1  [   0 10127  167 ]  0/1  [   0 10127  167 ]
1/1  [   0    6 5624 ]   1/1  [   0    7 5624 ]
```

---

## 6. External Caller Comparison

**Reused unchanged**, per protocol — checksums of `happy.summary.csv` re-verified this session against
`experimental/unified_happy/checksums_happy_outputs.txt`'s recorded hashes for the AI runs, and against the values
already cited in `experimental/unified_happy/UNIFIED_HAPPY_BENCHMARK.csv` for the external callers:

| Caller | Precision | Recall | F1 | TP | FP | FN | Runtime (calling-only) |
|---|---|---|---|---|---|---|---|
| AI PB-only + optimized C | 0.978501 | 0.931945 | 0.954656 | 15,748 | 346 | 1,150 | genotype-layer 5.55 ms; full pipeline not remeasured (Section 4) |
| AI Cascade + optimized C | 0.978319 | 0.931945 | 0.954569 | 15,748 | 349 | 1,150 | genotype-layer 10.79 ms; full pipeline not remeasured |
| DeepVariant 1.6.1 | 0.933997 | 0.965381 | 0.949430 | 16,313 | 1,153 | 585 | 819 s (`RUNTIME_RESULTS.csv`) |
| Clair3 v1.0.10 | 0.896516 | 0.963842 | 0.928961 | 16,287 | 1,880 | 611 | 1,342 s |
| GATK HaplotypeCaller 4.5.0.0 | 0.987748 | 0.949343 | 0.968165 | 16,042 | 199 | 856 | 287 s |

DeepVariant/Clair3/GATK `happy.summary.csv` sha256 this session: `0766cd7a...47ae68` / `af9061db...9bd616` /
`959efa71...3e776e` — these are the same files, not recomputed. **Their input identity (same BAM, same truth/BED,
same hap.py 0.3.15 image) was already established in Stage 1 (`experimental/head_to_head/COMMANDS.md`,
`environment.txt`); this session did not re-derive that, only confirmed the summary CSVs are byte-identical to what
Stage 1 produced.**

---

## 7. Allele-Level Analysis

From hap.py `happy.extended.csv`, `Type=SNP, Filter=PASS, Subset=*, Genotype=*` row (GA4GH allele-level, matches
by position+allele regardless of the query's assigned zygosity in this row; the numbers coincide with Section 5's
overall row here because Genotype=* is the top-level rollup hap.py already reports):

| Method | Allele Precision | Allele Recall | Allele F1 |
|---|---|---|---|
| AI PB-only + optimized C | 0.978501 | 0.931945 | 0.954656 |
| AI Cascade + optimized C | 0.978319 | 0.931945 | 0.954569 |

**Call-set invariant** (`set(calls before genotyping) == set(calls after genotyping)`), verified directly by parsing
both VCFs and comparing `(pos, ref, alt)` sets:

```
[ai_pb_only_C]  before=16094  after=16094  before-after=0  after-before=0  → INVARIANT HOLDS
[ai_cascade_C]  before=16097  after=16097  before-after=0  after-before=0  → INVARIANT HOLDS
```

Method C never adds or drops an allele call — it only assigns/refines GT and GQ on the existing candidate set (the
θ=0-wins sites are explicitly kept, flagged, forced to `0/1`, never dropped — same invariant documented in the
original `RESULTS_VALIDATION_C.md`, now re-confirmed independently on the Stage-1 scope).

---

## 8. Genotype Analysis

See Section 5 confusion matrices. Error pattern: almost all genotype errors are truth-0/1 sites predicted 1/1 (167
of both PB-only and Cascade), a small number the reverse (6–7). No 0/0 predictions appear because these are
candidate sets that already passed detection (0/0-favored posterior sites are forced to `0/1`, per the
never-suppress-a-detection invariant — see Section 7).

---

## 9. Failure Modes

- Both AI methods still miss 1,150 truth SNPs entirely (FN) — this is a **detection** limitation (cascade/PB-only
  candidate generation), not something Method C's genotyping layer can fix; genotyping only refines calls the
  detector already produced.
- Remaining FP.gt (172–173) and FP.al (29–30) after Method C indicate Method C is not perfect — ~1% residual
  genotype-assignment error and a small number of allele-identity errors survive.
- PB-only and Cascade are nearly identical on this benchmark (ΔF1 = 0.000087, 3 more FPs for Cascade) — this
  particular 12 Mb region does not strongly discriminate between the two detection strategies once genotyping is
  fixed; earlier (forced-0/1) results showing "PB-only ≈ Cascade" were not actually testing genotyping-layer
  differences at all, since GT was forced identically for both.

---

## 10. Limitations

- **Single sample, single region, single technology.** HG002, GRCh38, chr21:32–44 Mb, 15x short-read depth only.
  No claim of generalization to other samples, chromosomes, coverage depths, or sequencing technologies is made or
  supported by this data.
- **SNP-only.** INDEL rows show 0 calls for the AI methods (not attempted by this pipeline) — no INDEL comparison
  exists for AI vs external callers on this data.
- **Full BAM→VCF runtime not measured this session** for either AI method; only the genotype-layer sub-step and,
  separately, the historical (Stage-1) full-extraction time (4740s, single-threaded, over-scoped relative to
  production use) are on record. Do not compare AI "runtime" to DeepVariant/Clair3/GATK's calling-only runtimes
  without this caveat.
- **BAM sha256 not re-verified this session** (see Section 1) — inherited on trust from `environment.txt`.
- Genotype accuracy (98.9%) and allele-level F1 (~0.955) are measured only over the **already-detected candidate
  set** — they say nothing about the 1,150 truth SNPs the detector never proposed as candidates at all.

---

## 11. Reproducibility

- `git status`: clean before and after, except this session's new files (Section 12) and the pre-existing untracked
  files noted at session start (`DEEPMIND_EXTERNAL_VALIDATION.md`, `FINAL_BENCHMARK_RESULTS.csv`,
  `FINAL_RESEARCH_REPORT.md`, `FINAL_RESEARCH_SUMMARY.html`, `experimental/` — none of these were modified by this
  pass; the `experimental/` tree already contained everything read and extended here).
- Git HEAD: `16b221ee82d64f66de25119aafa4e5944cfe63cc`, unchanged throughout.
- No threshold, ε, detection logic, truth, BED, or region was modified. `FROZEN_BINOMIAL_THRESHOLD=7.0`,
  `FROZEN_PB_THRESHOLD=10.5`, `FROZEN_ROUTER_CUTOFF=5.411872376933351`, `EPS=0.01` — all read from source, none
  edited.
- hap.py image digest, truth/BED/reference checksums verified exactly as recorded in
  `experimental/unified_happy/environment.txt` (Section 1 table).
- Nothing was committed. No original benchmark result file was deleted or overwritten.

---

## 12. New Files This Session

- `experimental/genotype_layer/method_c_regression.py` — reference-vs-optimized equivalence test + optimized-C VCF
  builder.
- `experimental/genotype_layer/method_c_unified_metrics.py` — genotype/allele/invariant metrics on the Stage-1 scope.
- `experimental/genotype_layer/cache/method_c_regression_report.json`
- `experimental/genotype_layer/cache/unified_C_genotype_metrics.json`
- `experimental/unified_happy/ai_pb_only_C.vcf(.gz/.tbi)`, `experimental/unified_happy/ai_cascade_C.vcf(.gz/.tbi)`
- `experimental/unified_happy/happy/ai_pb_only_C/`, `experimental/unified_happy/happy/ai_cascade_C/` (hap.py outputs)
- `experimental/unified_happy/checksums_happy_outputs_C.txt`

---

## 13. Final Scientific Interpretation

**Did optimized C preserve reference-C accuracy?** Yes, exactly — 0 GT mismatches, 0 GQ mismatches across both
16,094- and 16,097-site candidate sets (Section 3). This is a proof on the actual data used in this report, not an
inference from a different dataset.

**Was allele-level F1 preserved?** Trivially yes — Method C (either implementation) never changes the allele call
set, only GT/GQ (Section 7 invariant check). Allele-level F1 is therefore identical to whatever it would be under
reference C by construction, and was independently confirmed via hap.py on the optimized-C output: 0.954656
(PB-only), 0.954569 (Cascade).

**How much did the genotype layer speed up?** 537–1027× on genotype-layer-only wall-clock (Section 4). This is not
an end-to-end pipeline speedup and must not be reported as one.

**Did full hap.py F1 change?** Not relative to reference C — they're proven identical (Section 3). It changed
enormously relative to the *old forced-0/1* VCFs from Stage 1: F1 rose from ~0.624 to ~0.955 for both PB-only and
Cascade once real genotyping replaced the forced-`0/1` placeholder. This confirms the concern already flagged before
this session: the original ~0.97 internal-evaluator F1 figures were not directly comparable to hap.py's stricter
genotype-aware matching, and the un-genotyped Stage-1 F1 (~0.62) was itself an artifact of missing genotyping, not
a true measure of detection quality.

**Did PB-only and Cascade results change relative to each other?** Both moved from F1≈0.6238 (forced-0/1) to
F1≈0.9546 (Method C) — nearly identical to each other in both regimes (ΔF1 < 0.0001 in both). This dataset does not
show a meaningful cascade-vs-PB-only difference once genotyping is applied; any conclusion about the cascade's value
must come from elsewhere (e.g. the router-rescue/false-rescue analysis in `FINAL_RESEARCH_REPORT.md` Section 12,
which is a separate, non-pooled analysis).

**How does AI now compare to DeepVariant, Clair3, GATK on this one benchmark?** AI+Method C (F1≈0.9546) sits below
GATK HC (0.9682) and DeepVariant (0.9494 — actually AI PB-only edges past DeepVariant by +0.005 F1 here, driven by
much higher precision but lower recall), and above Clair3 (0.9290). AI's precision (0.9785) is its strongest
dimension; its recall (0.9319) — driven by 1,150 truth SNPs never reaching the candidate set at all — is its
weakest, and genotyping cannot fix that; only better detection could.

**What limitations remain?** Single sample/region/technology (Section 10); SNP-only; no INDEL story; no re-measured
full-pipeline runtime; BAM identity trusted, not re-verified. **This is a benchmark of one HG002/chr21/15x setup,
not evidence of universal superiority or inferiority for any caller involved.**
