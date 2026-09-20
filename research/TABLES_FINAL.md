# Final Tables

All numbers are traceable through `CLAIM_EVIDENCE_MAP.csv`. "Internal" = the row-index boolean evaluator
(`robustness_benchmark.accuracy_block`); "hap.py" = hap.py 0.3.15 xcmp (SNP, PASS). The two are **not
comparable**. Tables 1–6 correspond to the main-text tables; S1–S5 are supplementary.

---

## Table 1. Datasets and validation experiments

| ID | Date | Sample (relation) | Platform / aligner | Region(s) | Depth cells | Loci scored | Truth SNPs (scope) | Evaluator | Purpose | Notes |
|---|---|---|---|---|---|---|---|---|---|---|
| Fit | 2026-08-13 | HG002 | Illumina 2×250, novoalign | chr21:32–40 Mb (train blocks 0,4; validation 1,5; test 2,3,6,7) | ≈15× (subsample 0.203, seed 20260811) | 7,611,136 | 10,658 | internal | selection of 7.0 / 10.5 / 5.4119 | validation 1.91 Mb / 2,410 SNP |
| v13 | 08-16 | HG002 | same | 13 chr21 regions | 15×/30×/full | 24,583,412 | 33,893 | internal | robustness | 3 regions inside the fit span ("tainted" for router-generalisation) |
| v14 | 08-17 | HG002 | same | chr20:33–34, chr19:1–2, chr4:101–102, chr1:120–121 Mb | full/30×/15× | 10,368,078 | 12,261 | internal | cross-chromosome | selected by script from annotation only |
| v15 | 08-17 | HG002 | same | chr16, chr15, chr7 (segdup), chr14 (neutral) | full/30×/15× | 9,498,705 | 10,863 | internal | held-out test of safety layer | validation = v14 cells |
| v16 | 08-18 | HG002 | same | chr13:70–73 Mb | full/30×/15× | 8,726,751 | 14,274 | internal | median-difficulty unseen contig | |
| v17 | 08-19 | HG002 | same | chr17:18–21 Mb | full/30×/15× | 7,400,175 | 8,484 | internal | hardest unseen contig (17p11.2 segdup) | 14.9 % of window unscorable |
| v18 | 09-04 | HG003 (father; same trio) | same | chr13:70–73 Mb | full/30×/15× | 8,699,208 | 13,266 | internal | cross-sample | same coordinates as v16 |
| v19 | 09-14 | HG004 (mother; same trio) | same | chr2:210–213, chr3:75–78, chr5:27–30 Mb | full/30×/15× | 25,336,365 | 38,340 | internal | new sample, 3 region roles | block bootstrap added |
| v20 | 09-20 | HG002 | Illumina 2×250 novoalign, downsampled to 15× (`-s 42.2114`) | **chr20 whole chromosome** (0–64,444,167) | 15× only | 56,242,693 (56,266,816 extracted) | 70,324 (locus frame) / 71,333 (hap.py) | internal + hap.py | chromosome-scale check; same sample and depth as constant selection | 1 Mb (chr20:33–34) previously in v14; see Table 3b |
| R3 | 09-15 | HG004 | same (re-read v19 bytes) | same 9 cells | same | 25,336,365 | 38,340 | internal | independent block bootstrap; GIAB v3.1 enrichment | unmerged worktree; no new sequencing evidence |
| U-H2 | 09-17 | HG002 | same, 15× | chr21:32–44 Mb (12 Mb) | 15× | 10,975,654 | 16,898 | hap.py | genotype layer; external callers | **contains the 8 Mb fit span** |
| U-H5 | 09-19 | HG005 (unrelated) | NHGRI 300× HiSeq, novoalign 3.02.07, ≤250-bp reads | chr1:1,000,001–4,000,000 | ≈30× (`samtools view -s 42.09`) | 2,402,190 (frame) / 2,403,008 rows | 4,090 (BED-scoped) | hap.py | independent sample | region-scoped truth; 3,948 in retained windows |

Truth: GIAB/NIST v4.2.1, GRCh38, SNP-only; strata: GIAB v3.1. Reference: Ensembl r110 per-chromosome FASTA.

---

## Table 2. Cascade versus PB-only (internal evaluator, depths pooled)

| Exp. | Sample / region | PB-only F1 | Cascade F1 | ΔF1 [95 % CI, paired locus bootstrap] | Routed | Disagreements | Lost true SNPs | Pre-registered verdict |
|---|---|---|---|---|---|---|---|---|
| v13 | HG002 chr21 pool† | 0.98162 | 0.98172 | +0.000103 [0.00000, +0.00022] | 0.1084 % | 15 | 0 | PRESERVED (12/13 regions; 1 IMPROVED) |
| v14 | HG002 four 1-Mb windows | 0.955726 | 0.956306 | +0.000581 [+0.000075, +0.001089] | 0.0704 % | 41 | 5 (4 loci) | **NEGATIVE** (chr20 full: ΔF1 −0.001262 DEGRADED); pooled IMPROVED |
| v15 | HG002 held-out (chr16/15/7/14) | 0.925079 | 0.927766 | +0.002686 [+0.001968, +0.003414] | 0.0581 % | 80 | 1 | cascade vs PB IMPROVED |
| v16 | HG002 chr13 | 0.994600 | 0.994565 | −0.000035 [−0.000171, +0.000070] | 0.0392 % | 3 | 0 | PRESERVED |
| v17 | HG002 chr17 segdup | 0.948070 | 0.949962 | +0.001892 [+0.001250, +0.002615] | 0.0582 % | 38 | 2 | IMPROVED |
| v18 | HG003 chr13 | 0.993751 | 0.994013 | +0.000262 [+0.000076, +0.000483] | 0.0409 % | 7 | 0 | IMPROVED |
| v19 | HG004 chr2/3/5 | 0.986405 | 0.987419 | +0.001014 [+0.000787, +0.001254]; block [+0.000520, +0.001601] | 0.0409 % | 87 | 1 | IMPROVED |
| **v20** | HG002 chr20, whole chromosome, 15× | 0.97809 | 0.97803 | **−0.0000564** [−0.000105, −0.0000141]; block (1,185 blocks) [−0.000108, −0.0000075] | 0.137 % | 10 | 1 | **DEGRADED** (strict rule; extremely small magnitude; McNemar 9:1 p = 0.021) |
| H2H | HG002 chr21:32–44 Mb 15× | 0.9742 | 0.9741 | −0.0000282 [−0.000122, +0.0000647] | 0.155 % | — | — | PRESERVED (internal evaluator; 16,597 denominator) |
| hap.py | HG002 chr21:32–44 Mb (Method C) | 0.9547 | 0.9546 | −0.00009 (no CI) | — | 3 of 10,975,654 | — | n/a |
| hap.py | HG005 chr1:1–4 Mb (Method C) | 0.9521 | 0.9519 | −0.00024 (no CI) | 0.0666 % | 2 of 2,402,190 | 0 | n/a |
| hap.py | HG002 chr20 (Method C) | 0.96285 | 0.96280 | −0.000055 (no CI) | 0.137 % | 10 (locus frame) | 1 | n/a (TP 67,700 / 67,699; FP 1,591 / 1,598; FN 3,633 / 3,634) |

† 3 of 13 regions lie inside the constant-fit span. Untouched-region pool: PB = cascade = 0.99233 (2,306,207 loci; 2 disagreements).
v15's "WEAK POSITIVE / NEGATIVE (mixed)" applies to the safety-layer arms, not to this row.

## Table 3. External-caller comparison (hap.py; HG002 chr21:32–44 Mb, 15×, one region/sample; TRUTH.TOTAL 16,898)

| Caller | Precision | Recall | F1 | TP | FP | FN | Runtime (calling stage) |
|---|---|---|---|---|---|---|---|
| AI PB-only + Method C | 0.9785 | 0.9319 | 0.9547 | 15,748 | 346 | 1,150 | n/a |
| AI cascade + Method C | 0.9783 | 0.9319 | 0.9546 | 15,748 | 349 | 1,150 | n/a |
| DeepVariant 1.6.1 | 0.9340 | 0.9654 | 0.9494 | 16,313 | 1,153 | 585 | 819 s, 8 threads |
| Clair3 v1.0.10 (ilmn) | 0.8965 | 0.9638 | 0.9290 | 16,287 | 1,880 | 611 | 1,342 s, 4 threads |
| GATK HC 4.5.0.0 (no BQSR) | 0.9877 | 0.9493 | 0.9682 | 16,042 | 199 | 856 | 287 s |
| Strelka2 2.9.10 | blocked (bundled htslib assertion vs host zlib) | | | | | | |
| FreeBayes | not run | | | | | | |
| *Forced 0/1 (superseded)* PB-only | 0.6394 | 0.6089 | 0.6238 | 10,290 | 5,804 | 6,608 | — |

Caveats: 8 of the 12 Mb overlap the AI-constant fit span; GATK without BQSR; novoalign alignments; single runs, unnormalised threads; the AI pipeline's runtime at this scope was never measured for a like-for-like comparison. Not a ranking.

## Table 3b. Denominators — HG002 chr20 (whole chromosome), from `results/bench_v20/chr20_results.json` (`m1_audit`) and the hap.py summaries

| Quantity | Count | Meaning |
|---|--:|---|
| Truth SNP bases in GIAB BED (`N_BED`) | 71,387 | truth restricted to the confident BED (73,421 truth SNP bases on chr20 overall) |
| Excluded by retained-window filtering | 1,030 | no candidate possible (1,026 window partly outside BED; 4 chunk-end remainders) → 1.44 % |
| Retained-window truth SNPs | 70,357 | 71,387 − 1,030 |
| Excluded by indel-label overlap | 33 | outside the SNP/Normal scoring frame |
| **Locus-evaluator denominator** | **70,324** | 70,357 − 33; all locus-level metrics |
| hap.py `TRUTH.TOTAL` (SNP) | 71,333 | hap.py normalisation; includes window-dropped SNPs; 54 fewer than `N_BED` |

Other M-1 denominators: HG005 chr1:1–4 Mb 4,090 / 3,948 / 142 (3.5 %; *preserved record, raw unavailable*); HG002 chr21:32–44 Mb 16,898 (hap.py) vs 16,597 (internal), 301 (1.8 %).

---

## Table 3c. Chromosome-scale validation (HG002 chr20, 15×, frozen constants)

| Metric | PB-only | Cascade | Δ (cascade − PB) |
|---|--:|--:|--:|
| Locus F1 (70,324 SNP) | 0.97809 | 0.97803 | −0.0000564 (95 % CI [−0.000105, −0.0000141]; block CI [−0.000108, −0.0000075]) |
| hap.py F1 (TRUTH.TOTAL 71,333) | 0.96285 | 0.96280 | −0.000055 (no CI) |
| Routed to PB | 100 % (control) | 0.137 % (77,248) | — |
| Disagreements | — | 10 of 56,242,693 (8 `router_fp`, 1 `router_avoids_fp`, 1 `router_fn`) | net −1 TP, +7 FP, +1 FN |
| Verdict | — | **DEGRADED** (pre-registered rule, unchanged) | ΔF1 ≈ −0.0056 percentage points; 18× inside the 0.001 margin |

---


## Table 4. Failure taxonomy

| Class | Experiments | Events | VAF | Mean BQ | Binomial LLR | PB LLR | Notes |
|---|---|---|---|---|---|---|---|
| Lost true SNP (`router_fn`) | v14 (5 events/4 loci), v15 (1), v17 (2), v19 (1); 0 in v13 (chr21, 24.6 M loci), v16, v18, HG005; **v20 (1: chr20:46,110,517)** | **9 historical consolidated (pre-chr20); 10 expanded incl. v20** | 0.110–0.130 (v20: 0.125) | 36.7–38.8 (v20: 39.2) | −5.30…+1.57 (v20: 0.67) | 11.8–24.1 (v20: 13.6) | events counted per validation cell; v20 event is `alldifficult`, not `lowmap_segdup`; covariates for 8 of 9 pre-chr20 events |
| Cascade-introduced FP (`router_fp`) | v14 (8), v16 (2), v19 (3), HG005 (2); **v20 (8)** | 15 pre-chr20; **23 incl. v20** | 0.15–0.60 | 9.6–26.5 (v20: 11–22) | 12.5–27.0 | −5.5…+10.45 | fixed ε=0.01 over-optimistic at low BQ; dominant class on chr20 (8/10 departures) |
| PB FP avoided (`router_avoids_fp`) | v14 (26), v16 (1), v17 (36), v18 (7), v19 (83); **v20 (1)** | 153 pre-chr20; **154 incl. v20** | 0.07–0.13 | 24.8–38.9 | −21…+1.5 | ≈10.7–29.4 | paralogous low-VAF reads in segdup; almost absent on chr20 |
| Window-drop artefact (M-1) | HG005 (record only); HG002 chr21; HG002 chr20 | 142/244 FN; 301/16,898; 1,030/71,387 (1.44 %) | — | — | — | — | no candidate generated; engineering artefact |
| Region enrichment (HG005) | FN lowmap_segdup 7.73×; FP 6.35×; FP tandemrepeats 7.46×; FN tandemrepeats 0.57× (n=6) | — | median 0.0 in segdup vs 0.19 elsewhere (102 FNs with candidate) | | | | no CIs |
| Genotype residuals (Method C) | HG002 hold-out: 169 het→hom, 9 hom→het; HG005: FP.gt 17, FP.al 3 | — | hets with VAF > 0.69 | | | | |

## Table 5. Performance

**5a. Extraction stage only (`load_counts`; function-level for the counts stage, not pipeline)**

| Region | Loci | pysam | native C/htslib | Speedup | Origin |
|---|---|---|---|---|---|
| HG002 chr21:32.0–32.2 Mb | 187,200 | 16.14 s | 0.54 s | 29.9× | this work, MEASURED |
| HG005 chr1:1.0–4.0 Mb (one call) | 2,403,648 | 373.11 s | 17.74 s | 21.03× | this work, MEASURED |
| HG002 chr21:32–44 Mb | 10,981,312 | 1,092.07 s | 16.76 s | 65.14× | original worktree (2026-09-18), not re-run |
| HG005 in-pipeline (500-kb chunks) | 2,403,008 | 344 s | 15 s | 22.9× | this work, MEASURED |

**5b. End-to-end, HISTORICAL counts-only configuration (only `load_counts` native; HG005 chr1:1.0–4.0 Mb, `extract_bench_v12.py --features none`; PB computed for every locus; superseded for current code by 5d/5e)**

| Stage | pysam | native | Share of native wall time |
|---|---|---|---|
| `load_counts` | 344 s | 15 s | 1.3 % |
| `load_reads` | 696 s | 775 s | 68.0 % |
| Poisson-binomial LLR | 351 s | 349 s | 30.6 % |
| **Total wall** | **1,392 s** | **1,140 s** | **1.22× speedup** |

**5c. Other performance measurements (kept distinct; different regions, scales and loads)**

| Mechanism | Measurement | Result | Status |
|---|---|---|---|
| htslib threading (`threads=1/2/4/8`) | `load_counts` 300 kb: 15.25/16.23/15.74/16.00 s; `load_reads` 47.88/47.26/46.39/48.98 s (n=3) | ≈ no speedup | rejected |
| Chunked multiprocessing, `load_counts` | 0.6 Mb: 34.90 → 23.47 / 15.16 / 11.25 s (1.49/2.30/3.10×) | exact | experimental |
| Chunked multiprocessing, `load_reads` | 0.3 Mb, 8 workers: 43.60 → 14.25 s (3.06×) | exact | experimental |
| Chunked multiprocessing, end-to-end | 2 Mb: 993.3 → 348.3 s (2.85×; median of 2–3 noisy reps) | exact | experimental |
| Chunked multiprocessing, independent re-implementation (round 3) | 1 Mb HG004: 800.07 → 149.72 s (5.34×, 8 workers; median of 3) | decision-identical | unmerged worktree |
| Cascade vs PB-only, measured BAM→calls (PB restricted to routed loci; both arms extract everything) | HG004 chr2:210.0–210.5 Mb, n=3: 469.96 s (sd 64.5) vs 284.66 s (sd 1.5) | 1.65× | measured |
| Projected caller-stage speedup (per-locus throughput × routed counts) | 171–396× (v13–v19) | projection | **not a pipeline result** |
| Whole-genome time | ≈ 6.1 h (`load_counts` only) / ≈ 399 h (full pipeline) | linear extrapolation | **not a measurement** |
| Genotype layer vectorisation | 5,701 → 5.5 ms (1,027×); 5,796 → 10.8 ms (537×) | exact | genotype layer only |

**5d. Function-level speedups — `load_reads` (HG005 30× 0.4 Mb; HG002 15× 0.4 Mb; medians of 3 unless noted; desktop load ≈ 1–2 cores)**

| Configuration | Time | Speedup vs serial pure Python |
|---|--:|--:|
| HG005 30×, serial pure Python | 79.95 s | 1.00× |
| Python threads ×4 | 109.14 s | 0.73× |
| Python multiprocessing ×4 / ×8 | 27.95 / 18.22 s | 2.86× / 4.39× |
| **Native, serial** | 2.98 s | **26.9×** |
| HG002 15×, serial pure Python (n = 2) | 89.63 s | 1.00× |
| **HG002 15× native, serial** | 2.56 s | **35.1×** |

Async prefetch was not benchmarked as a stand-alone `load_reads` change.

**5e. Pipeline-level (HG005; PB computed for every locus) — do not compare with 5d**

| Benchmark | Configuration | Result |
|---|---|---|
| chr1:1.0–1.8 Mb, 8 × 100 kb | serial pure Python | 271.7 s (1.00×) |
| same | serial native | 99.9 s (**2.72×**) |
| same | async prefetch, Python reads / native reads | 1.21× / 1.05× over native serial |
| same | native + 4 processes (multiprocessing) | 38.2 s (**7.11×** vs pure-Python serial; ≈ 2.6× vs native serial) |
| chr1:1–4 Mb | native, 1 worker / 4 workers | **309 s / 146 s** |
| Bottleneck | PB share | 287/309 s ≈ 93 % (HG005); 95.9 % of stage compute (chr20, 6 workers, contention-inflated) |
| Amdahl | eliminating `load_reads` | ceiling ≈ 3.06×; observed 2.72× |

Equivalence: 0 mismatching cells in > 6.7×10⁸ tensor cells; 8 fuzz BAMs; 5 real chunks; ≈ 574 k BED queries (0 mismatches); byte-for-byte identical for the validated output arrays. First single-pass design: 30 mismatching cells (rejected). Table 5b above is the *counts-only* historical configuration.


## Table 6. Important scientific and engineering corrections

| # | Issue | Discovered (date) | Evidence of what was affected | Correction | Residual |
|---|---|---|---|---|---|
| C1 | Coordinates written as `chunk_start + row_index` | documented benign 08-14/16; head-to-head join fails 09-16; HG005 F1 0.0085 09-18; fixed 09-19 | one script invalidated; `counts`/`labels` byte-identical pre/post on HG004; 0/2,403,008 REF mismatches post-fix | `load_counts(return_positions=True)` at 3 call sites; 11 tests | v16–v19 not re-run with fixed extractor; HG004 locus documented at wrong coordinate in two reports |
| C2 | Internal evaluator ≠ hap.py; forced GT 0/1 | 09-17 | F1 0.6238 vs 0.9742 | Method C; evaluator column | naming ambiguity in source (M-2) |
| C3 | hap.py truth 24 Mb vs 3-Mb query | 09-19 | F1 0.2372 (diagnostic) | region-scoped truth/BED (4,090) | standing gate needed |
| C4 | Native implementation "not found" in audit | 09-19 | erroneous claim that 65× was unsubstantiated | located, audited, rebuilt, integrated | stale statements remain in 5 documents (see open-questions file) |
| F10 | First native attempt not equivalent | 09-18 | up to 282,335 / 1,684,800 cells | `ignore_orphans`, `ignore_overlaps` | contract holds for the stated filter set only |
| C5 | Empty-region `np.concatenate([])` | 09-19 | crash | empty arrays + test | — |
| PB | Depth clipping (ALT count clipped to 48) | 08-14 (found), 08-16 (fixed) | PB F1 0.684 at 69×; 1,371 SNPs recovered | recount ALT | ≤ 48 reads scored |
| M-1 | Window-level vs SNP-level BED | 09-17/19 | HG005: 142 (3.5 %) truth SNPs ineligible | none (disclosed) | open |
| D1 | v14 verdict presented as IMPROVED | 2026-09-19 (this review) | pre-registered NEGATIVE | reported both | — |
| D2 | v15 safety-layer verdict attributed to cascade | 2026-09-19 (this review) | mislabel only | corrected | source docs uncorrected |
| D3 | "Never loses a true variant" | 08-17 | 9 events | claim withdrawn | — |
| D4 | HG005 recorded as blocked / no hap.py | 09-14/16 | superseded | scoped as historical | — |
| D5 | 65.1× "end-to-end" and "17–80× faster than callers" | 09-18 | frozen-mask scope | withdrawn | — |
| C7 | Native counts backend (whole-span pass) ≠ per-window on adversarial repeated-name synthetic BAMs | 09-20 | none on real validated data | disclosed; expected-failure tests; code unchanged | backend not claimed exact in general |
| C8 | First native `load_reads` design (one pass over run of windows) | 2026-09 | 30 mismatching cells (synthetic overlapping mates, repeated names) | per-window passes in C | — |
| C9 | Lost HG005 post-fix raw artefacts; two lost tests | 09-20 | HG005 numbers "preserved record" only; 102/102 not reproducible | marked; tests reconstructed | HG005 numbers not independently reproducible |

---

## Supplementary tables

### S1. v14 per-cell results (HG002; ΔF1 = cascade − PB; CI = paired locus bootstrap)

| Cell | Depth | SNPs | PB F1 | Cascade F1 | ΔF1 | 95 % CI | Routed | Disagree | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| chr20_neutral | full | 765 | 0.983290 | 0.982028 | −0.001262 | [−0.003849, +0.001228] | 0.0053 % | 4 | DEGRADED (rule clause) |
| chr20_neutral | 30× | 765 | 0.987718 | 0.987718 | 0 | [0, 0] | 0.0325 % | 0 | PRESERVED |
| chr20_neutral | 15× | 765 | 0.984375 | 0.984375 | 0 | [0, 0] | 0.1743 % | 0 | PRESERVED |
| chr19_gcrich | full | 1,083 | 0.972997 | 0.972122 | −0.000875 | [−0.002668, +0.000856] | 0.0255 % | 4 | PRESERVED |
| chr19_gcrich | 30× | 1,083 | 0.983097 | 0.983097 | 0 | [0, 0] | 0.0599 % | 0 | PRESERVED |
| chr19_gcrich | 15× | 1,083 | 0.974860 | 0.974407 | −0.000454 | [−0.001415, 0] | 0.2184 % | 1 | PRESERVED |
| chr4_atrich | full | 1,111 | 0.996413 | 0.995966 | −0.000447 | [−0.001393, 0] | 0.0030 % | 1 | PRESERVED |
| chr4_atrich | 30× | 1,111 | 0.998203 | 0.998203 | 0 | [0, 0] | 0.0133 % | 0 | PRESERVED |
| chr4_atrich | 15× | 1,111 | 0.987799 | 0.987799 | 0 | [0, 0] | 0.0936 % | 0 | PRESERVED |
| chr1_segdup | full | 1,128 | 0.882949 | 0.888167 | +0.005218 | [+0.001418, +0.008972] | 0.0507 % | 19 | IMPROVED |
| chr1_segdup | 30× | 1,128 | 0.875817 | 0.878873 | +0.003056 | [+0.000211, +0.006140] | 0.0597 % | 12 | IMPROVED |
| chr1_segdup | 15× | 1,128 | 0.858105 | 0.858105 | 0 | [0, 0] | 0.1233 % | 0 | PRESERVED |

### S2. v19 per-cell results (HG004)

| Cell | Loci | SNPs | PB F1 | Cascade F1 | ΔF1 | 95 % CI | Routed | Disagree | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| chr2 full | 2,909,002 | 4,026 | 0.996288 | 0.996781 | +0.000493 | [+0.000121, +0.001005] | 0.0030 % | 4 | IMPROVED |
| chr2 30× | 2,909,002 | 4,026 | 0.996902 | 0.996902 | 0 | [0, 0] | 0.0143 % | 0 | PRESERVED |
| chr2 15× | 2,909,002 | 4,026 | 0.984067 | 0.984067 | 0 | [0, 0] | 0.1070 % | 0 | PRESERVED |
| chr3 full | 2,643,651 | 3,877 | 0.977284 | 0.982866 | +0.005582 | [+0.004019, +0.007273] | 0.0042 % | 45 | IMPROVED |
| chr3 30× | 2,643,651 | 3,877 | 0.976744 | 0.979716 | +0.002972 | [+0.001853, +0.004211] | 0.0146 % | 24 | IMPROVED |
| chr3 15× | 2,643,651 | 3,877 | 0.970132 | 0.970256 | +0.000124 | [0, +0.000382] | 0.1032 % | 1 | PRESERVED |
| chr5 full | 2,892,802 | 4,877 | 0.994686 | 0.995194 | +0.000507 | [−0.000106, +0.001152] | 0.0036 % | 11 | PRESERVED |
| chr5 30× | 2,892,802 | 4,877 | 0.995495 | 0.995597 | +0.000102 | [0, +0.000311] | 0.0142 % | 1 | PRESERVED |
| chr5 15× | 2,892,802 | 4,877 | 0.982241 | 0.982139 | −0.000102 | [−0.000313, 0] | 0.1037 % | 1 | PRESERVED |
| **pooled** | 25,336,365 | 38,340 | 0.986405 | 0.987419 | +0.001014 | [+0.000787, +0.001254] | 0.0409 % | 87 | IMPROVED |

(The per-cell locus counts printed in `VALIDATION_ROUND_2.md` were copied from v18; these values are from `results/bench_v19/benchmark_results.json`.)

### S3. Genotype layer

| Set | Method | Calls | GT acc. | hap.py P / R / F1 | Errors |
|---|---|---|---|---|---|
| Dev chr21:30.0–30.47 Mb (555 calls) | A | | 0.2914 | 0.2883 / 0.2893 / 0.2888 | |
| | B (τ = 0.69 fit) | | 1.0000 | 0.9892 / 0.9928 / 0.9910 | 0 |
| | C | | 1.0000 | 0.9892 / 0.9928 / 0.9910 | 0 |
| Hold-out chr21:32–44 Mb (16,124 calls) | A | | 0.6464 | 0.6382 / 0.6089 / 0.6232 | 5,637 |
| | B | | 0.9250 | 0.9135 / 0.8716 / 0.8921 | 1,196 |
| | C | | 0.9888 | 0.9767 / 0.9319 / 0.9538 | 178 |
| Full unified scope (16,094 / 16,097 calls) | forced 0/1 | | — | 0.6394 / 0.6089 / 0.6238 (PB); 0.6393 / 0.6089 / 0.6237 (cascade) | |
| | C (vectorised) | | 0.989 (approx.) | 0.9785 / 0.9319 / 0.9547 (PB); 0.9783 / 0.9319 / 0.9546 (cascade) | |

Method C stratified GT accuracy (hold-out): het VAF > 0.69: 0.8586 (B 0.0000); depth < 10: 0.9642; depth 10–19: 0.9928; segdup/lowmap: 0.9856; tandem repeat: 0.9766.

### S4. Rescue analysis

| Dataset | Rescuable | Rescued | Missed | Routed | Unnecessary-PB | False rescue | Note |
|---|---|---|---|---|---|---|---|
| HG005 chr1:1–4 Mb | 142 (13 true-SNP rescues + 129 FP avoidances) | 140 (13 + 127) | 2 (both FP avoidance) | 1,600 | 1,410 (88.1 %) | 9 (0.56 %) | recomputed from cache |
| HG002 chr21:32–44 Mb | 714 | 712 | 2 | 17,007 | 15,955 (93.8 %) | — | internal index |
| HG004 (9 cells) | 504 | 500 | 4 | 10,517 | 90.5 % | — | earlier caches removed |

### S5. Router-boundary cutoff sweep (chr1_segdup, full depth, diagnostic only; cutoff never changed)

| Cutoff | Routed | F1 | ΔF1 vs PB | Departures | Projected caller-stage speedup |
|---|---|---|---|---|---|
| 2.5 | 0.0102 % | 0.887959 | +0.005010 | 23 | 322× |
| 4.0 | 0.0261 % | 0.887651 | | 20 | 306× |
| 5.0 | 0.0413 % | 0.888167 | | 19 | 292× |
| **5.412 (frozen)** | **0.0507 %** | **0.888167** | **+0.005218** | **19** | **284×** |
| 6.0 | 0.0516 % | 0.888167 | | | |
| 7.0 | 8.6716 % | 0.888683 | | 16 | 11× |
| 10 | 12.5987 % | 0.886627 | | 11 | 8× |
| 20 | 21.9849 % | 0.884686 | | 4 | 5× |
