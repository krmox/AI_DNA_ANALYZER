# chr20 holdout validation of the frozen AI_DNA_ANALYZER cascade

Status: scientific validation. No threshold, router cutoff, PB threshold, Method C parameter, truth set or benchmark protocol was changed. Nothing was tuned on chr20. Nothing is committed; the manuscript is untouched.

Provenance note. The first chr20 run was lost when the machine rebooted and wiped `/tmp` (where the uncommitted worktree lived). Everything here comes from a **complete re-run** in `.claude/worktrees/opt-chr20`. The re-run reproduced the lost run exactly (same 56,266,816 loci, 70,324 SNP, 69,291 / 69,297 calls, 2,826 rescuable loci), which is itself a determinism check.

## 1. Frozen model (B1)

| Item | Value |
|---|---|
| Git base commit | `a3d5761b9ad74a04d0344d94b4ec088c3440aff6` (branch `fix/hg005-extraction-integrity`), plus uncommitted working-tree changes, diff sha256 `e1865f4c72fc298e…` |
| Binomial threshold / router cutoff / PB threshold | 7.0 / 5.411872376933351 / 10.5 (asserted equal to the frozen values at run time) |
| Method C | eps = 0.01 fixed, GQ cap 99, 0/0 forced to 0/1 with GQ cap 5 |
| `cascade.py` sha256 | `b0ee9f4b24fe06dc…` (identical to the HG005 record) |
| `method_c_regression.py` sha256 | `b5cea1f1660eb472…` (identical to the HG005 record) |
| `native/pileup_native.c` (counts backend) | `85b08e2d8eb59ca0…` (source identical to the validated one; the `.so` was rebuilt) |
| PB / Binomial | `quality_error_model.py` `a69cb49591cbac13…`, `binomial_baseline.py` `f0408ecae0929aef…` |
| Optimised files (Part A; passed exact-equivalence first) | `read_level_pileup.py` `31bf5eb80072d609…`, `providers.py` `144ac4c276c03bdd…`, `pileup_counts.py` `53625bf4e79596a8…`, `extract_bench_v12.py` `dbc5af5b2ef4c8b2…`, `native/reads_native.c` `884a156794cabb93…` |
| Analysis scripts, hashed **before** the results were opened | `bench_v20_chr20.py` `74a0f3f249cf1a1e…`, `run_chr20_pipeline.py` `b85d64b3fb11d8d0…`, `run_happy_chr20.sh` `d94143a977952e32…`; full hashes in `experimental/chr20_validation/FROZEN_MANIFEST.json` and `ANALYSIS_FREEZE.json` |

Because the production extraction path now includes the Part A optimisations, it was gated against the original path on chr20 itself (section 11).

## 2. Dataset, defined before running (B2)

| Item | Value |
|---|---|
| Sample | HG002 / NA24385 (same sample as the constant-selection data; see limitations) |
| BAM | GIAB NIST Illumina 2x250 novoalign GRCh38, region `chr20`; native mean depth 70.96x over the high-confidence BED; downsampled with the project recipe `samtools view -s 42.2114` (seed 42) to **15x** (achieved 14.99x by `samtools depth`; 14.87x mean over the scored loci). 15x is the depth at which the cutoff was selected and at which the hap.py headline was measured, so only the chromosome changes. |
| Truth | GIAB v4.2.1 GRCh38 `HG002_..._benchmark.vcf.gz`, chr20: 85,951 records |
| Confident regions | `..._benchmark_noinconsistent.bed`, chr20 rows (10,192 intervals) |
| Reference | `data/reference/chr20_full.fa` (Ensembl r110, contig `20`); a `chr20`-renamed copy for hap.py |
| Region | **The entire chromosome**, `chr20:0-64,444,167`, in 500 kb chunks. Only 64-bp windows fully inside the BED are scored; 9 chunks (26.5-31.0 Mb, centromere) contain none. No region was chosen or filtered after seeing results. |
| Data checksums | in `FROZEN_MANIFEST.json` (`chr20_15x.bam` `0c76fd67…`, truth `527f2ece…`, BED `52d9c15a…`, reference `2d4a4c98…`) |

**Independence, stated explicitly.**
- The router cutoff was fitted on chr21:32-40 Mb. chr20 is a different chromosome, so this is **chromosome-level** separation.
- It is **not** sample-level separation: chr21 and chr20 are both HG002, same library, same depth.
- chr20 is not entirely unseen: `chr20:33-34 Mb` was already a bench_v14 cell (with frozen constants, not fitted). That is 1 Mb of 64 Mb. A sensitivity analysis excluding it gives the same result (section 5).

## 3. M-1 BED / window audit (B3)

Truth SNP bases (provider labelling: biallelic SNPs, MNPs decomposed):

| Quantity | Value |
|---|---|
| Truth SNP bases on chr20 | 73,421 |
| **N_BED** (inside the GIAB BED) | **71,387** |
| **N_retained** (inside a retained 64-bp window) | **70,357** |
| of which scored as SNP in the evaluation frame | 70,324 |
| of which label overwritten by an overlapping indel label (excluded from the frame) | 33 |
| **N_dropped** (no candidate locus possible) | **1,030** |
| **fraction_dropped** | **1.44 %** |
| drop reason | 1,026 window partly outside the BED; 4 tile remainder at a chunk end; 0 N-rich; 0 other |

Every locus-level recall below is therefore **conditional on the 70,324 evaluated SNPs**, i.e. on 98.5 % of BED-confident truth SNPs. The hap.py denominators (71,333) include the dropped SNPs. This is smaller than on HG005 (3.5 %) and comparable to HG002 chr21 (1.8 %); it is small enough to interpret the benchmark, but hap.py recall is bounded above by (71,333 - 1,030)/71,333 = 0.986 by construction.

Integrity gates on the whole population (56,266,816 loci): 0 reference-index mismatches against the FASTA; positions strictly increasing; every locus inside the BED; whole windows only.

## 4. Pipeline and routing (B4)

| Quantity | Value |
|---|---|
| Loci extracted / scored (labels 0 or SNP) | 56,266,816 / 56,242,693 |
| Truth SNP in the scored frame | 70,324 |
| Router-selected (PB-routed) loci | 77,248 = **0.137 %** of scored loci (HG005: 0.067 %; earlier HG002 cells 0.04-0.16 %) |
| PB-only calls / cascade calls (after Method C) | 69,291 / 69,297 |
| Method C 0/0 forced to 0/1 | 16 (PB-only) / 14 (cascade) |

## 5. PB-only control versus cascade (B5), project's existing statistical framework

Locus-level frozen evaluator (`evaluate_arms`, restricted to the scored frame):

| Arm | TP | FP | FN | Precision | Recall | F1 |
|---|--:|--:|--:|--:|--:|--:|
| Binomial only | 67,000 | 2,339 | 3,324 | 0.96627 | 0.95273 | 0.95945 |
| **PB only** | 68,278 | 1,013 | 2,046 | 0.98538 | 0.97091 | **0.97809** |
| **Cascade** | 68,277 | 1,020 | 2,047 | 0.98528 | 0.97089 | **0.97803** |

- ΔF1 (cascade - PB) = **-5.64e-5**; paired bootstrap (10,000 resamples) 95 % CI [-1.05e-4, -1.41e-5], two-sided p = 0.005.
- 50-kb block bootstrap (1,185 blocks, 2,000 resamples) 95 % CI [-1.08e-4, -7.5e-6].
- Disagreements: 10 loci of 56.2 M. McNemar: cascade wrong / PB right = 9, cascade right / PB wrong = 1, exact p = 0.021.
- Cascade-specific false negatives: 1. Cascade-specific false positives: 7 net (8 introduced, 1 avoided). Net: -1 TP, +7 FP, +1 FN.
- Sensitivity, excluding chr20:33-34 Mb: ΔF1 = -5.70e-5, CI [-1.06e-4, -1.42e-5]. Unchanged.
- Random-routing control at matched PB coverage (77,248 loci): ΔF1 = -0.0186, so the frozen router is doing real work; a random router of equal cost would lose 1.9 F1 points.

## 6. hap.py (B6): genotype-aware, both arms in identical conditions

Same image as all prior project runs (`hap.py:0.3.15`, digest `d63b963a…`), same truth VCF, same BED, same reference, default settings, Method C genotypes. SNP row, filter ALL:

| Arm | TRUTH.TOTAL | TP | FP | FN | Precision | Recall | F1 | FP.gt | FP.al |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| PB-only + Method C | 71,333 | 67,700 | 1,591 | 3,633 | 0.97704 | 0.94907 | **0.96285** | 537 | 158 |
| Cascade + Method C | 71,333 | 67,699 | 1,598 | 3,634 | 0.97694 | 0.94906 | **0.96280** | 537 | 162 |

ΔF1 = -5.5e-5 (no CI is available for hap.py). Allele-level view (TP + FP.gt): 68,237 versus 68,236. The hap.py conclusion is identical to the locus-level one. The old coordinate-only evaluator was not used as a final metric. hap.py F1 (0.963) is lower than the locus F1 (0.978) because its denominator includes the 1,030 window-dropped SNPs and because it scores genotypes.

## 7. Error analysis of every disagreement (B7)

Thresholds were fixed before opening results: low VAF < 0.15, low depth < 10, high depth >= 30, low BQ < 25, low MAPQ < 40; annotations from GIAB stratifications v3.1.

| Class | n | Signature |
|---|--:|---|
| `router_fp` (cascade false positive that PB rejects) | 8 (6 distinct sites + 2 adjacent pairs) | VAF 0.27-0.60 (**high**), depth 6-17 (**low**), mean BQ 11-22 (**all low BQ**), mean MAPQ 70; binomial LLR 12.5-13.6 (router margin 5.50-6.65, just **outside** the cutoff 5.41); PB LLR 9.0-10.4 (just under 10.5). One in tandem repeat, one in lowmap/segdup, six unannotated. |
| `router_avoids_fp` (cascade correct, PB false positive) | 1 | chr20:25,313,342: VAF 0.13, depth 31, BQ 35.9, in lowmap_segdup; binomial LLR 1.35, PB LLR 12.2 |
| `router_fn` (cascade loses a true SNP that PB calls) | 1 | chr20:46,110,517: VAF 0.125, depth 28, BQ 39.2, binomial LLR 0.67, PB LLR 13.6, router margin 6.33 (unrouted); annotated `alldifficult` **but not** lowmap/segdup |

Taxonomy against the existing record:
- The **previously observed lost-SNP mechanism** (low VAF 0.11-0.13, high BQ, depth 27-107, unrouted because the binomial LLR is near zero) **reproduces** in the single `router_fn` event (VAF 0.125, BQ 39.2, depth 28, margin outside the band). The segmental-duplication part of the signature does **not**: this locus is only in `alldifficult`, so by the pre-specified definition (low VAF AND lowmap_segdup AND unrouted) it does not match. One event is too few to say more.
- The **cascade-introduced FP class** (`router_fp`: fixed epsilon = 0.01 is over-optimistic for low-BQ data, PB rejects) is **the dominant departure here** (8 of 10). It matches the earlier signature (low BQ, PB LLR just under 10.5).
- The **class that made earlier cells IMPROVED** (PB false positives from low-VAF paralogous reads in segdup, which the cascade avoids: 153 events elsewhere) is **almost absent on chr20** (1 event). That is a new observation and it explains the sign of the result: on chr20 the cascade has almost nothing to gain, so its few losses are not offset.
- No new failure mode appears; no departure fits "other".
- Departure rate by depth: 7.9e-7 per locus at 1-9x, 6.2e-8 at 10-29x, 4.3e-6 at 30-47x (a single event; not interpretable).

## 8. Router rescue (B8), populations named

"Rescuable" = PB call correct where the binomial call is wrong. The two populations are reported separately.

| Quantity | chr20 (this run) | HG005 (record) |
|---|--:|--:|
| Total PB-rescuable loci | 2,826 | 142 |
| True-SNP rescuable (PB calls a real SNP the binomial misses) | 1,328 | 13 |
| False-positive avoidance (binomial FP that PB rejects) | 1,498 | 129 |
| True-SNP rescuable that the router sent to PB | 1,327 | - |
| False-positive avoidance that the router sent to PB | 1,490 | - |
| **True-SNP rescue recall** (1,327 / 1,328) | **0.9992** | - |
| **FP-avoidance capture** (1,490 / 1,498) | **0.9947** | - |
| All rescuable captured (2,817 / 2,826) | 0.9968 | 140/142 |
| PB errors where the binomial is right (PB harm) | 222, of which 221 routed (the cascade inherits them) | - |
| False-rescue rate (routed, binomial right, PB wrong / routed) | 0.29 % | 0.56 % |
| Unnecessary PB (routed, binomial already right / routed) | 94.4 % | 88.1 % |

The composition does **not** generalise from HG005: chr20 has 1,328 true-SNP opportunities (47 % of rescuable), not 13 (9 %). Rescue recall is high in both populations. The `false_rescue` figure is the rate of routed loci where PB is worse than the binomial.

## 9. Generalisation decision (B9), existing definitions

Frozen rule (`bench_v14_crosschrom.classify`, unchanged): PRESERVED = CI contains 0 and |ΔF1| < 0.001; DEGRADED = CI entirely below 0 or |ΔF1| >= 0.001; IMPROVED = CI entirely above 0; UNDERPOWERED = < 100 SNP or CI half-width > 0.01.

**Verdict: DEGRADED** (locus bootstrap CI [-1.05e-4, -1.4e-5] is entirely below zero; the block-bootstrap CI gives the same verdict; the exclusion sensitivity gives the same verdict). The magnitude is negligible: ΔF1 = -0.000056, 18x inside the 0.001 margin, 10 discordant loci out of 56 M, net -1 TP / +7 FP / +1 FN out of 70,324 SNP. The verdict is DEGRADED because the rule's first clause is a statistical one and this cell has the power (70,324 SNP, CI half-width 5e-5) to detect an effect that small. I report the letter of the rule and its magnitude; I do not relabel it PRESERVED.

1. Does the cascade preserve PB performance? In magnitude, yes (-0.006 F1 points). In the pre-registered statistical sense, no: a small loss is detectable.
2. Does the router behave similarly? Yes: 0.137 % routed, 99.7 % of rescuable loci captured, cutoff sweep, diagnostic only and not used to change anything: F1 0.97722 at -30 % to 0.97808 at +20 %, so the loss is not removed by a wider band.
3. Does the failure mechanism reproduce? Partly: low-BQ cascade FPs and the low-VAF/high-BQ lost SNP reproduce; the segdup-driven PB false positives that offset losses elsewhere are nearly absent here.
4. Is M-1 small enough to interpret the benchmark? Yes (1.44 %, and disclosed as conditioning).
5. Are results statistically powered? Yes, by a wide margin (70,324 SNP; 10 discordant loci is small in absolute terms, so the McNemar p = 0.021 rests on 10 events).
6. Any new failure mode? No.

No winner or ranking is claimed.

## 10. Performance (B10)

Measured on this machine (Intel i7-6700, 4 cores / 8 threads, with a desktop session using roughly 1-2 cores), extraction with 6 worker processes:

| Stage | Result |
|---|---|
| Whole extraction wall time | 2,640 s in-process (44.0 min); 45 min 34 s including start-up and writing; CPU 556 % |
| Native counts extraction (`load_counts`) | 194 s summed over workers |
| `load_reads` (native) | 444 s summed over workers |
| Binomial LLR | not timed separately in the extraction log (it is computed per chunk with the PB stage; its share of the total is not separately measured) |
| Poisson-binomial LLR | **15,036 s summed over workers = 95.9 % of stage compute** (inflated by 6-way contention on 4 physical cores; single-process HG005 PB rate is 0.145 ms per locus) |
| Method C genotype | 0.04 s per arm |
| VCF writing | 0.14 s per arm |
| Cascade/PB/Method C assembly script (incl. reading the 738 MB cache) | 98 s wall, peak RSS 9.9 GB |
| Locus-level scoring script (`bench_v20`, includes controls and bootstraps) | 7.8 min, peak RSS 24.2 GB |
| hap.py | about 11 s per arm |
| Peak RAM, extraction (whole process tree, sampled) | 12.1 GB (reached at the final concatenation in the parent) |

Reading these numbers:
- Extraction and `load_reads` are no longer the bottleneck at chromosome scale: together they are 4 % of stage compute. PB is 96 %.
- The router sends 0.137 % of loci to PB, but this pipeline **still computes PB for every locus** (it is needed for the PB-only control), so the routing fraction is not a CPU saving here.
- **Extrapolation, not a measurement:** scaling by scored-locus count (about 48x chr20 for the autosomes) suggests a whole-genome extraction of the order of a day and a half on this hardware with the same settings. This is a linear extrapolation from one chromosome, not a result.

## 11. Reproducibility (B11)

- Exact commands: `experimental/chr20_validation/run_extract_chr20.sh` (`extract_bench_v12.py --region 0 64444167 --contig chr20 --chunk-bp 500000 --features none --workers 6`), `run_chr20_pipeline.py`, `bench_v20_chr20.py`, `run_happy_chr20.sh`.
- Freeze records: `FROZEN_MANIFEST.json`, `ANALYSIS_FREEZE.json`. Cache (re-run) sha256 `18c1339113469846…`. The lost first run's cache had a different container hash (`c0ca994a…`); the .npz container is not bytewise comparable across runs, but locus count, SNP count, call counts and rescue counts matched exactly.
- Environment: Python 3.14.3, numpy 2.4.4, scipy 1.17.1, pysam 0.24.0, samtools/bcftools 1.23.1, gcc 15.3.1, hap.py 0.3.15 (docker).
- Result files: `results/bench_v20/chr20_results.json`, `chr20_disagreements.json`, `experimental/chr20_validation/happy/*`, `ai_*_chr20.vcf(.gz)`, `pipeline_meta_chr20.json`, `research/CHR20_RESULTS.csv`.
- **Gate of the optimised path against the original path, on chr20:** two full chunks (chr20:40.0-40.5 Mb and 46.0-46.5 Mb; 980,928 loci, the second contains the lost-SNP locus) were recomputed with native code disabled, pure-Python counts, pure-Python reads and the original BED scan. Every array (`counts`, `labels`, `positions`, `pb_llr`, `binomial_llr`, `depth`, `k_full`, `k_retained`, `n_counted`) is bytewise identical to the production cache (`experimental/chr20_validation/gate/`).

## 12. Limitations

- Same sample as the constant-selection data; chromosome-level, not sample-level, separation. Depth is a single value (15x); 30x and native-depth chr20 were not run.
- 1 Mb of chr20 was previously scored in bench_v14 (frozen).
- The conclusion rests on 10 discordant loci; the direction (cascade slightly worse) is statistically supported but the absolute effect is tiny and the taxonomy classes have single-digit counts.
- Locus-level metrics are conditional on evaluated windows (M-1, 1.44 %).
- No SV/indel scoring: the pipeline is SNP-only; hap.py INDEL rows are 0 by construction.
- The i7-6700 has a desktop workload running; absolute times carry that noise and PB time under 6 workers is contention-inflated.
