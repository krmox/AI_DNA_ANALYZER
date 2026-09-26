# Limitations and Open Questions

Written 2026-09-19 while assembling `PAPER_MANUSCRIPT.md`; **updated 2026-09-20** (chr20, `load_reads`, lost artefacts, reconstructed tests). Sections A–B restate scientific limitations;
C lists unresolved questions; D lists statements in earlier project documents that this review found
superseded or wrong (they were **not edited** — original files are left as-is for the historical record);
E lists provenance/hygiene risks.

## A. Limits on what the evidence supports

1. SNP-only; no indel, MNP or SV calling or scoring (indel truth appears as 0 calls in hap.py output).
2. Evidence: one whole chromosome (HG002 chr20, 15× only, 64.4 Mb) plus regional runs (largest other span 12 Mb, largest pooled 25 Mb; largest HG005 full-pipeline run 3 Mb). No genome-scale run for any sample.
3. GIAB high-confidence regions only; truth is least reliable in segdup where the cascade differs from PB most.
4. HG002/HG003/HG004 are a related trio (one batch, one pipeline; HG003/HG004 are weak evidence of independence); one unrelated sample (HG005) on one 3-Mb region; depth cells are downsamples of one library. No population/ancestry diversity.
5. Illumina short reads, novoalign only (2×250; HG005 ≤250-bp reads from a 300× product); no HiFi/ONT; no other aligner.
6. The three frozen constants were fit at ≈15× on chr21:32–40 Mb and evaluated off-optimum elsewhere; the HG002 hap.py / external-caller region contains that fit span.
7. Statistics: locus-level bootstraps treat loci as independent and may not fully account for spatial correlation; block CIs only for v19, v20 and round-3; no CI for hap.py numbers or strata; no multiplicity correction; the pre-registered rule has no equivalence test; the v20 DEGRADED verdict rests on 10 discordant loci.
8. External callers: one region/sample/depth; GATK without BQSR; Strelka2 not evaluated; FreeBayes not run; no matched-scope AI runtime.
9. Performance: one machine (desktop load); mostly single runs; function-level (21–65×, 26.9×, 35.1×) ≠ pipeline speedups (1.22×, 2.72×, ≤ 7.1×); PB is now ≈ 93–96 % of compute, so the native optimisation improved efficiency without eliminating pipeline cost; router-gated *extraction/compute* not built.
10. No clinical validity claim of any kind.
11. **chr20 is not sample-independent:** same HG002 sample, library and 15× depth as the constant-selection data; 1 Mb (chr20:33–34 Mb) was already scored in v14 (exclusion sensitivity: verdict unchanged); 15× only.
12. **True-SNP losses are counted per validation cell:** 9 historical consolidated (v14, v15, v17, v19) + 1 on chr20 = 10 expanded.
13. **Lost artefacts:** the raw post-fix HG005 cache, region-scoped truth/BED, post-fix hap.py output, error/rescue CSVs and `experimental/performance/` reports are not in the repository (lost with a wiped `/tmp` worktree, 2026-09-20). HG005 numbers (hap.py 0.9521/0.9519, M-1 4,090/3,948/142, rescue 13/129, enrichment) stand as "reported in the preserved experiment record; raw artefact unavailable" and are flagged `RECORD_ONLY_RAW_UNAVAILABLE` in `CLAIM_EVIDENCE_MAP.csv`.
14. **Reconstructed tests** (`test_native_pileup_equivalence.py`, `test_hg005_extraction_integrity.py`, `test_counts_backend_fuzz.py`) are new, not the lost originals; the historical "102/102" count is not reproducible; `test_providers.py` does not collect (missing `bimamba_variant_caller`).
15. **Native counts backend** uses a whole-span htslib pass and differs from per-window passes on adversarial synthetic BAMs (repeated read names + inconsistent mate fields); equal on realistic synthetic BAMs and all validated real regions; code unchanged (frozen).

## B. Known, unaddressed artefact — M-1 (window-level vs SNP-level BED filtering)

Current state per the latest audit (`PIPELINE_INTEGRITY_AUDIT.md`, MEDIUM, "not fixed in this task") and confirmed by direct recomputation: **unresolved**.
* HG005 chr1:1.0–4.0 Mb: 4,090 truth SNPs in the BED; 3,948 in retained windows; **142 (3.5 %) cannot be called at all**; 142/244 (58.2 %) of FNs have no candidate. hap.py recall 0.9403 → 0.9742 conditional on retained windows; ≤ 0.9751 if all dropped SNPs were called.
* HG002 chr21:32–44 Mb: 16,898 in BED vs 16,597 in retained windows (301; 1.8 %).
* HG002 chr20 (whole): 71,387 truth SNPs in BED; 70,357 in retained windows; **1,030 (1.44 %) dropped**; 33 more excluded by an overlapping indel label → **70,324** in the locus frame; hap.py `TRUTH.TOTAL` 71,333 (Table 3b of `TABLES_FINAL.md`).
* Not quantified for v13–v19 (internal evaluator uses window-level denominators throughout).
* The HG005 report's "~35 % of all region-scoped truth SNPs" is wrong (3.5 %).
* Open: implement per-SNP/per-locus BED handling and re-measure; decide whether to rescore historical experiments.

## C. Open scientific and engineering questions

1. **How many cascade-specific true-SNP losses exist?** This review counts 9 cell-level events (v14 5 [4 distinct loci], v15 1, v17 2, v19 1). The consolidated report counts 8 (omits v15); the brief counts 7 (v14+v17). Identity and covariates of the v15 (chr16_segdup, full depth) event were not extracted `[DATA VERIFICATION NEEDED]`; whether any events coincide physically across depth cells other than v14's frame_index 604015 is unchecked.
2. **Where were the HG004 lost-SNP coordinates wrong?** Resolved for the loss (29,487,976 is correct; 29,464,488 is a naive-coordinate value). The three cascade FPs are recorded at 28,377,427 / 29,787,338 / 27,696,940 (FRR, round 3) versus 28,391,763 / 29,796,874 / 27,702,444 (`disagreements.json`); the `disagreements.json` values are probably the true coordinates (that file's lost-SNP position matched a truth SNP exactly), but no per-FP truth-context check beyond a nearby-record lookup was done. Note that 29,796,874 coincides with a *truth indel* record (CT>C), so its "false positive" label is a scoring-frame effect in a SNP-only evaluation.
3. **Mechanism of the low-VAF loss.** Interpretation (fixed-ε binomial ≈ 0 at VAF ≈ 0.12; router sees uncertainty not error) is supported by ε-shift statistics but no read-level analysis (PSV/phasing) was done. Testing whether low-VAF segdup "SNPs" in truth are real paralogous variants was proposed and not done.
4. **Is the v15 69:5 damage ratio a truth-set artefact in segdup?** Proposed in devlog 15 §21; not done.
5. **Thresholds across depths/samples.** Router and thresholds are 15×-specific; no calibration study at 30×/60× or other libraries.
6. **Does the fixed extractor change v16–v19?** Audit argues no (row-index metrics; coordinates rebuilt client-side) but the fixed `extract_bench_v12.py` was not re-run on those cells; recorded frozen-file hashes predate the fix and the native backend.
7. **Runtime of the AI pipeline at the external-caller scope** (HG002 chr21:32–44 Mb) has not been measured; the earlier 4,740 s number is pysam-era evaluation-cache construction and is not comparable.
8. **Router-gated extraction** (reading reads only for routed loci) has not been designed or timed; it is the only design in which selective escalation would avoid I/O.
9. **`load_reads` and PB optimisation**: not attempted; 48-read cap unresolved; `load_reads` run-to-run variation (696 vs 775 s) not characterised.
10. **Block bootstrap** for v13–v18 and hap.py results not computed; the v18 "IMPROVED" verdict rests largely on one 200-bp cluster.
11. **HG005 library details** (instrument/prep vs HG002–4) not verified beyond read length and aligner.
12. **Was any AI system other than Claude used** (the brief mentions Gemini)? No repository document says so `[DATA VERIFICATION NEEDED]`. The identity of the "third-party review" in the Aug 7 devlog and of the source of commit `1dfa759`'s message is not recorded.
13. **Exact genomic coordinates of the constant-selection train/validation/test blocks** (only block indices are recorded; 1-Mb blocks of chr21:32–40 Mb, index mod 4) `[DATA VERIFICATION NEEDED]`.
14. **Provenance of `FROZEN_PB_THRESHOLD = 10.5`**: recorded as validation F1-maximising in the router-freeze worktree; it numerically equals an unrelated binomial-v1 threshold in another region — no evidence of confusion found, but not independently re-derived.
15. **Unverified citations:** Zook 2016 journal record; Novoalign primary reference; pysam formal citation; Efron & Tibshirani co-author metadata; Dwarshuis 2024 article number; prior cascade-classifier literature *for variant calling* (no search done).
16. **Unmerged evidence:** round-3 (HG004 re-score, 5.34× multiprocessing, GIAB-strata enrichment), the threading/profiling/multiprocessing reports, and the original native-backend report live only in untracked/unmerged worktrees (`.claude/worktrees/agent-a5e7…`, `agent-afd5…`, `agent-a4ce…`).
17. **HG005 `false_rescue_count` (9) vs. manuscript Table 5 `router_fp` HG005 contribution (2, RO) do not reconcile.** `RESCUE_ANALYSIS.csv` defines `false_rescue_count = 9` as "routed & binomial_correct & pb_wrong" for HG005 chr1:1–4 Mb — conceptually the same event class as `router_fp` (a routed locus where the cascade, which adopts PB's call once routed, is wrong while the binomial screen was right). Table 5 of `PAPER_MANUSCRIPT_V2.md` records only 2 (RO) for HG005 in that class. Neither number is demonstrably wrong; they may use different locus-level definitions (event vs. distinct-locus counting, or a broader "any wrong call" vs. an FP-specific criterion), but no document states which, or even records that both numbers exist. `[DATA VERIFICATION NEEDED]` — raised by `FINAL_SCIENTIFIC_AUDIT.md` Finding 4. Resolving this from first principles would need the original per-locus HG005 classification re-run against the surviving cache; `NEW EXPERIMENT REQUIRED` if that cache is no longer sufficient to reproduce the distinction, since the raw post-fix HG005 artefacts behind both numbers are themselves record-only (item A.13).

## D. Statements in existing project documents that are superseded or incorrect

(Left unedited; do not quote these without this correction.)

| File | Statement | Problem |
|---|---|---|
| `research/PAPER_DRAFT.md` §8, §14, abstract | "no compiled/native extension exists in this codebase (I-2)"; "6/7 internal experiments…one mixed"; blocked-HG005 story | native code is integrated (working tree); the "mixed" v15 label belongs to the safety layer; v14 pre-registered verdict is NEGATIVE |
| `research/PIPELINE_INTEGRITY_AUDIT.md` I-2 and the preserved original H-1 text | "No compiled/native extension exists anywhere in the tracked repository"; "no file substantiates 65×" | true only of tracked files/one worktree; superseded by H-1 RESOLVED in the same file |
| `research/REPRODUCIBILITY.md` (software table, HTSlib row) | "no C extension in this repo" | native module now exists in the working tree; also the branch is not committed |
| `experimental/stress_test/HG005_POSTFIX_STRESS_TEST_REPORT.md` §9, §10, §12 | "~35 % of ALL region-scoped truth SNPs"; "No C/htslib-accelerated extractor is integrated" | 142/4,090 = 3.5 %; written before the native integration |
| `research/RESULTS.md` §1 table; `research/BENCHMARK_MASTER.csv` | v14 "IMPROVED"; v15 "WEAK POSITIVE / NEGATIVE (mixed)" attributed to the cascade; v13 row pools 3 in-span regions unmarked | see Section 3.1 of the manuscript |
| `research/PERFORMANCE_RESULTS.csv` | thread-scaling region "chr21:0-2000000, 2 Mb" | the sweep used a 300-kb region; 2.85× end-to-end multiprocessing and 5.34× round-3 figures are absent |
| `experimental/performance/PERFORMANCE_OPTIMIZATION_REPORT.md` (worktree) §8, §10, "Final decision" | "65.1× end-to-end"; "faster than DeepVariant/Clair3/GATK (48.9×/80.1×/17.1×)" | scope = frozen masks; excludes `load_reads`/PB; thread counts unnormalised — not a pipeline claim |
| `VALIDATION_ROUND_2.md` | per-cell loci table copied from v18; "2000 resamples" for the paired bootstrap; "first cascade-caused true-SNP loss"; low-BQ mechanism for the lost SNP; "HG005 untested" | JSON has different loci counts, 10,000 resamples; v14 had 5 and v17 had 2 losses; the lost SNP has BQ 37.3; HG005 later analysed |
| `FINAL_RESEARCH_REPORT.md`, `AI_DNA_ANALYZER_Scientific_Report.pdf` | HG005 "blocked"; "no hap.py/vcfeval"; HG004 loss coordinate 29,464,488 (FRR) | superseded by later work; coordinate is naive |
| `History/11_DEVLOG.md` §2 (main) | cutoff chosen as F-beta optimum on a training split | router-freeze worktree devlog 13 §7: validation 0.1 % quantile, F1-tie rule |
| `History/13_DEVLOG.md` §17 | "the cascade never lost a true SNP that PB found" | refuted by devlog 14 (5), 17 (2), 19 (1), 15 (1) |
| `History/15_DEVLOG.md` / `safety_layer.py` docstring | "36-candidate grid"; "~1/175 of PB per locus" | 50 candidates; measured 492× per-locus cost ratio |
| `History/17_DEVLOG.md` §"35 of 36 lowmap_segdup" | | JSON recount: 34 of 36 |
| `test_providers.py`, `test_bimamba.py` | importable tests | fail collection (`ModuleNotFoundError: bimamba_variant_caller`); excluded from the 102-test count |

## E. Provenance and hygiene risks

* The work described in the coordinate-fix, native-integration, and this manuscript exists **only** as uncommitted/untracked files in a git worktree under `/tmp` (branch `fix/hg005-extraction-integrity`, based on `main` a3d5761). A reboot or temp-directory cleanup could destroy it. A copy of the `research/` folder was made in the main checkout (see the manuscript's availability statement); code changes (`pileup_counts.py`, `extract_bench_v12.py`, `native/`, tests) were **not** copied or committed.
* The compiled `.so` is not portable; `native/build.sh` needs Python headers (`python3-devel`) or downloads them.
* The frozen constant `5.411872376933351` is hard-coded in several files rather than imported (audit H-2).
* No lockfile / environment file exists.
* Large caches (`cache/hg005_stress`, `cache/bench_v19`, …) are not tracked; earlier bench_v12–v18 caches were deleted for disk space (pre-approved as recomputable), which is why rescue analysis exists only for HG004 among v13–v19.
