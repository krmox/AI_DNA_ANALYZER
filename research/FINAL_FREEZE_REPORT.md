# Final Freeze Report

Date: 2026-09-20. Scope: closure of the MUST FIX / SHOULD FIX items of `FINAL_FORENSIC_AUDIT.md`. No new experiments were run; no commits were made.
Backup of all pre-edit research files, the quarantined outputs and the pre-edit scripts: `/home/mark/Documents/Projects/AI_DNA_ANALYZER_backup_pre_freeze_20260920/` (outside the repository, not in `/tmp`).

## 1. Changes made

**Manuscript (`PAPER_MANUSCRIPT.md`)**
- Added the HG002 whole-chromosome chr20 validation (v20): Table 1/2 rows, §3.1 results block, Table 3b denominators, failure taxonomy (8 / 1 / 1), rescue composition, sensitivity, limitations.
- Replaced the unattainable global "CONFIRMED" rule with a per-cell reporting rule and an explicit summary statement (predominantly PRESERVED/IMPROVED; DEGRADED in the v14 chr20:33–34 Mb cell and in v20; v20 magnitude extremely small).
- Fixed the true-SNP loss count definition: **9 historical consolidated (pre-chr20) / 10 expanded including chr20** (raw records: v14 5, v17 2, v19 1 in `results/bench_v*/disagreements.json`; v15 1 in devlog 15 §12.3; chr20 1 in `chr20_disagreements.json`).
- Added §3.8.2–3.8.4: native `load_reads` (function-level 26.9× / 35.1×), pipeline-level results kept apart (2.72× serial; 7.11× native + 4 processes vs pure-Python serial; 309 s / 146 s HG005 3 Mb), PB now ≈ 93–96 % of compute, Amdahl ceiling ≈ 3.06× vs 2.72× observed, threads 0.73×, async only pipeline-level (1.21× / 1.05×), and the first single-pass design failure (30 mismatching cells) with its per-window fix.
- Removed/qualified "102/102", "12-category suite", "11 regression tests", "bit-exact"; equivalence is stated as "identical on the listed validated inputs".
- Updated Abstract, Conclusion, §1 contributions, §3.4, §3.9 (F5 disambiguation: v14 chr20_neutral vs v20), §4.1/4.2/4.5/4.6/4.8/4.9, §5 Limitations (all ten requested items), figure legends, Data and code availability, Table 6 (C7–C9).
- HG005 statements that rest on lost artefacts now say "reported in the preserved experiment record; raw post-fix cache/report is not currently present in the repository".

**Research package**: `TABLES_FINAL.md` (Tables 1, 2, 3b, 3c, 4, 5b–5e, 6), `CLAIM_EVIDENCE_MAP.csv` (107 rows; 18 rows flagged `RECORD_ONLY_RAW_UNAVAILABLE`, 17 new rows for chr20 / load_reads / tests / verdicts), `LIMITATIONS_AND_OPEN_QUESTIONS.md`, `RESEARCH_TIMELINE.md` (Phase 6), `FIGURE_PLAN.md` (no figure regenerated; status noted), and "SUPERSEDED IN PART" banners on eleven older documents.

## 2. Scientific results preserved unchanged
No threshold, frozen constant, router rule, PB algorithm, Method C, verdict rule, raw result file, CI, error count, VCF/BAM/BED, hap.py output or external-caller number was changed. SHA-256 re-checked after all edits: `cascade.py` `b0ee9f4b24fe06dc…`, `method_c_regression.py` `b5cea1f1660eb472…`, `native/pileup_native.c` `85b08e2d8eb59ca0…`, `native/reads_native.c` `884a156794cabb93…` (all equal to the freeze records). Every chr20 number in the manuscript was re-read from `results/bench_v20/*.json` and the hap.py `happy.summary.csv` files; denominators checked arithmetically (71,387 − 1,030 = 70,357; 70,357 − 33 = 70,324; hap.py 71,333). The lost SNP is `chr20:46,110,517` (1-based; the JSON stores 0-based 46110516). **No scientific number changed.**

## 3. Manuscript claims corrected
See the list in §1; principal corrections: chr20 added (verdict DEGRADED kept, magnitude −0.0056 percentage points stated alongside); "no chromosome-scale run / largest span 12 Mb" removed; "`load_reads` no speedup", "98.7 %", "1.22× end-to-end" re-scoped as the historical counts-only configuration; "seven experiments / two hap.py" → eight / three; routed-fraction range widened to include 0.137 %; lost-SNP signature widened (BQ 39.2; `alldifficult` without segdup); rescue recall split into true-SNP and FP-avoidance populations.

## 4. Reconstructed tests
`test_native_pileup_equivalence.py`, `test_hg005_extraction_integrity.py`, `test_counts_backend_fuzz.py` were written from scratch; they are **reconstructions, not the lost originals**, and no expected output was copied from historical results (the original HG005 cache does not exist). Coverage: positions, counts, labels, ordering, overlapping mates, duplicate read names, orphans, duplicate/QC-fail/secondary/supplementary/low-MAPQ flags, BQ filtering, BED boundaries, empty/partial windows, dropped-window coordinates (the C1 bug), reference_index vs FASTA, truth-SNP REF vs FASTA, determinism, native vs pysam, and whole-span-vs-per-window native passes (the analogue of the failed first `load_reads` design).

**New finding while writing them (disclosed, not hidden):** the native *counts* backend equals pysam on every validated real region and on synthetic BAMs with unique pair names and consistent mate fields (20/20 pass), but differs from per-window passes on adversarial synthetic BAMs (12-name pool reused across hundreds of reads, inconsistent mate fields; thousands of cells differ). These cases are non-strict expected failures (xfail, with reason). The production code was **not** changed (the frozen chr20 extraction used the whole-span design); the manuscript scopes all equivalence claims accordingly (§2.7, §3.7, §4.5, §5).

## 5. Reproducibility fixes
- Hard-coded `/home/mark/…`, `.claude/worktrees/opt-chr20` and `/tmp/head_to_head_work` removed from 21 Python and 6 shell scripts (project root now derived from the script location). Path-only edits of the four frozen chr20 scripts are recorded, with old/new hashes and the diff, in `experimental/chr20_validation/PATH_ONLY_EDIT_2026-09-20.md`. Documented exceptions: `experimental/genotype_layer/method_c_regression.py` (frozen model file; its hash is in the freeze record) and a generated `.dict` header line in `experimental/head_to_head/shared_ref/`.
- 79 tracked `.pyc` removed from the git index (`git rm --cached`; files left on disk); `.gitignore` now ignores `__pycache__/`, `*.pyc`, `native/*.so`, `.pytest_cache/`.
- Invalidated pre-fix HG005 outputs (F1 0.0085) moved, unmodified, to `experimental/stress_test/INVALIDATED_prefix_coordinate_bug/` with a README.
- `CLAIM_EVIDENCE_MAP.csv`: every remaining evidence path exists (except this report's own path, created now); missing sources are labelled instead of silently cited.

## 6. Remaining limitations
Raw post-fix HG005 artefacts (cache, region-scoped truth, post-fix hap.py, error/rescue CSVs, `experimental/performance/` reports) are unavailable — HG005 numbers are "preserved record", not independently reproducible; chr20 is same-sample / same-depth with 1 Mb previously scored; one chromosome only; SNP-only, GRCh38, GIAB; related trio; locus bootstrap ignores spatial correlation (block CIs exist for v19/v20); the counts-backend adversarial-input limitation; all work is uncommitted; no lockfile; `test_providers.py` / `test_bimamba.py` do not collect (missing `bimamba_variant_caller`, pre-existing); figures 2, 3, 4, 6 were not regenerated (legends state what they omit).

## 7. Test results (full suite, `python -m pytest --continue-on-collection-errors`, 2026-09-20)
```
PASS  = 490
FAIL  = 0
SKIP  = 1
XFAIL = 21   (documented adversarial repeated-name cases)
XPASS = 5    (same category, non-strict)
ERROR = 2    (test_bimamba.py, test_providers.py: ModuleNotFoundError: bimamba_variant_caller — pre-existing, unrelated)
```
The four native-backend files alone: 63 passed, 21 xfailed, 5 xpassed (includes the 17 `load_reads` tests).

## 8. Git cleanliness
`git diff --check` clean. 0 tracked `.pyc`. Tracked changes: `.gitignore`, path-only edits in `experimental/{genotype_layer,head_to_head,unified_happy}`, and the pre-existing uncommitted production changes (`extract_bench_v12.py`, `pileup_counts.py`, `providers.py`, `read_level_pileup.py`). Everything else is untracked (native/, chr20, tests, research/). No commit was made. No `/tmp`, `opt-chr20` or `head_to_head_work` dependency remains in scripts.

## 9. Final status
**READY FOR MANUSCRIPT FREEZE** — with the following stated exceptions, none of which affects a reported number: raw HG005 post-fix artefacts are unavailable (labelled), the counts-backend adversarial-input limitation is disclosed (not fixed), and the work remains uncommitted.
