# Rewrite pass 1 — research map, source conflicts, claim audit, manuscript design

Status: editorial audit for the manuscript rewrite. No experiment, threshold, verdict or raw file was touched. `PAPER_MANUSCRIPT.md` is unmodified; the rewrite is `PAPER_MANUSCRIPT_V2.md`, the hostile review is `REWRITE_03_PEER_REVIEW.md`.

Evidence hierarchy applied: (1) `FINAL_FORENSIC_AUDIT.md`, (2) `CLAIM_EVIDENCE_MAP.csv`, (3) final tables/reports, (4) `PAPER_MANUSCRIPT.md`, (5) historical reports.

What I recomputed myself in this pass (everything else is carried from the sources and labelled by their status):

- chr20 hap.py P/R/F1 and TP/FP/FN for both arms, from `experimental/chr20_validation/happy/*/happy.summary.csv`: match the manuscript (PB 0.962851, cascade 0.962796, ΔF1 −5.5×10⁻⁵).
- chr20 locus-level arithmetic from `CHR20_RESULTS.csv`/report: TP+FN = 70,324 in both arms; ΔFP = +7; routed 77,248/56,242,693 = 0.13735 %; rescue 1,328+1,498 = 2,826, 1,327+1,490 = 2,817 (0.99682); McNemar 9:1 exact p = 2·11/1024 = 0.0215.
- Performance ratios: 271.7/99.9 = 2.72; 271.7/38.2 = 7.11; 99.9/38.2 = 2.61; 309/146 = 2.12; 1392/1140 = 1.22; 287/309 = 0.929; 183/272 → 1/(1−p) = 3.06.
- Repository state: `cache/hg005_stress_postfix` and `cache/hg005_stress_native` do not exist. `cascade.py` sha256 prefix `b0ee9f4b24fe06dc` unchanged. Commit `79155df` (2026-09-20 17:20) now tracks `native/`, `research/`, `results/bench_v20/`, `experimental/chr20_validation/` and the tests; it is **not pushed** (`origin/main` = `16b221e`).

---

## 1. Research map

| Item | What the record shows | Status |
|---|---|---|
| Research question | Can a cheap fixed-ε binomial screen plus a fixed-cutoff confidence router, escalating only ambiguous loci to a per-read Poisson-binomial (PB) model, keep the accuracy of PB-everywhere for short-read germline SNP calling? | clear |
| Hypothesis (as pre-specified) | Accuracy is not lost: pooled ΔF1 = F1(cascade) − F1(PB-only) is not below 0 by the rule below. Compute saving was a motivation, not a tested hypothesis. | clear |
| Baseline | PB-only (PB evaluated at every locus, threshold 10.5); binomial-only reported as a lower reference. | clear |
| Cascade | (1) binomial LLR, ε = 0.01, call if ≥ 7.0; (2) route locus to PB iff \|LLR_bin − 7.0\| ≤ 5.411872376933351 (band [1.588, 12.412]); (3) routed loci: call if PB LLR ≥ 10.5. Genotype layer (Method C) assigns GT/GQ afterwards without changing the allele set. | clear |
| Frozen parameters | 7.0 / 10.5 / 5.411872376933351; Method C ε = 0.01; PB read cap 48; MAPQ ≥ 20, BQ ≥ 13. Selected on HG002 chr21:32–40 Mb validation blocks (1.91 Mb, 2,410 SNPs) at ≈15×. `cascade.py` sha256 `b0ee9f4b…`. | clear; cutoff derivation differs between main-branch devlog 11 and router-freeze worktree (worktree used) |
| Datasets | HG002/3/4 Illumina 2×250 novoalign GRCh38 (related trio); HG005 300× HiSeq novoalign, downsampled ≈30× (unrelated). Depth cells are downsamples of one library. | clear |
| Truth | GIAB v4.2.1 GRCh38, high-confidence BEDs; SNP-only. Strata: GIAB v3.1. | clear |
| Evaluation protocol | Internal locus evaluator (v13–v20) and hap.py 0.3.15 xcmp (U-H2, U-H5, v20). Paired bootstrap on ΔF1 (10,000 resamples); 50-kb block bootstrap for v19, v20, round 3; McNemar from v19. No CI for hap.py. | clear |
| Primary endpoint | Verdict per cell/experiment from `bench_v14_crosschrom.classify` (PRESERVED / IMPROVED / DEGRADED / UNDERPOWERED). | clear |
| Secondary endpoints | Routed fraction; disagreement taxonomy; rescue capture; genotype accuracy; hap.py P/R/F1; external callers; stage/pipeline timing. | clear |
| Validation experiments | v13–v20 (8 internal-evaluator experiments), U-H2, U-H5, v20 hap.py = 3 hap.py evaluations. | clear |
| Negative results | v14 experiment-level NEGATIVE; v20 DEGRADED; 10 true-SNP losses (per cell); safety layer rejected (fixed 5, broke 69); Method B genotype failed hold-out; neural stages null/negative; first native counts attempt failed (≤282,335/1,684,800 cells); first native `load_reads` design failed (30 cells); PB depth-clipping defect; coordinate bug (HG005 F1 0.0085, invalidated). | clear |
| Failure modes | Low-VAF (0.11–0.13), high-BQ true SNPs at binomial LLR just outside the band (lost); low-BQ high-VAF binomial over-calls (cascade-introduced FP); PB false positives in paralogous sequence (avoided by the cascade). | clear (mechanism is interpretation, not read-level tested) |
| Performance | Function-level, extraction-stage, pipeline-level and projected numbers all exist; kept apart in the rewrite. | clear |
| External comparison | One region (HG002 chr21:32–44 Mb, 15×) overlapping the constant-selection span; DeepVariant 1.6.1, Clair3 1.0.10, GATK 4.5.0.0 (no BQSR). Strelka2 blocked, FreeBayes not run. | clear |
| Limitations | See manuscript §5. | clear |
| Reproducibility | Hashes, image digests, seeds, commands recorded; chr20 re-run reproduced the lost run exactly; chr20 optimised path gated bytewise on two chunks. Raw post-fix HG005 artefacts lost; two original tests lost; no lockfile. | clear |

Unclear, left unclear in the rewrite:

- Distinct-locus count of true-SNP losses across cells (only v14's 5 events = 4 loci is stated; whether v15/v17/v19/v20 events coincide with any other locus is unchecked).
- Covariates of the v15 lost SNP (not extracted).
- Whether the pre-specified verdict rule was written before or after the v13 numbers were seen (it is in devlog 13 with them).
- HG005 library/instrument differences beyond read length and aligner.
- Exact genomic coordinates of the constant-selection train/validation/test blocks.
- Whether any AI system other than Claude was used (the brief once mentioned another; no repository document names one).

---

## 2. Source conflicts found (not silently resolved)

| # | Conflict | Authoritative resolution | Consequence for the text |
|---|---|---|---|
| S1 | `CLAIM_EVIDENCE_MAP.csv` rows 43, 45, 46, 48 (routed 1,600; rescue composition 13/129; M-1 4,090/3,948/142; recall 0.9742/0.9751) and row 63 (native vs pysam cache identity) are `VERIFIED_RAW` and cite `cache/hg005_stress_postfix/…npz`. Rows 38, 42, 44 cite the same cache as `RECORD_ONLY_RAW_UNAVAILABLE`. The audit (CF-2) and my directory check say the cache does not exist. | Forensic audit outranks the map. The `VERIFIED_RAW` flag on these rows can only refer to a check made before the 2026-09-20 wipe. | All of these are RECORD-ONLY in the rewrite. Manuscript §3.7 says the HG005 cache identity was "re-verified for this paper"; that cannot be true of a cache that is absent, so the rewrite drops it. **Map rows 43, 45, 46, 48, 63 re-labelled `RECORD_ONLY_RAW_UNAVAILABLE` (values unchanged).** |
| S2 | Manuscript Data-availability says everything after `a3d5761` is uncommitted and no public remote is configured. Repository: commit `79155df` tracks it all; `origin` = github.com/Kramkost/AI_DNA_ANALYZER; `origin/main` is at `16b221e`, so nothing after `16b221e` is pushed. | Repository state. | Rewrite states commit `79155df`, local, not pushed. Public visibility of the remote is not established. |
| S3 | Earlier consolidated loss count: manuscript says "7 (v14+v17)", limitations file says consolidated report = 8 (omits v15) and the brief = 7. | Claim map row 56 (8 raw records in v14/v17/v19 + v15 in devlog). | Rewrite: 9 = 5+1+2+1, of which 8 have raw records and covariates. |
| S4 | v14: 41 departures; class counts in the manuscript give 8 + 26 + 5 = 39. | **Resolved from raw data** (`results/bench_v14/disagreements.json`, all 12 cells, none truncated): 26 `router_avoids_fp`, 8 `router_fp`, 5 `router_fn`, 2 `router_recovers_fn` = 41. | Rewrite reports the four-class split. The source documents never define `router_recovers_fn`; the rewrite says so. |
| S5 | Manuscript/map: the U-H2 region "contains the 8 Mb used to fit constants". Method record: only the validation blocks (1.91 Mb) were used to select thresholds; train blocks trained the (rejected) logistic router; test blocks were scored once. | Method record. | Rewrite: the region contains the 8-Mb selection span; we treat the whole span as non-independent, stricter than "validation blocks only". |
| S6 | PB/binomial per-locus cost ratio: 240–490× (manuscript §4.2), 260–520× (audit, from 2.6–3.1×10⁶ vs 6–10×10³ loci/s), 492× (measured, limitations file §D). | Not resolvable. | Rewrite gives the two throughputs and says "several hundred-fold". |
| S7 | 65.14× (map) vs 1,092.07/16.76 = 65.16×; 26.9× vs 79.95/2.98 = 26.8×; 35.1× vs 89.63/2.56 = 35.0×. | Rounding of unrounded stage times. | Rewrite quotes ratios as recorded and tags them "ratio of medians of unrounded times". A reader recomputing from the printed seconds will differ in the last digit. |
| S8 | Forensic audit (14:52) predates the final manuscript (15:53) and freeze report. Audit MUST-FIX 1, 2, 3, 5, 6, 7, 8 are addressed in the manuscript; MUST-FIX 4 (recover/recreate HG005 raw) is not. | Freeze report §9 accepts it as a labelled exception. | Rewrite keeps HG005 as RECORD-ONLY throughout. |
| S9 | `RESULTS.md`, `RESEARCH_SUMMARY_1PAGE.md`, `LIMITATIONS.md`, `SCIENTIFIC_CLAIMS.md`, `PAPER_DRAFT.md`: "preserves PB-only accuracy in every tested experiment", rescue recall 0.986 unqualified, "generalizes across trio members", ~65×. | Audit C64, C71–C73. | Not used as sources. |

---

## 3. Claim audit

Labels: **S** supported · **SQ** supported with qualification · **RO** record-only (raw artefact unavailable) · **OD** outdated · **U** unsupported · **X** contradicted.

| ID | Claim (as it circulates) | Label | Basis | Defensible wording used in the rewrite |
|---|---|---|---|---|
| K1 | Constants 7.0/10.5/5.4119 never refit; `cascade.py` unchanged | S | sha256 recomputed by audit and by me | as stated |
| K2 | Constants were frozen before all later validation | SQ | v13 has 3 regions inside 32–40 Mb; cutoff derivation conflict (S5-type) | "frozen before v14–v20; v13 partly overlaps the selection span" |
| K3 | v13–v19 pooled ΔF1/CI/verdicts (Table 2) | SQ | claim map `VERIFIED_REPORT`; audit did not recompute bootstraps | cite as recorded results; say bootstraps were not recomputed in the freeze audit |
| K4 | v14 pooled IMPROVED | SQ | true for pooled only | report NEGATIVE (experiment) and pooled IMPROVED together |
| K5 | v15 "WEAK POSITIVE/NEGATIVE" | SQ | belongs to the safety-layer arms | not attributed to cascade |
| K6 | chr20 verdict | S / X | DEGRADED under the rule; "no significant degradation" or "PRESERVED" is X | "DEGRADED under the pre-specified rule; ΔF1 −5.6×10⁻⁵" |
| K7 | chr20 hap.py ΔF1 −5.5×10⁻⁵ | S | recomputed | as stated, no CI |
| K8 | chr20 as independent validation | U | same sample, library, depth; 1 Mb in v14 | "chromosome-level separation from the selection region; not sample-independent" |
| K9 | True-SNP losses: 9 historical / 10 with chr20 | SQ | counted per validation cell | both numbers, definition, no distinct-locus claim |
| K10 | "Never loses a true SNP", "zero losses" | X | v14, v15, v17, v19, v20 | not used |
| K11 | Lost-SNP signature is low-VAF, high-BQ, segdup | SQ | chr20 event is `alldifficult`, not `lowmap_segdup` | "segdup in most annotated events, not all" |
| K12 | v14 `chr20_neutral` = v20 chr20 | X | 1-Mb cell vs whole chromosome containing it | always "v14 chr20:33–34 Mb cell" |
| K13 | "Cascade preserves PB-only accuracy in every experiment" | X | v14 NEGATIVE, v20 DEGRADED | not used |
| K14 | "Accuracy preserved" without scope | U | | "within small margins, in the evaluated setting" |
| K15 | Overall CONFIRMED | U | rule unsatisfiable with five IMPROVED pooled results | no global verdict; per-experiment verdicts |
| K16 | HG005 hap.py 0.9521 / 0.9519, 2 discordant loci | RO | raw absent; the pre-fix run left in the tree is invalid (quarantined) | "reported in the project record; raw output not available" |
| K17 | HG005 M-1 4,090 / 3,948 / 142 (3.5 %) | RO | arithmetic checks, raw absent | same |
| K18 | "~35 % of truth SNPs" | X | 142/4,090 = 3.47 % | not used |
| K19 | HG005 rescue 142 = 13 + 129 | RO | S1 | record-only; not carried to other datasets |
| K20 | HG005 "rescue recall 0.986" as recall of true SNPs | U | 129 of 142 are FP avoidance | not used; two populations reported separately |
| K21 | chr20 "true-SNP rescue recall 0.9992" | SQ | 1,327/1,328 is capture of PB-rescuable true SNPs (PB-correct, binomial-wrong), not recall against truth | "router capture of PB-rescuable true SNPs" |
| K22 | HG005 shows generalisation across samples | U | one 3-Mb region, no CI, raw absent | "one region of one unrelated sample; no population claim" |
| K23 | External callers: AI between DeepVariant and GATK | SQ | F1 order true (0.9494 < 0.9547 < 0.9682) | with overlap and single-region caveats; no ranking |
| K24 | External region = blind holdout | X | contains selection span | never used |
| K25 | AI P 0.9893 / R 0.9228 | X | not found in any hap.py output | not used |
| K26 | GT accuracy 0.989 and F1 0.955 together | SQ | scopes 16,124 vs 16,094 calls | one scope per number |
| K27 | Genotype layer leaves allele calls unchanged | S | symmetric difference 0 | as stated |
| K28 | Extraction speedup 21.0× / 65.1× | SQ | stage-level; 65.1× reused single run | "extraction stage only" |
| K29 | "65× end-to-end", "17–80× faster than callers" | X | frozen-mask scope; withdrawn | stated as withdrawn |
| K30 | 1.22× end-to-end | SQ / OD / RO | historical counts-only configuration; source report lost | labelled historical, record-only |
| K31 | `load_reads` no speedup; PB+`load_reads` = 98.7 % | OD | native `load_reads` exists | historical configuration only |
| K32 | `load_reads` 26.9× / 35.1× | SQ | function-level; medians of 3 (baseline n = 2 on HG002); loaded desktop | function-level only |
| K33 | Pipeline 2.72× serial; 7.11× with 4 processes | SQ | 7.11× baseline is pure-Python serial; process gain ≈ 2.6× | decomposed |
| K34 | PB ≈ 93–96 % of compute | SQ | 93 % (HG005, 1 worker), 95.9 % (chr20, 6 workers on 4 cores, contention-inflated) | both figures with conditions |
| K35 | Projected caller-stage speedup 171–396× | SQ | projection from throughput × routed counts | labelled projection; not used as a result |
| K36 | Measured BAM→calls 1.65× cascade vs PB-only | S | HG004 0.5 Mb, n = 3, sd 1.5 / 64.5 s | as stated with n |
| K37 | Whole-genome time estimates | U | linear extrapolation | omitted from results; one line in limitations |
| K38 | Native counts backend "bit-exact" | U | whole-span pass differs from per-window on adversarial repeated-name BAMs | "identical on listed validated inputs"; limitation stated |
| K39 | Native `load_reads` 0 mismatches in 6.7×10⁸ cells | SQ | report-level; 17 tests re-run | "on the listed inputs" |
| K40 | 102/102 tests, 12-category suite, 11 regression tests | OD | originals lost | not cited; current status only (490 pass, 21 xfail, 5 xpass, 2 collection errors) |
| K41 | Locus bootstrap CIs | SQ | adjacent loci treated as independent; block CI only v19/v20/R3 | said where CIs are used |
| K42 | Coordinate bug affected one script only | RO | historical audit file lost; v16–v19 not re-run with the fixed extractor | stated with both caveats |
| K43 | "Band captured 98.6–99.7 % of rescuable loci" | SQ | mixes populations | separate populations |
| K44 | "IMPROVED" means cascade is a better caller | U | driven by PB false positives in paralogous sequence | interpretation, flagged as such |
| K45 | "Pre-registered" | SQ | rule is in project devlog; no external registry; timing vs v13 unclear | "pre-specified in the project devlog" |
| K46 | Constants are depth-specific (30× optimum 11.5/14.0) | SQ | worktree devlog record | as stated, record-level |
| K47 | Novelty of selective escalation for SNP calling | not established | no literature search done | no novelty claim; `[CITATION REQUIRED]` |
| K48 | "No chromosome-scale run"; "largest span 12 Mb"; "seven experiments" | OD / X | v20 exists | not used |

---

## 4. Manuscript design decisions

- Length: the working manuscript is ~110 KB. The rewrite targets roughly a third of that. The neural-phase history, per-cell tables, genotype method comparison, correction ledger and cutoff sweep go to supplementary tables that already exist in `TABLES_FINAL.md`.
- Results order follows the requested structure. Rescue analysis sits under 3.4 (failure analysis) because it is about what the router does and does not catch.
- Every DEGRADED/NEGATIVE result appears in the Results text, not only in limitations.
- HG005 numbers carry an "RO" marker in every table row.
- Speedups: four labelled levels (function, extraction stage, pipeline, projection). Only the first three are results.

### Tables

| # | Content | Source |
|---|---|---|
| T1 | Frozen parameters | `cascade.py`, freeze manifests |
| T2 | Datasets and experiments | `TABLES_FINAL.md` T1 |
| T3 | Primary results, cascade vs PB-only, with verdicts | `TABLES_FINAL.md` T2 |
| T4 | Denominators (chr20, chr21, HG005) | `chr20_results.json` `m1_audit`; HG005 RO |
| T5 | Failure taxonomy | `results/bench_v14/17/19/20/disagreements.json`; claim map row 56 |
| T6 | Router capture / rescue populations | `chr20_results.json`; HG005 RO |
| T7 | External callers | `HEAD_TO_HEAD_RESULTS.csv`; `unified_happy` |
| T8 | Performance, four levels | `LOAD_READS_OPTIMIZATION_REPORT.md`, `PERFORMANCE_RESULTS.csv`, `CHR20_VALIDATION_REPORT.md` §10 |
| T9 | Limitations register | manuscript §5 |
| S1–S4 | v14 cells, v19 cells, genotype methods, correction ledger | `TABLES_FINAL.md` |

### Figures (six, main text unless noted)

| Fig | Purpose | x | y | Grouping / uncertainty | Source | Exact numbers | Main/Supp |
|---|---|---|---|---|---|---|---|
| 1 | Architecture and data flow, with the two evaluators and the two backends | — | — | boxes; constants annotated | `cascade.py`; existing `fig1_architecture_v2.png` (no numeric content; can be reused) | 7.0, 5.4119 (band 1.588–12.412), 10.5 | main |
| 2 | Accuracy difference across experiments | ΔF1 (×10⁻³), grey band ±1 | experiment (v13…v20, three hap.py points) | colour = sample; bars = 95 % paired locus bootstrap; second bar = block CI where it exists (v19, v20); hap.py points hollow, no CI; v20 labelled DEGRADED, v14 cell −1.262 marked | Table 3, `TABLES_FINAL.md` T2 | all ΔF1 and CIs in T3 | main. Must be **regenerated**: existing Fig. 3 omits v20 and the third hap.py point. A depth-resolved version is possible only for v14 and v19 (S1, S2); per-cell v15–v18 are not in the tables → supplement only for those two |
| 3 | Router behaviour: how much goes to PB and what it captures | (a) experiment; (b) rescue class | (a) routed fraction (%, log); (b) loci count | (a) pooled and per-depth range; (b) chr20 stacked: true-SNP 1,327/1,328, FP-avoid 1,490/1,498; random-router control ΔF1 −0.0186 as annotation | `chr20_results.json`; Table 3 | 0.137 %, 77,248; 2,826; 1,328/1,498 | main. HG005 13/129 as hatched RO bars, or omit |
| 4 | Cascade-specific failure analysis | (a) VAF; (b) experiment | (a) binomial LLR with band [1.588, 12.412] and dashed 7.0; (b) events | (a) 8 pre-chr20 + chr20 event = 9 points (v15 event has no covariates, stated); (b) stacked classes 15/153 pre-chr20 plus v20 8/1/1 | `results/bench_v14/17/19/20/disagreements.json` | per-event VAF and binomial LLR | main |
| 5 | External-caller comparison | recall | precision | F1 iso-lines; single points, no CI (none exist); AI arms marked "region overlaps selection span" | `HEAD_TO_HEAD_RESULTS.csv`; `unified_happy` | T7 | main. No runtime axis (no matched runtime) |
| 6 | Where the time goes | (a) configuration; (b) stage | (a) speedup, four labelled levels; (b) share of wall time | (a) function 26.9/35.1; pipeline 2.72 serial, 7.11 = 2.72 × 2.6; historical 1.22; (b) PB 287/309 s | `LOAD_READS_OPTIMIZATION_REPORT.md`; `PERFORMANCE_RESULTS.csv` | as listed; medians of 3 with n noted | main. Projections not plotted |

Not proposed: forest of hap.py CIs (none computed); cost-vs-genome-size (no measurements beyond 3 Mb pipeline / 64 Mb chr20 extraction); runtime comparison with callers (none matched); genotype confusion matrices (table suffices).


Note (figure numbering, final): in `PAPER_MANUSCRIPT_V2.md` the failure-analysis figure is Figure 3 and the router-behaviour figure is Figure 4, so that figures are numbered in order of appearance.
