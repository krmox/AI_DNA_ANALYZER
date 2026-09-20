# Rewrite pass 3 — hostile but fair peer review of `PAPER_MANUSCRIPT_V2.md`

Reviewer stance: try to break the paper. The manuscript file is not edited by this review. Items I fixed **before** the review, while re-reading my own draft against the sources, are listed first so nobody has to discover them.

## 0. Errors in my own draft, caught and corrected before this review

| # | Error | Fix |
|---|---|---|
| E1 | Discussion said the F1 = 0.624 artefact "cost a day to diagnose". No source says that. Invented. | Removed. |
| E2 | Conclusion said accuracy was kept within small margins "in seven of eight experiments", which skipped v14's experiment-level NEGATIVE. | Rewritten: all eight pooled |ΔF1| ≤ 0.0027; 2 PRESERVED, 5 IMPROVED, 1 DEGRADED; v14 NEGATIVE stated. |
| E3 | Abstract put genotype accuracy 0.989 (16,124-call scope) and hap.py F1 0.955 (16,094-call scope) in one clause, the exact scope mixing the audit flagged (C31). | Scopes now separated. |
| E4 | Introduction said every neural stage was null or negative; the record has one weak positive (2/3 seeds) on a side branch. | Reworded. |
| E5 | Discussion said every methodological defect was found because a number was "too bad to accept"; the window filter was found as a denominator mismatch. | "The first three" wording. |
| E6 | Methods said "two DEGRADED or NEGATIVE" pooled results, double-counting v14. | Reworded. |
| E7 | I first wrote that 2 v14 departures were unaccounted for. The raw JSON resolves them (`router_recovers_fn`, 2). | Corrected in manuscript and audit. |
| E8 | The chr20 gate was called "reproducible from the repository"; the repository holds the recorded comparison (`gate/*.json`, `all_bytewise_equal: true`), not the arrays. | "Recorded in the repository". |

Also verified: no banned filler phrases in the text (grep for "importantly", "furthermore", "taken together", "significantly", "robust", "novel", "seamless" etc. returns only the sentence disclaiming novelty).

---

## CRITICAL ISSUES

**CR-1. The evidence map asserts raw verification of files that do not exist.** `CLAIM_EVIDENCE_MAP.csv` rows 43, 45, 46, 48 and 63 are `VERIFIED_RAW` and cite `cache/hg005_stress_postfix/…npz` (row 63 also `cache/hg005_stress_native`). Neither cache exists (checked). Rows 38, 42, 44 cite the same cache as `RECORD_ONLY_RAW_UNAVAILABLE`. The manuscript's Supplement points readers at this map as the claim→source index. A reviewer who opens it will read "VERIFIED_RAW" next to numbers (routed 1,600; rescue 13/129; M-1 4,090/3,948/142; native-vs-pysam identity) that cannot be reproduced. The V2 text treats all of them as record-only; the map does not. The old manuscript also said the HG005 cache identity was "re-verified for this paper". **Action (author):** relabel those five rows before the map is shared. Do not describe the pre-wipe check as a current verification. **RESOLVED:** the five rows are now `RECORD_ONLY_RAW_UNAVAILABLE`, values unchanged; `TABLES_FINAL.md` S4 wording changed from "recomputed from cache" to "record-only".

## MAJOR ISSUES

**MA-1. The chr20 result rests on a backend validated, inside the repository, on 1.7 % of its loci.** The native-vs-reference gate on chr20 covers two chunks, 980,928 of 56,266,816 loci. The 12-Mb and HG005 counts-equivalence evidence is record-only. The counts backend is known to differ from per-window passes on adversarial inputs. PB-only and cascade share the same extraction, so the paired ΔF1 is less exposed than absolute F1 or hap.py F1, but a difference in extraction that changes which loci are routed is not covered by that argument. The V2 text states the facts but not this consequence. **Options:** (a) state it in §3.6 or §5 as written here; (b) extend the gate to more chr20 chunks with the pure-Python path. **Decision: option (b) not taken; scope stated in the manuscript as tested chunks only.** Option (b) would have been a new validation run on frozen code, not a change to any result; whether to do it is the author's decision.

**MA-2. Code and results are not available to a reviewer.** Commit `79155df` is local; `origin/main` is `16b221e`. The remote's visibility is unknown. The V2 Data availability says this. As it stands the paper cannot be reproduced by anyone but the author. **Action:** push, and archive with a DOI, or state that the code is available on request.

**MA-3. The decision rule cannot support "preserved within small margins" as a statistical statement.** It has no equivalence test. The 0.001 margin has no stated derivation (about 1.5 false calls at 765 SNPs). DEGRADED can be triggered by a tiny difference in a well-powered cell (chr20), while PRESERVED in a small cell only means the interval contained zero. The V2 text uses "small margins" as a description of |ΔF1| and says so, but a reviewer will still ask why 0.001. Origin of the margin: **not established** in the sources. **Action:** say who chose it and on what basis, or say it was an a-priori convention.

**MA-4. Multiplicity and independence of the negative result.** Eight experiments and about 60 cells were examined with no correction; the v20 McNemar p = 0.021 rests on ten events; the locus bootstrap treats adjacent loci as independent, and the block interval exists only for v19, v20 and round 3. The v20 verdict survives the block bootstrap, which helps, but no block interval exists for the experiments that produced the five IMPROVED verdicts (v14–v18). The V2 text discloses all of this; the reader still cannot tell how many of the IMPROVED verdicts survive block resampling. **Action:** optional; block-bootstrap re-scoring of v13–v18 from existing caches is not possible for v13–v18 because "earlier bench_v12–v18 caches were deleted" (limitations file §E).

**MA-5. The external comparison sits on a region that overlaps the selection span.** V2 says so three times. A reviewer may argue it should be dropped. A cheaper alternative that needs no new run: stratify the existing U-H2 calls into chr21:32–40 Mb (inside the span) and chr21:40–44 Mb (outside). That is a new post-hoc analysis and would need to be labelled as such; it was not done here.

**MA-6. The cascade-introduced false-positive mechanism (fixed ε too optimistic at low BQ) is an interpretation on 23 events**, and the lost-SNP mechanism is an interpretation on 10 events with no read-level test. V2 labels both as interpretation. The Discussion's sentence "a different second opinion, not a wider band, is the direction the data point to" goes a step further than the data: the safety-layer failure and the cutoff sweep support "a wider band does not help", not "a different second stage would". **RESOLVED:** the sentence now says a different second stage is untested.

## MINOR ISSUES

- **MI-1.** The cutoff derivation conflict (main-branch devlog 11 vs router-freeze worktree) is disclosed but unresolved; a reviewer with repo access will find both.
- **MI-2.** "Pre-specified" is defensible; "pre-registered" (used in the sources) is not, because no external registry exists. V2 uses "pre-specified". Whether the rule preceded the v13 numbers is not established.
- **MI-3.** Speedup ratios are ratios of unrounded times; printed seconds give 26.8× and 35.0× (recomputed 79.95/2.98, 89.63/2.56) against 26.9× and 35.1×. Stated in Table 8. The 65× is quoted as "≈ 65×" because the map says 65.14× and the printed times give 65.16×.
- **MI-4.** The async 1.21× baseline is ambiguous in the sources (Python reads vs native serial). Flagged in Table 8.
- **MI-5.** Loss rates ("1 per 10⁷", "1 per 5.6×10⁷") use per-cell events and different denominators; treat as descriptive.
- **MI-6.** Distinct-locus count of the 10 lost SNPs is not established (only v14's 5 events = 4 loci). Checking `results/bench_v*/disagreements.json` for coinciding positions would settle 8 of 10; the v15 event has no stored position.
- **MI-7.** PB/binomial cost ratio: three different ranges in the sources; V2 gives the two throughputs and "several hundred-fold".
- **MI-8.** Reference-list entries marked ‡ were not re-verified; three citations are placeholders.
- **MI-9.** The 2×250 novoalign HG005 product was unavailable; the 300× product with ≤ 250-bp reads is not platform-identical. Stated in §2.2; the library difference is not established.
- **MI-10.** Author Contributions and Conflict of Interest are placeholders; only recorded facts are filled in.
- **MI-11.** Figures 2, 3, 4 and 6 need to be regenerated to include v20 and the current bottleneck picture; nothing was generated in this pass.

## NO ACTION REQUIRED (checked)

1. **Every major claim has evidence.** Each headline number traces to a raw file or a labelled record (audit and map). chr20 numbers were recomputed from raw hap.py summaries and the report tables in this pass.
2. **Numeric consistency.** chr20: 71,387 − 1,030 = 70,357; − 33 = 70,324; TP+FN = 70,324 in both arms; ΔFP +7 = 8 introduced − 1 avoided; 1,328 + 1,498 = 2,826; 1,327 + 1,490 = 2,817; routed 0.13735 %; McNemar p = 0.0215. v14 taxonomy 26 + 8 + 5 + 2 = 41. v19 83 + 3 + 1 = 87; v17 36 + 2 = 38. Pooled verdict count 2 + 5 + 1 = 8.
3. **Historical numbers not presented as final.** 1.22×, 98.7 %, "load_reads: none", 65× and 17–80× are labelled historical, record-only or withdrawn. "102/102" is not cited.
4. **Negative results disclosed in Results, not only in Limitations:** v14 NEGATIVE, v20 DEGRADED, ten lost SNPs, safety layer rejected, both native-design failures, F1 0.624, HG005 raw loss.
5. **Speedup definitions.** Function, stage, pipeline and projection are separated in Table 8; 7.11× is decomposed (2.72× × ≈ 2.6×); routing is not equated with CPU saving; projections are not plotted or used.
6. **Statistical conclusions match the frozen rule.** v20 = DEGRADED by clause 1, block CI and sensitivity agree; magnitude stated (−5.6×10⁻⁵, 18× inside the margin) without "substantial" or "no significant".
7. **chr20 is not called independent.** Same sample, library, depth; 1 Mb previously in v14; one depth.
8. **HG005 raw loss disclosed** in Abstract, §3.2, §3.5, Table 4, 6, 10, Data availability.
9. **Native counts limitation disclosed** in Abstract, §3.6, §5. "Bit-exact" never claimed.
10. **External comparison caveated** (overlap with selection span, GATK without BQSR, single region, unnormalised threads, no matched runtime). No ranking.
11. **Novelty:** none claimed.
12. **Reproducibility of calculations:** the Methods give the formulas, constants, seeds, hashes, image digests and scripts; the HG005 rows are the exception.
13. **AI-sounding filler:** grep-clean of the listed phrases; sentences that only restated a table were removed. Whether the prose reads as the author's own voice cannot be judged from the prompt alone.
14. **Length:** ≈ 9,100 words including tables and markup, against ≈ 16,400 for the working manuscript. Still long for a regional-benchmark paper; the Discussion and §3.6 are the candidates for further cutting.
