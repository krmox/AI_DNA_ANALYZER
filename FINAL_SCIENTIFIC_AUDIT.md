# Final Scientific Audit — AI_DNA_ANALYZER

Independent adversarial audit of `research/PAPER_MANUSCRIPT_V2.md` (the manuscript that feeds
`research/latex/paper.tex` / `paper.pdf` and is the one linked from `README.md` §Publication) and
its supporting evidence, code, figures and repository metadata, run against the question: **can this
be shown to a scientific mentor as a serious research project without hiding methodological
weaknesses?**

Scope read in full: `PAPER_MANUSCRIPT_V2.md`, `PAPER_MANUSCRIPT_V2.html`, `latex/paper.tex`,
`latex/preamble.tex`, `latex/md2tex.py`, `FINAL_FORENSIC_AUDIT.md`, `FINAL_FREEZE_REPORT.md`,
`CLAIM_EVIDENCE_MAP.csv`, `TABLES_FINAL.md`, `FIGURE_PLAN.md`, `RESEARCH_TIMELINE.md`,
`LIMITATIONS_AND_OPEN_QUESTIONS.md`, `CHR20_VALIDATION_REPORT.md`, `CHR20_RESULTS.csv`,
`AI_ASSISTANCE_DISCLOSURE.md`, `RESCUE_ANALYSIS.csv`, `PERFORMANCE_RESULTS.csv`,
`BENCHMARK_MASTER.csv`, `ERROR_ANALYSIS.csv`, all six live figures and their generators, `README.md`,
`README_RESEARCH.md`, `CITATION.cff`, `LICENSE`, `REPRODUCIBILITY.md`, and the code for the router,
Binomial/PB models, Method C genotyping, the native C/htslib backend and its equivalence tests, the
`classify()` verdict rule, and the bootstrap implementations. Method: direct reading by the auditor
plus five independent adversarial sub-investigations (threshold-fitting independence, HG005
record-only evidence, native-backend/Method C code-vs-paper, performance-claim conflation +
statistics + benchmark fairness, reproducibility + figure integrity), cross-checked against each
other and against the manuscript text quoted above.

**Headline result of this audit: no FATAL issue was found.** The manuscript is unusually — almost
aggressively — self-critical: it already discloses a NEGATIVE experiment (v14), a DEGRADED
chromosome-scale result (v20/chr20), ten lost true SNPs, a withdrawn "65.1× end-to-end" and
"17–80× faster than GATK/DeepVariant/Clair3" claim, an explicit refusal to claim "bit-exact" for the
native counts backend, and a region that overlaps its own constant-fitting span used for the external
caller comparison. This is not the common failure mode of a manuscript overselling itself. The
findings below are real, but they are almost entirely about **provenance/navigation hygiene around
the manuscript**, not about the manuscript's own scientific content, which held up under adversarial
cross-checking.

---

## Pass A — Scientific validity (leakage, threshold independence)

No FATAL or MAJOR issue. Threshold-fitting overlap (HG002 chr21:32–40 Mb, ≈15×) is disclosed
directly in the abstract, in Table 1's "Relation to selection data" column (per-experiment), and
repeated in bold at the external-caller comparison ("a region that contains the 8-Mb span from which
the AI constants were selected," §3.7 and Figure 6 caption). chr20 (v20) is correctly scoped as
chromosome-separated but *not* sample/library/depth-independent, and this is stated in §3.3
*Independence* and Table 10. The one candidate leakage-disclosure gap raised by a sub-investigation —
an internal-record conflict about how the router cutoff was originally derived (F-beta optimum on a
training split, per one devlog, vs. the 0.1% validation-quantile rule the project actually froze) —
**was checked directly against the manuscript text and is already disclosed**, verbatim, in §2.3:
*"One internal project record describes the cutoff as an F-beta optimum on a training split; that
contradicts the record kept at the time the router was frozen, which we follow."* No action needed.

## Pass B — Statistical validity

No FATAL or MAJOR issue. The locus bootstrap is a genuine paired multinomial resample
(`evaluate_quality_error.py`); the block bootstrap (v19, v20) is a real 50-kb block resample that
correctly produced a *wider* CI than the naive locus bootstrap (2.3× on v19) — this is reported, not
hidden. The i.i.d.-loci assumption of the locus bootstrap is disclosed as a limitation in §2.4 and
Table 10, not silently assumed. `classify()` in `bench_v14_crosschrom.py` was read directly: it
implements PRESERVED → IMPROVED → DEGRADED → UNDERPOWERED exactly as the manuscript describes,
including the disclosed clause defect (a zero-crossing CI with `|ΔF1| ≥ 0.001` falls through to
DEGRADED regardless of sign) — code and text agree. The chr20 "DEGRADED but ΔF1 ≈ −0.0056 percentage
points" case is presented exactly as intended: a frozen categorical rule that fired, stated alongside
its own effect size, not as a contradiction and not softened into "no significant degradation" (§3.3
explicitly rejects that wording). No confusion found between the ε-sensitivity "safety layer"
sensitivity analysis (§3.4.3, explicitly "nothing was changed on the basis of these sweeps") and
production threshold selection.

## Pass C — Reproducibility

No FATAL issue. Two MODERATE gaps (see findings table): `REPRODUCIBILITY.md`'s "Exact commands
(HG005 post-fix)" section cites two scripts that no longer exist in the repository
(`run_hg005_pipeline_postfix.py`, `run_happy_postfix_regionscoped.sh` — casualties of the same
`/tmp` wipe that took the raw HG005 artefacts), presented as if runnable; and the chr20 full
extract→pipeline→hap.py→verdict chain is traceable to frozen, checksummed outputs but was not
re-executed end-to-end in this or the prior audit session (this is disclosed in README.md itself:
*"This chain was not re-run end to end for this README"*). A live re-run of the documented
"Verify the checkout" test command was performed by one sub-investigation and reproduced the exact
counts README claims (`117 passed, 21 xfailed, 5 xpassed`, contingent on running `native/build.sh`
first — also documented). Data/code availability is concrete (GIAB accession paths, SHA-256 checksums,
`FROZEN_MANIFEST.json`, pinned Docker image digests for hap.py/DeepVariant/Clair3), not a vague
"available on GitHub."

## Pass D — Claim discipline

No FATAL or MAJOR issue in the manuscript body. Terminology was checked exhaustively: no instance of
"lossless," "never loses a true SNP," or unqualified "no accuracy degradation" exists in the current
manuscript (the one historical use of that exact claim, in an old devlog, is explicitly named and
refuted, §3.4.1). "AI" in the project name is directly disclaimed in README.md ("statistical cascade
with no neural network... not a deep-learning caller"). Performance numbers carry an explicit
"Level" column (function / pipeline / measured / projection) and the manuscript actively withdraws two
historical conflated claims rather than repeating them. Genotype accuracy, allele-level F1 and hap.py
F1 are kept explicitly distinct throughout §3.5, with a direct warning against pooling two differently
-scoped genotype accuracies. The **one live claim-discipline problem found is external to the
manuscript body**: the citation-facing title in `CITATION.cff` and the H1 subtitle of `README.md`
call the method "**Adaptive** statistical SNP variant calling with **confidence-aware** compute
routing," while the manuscript's central, repeated methodological claim is the opposite — three
constants "selected once... and never changed" (Abstract), frozen, non-adapting. See Finding 2.

## Pass E — Internal consistency

**This is where the real, fixable problems are.** Numeric spot-checks (chr20 P/R/F1, ΔF1 and CI,
McNemar statistic, HG005 rescue decomposition, true-SNP loss reconciliation 9→10, test-suite counts)
all matched exactly across the manuscript, `CHR20_RESULTS.csv`, `RESCUE_ANALYSIS.csv` and
`CLAIM_EVIDENCE_MAP.csv`. But the repository carries **two full manuscripts** — `PAPER_MANUSCRIPT.md`
(V1) and `PAPER_MANUSCRIPT_V2.md` (a further rewrite) — and the cross-reference network inside
`research/` almost uniformly points to the wrong one. See Finding 1. Separately, one derived-metric
discrepancy was found and not reconciled anywhere in the evidence chain: see Finding 4.

## Pass F — Figure integrity

No issue with the six figures actually cited by the manuscript. All six were viewed directly. Every
one prints its own scope caveat on the plot itself (open/hatched markers or explicit text for
record-only HG005 points, an explicit "no ranking implied" on the external-caller figure, an explicit
"projections, not results, not plotted" note for the 171–396× projected speedups, per-experiment
verdict labels rather than one pooled arrow). Generator provenance (`make_v2_figures.py`,
`make_v2_figures_more.py`) was read and traces every plotted number to `TABLES_FINAL.md`,
`CHR20_RESULTS.csv` or the raw `disagreements.json` files — no stale or alternate data source feeds a
live figure. One MINOR housekeeping issue: 12 of 18 PNGs in `research/FIGURES/` are orphaned v1
duplicates not referenced by the current manuscript or `paper.tex` (see Finding 5).

## Pass G — Code-paper consistency

No FATAL or MAJOR issue. `native/pileup_native.c` and `native/reads_native.c` were read directly
against the manuscript's equivalence claims: the manuscript states, verbatim, **"'Bit-exact' is not
claimed for this backend"** for the counts backend, gives the exact tested fraction of chr20
(980,928 / 56,266,816 loci, 1.7%), and separately scopes the read-tensor backend's much larger
equivalence gate (>6.7×10⁸ tensor cells). Method C's code (`method_c_regression.py`) matches the
manuscript's posterior formula line-for-line (three binomial hypotheses θ∈{ε, 0.5, 1−ε}, ε=0.01
fixed, GQ from the best/second-best likelihood ratio, the `0/0→written 0/1` special case), and the
"0 GT mismatches / 0 GQ mismatches" claim is backed by an actual reference-vs-vectorised comparison
with a stated population size. One MINOR finding: several pre-V2 research documents still carry the
old, unscoped "proven bit-exact" wording that the current manuscript deliberately walked back (Finding
6) — this doesn't touch the manuscript itself but sits in the same directory a reader might browse.

## Pass H — Publication readiness

Conditional on Findings 1–4 below being fixed (which this audit did, see "What was changed"), the
manuscript is ready to show a mentor. Its own candor about its NEGATIVE experiment, its DEGRADED
chromosome result, and its lost raw artefacts is exactly the kind of disclosure the audit brief asked
whether a reviewer could honestly say exists. Before this audit's fixes, it could not be called ready
for that reason: the repository's own index and citation metadata actively pointed a reader at the
wrong manuscript and used the one word ("adaptive") that most directly contradicts the paper's central
methodological identity. These are now fixed (see below).

---

## Findings

| Severity | Finding | Evidence | Blocking? | Action |
|---|---|---|---|---|
| **MAJOR** | **(1) Provenance fork: two "current" manuscripts, wrong one privileged repo-wide.** `research/PAPER_MANUSCRIPT.md` (V1) predates the chr20/load_reads rewrite that became `PAPER_MANUSCRIPT_V2.md`, but V1 carries no pointer to V2's existence, and *every* superseded-banner in `research/` (`METHODS.md`, `RESULTS.md`, `SCIENTIFIC_CLAIMS.md`, `OUTREACH_PACKAGE.md`, `PAPER_DRAFT.md`, `RESEARCH_SUMMARY_1PAGE.md`) plus `README_RESEARCH.md`'s own "how to read this dossier" index all point to V1 as authoritative. Only the top-level `README.md` §Publication and the LaTeX/PDF build (`md2tex.py:13`, hardcoded) correctly use V2. | `research/PAPER_MANUSCRIPT.md:1-2`; `research/README_RESEARCH.md:1,11`; `research/latex/md2tex.py:13`; six banner files checked directly | Yes, until fixed | **FIXED** — added a redirect banner to the top of `PAPER_MANUSCRIPT.md` and corrected `README_RESEARCH.md`'s reading-order index to name V2 as canonical. Did not touch the six individual historical banners (they were accurate when written; the one-file fix at the root of the reference chain resolves the ambiguity without rewriting history). |
| **MAJOR** | **(2) Citation-facing title contradicts the manuscript's central claim.** `CITATION.cff:3` and `README.md:3` call the method "**Adaptive** statistical SNP variant calling with **confidence-aware** compute routing." The manuscript's repeated, load-bearing claim is the opposite: three constants "selected once on HG002 chr21:32–40 Mb... and never changed" (Abstract); "frozen," never refit, is used throughout as the entire point of the freeze-and-validate design. | `CITATION.cff:3`; `README.md:3`; contrast with `PAPER_MANUSCRIPT_V2.md:9` ("never changed"), Table 2 ("frozen parameters") | Yes, until fixed | **FIXED locally** — retitled to "frozen-constant statistical SNP variant calling with confidence routing" (reusing the manuscript's own keyword, "confidence routing," rather than inventing new language) in both files. **Cannot fix retroactively**: `CITATION.cff`'s DOI (`10.5281/zenodo.22895746`) is already minted and archived on Zenodo with the old title baked into that permanent record; per this audit's instructions, Zenodo itself was not touched. The user would need to submit a Zenodo metadata correction separately if the archived title should also change. |
| **MAJOR** | **(3) Manuscript's own Data-and-Code-Availability section is stale.** `PAPER_MANUSCRIPT_V2.md`, `paper.tex` and `PAPER_MANUSCRIPT_V2.html` all still say *"The repository is currently private `[AUTHOR TO CONFIRM: make it public before submission, and archive a release with a DOI]`."* Verified directly via `gh api repos/Kramkost/AI_DNA_ANALYZER`: the repository is **public**, `origin/main` is fully synced with local `main` (top commit `ca60e12`), and a Zenodo DOI (`10.5281/zenodo.22895746`, release `v1.0.0`) is already archived and referenced from `README.md`/`CITATION.cff`. The bracketed author placeholder in the actual manuscript text was simply never removed after this happened. | `gh api` output (`"private": false, "visibility": "public"`); `git log origin/main` matches local `main`; `README.md:5,220` DOI badge/section | Yes, until fixed | **FIXED** — replaced the placeholder in `PAPER_MANUSCRIPT_V2.md` with the actual current state (public repo, DOI, release tag), then regenerated `paper.tex` and `PAPER_MANUSCRIPT_V2.html` from the corrected Markdown via the existing pipeline. No scientific content touched. |
| MODERATE | **(4) Unreconciled HG005 discrepancy: `false_rescue_count` (9) vs. Table 5 `router_fp` HG005 contribution (2, record-only).** `RESCUE_ANALYSIS.csv` defines `false_rescue_count = 9` as "routed & binomial_correct & pb_wrong" for HG005 chr1:1–4 Mb — conceptually the same event class as the manuscript's `router_fp` taxonomy (a routed locus where the binomial screen was right and the cascade, which adopts the PB call once routed, was wrong). Table 5 records only 2 (RO) for HG005 in that class. Neither number is wrong on its face — they may use different locus-level definitions (e.g. event vs. distinct-locus counting, or a broader "any wrong call" vs. FP-specific criterion) — but nothing in the manuscript, `CLAIM_EVIDENCE_MAP.csv` or any audit document states why they differ, and no current document even acknowledges both numbers exist side by side. | `research/RESCUE_ANALYSIS.csv` row `false_rescue_count`; `PAPER_MANUSCRIPT_V2.md` Table 5, HG005 row of `router_fp` (pre-chr20 column, "15 (v14 8, v16 2, v19 3, HG005 2 RO)") | No — both figures concern record-only HG005 data already extensively caveated, and neither feeds the chr20/pooled headline conclusion | **Disclosed, not resolved.** Added as a new open question to `research/LIMITATIONS_AND_OPEN_QUESTIONS.md` §C (matching that document's existing style for this exact kind of unreconciled-count item, e.g. its existing item 1 on true-SNP loss counting). **`NEW EXPERIMENT REQUIRED`** if a resolution needs the original per-locus classification script re-run against the surviving HG005 cache — did not attempt this, since re-deriving it risks producing a number that looks resolved but wasn't independently checked against lost raw artefacts. |
| MODERATE | **(5) `REPRODUCIBILITY.md` "Exact commands (HG005 post-fix)" cites two scripts that no longer exist.** `experimental/stress_test/run_hg005_pipeline_postfix.py` and `experimental/stress_test/run_happy_postfix_regionscoped.sh` are both absent from the working tree (casualties of the same `/tmp` wipe that took the raw HG005 artefacts) but are presented as a runnable four-step command block. | Verified: `ls experimental/stress_test/run_hg005_pipeline_postfix.py` → no such file; same for the hap.py script; `research/REPRODUCIBILITY.md:76,79` | No — the document already carries a document-level "superseded in part" banner, and the surrounding text elsewhere discloses the HG005 raw-artefact loss | **FIXED** — added an inline note next to the command block stating both scripts were lost in the 2026-09-20 wipe and are not currently runnable, rather than leaving a silently-dead command block. |
| MINOR | **(6) 12 of 18 PNGs in `research/FIGURES/` are orphaned v1 duplicates**, not referenced by the current manuscript or `paper.tex` (only the six `*_v2_*` files are live). Confirmed by grep of both `PAPER_MANUSCRIPT_V2.md` and `research/latex/paper.tex` `\includegraphics` lines. | File listing vs. grep of both citation paths | No | Not fixed — noted here only. Deleting or moving tracked files is a repo-structure change beyond what this audit's brief authorizes ("do not delete inconvenient results"; these aren't results, but the same caution applies to unrequested deletions). Recommendation for the author: move the 12 orphans to `research/FIGURES/archive/` in a future, deliberate pass. |
| MINOR | **(7) Stale unscoped "bit-exact"/"proven bit-exact" wording survives in pre-V2 documents** (`METHODS.md`, `OUTREACH_PACKAGE.md`, `SCIENTIFIC_CLAIMS.md`, `PAPER_DRAFT.md`, `RESULTS.md`, `RESEARCH_SUMMARY_1PAGE.md`) that the current manuscript deliberately walked back to "identical on the listed validated inputs." Each of these six files already carries a document-level "SUPERSEDED IN PART" banner (verified directly), so this is a residual browsing risk, not a live overclaim. | Banner text confirmed present at the top of all six files; `PAPER_MANUSCRIPT_V2.md:336` "'Bit-exact' is not claimed for this backend" | No | Not fixed — already adequately covered by the existing banner convention; re-editing six historical documents' body text was judged out of proportion to the residual risk. |
| OBSERVATION | Router-cutoff derivation conflict (F-beta optimum vs. validation-quantile rule) — checked directly against manuscript text and found **already disclosed** verbatim in §2.3. One sub-investigation flagged this as a gap; it is not one. | `PAPER_MANUSCRIPT_V2.md:114` | No | None needed. |
| OBSERVATION | Test-suite reproducibility claim in README (`117 passed, 21 xfailed, 5 xpassed`) was independently re-run live by a sub-investigation and matched exactly, contingent on the documented `native/build.sh` step. | Live re-run, reported by sub-investigation | No | None needed. |
| OBSERVATION | The "AI" in the project name, the "9 vs 10" true-SNP loss count, the "3.5% vs 35%" historical arithmetic error, the "rescue recall 0.986" mixed-population figure, and the forced-`0/1` F1=0.624 genotype artifact are all explicitly named, corrected and never presented as live results anywhere in the current manuscript or its cited tables. | Multiple, see Pass D/G above | No | None needed. |

---

## Addendum — post-audit revision (author-requested, 2026-09-22)

After this audit closed, the author asked to shorten and reword the manuscript's **Author
Contributions** and **AI assistance disclosure** sections. The author's initial proposed text stated
the author "implemented the pipeline, conducted all experiments... and wrote the final manuscript"
and characterised Claude's role as limited to "code refactoring, assistance with test automation, and
preliminary drafting/language editing." **That specific wording was declined**: it contradicts
`research/AI_ASSISTANCE_DISCLOSURE.md` and this audit's own direct observation (Claude wrote the
coordinate fix, ran extraction/cascade/Method C/hap.py, computed the statistics behind the paper's
own numbers, and drafted this manuscript's V2 rewrite) — using it would have reintroduced exactly the
kind of claim-stronger-than-evidence problem this audit was built to catch, in the one section a
reader would use to judge how much of the science to attribute to the author versus the tool.

A shorter, accurate alternative was proposed instead and the author agreed to it. Both sections in
`PAPER_MANUSCRIPT_V2.md` were replaced with wording that keeps the same underlying facts as the
pre-audit text — Claude wrote code, ran the experiments/analyses, and drafted/edited the manuscript
under the author's direction; the author set the research question, design, frozen constants and
acceptance criteria, and reviewed/verified/is responsible for the final content — stated more
concisely. `paper.tex`, `PAPER_MANUSCRIPT_V2.html` and `paper.pdf` were regenerated via the existing
pipeline and the changed page was visually re-checked (24 pages, no new errors, no undefined
references). No scientific number, threshold, or result was touched by this change.

---

## FINAL VERDICT

### SCIENTIFICALLY READY FOR MENTOR REVIEW

No FATAL issue was found anywhere in the science, statistics, code, or evidence chain. Three MAJOR
issues were found — all of them about *pointing a reader at the right document and describing the
project's own current state accurately*, none about the science itself — and all three are fixed as
of this audit (with one residual caveat: the already-minted Zenodo DOI record cannot be retroactively
retitled from here). One MODERATE finding (the HG005 false-rescue-count discrepancy) remains
genuinely unresolved and is now honestly flagged as an open question rather than fixed or hidden,
consistent with how this project has handled every other unresolved discrepancy it has found in
itself. Two MINOR findings (orphaned figure files, residual stale wording in historical documents)
are noted and left for the author's discretion.

`SUBMISSION-READY` is not used: no specific journal or conference venue was given or checked against,
per the audit brief's own instruction that this status requires a named, verified venue.

---

## What remains unresolved (stated plainly, not minimized)

- HG005 `false_rescue_count` (9) vs. `router_fp` (2, RO) is not reconciled (Finding 4). Both are
  record-only.
- 5 of the `CLAIM_EVIDENCE_MAP.csv` rows behind HG005 numbers are `RECORD_ONLY_RAW_UNAVAILABLE`; the
  raw post-fix HG005 cache, error/rescue CSVs and native-validation reports were permanently lost in
  the 2026-09-20 `/tmp` wipe and are not independently reproducible from this repository, ever.
- The native counts backend is not, and does not claim to be, bit-exact on adversarial synthetic BAMs
  with repeated read names and inconsistent mate fields; on chr20 the byte-identical equivalence gate
  covers only 2 of the ~108 500-kb chunks (1.7% of loci).
- The chr20 DEGRADED verdict rests on 10 discordant loci out of 56.2 million; it is statistically real
  (the interval excludes zero) but the manuscript itself declines to say whether a magnitude this small
  "matters" for any particular use — correctly, since that is a question about use, not evidence.
- No genome-scale run exists for any sample; all conclusions are scoped to the tested regions/depths.
- No lockfile exists; environment reproducibility depends on the version table in
  `REPRODUCIBILITY.md`/README rather than a pinned manifest.
- The external-caller comparison (Table 9/Figure 6) sits on a region that contains the constant
  -fitting span; the manuscript already declines to call it a ranking, and this audit agrees that is
  the correct posture, not a fixable defect.

None of the above was hidden by the project before this audit. This audit's job was to check whether
that was true, not to discover it for the first time — and it was true.
