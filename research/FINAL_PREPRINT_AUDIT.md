# Final Preprint Audit — AI_DNA_ANALYZER

Second-pass adversarial audit, run after `FINAL_SCIENTIFIC_AUDIT.md` (verdict: SCIENTIFICALLY READY
FOR MENTOR REVIEW) to independently check a list of concerns raised by an external multimodal reviewer
(Gemini) against `research/PAPER_MANUSCRIPT_V2.md`, its LaTeX/PDF build, its figures, and the
supporting evidence in `research/`. Method: direct reading of the full manuscript (all 491 lines),
the six live figures (viewed as images), `research/latex/paper.log`/`paper.tex`, `README.md`,
`CITATION.cff`, `AI_ASSISTANCE_DISCLOSURE.md`, `REPRODUCIBILITY.md`, `LIMITATIONS_AND_OPEN_QUESTIONS.md`,
and exhaustive text search for unresolved placeholders across the manuscript, its generated HTML/TeX,
and top-level metadata files. No new experiment, threshold, dataset or benchmark was run. No number was
changed. No git commit was made.

**Headline finding of this pass: the manuscript already implements almost every substantive point
Gemini raised.** This is not a case of accepting Gemini's list uncritically — each item was checked
independently against the actual text, tables, figures and code, per the audit brief's instruction not
to fix what isn't broken. Of Gemini's 12 points, 9 were found ALREADY FIXED (the manuscript already
does exactly what the point asks for, in most cases with more precision than the point requests), 1 was
found INCORRECT as a change request (repositioning the paper — it is already positioned that way), and
2 were PARTIALLY CORRECT, leading to one genuine fix (a stale epigraph line) and reconfirmation of a
gap the prior audit had already disclosed rather than resolved.

---

## Gemini review reconciliation

| # | Gemini claim | Evidence in repository | Correct? | Severity | Action |
|---|---|---|---|---|---|
| 1 | Presentation/readiness issues | `paper.log` after rebuild: 24 pages, 0 undefined references, only cosmetic `pdfTeX` duplicate-destination notices (an `xltabular` artifact, not a rendering defect) and two overfull/underfull `\hbox` warnings (table cell wrapping, not content loss). All 6 figures viewed directly: axis labels, units, legends, per-plot caveats (record-only markers, "no ranking implied," projection exclusions) all present and accurate. One stale sentence found (see #8). | PARTIALLY CORRECT | MINOR | **FIXED** — see #8. Everything else checked and found in order; no other presentation defect found. |
| 2 | Computational-savings claim may be weaker than implied because PB still dominates runtime | Abstract itself states "PB is now 93–96 % of compute." §3.6 bold-states **"Routing is not a CPU saving here."** Table 8 has an explicit `Level` column separating function/pipeline/measured/projection. Discussion: "On compute, the results are mostly negative." Fig. 5b plots PB's 93–96% share directly. | ALREADY FIXED | — | None needed. The manuscript already leads with this as a negative result, not a footnote. |
| 3 | Real E2E routing speedup is much smaller than function-level speedups | Table 8 lists 21–65× (function) beside 1.65× (the one measured BAM→calls figure) in adjacent rows with a `Level` column making the distinction explicit. §3.6: "Two statements from the earlier native-backend report are withdrawn: that 65.1× was end-to-end, and that the pipeline was '17–80× faster' than GATK, DeepVariant and Clair3." Fig. 5a plots them on the same log axis in visually distinct colour groups with a legend. | ALREADY FIXED | — | None needed. |
| 4 | HG005 has historical record-only raw artifacts | Tagged `(RO)`/`record-only` at every point of use (Tables 1, 3, 4, 6; Figs. 2, 4, 5); Data & Code Availability section names the exact incident (2026-09-20 `/tmp` wipe) and which evidence-map rows it affects; `REPRODUCIBILITY.md` now states inline that the two HG005 post-fix scripts are not runnable (prior-audit fix, verified still present). | ALREADY FIXED | — | None needed. |
| 5 | Scope is SNP-only | §2.1: "Indels, MNPs and structural variants are neither called nor scored." Table 10 row 1. Title, abstract, README badge (`scope: SNP-only`) all state it. hap.py indel rows are noted as "0 by construction." | ALREADY FIXED | — | None needed. |
| 6 | Validation is not whole-genome independent validation | §3.3 *Independence*: "chr20 is a different chromosome from the selection region. It is the same sample, library and depth. 1 Mb of it was scored in v14... Chromosome-level separation is all this run provides." Table 10: "One chromosome-scale run... no genome-scale run for any sample." §5 Limitations restates it a third time. | ALREADY FIXED | — | None needed — this is stated more precisely than the Gemini point itself. |
| 7 | Native backend has an adversarial repeated-read-name limitation | §3.6 "Equivalence gates and their failures" gives exact figures: 20/20 pass on realistic synthetic BAMs, adversarial repeated-name BAMs diverge in "thousands of cells per BAM" (20/20 fuzz cases fail as expected; 1/6 in a separate regression suite), explicit: **"'Bit-exact' is not claimed for this backend."** Table 10 restates the scope restriction. | ALREADY FIXED | — | None needed. |
| 8 | Possible unresolved `[placeholder]` text | Exhaustive search (`TODO`, `FIXME`, `PLACEHOLDER`, `TBD`, `AUTHOR TO CONFIRM`, plus a full manual bracket-content audit of every `[...]` span in the manuscript body) found **zero genuine unresolved placeholders** in the manuscript, its generated HTML, or `paper.tex`. But line 3 of `PAPER_MANUSCRIPT_V2.md` (the epigraph) still asserted *"Placeholders in square brackets are unresolved on purpose"* — true of an earlier draft state, false now that the one real placeholder (`[AUTHOR TO CONFIRM...]`, fixed by the prior audit) is gone. A reader could search for placeholders on the strength of that sentence and be confused not to find any, or worse, treat any bracketed citation number as suspect. | PARTIALLY CORRECT | MINOR | **FIXED** — removed the stale clause from the epigraph in `PAPER_MANUSCRIPT_V2.md`; regenerated `paper.tex`, `PAPER_MANUSCRIPT_V2.html` and `paper.pdf`. All other bracket content in the manuscript (citation numbers, CI intervals, figure refs, math notation, the DOI, `[software]` reference tags) was individually inspected and is legitimate. |
| 9 | Wording around AI assistance | `AI_ASSISTANCE_DISCLOSURE.md`, the manuscript's Author Contributions and AI assistance disclosure sections, and `FINAL_SCIENTIFIC_AUDIT.md`'s addendum were all read directly. Wording is specific, professional, does not overstate or understate AI's role (it names exact files/tests Claude wrote, and separately names a case where an AI-generated audit conclusion was wrong and human-corrected), and is internally consistent across all three documents. | ALREADY FIXED | — | None needed — this was already tightened in the prior audit's post-audit addendum and is consistent. |
| 10 | Reposition the paper toward diagnostic/failure-mode analysis | Checked against the audit brief's own test: (1) primary research question, stated in §1: "does routing keep the accuracy of PB-everywhere, and where does the compute go?"; (2) the experiments (Table 1, eight cells plus chr20) test exactly that; (3) main contribution is the frozen-constant validation plus the failure catalogue found while running it; (4) the title is **already** "...frozen-constant validation on GIAB samples, **failure modes, and the compute bottleneck**" — failure modes and the compute bottleneck are already in the title, not buried. | INCORRECT (as a change request) | — | **No change made.** The paper is already positioned exactly as Gemini suggests it should be; repositioning it further would not change anything real, only reword an already-accurate title/abstract. |
| 11 | Tighten Limitations | `§5` Table 10 has 17 rows, each with a `Restricts` column naming exactly which claims it bounds (scope, independence, generalisation, statistics, timing, external comparison, reproducibility, clinical use). Cross-checked against the Discussion and Conclusion for redundant or missing items — none found missing, and it is already the single most exhaustive limitations table in the document set (denser than typical Bioinformatics limitations sections). | ALREADY FIXED | — | None needed. Further "tightening" would mean shortening it, which would reduce disclosure, not improve it. |
| 12 | Make performance claims more precise | Table 8 already carries `Level` (function/pipeline/measured/projection) and `Conditions` (exact stage, region, depth, hardware, n) columns for every one of its 17 rows. §2.7 defines the four levels explicitly before any number is given. | ALREADY FIXED | — | None needed. |

**Net effect of this pass:** one stale sentence removed; zero scientific numbers, thresholds, or
verdicts touched; zero new experiments run.

---

## Changes made

1. `research/PAPER_MANUSCRIPT_V2.md` line 3: removed the clause *"Placeholders in square brackets are
   unresolved on purpose"* from the epigraph (no longer true — the one real placeholder was closed by
   the prior audit, and an exhaustive re-scan of this pass found none remaining). Kept the adjacent RO
   ("record-only") legend sentence, which is still accurate and load-bearing for reading the tables.
2. Regenerated `research/latex/paper.tex` from the corrected Markdown via `md2tex.py` (0 unmapped
   citations).
3. Regenerated `research/PAPER_MANUSCRIPT_V2.html` via `md2html.py`.
4. Rebuilt `research/latex/paper.pdf` with two `pdflatex` passes (24 pages, 0 undefined references,
   only cosmetic duplicate-destination/hbox warnings that predate this pass and are unrelated to the
   edit).
5. Published `research/latex/paper_final_preprint.pdf` as a copy of the rebuilt PDF. The prior
   `research/latex/paper_final_audited.pdf` and `research/latex/paper_pre_audit_original.pdf` were
   left untouched, per instructions.

No experiment was re-run. No threshold, verdict, table number, or figure value was changed. No file
outside the manuscript/build pipeline was edited in this pass.

---

## Remaining limitations (stated plainly)

These are unchanged from `FINAL_SCIENTIFIC_AUDIT.md` and are restated here because this pass did not
resolve them — it only checked whether they were still accurately disclosed, and they are:

- **HG005 `false_rescue_count` (9) vs. Table 5 `router_fp` HG005 contribution (2, record-only) is still
  unreconciled.** Both numbers stand, both are disclosed as existing side by side (in
  `LIMITATIONS_AND_OPEN_QUESTIONS.md` item 17, added by the prior audit), and neither is asserted to be
  wrong. Resolving which definition each uses would require re-deriving the classification from the
  surviving HG005 cache — not attempted here, consistent with "do not run new experiments to make the
  paper look better."
- chr20 (v20) is chromosome-separated from the constant-fitting data but not sample/library/depth
  -independent; this is stated three times in the manuscript (Abstract, §3.3, Table 10) and once more
  in this document, deliberately, because it is the single limitation most likely to be probed first by
  a mentor or reviewer.
- The native counts backend's byte-identical equivalence gate on chr20 covers 980,928 of 56,266,816
  loci (1.7 %); it is not claimed bit-exact on adversarial synthetic BAMs with repeated read names.
- The chr20 DEGRADED verdict rests on 10 discordant loci out of 56.2 million — statistically real
  (CI excludes zero, corroborated by block bootstrap and McNemar), but the manuscript correctly declines
  to say whether an effect this small "matters" for any given use.
- No genome-scale run exists for any sample. All conclusions are scoped to the tested regions/depths.
- HG005 is the project's only unrelated (non-trio) sample, and it is a single 3-Mb region with lost raw
  artefacts. Generalisation across ancestry/samples is not established.
- No lockfile; environment reproducibility rests on the version table in `REPRODUCIBILITY.md`/README.

None of these were hidden, softened, or newly discovered by this pass — they were already disclosed by
the manuscript and the prior audit. This pass's job was to check that the disclosure was still accurate
and complete, which it was.

---

## New experiments required

**None for mentor-review readiness.** The one unresolved item that *could* be closed with new work
— the HG005 `false_rescue_count` vs. `router_fp` discrepancy — would require re-deriving a
classification from a surviving cache and does not affect any pooled or chr20 conclusion; it remains an
open question in `LIMITATIONS_AND_OPEN_QUESTIONS.md`, not a blocker.

---

## Scientific claims audited

- Router routes 0.039–0.108 % of loci in regional experiments, 0.137 % on chr20 — spot-checked against
  Table 3/Fig. 4 source values, consistent.
- Pooled ΔF1 range (−5.6×10⁻⁵ to +2.7×10⁻³), verdict counts (2 PRESERVED, 5 IMPROVED, 1 DEGRADED, 1
  experiment-level NEGATIVE) — consistent across Table 3, §3.1, Discussion, Conclusion, and Fig. 2.
- chr20 P/R/F1 table, ΔF1, both CIs (locus and 50-kb block), McNemar p = 0.021 on 10 events — consistent
  between §3.3 text, Table 3/§3.3 table, and Fig. 2a's DEGRADED marker.
- True-SNP loss count (9 pre-chr20, 10 with chr20) — consistent across §3.4.1, Table 5, Fig. 3, Abstract,
  Conclusion; the "9 vs 10" reconciliation is stated, not silently resolved into one number.
- Genotype layer: Method C 0.9888 GT accuracy on hold-out, 0 GT/GQ mismatches on the vectorised
  reimplementation, hap.py F1 lift 0.624 → 0.955 — consistent between §3.5, Table 7, Abstract.
- Performance: function-level (21–65×, 26.9–35.1×) vs. pipeline-level (2.72×, 7.11×) vs. the one measured
  BAM→calls figure (1.65×) vs. projections (171–396×, explicitly excluded from the figure) — consistent
  across Table 8, §3.6, Fig. 5, Discussion, README.
- External caller comparison (Table 9/Fig. 6): AI F1 0.9546–0.9547 vs. DeepVariant 0.9494, Clair3 0.9290,
  GATK 0.9682, with the constant-fitting-span overlap and "no ranking implied" stated in the table, the
  figure caption, and the running text — consistent, and present in all three places independently
  (not just inherited from one).

No numeric discrepancy was found between the manuscript, its cited CSVs (`CHR20_RESULTS.csv`,
`RESCUE_ANALYSIS.csv`, `CLAIM_EVIDENCE_MAP.csv`), and the figures, beyond the already-disclosed HG005
reconciliation gap above.

---

## Reproducibility status

- **Directly reproducible from this repository:** the full test suite (`native/build.sh` then pytest;
  previously re-run live and matched exactly per the prior audit); the frozen thresholds and their
  SHA-256-verified source files; the chr20 run's inputs (GIAB checksums, `FROZEN_MANIFEST.json`); the
  regional experiments' caches.
- **Record-only, not independently reproducible from this repository:** all HG005 post-fix numbers
  (hap.py output, rescue/error CSVs, the counts-only 1.22× historical pipeline figure), two of the
  original native-backend validation reports, and the two HG005 pipeline scripts referenced by
  `REPRODUCIBILITY.md` — all lost in the disclosed 2026-09-20 `/tmp` wipe, all explicitly tagged `(RO)`
  wherever they are used.
- **Not re-executed in this pass:** the full chr20 extract→pipeline→hap.py→verdict chain (traceable to
  frozen, checksummed outputs; README itself discloses it "was not re-run end to end for this README").
  This pass rebuilt only the manuscript→LaTeX→PDF pipeline, which is fully deterministic from
  `PAPER_MANUSCRIPT_V2.md`.

---

## Final reviewer test (strict Bioinformatics reviewer posture)

**A. Central contribution.** An honest, pre-registered empirical test of whether confidence-based
routing between a cheap binomial SNP screen and an expensive Poisson-binomial model preserves accuracy
— together with a disclosed catalogue of every failure mode found while testing it (a NEGATIVE
experiment, a DEGRADED chromosome, 10 lost true SNPs, an evaluator artefact, and a compute story that
turned out mostly negative).

**B. Strongest evidence supporting it.** The chr20 whole-chromosome run: 56.2 M loci, matched verdicts
from two independent evaluators (internal, hap.py), corroborated by both a locus bootstrap and a 50-kb
block bootstrap and an exact McNemar test, reported as DEGRADED with its exact effect size rather than
rounded away.

**C. Strongest evidence weakening it.** chr20 is not sample/library/depth-independent from the
constant-fitting data; the cascade demonstrably loses real true SNPs (10, at VAF≈0.12, a mechanistically
coherent and disclosed pattern); the only unrelated sample (HG005) is a single 3-Mb record-only region;
and the entire compute-saving premise is undercut by PB now being 93–96 % of runtime with only one
1.65×, 0.5-Mb measurement of any end-to-end gain.

**D. Fatal methodological flaw?** None found. The paper's willingness to report v14 NEGATIVE and v20
DEGRADED against its own pre-registered rule, rather than reframe or drop them, is the opposite of the
usual failure mode this kind of review looks for.

**E. Unresolved major issue?** The HG005 `false_rescue_count`-vs-`router_fp` discrepancy remains open
(disclosed, not resolved); the native counts backend's byte-identical gate on chr20 covers only 1.7 % of
its loci, leaving most of the chromosome's extraction correctness argued by the paired-comparison logic
in §3.6 rather than independently checked cell-by-cell.

**F/G. Riskiest wording / claim most likely to be asked to soften.** Not a factual overclaim, but an
ordering risk: the Abstract states the 21–65× and 26.9–35.1× native speedups in the same paragraph as the
93–96 % PB-dominance and 1.65× figures. A reader who stops at the first sentence could walk away with
the wrong headline. The text immediately corrects this, and Table 8's `Level` column and Fig. 5 make it
unambiguous on a second read — but a reviewer may still ask for the negative compute finding to be
foregrounded ahead of the multiplier in the Abstract's sentence order.

**H. What a reviewer will most likely ask to verify further.** (1) the HG005 reconciliation; (2) whether
native-backend correctness holds beyond the tested 1.7 % of chr20; (3) a genome-scale run before any
claim broader than chromosome-scale; (4) a second unrelated, non-batch sample beyond the 3-Mb HG005
region before any generalisation claim.

No FATAL or MAJOR issue survives this pass that was not already disclosed by the manuscript itself.

---

## FINAL VERDICT

### SCIENTIFICALLY READY FOR MENTOR REVIEW

Unchanged from the prior audit's verdict. This pass found the manuscript already implements nearly
everything the external review asked for — in most cases with more precision than the request itself
— made the one genuine fix it found (a stale placeholder-claiming sentence), reconfirmed the one
genuinely open item (HG005 rescue-count reconciliation) is still accurately disclosed as open, and
introduced no new claim, number, or scope change.

`SUBMISSION-READY` is not used: no specific journal or conference submission checklist was checked
against.
