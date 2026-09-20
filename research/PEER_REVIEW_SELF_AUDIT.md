> **SUPERSEDED IN PART (2026-09-20).** Historical self-audit (pre-chr20). Current status: `FINAL_FORENSIC_AUDIT.md` and `FINAL_FREEZE_REPORT.md`.

# Peer Review Self-Audit

Eighteen adversarial questions, answered against the actual evidence in this
repository, not against the project's own framing of itself.

---

**1. Could the benchmark be leaking information?**
Issue: labels in `providers.py::_fetch_labels` come from the truth VCF; if
detection thresholds were ever fit using these labels, the benchmark would
be circular.
Evidence: `cascade.py`'s three frozen constants and Method C's formula are
verified byte-identical (sha256) across every experiment in this dossier,
including before and after this session's HG005 fix. `test_bench_v12.py`
contains `test_no_routing_function_sees_labels`, asserting `"truth" not in
inspect.signature(fn).parameters` for the routing functions.
Severity: LOW (checked, not merely assumed).
Resolution: verified by source inspection and an existing test.
Remaining limitation: this check confirms the *functions* don't take truth
as an argument; it does not prove a human never looked at an intermediate
result and informally adjusted a *different* file's threshold before this
dossier's session began. No such adjustment was found in the audited
history, but absence of evidence is not proof of absence for work predating
this session's direct observation.

**2. Were thresholds trained on test data?**
Issue: same class of concern as #1, specifically for train/test split
integrity.
Evidence: `FROZEN_BINOMIAL_THRESHOLD`/`FROZEN_PB_THRESHOLD`/`FROZEN_ROUTER_CUTOFF`
are named "frozen" throughout the codebase and are literal constants, not
computed from any data-loading call in `cascade.py`.
Severity: LOW.
Resolution: source inspection; no fitting code found anywhere near these
constants' definitions.
Remaining limitation: no formal train/validation/test region split registry
exists outside individual devlog documents; this dossier trusts each
experiment's own `independence_audit.json`/`region_selection.json`
provenance rather than an independently-audited master split file.

**3. Were any thresholds changed after seeing HG005?**
Issue: this session specifically investigated and fixed an extraction bug
after seeing a bad HG005 result — did that investigation touch any
threshold?
Evidence: `cascade.py` sha256 `b0ee9f4b24fe...` is identical before and
after this session's changes (verified directly in this session, not
assumed). Only `extract_bench_v12.py` and `pileup_counts.py` were modified,
and the diff (61 lines total) is coordinate-bookkeeping only — no numeric
threshold appears in either diff.
Severity: LOW.
Resolution: verified directly via `git diff` and checksum comparison in this
session.
Remaining limitation: none identified.

**4. Are samples independent?**
Issue: HG002/HG003/HG004 are a trio (son/father/mother) and share
substantial genomic identity by descent; treating them as 3 independent
generalization tests overstates independence.
Evidence: stated explicitly in `LIMITATIONS.md`.
Severity: MEDIUM.
Resolution: framed correctly in this dossier — "cross-sample within a trio,"
not "3 independent samples" — and HG005 (unrelated) is called out as the
only genuinely independent additional sample.
Remaining limitation: only one genuinely unrelated sample (HG005), and only
one 3 Mb region of it, has been evaluated with the genotype-aware evaluator.

**5. Are regions representative?**
Issue: hand- or script-selected regions could be unrepresentatively easy.
Evidence: regions from `bench_v14` onward were selected by a dedicated
script from public annotation *before* any caller ran (`independence_audit.json`
per experiment); v14 deliberately included a segdup region (`chr1_segdup`)
and v15/v17 explicitly targeted "hardest unseen contig" / difficult strata.
Severity: LOW-MEDIUM.
Resolution: selection protocol is pre-registered and includes deliberately
hard regions, not only easy ones; v15's mixed/weak result and v17's harder
performance profile are evidence the selection was not cherry-picked toward
favorable outcomes.
Remaining limitation: 12 Mb (chr21 pool) is the largest single tested span;
representativeness at genome scale remains untested.

**6. Are truth sets appropriate?**
Issue: is NIST v4.2.1 the right, current truth standard for each sample?
Evidence: v4.2.1 used consistently across all 4 samples; the same version
used for GIAB HG002/3/4/5 releases at the time this work was done.
Severity: LOW.
Resolution: consistent truth-set version across every experiment in this
dossier.
Remaining limitation: none identified within the version's own scope; a
newer truth-set release (if one now exists) has not been evaluated.

**7. Is hap.py used consistently?**
Issue: mixing hap.py configurations (flags, region scope) across
experiments would make F1 numbers non-comparable.
Evidence: identical hap.py image digest
(`sha256:d63b963a6cb01b4830393b22369e7b91d298e4156dde353739e74e4cfa4f96d0`)
and flag set (`xcmp`, `preprocessing_leftshift=true`,
`preprocessing_truth_confregions=true`, `window=50`) across HG002 chr21 and
HG005 runs. **However**, this session discovered a real scope inconsistency
mid-task: the first HG005 hap.py run reused the 24Mb-scope truth/BED (matching
the pre-fix run's config) against a 3Mb-scope query VCF, producing an
artificially deflated TRUTH.TOTAL=28,435/recall — caught, diagnosed, and
corrected with a region-scoped rerun before being reported as the primary
number (`HG005_POSTFIX_STRESS_TEST_REPORT.md` Sec. 4).
Severity: MEDIUM (a real inconsistency existed and was self-caught within
this same session, not by an external reviewer).
Resolution: region-scoped truth/BED rerun performed; both configurations
reported, with the flawed one explicitly labeled `unrestricted...DIAGNOSTIC
ONLY, DO NOT use as detection quality`.
Remaining limitation: this specific class of mistake (truth/query region
mismatch) should be checked as a standing gate on any future regional
benchmark, not just caught ad hoc.

**8. Are runtime measurements fair?**
Issue: comparing runtimes across different hardware, configurations, or
partial pipelines would be misleading.
Evidence: `PERFORMANCE_RESULTS.csv` labels every row `MEASURED` or
`EXTRAPOLATION`; the external-caller runtimes (DeepVariant 819s, Clair3
1342s, GATK 287s) come from `experimental/head_to_head/HEAD_TO_HEAD_RESULTS.csv`,
measured on the same region as their accuracy figures but the AI pipeline's
own runtime at that exact scope was never captured.
Severity: MEDIUM.
Resolution: `SCIENTIFIC_CLAIMS.md` explicitly states no runtime comparison
against any external caller currently exists, rather than implying one from
mismatched-scope numbers.
Remaining limitation: a same-region, same-hardware runtime comparison
against DeepVariant/Clair3/GATK does not exist and would need a new
experiment.

**9. Are genotype and allele metrics being confused?**
Issue: conflating GT-blind allele-level accuracy with genotype-aware hap.py
accuracy was the exact bug that originally produced "F1 0.974 -> 0.624"
confusion (`experimental/unified_happy/ROOT_CAUSE_ANALYSIS.md`).
Evidence: that confusion is now the subject of an entire prior investigation
whose conclusion (two evaluators, never interchangeable) is enforced by an
explicit `evaluator` column in every table in this dossier.
Severity: LOW (previously HIGH, resolved with a permanent labeling
convention).
Resolution: `evaluator` column in `BENCHMARK_MASTER.csv`;
`PIPELINE_INTEGRITY_AUDIT.md` M-2 documents the underlying source-code risk
that this labeling convention compensates for without fixing at the API
level.
Remaining limitation: the underlying source-code ambiguity (M-2) is not
fixed; a future contributor working directly from `robustness_benchmark.py`
without this dossier could still make the same mistake.

**10. Is any performance number extrapolated improperly?**
Issue: presenting a projection as a measurement.
Evidence: whole-genome time figures in this dossier are explicitly labeled
`EXTRAPOLATION`. The "~65x C/htslib" claim was investigated this session
specifically because it looked like exactly this failure mode — and this
investigation initially went wrong in the opposite direction: an incomplete
search concluded the claim was unsubstantiated, when a real, measured 65.1x
result existed (in a worktree not yet checked) and is now independently
re-verified and integrated. The corrected, final state: 21.0x-65.1x are both
real, measured, extraction-stage-only numbers (not whole-pipeline, not
extrapolated) — see `PIPELINE_INTEGRITY_AUDIT.md` H-1 RESOLVED and
`SCIENTIFIC_CLAIMS.md`.
Severity: HIGH (a wrong conclusion was nearly published in both directions —
first accepting an unverified number, then wrongly rejecting a real one).
Resolution: independently re-measured and re-verified in this session rather
than trusted either way secondhand.
Remaining limitation: the full-pipeline (not extraction-stage-only) speedup
is smaller than 65x because `load_reads`/PB are unaccelerated — this
distinction must always travel with the 65x figure, or it becomes exactly
the "projection presented as measurement" failure mode this question warns
against, just from the opposite direction.

**11. Are there hidden extraction artifacts?**
Issue: this is the central question this whole session investigated.
Evidence: the CRITICAL finding (C-1, fabricated coordinates) was found, root-
caused, fixed, and independently verified (full-population 0-mismatch
reference check). A second, smaller artifact (M-1, window-vs-SNP BED
granularity) was found *while investigating the first* and is reported as a
distinct, unfixed limitation rather than folded into "recall" silently.
Severity: CRITICAL (C-1, now fixed) / MEDIUM (M-1, open).
Resolution: see `PIPELINE_INTEGRITY_AUDIT.md`.
Remaining limitation: M-1 is unresolved and depresses measured recall by an
unknown, region-dependent amount across every experiment in this dossier
that uses window-tiled BED filtering (i.e. all of them).

**12. Are failures being omitted?**
Issue: selective reporting of only favorable results.
Evidence: v15's mixed/weak verdict is reported in `BENCHMARK_MASTER.csv` and
`RESULTS.md` exactly as its own devlog recorded it (`WEAK POSITIVE` /
`NEGATIVE` depending on acceptance rule), not smoothed into "IMPROVED." The
original invalidated HG005 F1=0.0085 result is preserved, not deleted, and
explicitly cross-referenced as superseded rather than erased.
Severity: LOW.
Resolution: both examples above are direct counter-evidence against omission.
Remaining limitation: this audit cannot rule out omission of an experiment
that was run and never written up anywhere in the repository at all — it can
only confirm that every experiment *found* in the repository is reported
faithfully.

**13. Is the C implementation actually equivalent?**
Issue: the original task prompt asked for exactly this proof.
Evidence: **initially answered wrong.** An earlier pass of this self-audit
claimed no C implementation existed, based on checking only one of several
relevant worktrees; the user identified the correct one
(`.claude/worktrees/agent-a4cedfdeffb76e3dd/experimental/performance/`),
which contains a real C extension (`native/pileup_native.c`) with its own
prior full-region equivalence proof (0/109,813,120 cell mismatches). This
session independently re-verified that proof rather than trusting it
secondhand: rebuilt the extension from source in a fresh worktree, and
re-ran a bit-exact comparison via the actual production `load_counts()`
entry point (native ON vs OFF) on both HG002 and HG005 data — 0 mismatches
across every check, including real paired-end-read edge cases (overlapping
mates, orphans, base-quality filtering) the original proof's synthetic-data
predecessor could not have exercised.
Severity: was HIGH (a wrong answer to a direct equivalence question);
resolved.
Resolution: `PIPELINE_INTEGRITY_AUDIT.md` H-1 RESOLVED;
`experimental/performance/CHTSLIB_INTEGRATION_REPORT.md`;
`test_native_pileup_equivalence.py` (permanent regression suite).
Remaining limitation: the equivalence contract is proven for the exact
filter configuration currently in use (`min_mapping_quality=20`,
`min_base_quality=13`, pysam's `stepper="samtools"` defaults); a future
caller changing these filter values would need to re-verify equivalence
before trusting the native path, since the C code hard-codes this specific
configuration rather than generically parameterizing every pysam `pileup()`
option (documented in the source worktree's own report, Section 13,
inherited here unchanged).

**14. Could the router result be explained by chance?**
Issue: small ΔF1 deltas could be noise.
Evidence: v14 reports a paired bootstrap (10,000 resamples) with
`delta_f1_ci95=[0.000075, 0.00109]`, excluding 0, p=0.0228 — a real, if
small, effect. Other cells (v16: -0.00003) are consistent with pure noise
around zero and are reported as `PRESERVED`, not misrepresented as
`IMPROVED`.
Severity: LOW.
Resolution: bootstrap CIs are used and reported where computed; the language
used per-cell (`IMPROVED` vs `PRESERVED`) tracks the CI, not just the point
estimate's sign.
Remaining limitation: no bootstrap CI was computed for the HG005 hap.py
result (single point estimate on one region) — flagged explicitly in
`LIMITATIONS.md`.

**15. Are confidence intervals appropriate?**
Issue: bootstrap methodology validity (resampling unit, independence
assumption).
Evidence: paired bootstrap over scored loci, fixed seed, 10,000 resamples —
standard for this kind of paired-classifier comparison.
Severity: LOW-MEDIUM.
Resolution: methodology documented in `METHODS.md` Sec. 5.
Remaining limitation: resampling individual loci as independent units is a
simplification (adjacent loci in a genome are not independent due to
linkage/local coverage correlation); no block-bootstrap or other
correlation-aware method was used. This could make reported CIs somewhat
too narrow.

**16. Is the sample size sufficient?**
Issue: statistical power, especially for stratified sub-analyses.
Evidence: `tandemrepeats` FN stratum on HG005 has n=6 and is explicitly
labeled "underpowered — do not over-interpret" in both
`HG005_POSTFIX_STRESS_TEST_REPORT.md` and this dossier's `ERROR_ANALYSIS.csv`.
Severity: LOW (correctly flagged, not hidden).
Resolution: explicit underpowered-stratum labeling throughout.
Remaining limitation: several other strata (e.g. `lowmap_segdup` FP n=60) are
adequately powered for a point estimate but not large enough for a tight CI;
no CI was computed for any stratification cell in this dossier.

**17. Are external caller comparisons fair?**
Issue: version pinning, model choice, region restriction, BQSR omission.
Evidence: `HEAD_TO_HEAD_RESULTS.csv`'s own notes column documents exactly
this kind of caveat per caller (GATK ran without BQSR "no known-sites
resource available for this region-only subset — documented protocol
deviation"; Clair3's `:latest` tag silently fell back to pileup-only and had
to be retried with a pinned version).
Severity: MEDIUM (real caveats exist and are disclosed, not hidden).
Resolution: caveats preserved verbatim in this dossier's
`BENCHMARK_MASTER.csv` notes column.
Remaining limitation: the GATK BQSR omission specifically could
systematically favor or disfavor GATK's reported numbers relative to a
full-pipeline GATK run; this was not corrected or quantified, only
disclosed.

**18. Are claims stronger than evidence?**
Issue: the central purpose of `SCIENTIFIC_CLAIMS.md`.
Evidence: every claim in that file is explicitly scoped and classified;
"universally faster," "clinically ready," "genome-wide validated," and
"universally superior to external callers" are all marked `NOT SUPPORTED`
or `NOT YET TESTED` rather than implied by omission.
Severity: LOW (actively managed).
Resolution: `SCIENTIFIC_CLAIMS.md`.
Remaining limitation: this self-audit cannot fully substitute for an
external reviewer independently trying to break each claim; it is a
structured self-check, not independent peer review (see `SCIENTIFIC_CLAIMS.md`
final entry on this exact point).

---

## Summary judgment

Of 18 questions: **1 CRITICAL-class issue found and fixed this session**
(extraction coordinate bug, question 11), **1 real methodological mistake
found and corrected within this same session** (question 7, HG005 truth-scope
mismatch), and **1 performance-claim investigation that itself went wrong
before being corrected** (questions 10 and 13: a real C/htslib
implementation and its genuine 65x result were first wrongly dismissed as
unsubstantiated, based on an incomplete search, then located, independently
re-verified bit-exact, and integrated into production once the error was
caught) — plus several genuine, unresolved, explicitly-disclosed limitations
(questions 4, 8, 11 [M-1], 15, 16, 17) that constrain how strongly any
result in this dossier should be read. No fabricated equivalence, citation,
or benchmark data was produced at any point; where this self-audit was
itself wrong, that error is disclosed rather than quietly patched over.
