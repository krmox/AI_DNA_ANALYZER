> **SUPERSEDED IN PART (2026-09-20).** Superseded in part: the ~65× figure is a `load_counts` function-level speedup on 12 Mb (reused measurement), not a pipeline speedup; `load_reads` is now also native (26.9×/35.1× function-level; pipeline 2.72×, ≤ 7.1× with 4 processes); chr20 (v20) added a strict DEGRADED verdict of extremely small magnitude. See `PAPER_MANUSCRIPT.md` §3.1, §3.8.

# Scientific Claims — Classified

Every claim this project could plausibly make, classified as
**SUPPORTED**, **PARTIALLY SUPPORTED**, **NOT SUPPORTED**, or **NOT YET
TESTED**, with the evidence it rests on. No claim here should be quoted
without its classification and scope qualifier.

## Architecture and detection accuracy

**"The frozen Binomial -> Router -> Poisson-Binomial cascade preserves or
improves PB-only accuracy while routing a small fraction of loci."**
→ **SUPPORTED**, within tested scope. 6/7 internal-evaluator experiments
show cascade F1 >= PB-only F1 (5 with a statistically real positive delta,
1 statistically indistinguishable); 1/7 (v15) shows a mixed/weak result on
its own held-out test set. Two hap.py-evaluated experiments (HG002 chr21,
HG005 chr1 3Mb) both show cascade and PB-only statistically indistinguishable.
Scope: GIAB HG002/HG003/HG004/HG005, Illumina short-read, GRCh38, SNP-only,
the specific regions tested.

**"Method C's posterior genotype layer substantially improves hap.py's
genotype-aware F1 over a naive forced-heterozygous VCF writer."**
→ **SUPPORTED**. Direct A/B measurement on identical HG002 chr21:32-44M
detection calls: F1 0.6238 (forced 0/1) -> 0.9547 (Method C), ΔF1=+0.331,
with the mechanism (truth-homalt genotype mismatch) independently confirmed
per-record from hap.py's own BD/BK tags.

**"The AI cascade generalizes to a new, non-trio GIAB sample (HG005)."**
→ **PARTIALLY SUPPORTED**. F1≈0.952 on one 3 Mb HG005 region is a reasonable
result, close to (slightly below) the HG002 chr21 Method-C figure of
≈0.9546 — but this is one region on one additional sample, not a systematic
cross-sample study at the scale of the HG002 experiments. Do not extrapolate
to "the model generalizes to arbitrary new samples."

**"AI_DNA_ANALYZER's extraction pipeline had a data-integrity bug that
invalidated the original HG005 result, which has now been fixed and
independently verified."**
→ **SUPPORTED**. Full-population 0-mismatch reference check, byte-identical
determinism, 0 VCF REF mismatches, and a historical audit proving no other
published result was affected. See
`experimental/stress_test/HG005_EXTRACTION_BUG_ROOT_CAUSE.md`.

## Performance

**"AI_DNA_ANALYZER's extraction achieves ~65x speedup via a C/htslib
implementation."**
→ **SUPPORTED, for the extraction stage specifically, now integrated into
production.** Corrected after an earlier version of this file wrongly
marked it NOT SUPPORTED (an audit error — the C extension existed in a
worktree that was not checked at the time; see
`PIPELINE_INTEGRITY_AUDIT.md`'s "Correction" section for the full account).
The C extension (`native/pileup_native.c`) is proven bit-exact against the
pure-Python reference (0 mismatches, both the original 12Mb HG002 proof and
this session's independent re-verification on HG002 and HG005 data), is now
the default path in `pileup_counts.load_counts`, and the ~65x figure is a
real, measured, full-region result — see `PERFORMANCE_RESULTS.csv` for the
exact numbers at both the original 12Mb benchmark scale and the HG005 3Mb
scale. **Scope qualifier that must travel with this claim:** 65x describes
the *pileup-counting extraction stage* specifically (`load_counts`), not
automatically the full BAM-to-VCF pipeline — `read_level_pileup.load_reads`
(feeds the Poisson-Binomial caller) is a separate extractor, not touched or
accelerated by this work, and dominates wall time in a full run that
includes PB scoring (see the next claim below). Do not quote "65x" as an
end-to-end pipeline number without checking which pipeline configuration is
being described.

**"The AI pipeline's full BAM-to-VCF wall time is 65x faster than before."**
→ **NOT SUPPORTED as a blanket claim; PARTIALLY SUPPORTED for one specific
configuration.** The original C/htslib report's own "65.1x end-to-end"
figure was measured for a pipeline that reused pre-computed, frozen
`pb_calls`/`cascade_calls` masks rather than recomputing them from
`read_level_pileup.load_reads` — i.e. it is end-to-end for "counts →
REF/ALT/VCF given already-known calls," not for a run that also performs
PB-caller read extraction. In a run that includes both stages (e.g. this
session's HG005 extraction, `counts_seconds=343.65s` vs
`reads_seconds=696.05s`), `load_reads` is the larger cost and is
unaccelerated, so the true full-pipeline speedup is smaller than 65x. See
`PERFORMANCE_RESULTS.csv` for both figures, clearly separated.

**"Process-level parallel chunking accelerates extraction with proven
correctness equivalence."**
→ **SUPPORTED**, narrowly. Measured up to 3.1x at 8 workers on a 0.6 Mb
region, with an explicit correctness gate (counts/labels/LLR/routing/calls
byte-identical parallel vs serial) that passed. Not integrated into the
active production pipeline; not measured beyond a 2 Mb region; not a claim
about any other acceleration mechanism.

**"AI_DNA_ANALYZER is faster than DeepVariant/Clair3/GATK."**
→ **NOT TESTED as a general claim; PARTIALLY reportable for one configuration.**
No timing was captured for the AI pipeline at the HG002 chr21:32-44M scope
where DeepVariant/Clair3/GATK runtimes ARE recorded (819s/1342s/287s). The
AI pipeline's own extraction+calling wall time was only measured at the
HG005 3 Mb scale (1392s) under different region/sample/scope conditions.
**No apples-to-apples runtime comparison against any external caller
currently exists.** Do not claim a speed advantage or disadvantage versus
any external caller.

**"AI_DNA_ANALYZER works genome-wide."**
→ **NOT YET TESTED.** No whole-genome run has ever been attempted. The
whole-genome time figure in `PERFORMANCE_RESULTS.csv` is a labeled linear
extrapolation, not a measurement, and should not be cited as a performance
claim.

## External comparison

**"AI Cascade + Method C's hap.py F1 exceeds DeepVariant's at the tested
HG002 chr21:32-44M region."**
→ **SUPPORTED, narrowly** (0.9546 vs 0.9494, one region/sample/depth,
driven by higher precision / lower recall — a trade-off, not a strict
dominance).

**"AI Cascade + Method C is universally superior, or superior to GATK."**
→ **NOT SUPPORTED.** GATK's F1 (0.9682) exceeds the AI cascade's (0.9546)
at the same region. Do not claim overall superiority over any external
caller in any general sense.

## Readiness claims

**"AI_DNA_ANALYZER is clinically validated / clinically ready."**
→ **NOT SUPPORTED.** No claim of clinical validity is made anywhere; this
is a research-stage variant-calling architecture study on GIAB reference
materials only.

**"AI_DNA_ANALYZER has been independently peer-reviewed."**
→ **NOT SUPPORTED.** `PEER_REVIEW_SELF_AUDIT.md` in this dossier is a
self-conducted adversarial review, explicitly not a substitute for external
peer review.

**"The router/threshold values were derived without seeing HG005 or any
other reported test result."**
→ **SUPPORTED.** The three frozen constants (`FROZEN_BINOMIAL_THRESHOLD=7.0`,
`FROZEN_PB_THRESHOLD=10.5`, `FROZEN_ROUTER_CUTOFF=5.411872376933351`) and
Method C's formula are verified byte-identical (sha256) across every
experiment in this dossier, including before and after the HG005 extraction
fix session, and no script in the audited scope recomputes or refits them
from any test-set outcome (`PIPELINE_INTEGRITY_AUDIT.md` I-3).
