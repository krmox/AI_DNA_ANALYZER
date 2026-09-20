> **SUPERSEDED IN PART (2026-09-20).** Written before chr20 and `load_reads`; do not use its accuracy or performance statements without checking `PAPER_MANUSCRIPT.md` (no global "preserved" claim; chr20 DEGRADED at extremely small magnitude).

# Outreach Package

Tone note: nothing below uses superlatives ("revolutionary", "state-of-the-art",
"best") that are not independently demonstrated. Every number quoted matches
`BENCHMARK_MASTER.csv` / `SCIENTIFIC_CLAIMS.md` exactly.

## A. Research profile paragraph

I am an independent researcher investigating cost/accuracy trade-offs in
short-read variant calling. Over this project I built a frozen Binomial ->
confidence-router -> Poisson-Binomial cascade with a closed-form genotype
posterior layer, and validated it against GIAB reference materials
(HG002/HG003/HG004/HG005) using both an internal evaluator (9 pre-registered
region/depth/sample cells) and the GA4GH-standard hap.py comparison engine
(2 regions). Along the way I found and fixed a coordinate-bookkeeping defect
in my own extraction pipeline that had silently invalidated an earlier
benchmark, and located and integrated a previously-validated but unmerged
C/htslib extraction accelerator (21-65x speedup on the extraction stage,
re-verified bit-exact against the Python reference before trusting it). I
used AI coding assistance
(Claude) extensively during implementation and documentation, disclosed in
full in `AI_ASSISTANCE_DISCLOSURE.md`; the research questions, experimental
design, and validity judgments are my own.

## B. Professor email template

> Subject: Independent research on confidence-routed variant calling — GIAB-validated, seeking feedback
>
> Dear Professor [NAME],
>
> I've been running an independent research project on short-read SNP
> calling: a frozen Binomial screening stage gated by a confidence router
> to a more expensive Poisson-Binomial caller, plus a closed-form genotype
> posterior layer, validated against GIAB HG002/HG003/HG004/HG005 with both
> an internal evaluator and hap.py (GA4GH-standard). The cascade preserves
> detection accuracy relative to running the expensive caller everywhere
> (ΔF1 within noise across the tested regions/samples) while routing under
> 0.1% of loci to it; the genotype layer recovers +0.33 F1 over a naive
> baseline. I also found and fixed a coordinate bug that had invalidated an
> earlier version of one benchmark, and integrated a validated C/htslib
> extraction accelerator (21-65x on the extraction stage, re-verified
> bit-exact before trusting it) into the production pipeline — documented in
> full in the project's audit and limitations files, alongside what the
> results do *not* support (no whole-genome or clinical claim is made).
>
> I'm reaching out because your work on [SPECIFIC PROFESSOR RESEARCH AREA —
> fill in] overlaps with [SPECIFIC OVERLAP — fill in, e.g. confidence
> calibration in genomic pipelines / GIAB benchmarking methodology /
> computational cost trade-offs in bioinformatics]. I'd value a short
> conversation about whether this approach or its limitations connect to
> problems your group is working on, and whether there's an opportunity
> for mentorship or a research discussion.
>
> Full write-up, code, and reproducibility details:
> [PROJECT LINK — fill in]
>
> Thank you for your time,
> [NAME]

## C. Research project description — four lengths

**50 words:**
> Independent research on short-read SNP calling: a frozen Binomial screen
> gated by a confidence router to a Poisson-Binomial caller, plus a
> closed-form genotype posterior. Validated on GIAB HG002-HG005 with an
> internal evaluator and hap.py. Detection accuracy preserved; genotype
> accuracy improved +0.33 F1 over a naive baseline.

**100 words:**
> I designed and validated a variant-calling cascade that routes only
> confidence-ambiguous loci from a cheap Binomial screen to a more expensive
> Poisson-Binomial caller, holding all thresholds frozen across every test.
> Across 9 pre-registered GIAB region/depth/sample cells (HG002, HG003,
> HG004) and 2 hap.py-evaluated regions (HG002, HG005), the cascade's F1 is
> statistically indistinguishable from running the expensive caller
> everywhere, while routing under 0.1% of loci. A closed-form genotype
> posterior layer closed a 0.33 F1-point gap versus a naive
> forced-heterozygous baseline. I found and fixed a data-integrity bug in my
> own extraction pipeline mid-project and documented it fully rather than
> discarding the invalidated result.

**250 words:**
> Poisson-Binomial per-read error models improve short-read SNP-calling
> accuracy over simpler binomial screens, but are too expensive to run
> uniformly across a genome. I investigated whether routing only a small,
> confidence-selected subset of loci to the expensive model — via a fixed
> cutoff around a cheap Binomial screen's decision boundary — preserves
> accuracy relative to running the expensive model everywhere, under a
> genuinely frozen, never-retuned configuration.
>
> I validated this across GIAB reference materials HG002, HG003, and HG004
> (an Ashkenazi trio, testing depth, chromosome, and cross-sample
> generalization: 9 pre-registered region/depth/sample cells with an
> internal row-index evaluator) and HG005 (an unrelated sample, evaluated
> with the GA4GH-standard hap.py comparison engine). The cascade's F1 is
> statistically indistinguishable from the expensive-everywhere baseline in
> every tested configuration but one (a held-out test set showing a mixed
> result, reported as such rather than omitted), while routing fewer than
> 0.1% of loci. A closed-form three-hypothesis genotype posterior layer
> (homozygous-ref / het / homozygous-alt) recovered a 0.33 F1-point gap
> against a naive forced-genotype VCF writer, verified against hap.py's own
> per-record genotype-decision tags.
>
> Midway through the project I discovered the extraction pipeline had a
> coordinate-bookkeeping defect that silently invalidated an earlier
> benchmark on a new sample; I traced it to its root cause, fixed it, wrote
> a permanent regression test suite, and re-ran the full validation rather
> than reporting the flawed number. I also located a previously-built but
> unmerged C/htslib extraction accelerator, independently re-verified it was
> bit-exact against the pure-Python reference before trusting its 21-65x
> speedup figure, and integrated it into the production pipeline. Full
> methodology, limitations, and reproducibility details are documented
> alongside the code.

**500 words:**
> **Problem.** Short-read SNP calling faces a persistent cost/accuracy
> trade-off: richer per-read statistical models (e.g. Poisson-Binomial
> mixtures over base-quality distributions) are more accurate than simple
> fixed-error-rate binomial tests, but too computationally expensive to
> apply uniformly across a whole genome. I investigated a routing-based
> architecture — screen every locus cheaply, then send only the
> confidence-ambiguous subset to the expensive model — as an alternative to
> choosing one model for the entire genome.
>
> **Method.** The pipeline is: pileup extraction from an aligned BAM, a
> Binomial log-likelihood-ratio screen against a fixed threshold, a
> confidence router that sends a locus to the Poisson-Binomial arm only if
> its Binomial LLR falls within a fixed band around the screening threshold,
> and — for genotyping — a closed-form binomial-posterior layer evaluating
> three fixed hypotheses (homozygous-reference, heterozygous,
> homozygous-alternate) rather than a learned or per-sample-calibrated
> model. Every threshold and the genotype formula were frozen before
> validation began and never adjusted after seeing any test result; I
> verified this by checksumming the relevant source files across every
> experiment.
>
> **Validation.** I ran 9 pre-registered region/depth/sample cells spanning
> GIAB HG002 (son), HG003 (father, cross-sample generalization test), and
> HG004 (mother, 3 new contigs), using a fast internal row-index evaluator,
> plus 2 experiments evaluated with hap.py, the GA4GH-standard genotype-aware
> comparison engine: HG002 chr21:32-44Mb and HG005 (an unrelated,
> independently-validated sample) chr1:1,000,001-4,000,000. I also ran a
> direct, apples-to-apples comparison against DeepVariant, Clair3, and GATK
> HaplotypeCaller on the HG002 chr21 region.
>
> **Findings.** The cascade's detection F1 is statistically indistinguishable
> from running the expensive Poisson-Binomial caller on every locus, across
> every tested configuration but one held-out test set (which showed a
> mixed result — I report this as-is rather than omitting it), while routing
> under 0.1% of loci to the expensive arm. The genotype posterior layer
> closed a 0.33 F1-point gap versus an earlier, naive forced-genotype VCF
> writer, a result I independently confirmed against hap.py's own
> per-record genotype-decision annotations rather than trusting the
> aggregate number alone. The dominant, reproducible weakness across both
> independently-tested samples is low-mappability/segmental-duplication
> regions, with a genuine low-allele-fraction component.
>
> **Failure discovery and correction.** Midway through validating a new
> sample (HG005), an accuracy result came back near zero. Rather than
> accepting or explaining it away, I traced it to a coordinate-bookkeeping
> defect in my extraction code: a genomic position was being computed by
> simple arithmetic that silently broke whenever a filtering step dropped a
> window of loci. I localized the exact mechanism, fixed it, wrote a
> permanent regression-test suite (11 new tests) to prevent recurrence, and
> re-ran the full validation, independently verifying zero remaining
> coordinate errors across the entire re-extracted dataset before trusting
> any accuracy number computed from it. Separately, I located a
> previously-built but unmerged C/htslib extraction accelerator in an
> isolated development branch, independently re-verified it was bit-exact
> against the pure-Python reference on real data before trusting its
> reported speedup, and integrated it into the production pipeline as the
> default path with an automatic fallback.
>
> All code, data checksums, exact commands, and a full limitations
> discussion are documented for independent reproduction.

## D. CV project entry

- Designed and validated a frozen confidence-router SNP-calling cascade
  (Binomial screen -> router -> Poisson-Binomial) plus a closed-form
  genotype posterior layer against GIAB HG002-HG005 (9 internal-evaluator
  cells + 2 hap.py-evaluated regions); cascade preserved detection F1 while
  routing <0.1% of loci to the expensive arm.
- Ran a direct, apples-to-apples benchmark against DeepVariant, Clair3, and
  GATK HaplotypeCaller on a shared GIAB region using the GA4GH hap.py
  standard.
- Discovered, root-caused, fixed, and regression-tested a coordinate-
  bookkeeping defect in the extraction pipeline that had silently
  invalidated an earlier benchmark; independently re-verified zero
  remaining errors across 2.4M+ genomic loci before re-reporting results.
- Located, independently re-verified (bit-exact, 0 mismatches), and
  integrated a C/htslib extraction accelerator into production, giving a
  measured 21-65x extraction-stage speedup with a 12-category permanent
  regression suite.
- Audited the pipeline for hidden tuning, duplicate implementations, and
  unverified performance claims; caught my own initial audit error (a real
  C/htslib implementation existed in a branch I had not yet checked) and
  corrected course before publishing a false conclusion.

## E. Research portfolio version

**Problem.** Poisson-Binomial per-read error models are more accurate than
simple binomial screens for short-read SNP calling but too expensive to run
genome-wide.

**Method.** A frozen cascade — Binomial screen, confidence router, Poisson-
Binomial arm for routed loci only, closed-form genotype posterior — with
every threshold fixed before validation and never retuned.

**Experiments.** 9 pre-registered GIAB region/depth/sample cells (HG002
baseline, cross-chromosome, held-out test, cross-sample to HG003/HG004) with
an internal evaluator; 2 regions evaluated with the GA4GH-standard hap.py
engine (HG002, HG005); a direct external-caller comparison against
DeepVariant/Clair3/GATK on one shared region.

**Failure discovery and correction.** A near-zero accuracy result on a new
sample (HG005) was not accepted at face value. I traced it to a
coordinate-bookkeeping bug in the extraction code, fixed it, wrote
permanent regression tests, and independently re-verified reference/
coordinate integrity across the full re-extracted dataset (zero mismatches)
before trusting any downstream accuracy number. Separately, I located a
previously-built but unmerged C/htslib extraction accelerator, independently
re-verified it was bit-exact against the pure-Python reference on real data,
and integrated it into production, giving a measured 21-65x extraction-stage
speedup with a permanent regression suite.

**Final findings.** The cascade preserves detection accuracy relative to
the expensive-everywhere baseline (statistically indistinguishable ΔF1 in
every configuration but one honestly-reported mixed result) while cutting
expensive-arm usage to under 0.1% of loci; the genotype layer closes a 0.33
F1-point gap versus a naive baseline; low-mappability/segmental-duplication
regions remain the dominant, reproducible weakness. No whole-genome,
cross-platform, or clinical claim is made — full limitations are documented
alongside the results.
