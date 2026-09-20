> **SUPERSEDED IN PART (2026-09-20).** Earlier draft, superseded by `PAPER_MANUSCRIPT.md`; several statements (native implementation "not found", no chr20, 1.22× end-to-end) are outdated.

# PAPER DRAFT (working draft — not submission-ready)

## 1. Title

A Frozen Confidence-Router Cascade for Short-Read SNP Calling: Accuracy
Preservation and a Closed-Form Genotype Posterior on GIAB Benchmarks

## 2. Abstract

Poisson-Binomial variant callers model per-read quality more precisely than
simple binomial screens but are computationally expensive to run at every
genomic locus. We investigate whether a cheap Binomial LLR screen, combined
with a fixed-cutoff confidence router that sends only a small, uncertain
subset of loci to the expensive Poisson-Binomial caller, can preserve
detection accuracy relative to running the expensive caller everywhere. We
further investigate whether a closed-form binomial-posterior genotype layer
(three fixed hypotheses: homozygous-reference, heterozygous, homozygous-alt)
recovers standards-compliant, genotype-aware accuracy without any retraining.
Using GIAB reference materials HG002, HG003, HG004, and HG005 across 9
pre-registered region/depth/sample cells (internal, row-index evaluation)
and 2 hap.py-evaluated regions (HG002 chr21:32-44Mb, HG005
chr1:1,000,001-4,000,000), we find the cascade preserves or slightly improves
PB-only accuracy in 6/7 internal experiments (one showing a mixed result) and
is statistically indistinguishable from PB-only under hap.py in both tested
regions (ΔF1 within measurement noise), while routing fewer than 0.1% of
loci to the expensive arm. The genotype layer closes a 0.33 F1-point gap
versus a naive forced-genotype baseline. A reproducible weakness persists in
low-mappability/segmental-duplication regions (7-8x error enrichment). We
also report and correct an internal, previously undetected coordinate-
bookkeeping defect that had invalidated an earlier HG005 benchmark. A
validated C/htslib extraction extension, proven bit-exact against the
Python reference implementation, was integrated into the production
extraction path in this work, giving a measured speedup (pileup-counting
stage only; see Performance) consistent with a previously-recorded ~65x
figure that an earlier pass of this document's own audit had (incorrectly)
flagged as unsubstantiated before the implementation was located and
re-verified. All results are regional (largest evaluated span 12 Mb) and
platform-specific (Illumina short-read, GRCh38); no whole-genome, long-read,
or cross-population claim is made.

## 3. Introduction

Variant calling from short-read sequencing data typically trades off
statistical fidelity against computational cost: richer per-read error
models (e.g. Poisson-Binomial mixtures over base qualities) improve accuracy
but are too expensive to apply uniformly across a genome, while cheap
screens (fixed-error-rate binomial tests) are fast but less accurate at
marginal loci. [REFERENCE NEEDED — prior work on cost/accuracy trade-offs in
variant calling]. This project investigates a simple architectural pattern —
route only ambiguous loci to the expensive model — as an alternative to
either uniformly cheap or uniformly expensive calling, and asks whether the
pattern holds up under a genuinely frozen, unretuned configuration across
multiple independent GIAB samples.

## 4. Related Work

[REFERENCE NEEDED — DeepVariant, Clair3, GATK HaplotypeCaller methodology
papers; GIAB benchmark methodology papers (Zook et al.); hap.py/vcfeval
comparison-engine literature; confidence-routing / cascade architectures in
other domains (e.g. cascaded classifiers in computer vision)]. This draft
does not fabricate citations; every claim requiring literature support is
marked `[REFERENCE NEEDED]` rather than invented.

## 5. Methods

See `METHODS.md` for full detail. Summary: pileup extraction (Python/pysam)
-> Binomial LLR screen (frozen threshold 7.0) -> confidence router (frozen
cutoff 5.411872376933351 around the binomial threshold) -> Poisson-Binomial
LLR (frozen threshold 10.5) for routed loci -> cascade decision -> Method C
genotype posterior -> VCF.

## 6. Cascade Architecture

```
BAM -> pileup extraction -> Binomial LLR -> router -> [PB LLR if routed]
    -> cascade decision -> Method C genotype -> VCF
```
Router: routes locus i to PB iff `|binomial_LLR_i - 7.0| <= 5.411872376933351`.
This is a symmetric confidence band around the binomial decision boundary,
not a learned or per-sample-tuned quantity.

## 7. Genotype Layer

Method C: for observed alt-allele count k out of n called reads, compute the
binomial log-likelihood under three fixed hypotheses θ∈{0, 0.5, 1}
(eps-adjusted: {0.01, 0.5, 0.99}), take the posterior argmax as the genotype
call, and report GQ = 10·log10(P_best/P_second), capped at 99. A 0/0-call
safety rule forces a 0/0 posterior-argmax result to 0/1 with GQ capped at 5
(a conservative choice favoring sensitivity over specificity at the genotype
layer specifically). [REFERENCE NEEDED — standard GATK-style genotype
likelihood formulation this resembles].

## 8. Implementation

Pure Python + pysam/numpy/scipy; no compiled/native extension exists in this
codebase (`PIPELINE_INTEGRITY_AUDIT.md` I-2). Read admission filters:
MAPQ>=20, base quality>=13, duplicate/QC-fail/unmapped/secondary/supplementary
excluded. Coordinates 0-based half-open internally throughout.

## 9. Experimental Design

Regions were selected mechanically from public annotation (mappability, GC
content, GIAB difficulty strata) by a dedicated selection script *before* any
caller ran on them (`select_v1N_region*.py`, independence verified
programmatically per-run). No threshold or genotype-formula parameter was
ever fit to, or adjusted after seeing, any test-set result (verified by
checksum across every experiment in this dossier).

## 10. Datasets

GIAB HG002 (son), HG003 (father), HG004 (mother) — an Ashkenazi trio — and
HG005 (Chinese trio son, unrelated to HG002-4). GRCh38 reference. NIST
v4.2.1 high-confidence truth VCF/BED for all samples. Illumina short-read,
novoalign-aligned. Regions: chr21 (multiple sub-windows, HG002 baseline);
chr20/19/4/1 (HG002 cross-chromosome); chr16/15/7/14 (HG002 held-out test);
chr13 (HG002, then repeated on HG003 for cross-sample); chr2/3/5 (HG004,
new contigs); chr1:1,000,001-4,000,000 (HG005).

## 11. Evaluation

Two evaluators, never conflated (`METHODS.md` Sec. 3): an internal row-index
boolean evaluator (fast, used for the 9-cell sweep) and hap.py 0.3.15
(`xcmp` engine, GA4GH-standard, genotype-aware, used for HG002 chr21 and
HG005).

## 12. Results

See `RESULTS.md` and `BENCHMARK_MASTER.csv` for full tables. Headline: 6/7
internal-evaluator experiments show cascade F1 >= PB-only F1; both
hap.py-evaluated experiments show cascade and PB-only statistically
indistinguishable (ΔF1 = -0.0001 and -0.0002 respectively). Method C
recovers +0.331 F1 over a forced-genotype baseline on HG002 chr21. External
comparison (HG002 chr21 only): AI Cascade+MethodC F1=0.9546, between
Clair3 (0.9290) and GATK (0.9682), ahead of DeepVariant (0.9494) on this one
region under a precision/recall trade-off, not a strict dominance.

## 13. Error Analysis

See `ERROR_ANALYSIS.csv`. Dominant, reproducible failure mode:
low-mappability/segmental-duplication regions (FN enrichment 7.7x, FP
enrichment 6.4x on HG005), with a genuine low-VAF component among false
negatives that have any extracted candidate at all (median VAF 0.0 in-segdup
vs 0.19 elsewhere). A separate, non-biological engineering artifact (window-
level vs SNP-level BED-confidence granularity) accounts for 58% of HG005's
false negatives having no candidate locus generated at all — explicitly
distinguished from a model/detection failure (`PIPELINE_INTEGRITY_AUDIT.md`
M-1).

## 14. Performance

A C/htslib extraction extension is integrated into production
`pileup_counts.load_counts`, proven bit-exact against the pure-Python
reference (0 mismatches, every check run). Measured extraction-stage
speedup: 21.0x (HG005 3 Mb region, this session, apples-to-apples
same-loci-count comparison) and 65.1x (original 12 Mb HG002 benchmark scale,
reused under a full-region equivalence proof). This figure is scoped to the
pileup-counting stage; the Poisson-Binomial caller's separate read
extraction (`read_level_pileup.load_reads`) is unaccelerated and dominates
full-pipeline wall time (HG005 3Mb: 696s vs 344s pre-acceleration for
`load_counts`), so full end-to-end speedup is smaller — measured, not
estimated, and reported separately in `PERFORMANCE_RESULTS.csv`. A separate,
smaller, unintegrated process-level parallelism mechanism (~3.1x at 8
workers, sub-linear) and a genuine null result from htslib-internal thread
scaling (~0x) are also documented, kept distinct from the C-extension
result above.

## 15. Limitations

See `LIMITATIONS.md` in full. SNP-only; short-read/Illumina/GRCh38 only;
regional (largest span tested 12 Mb); HG002-4 are a related trio; external-
caller and hap.py coverage is limited to one region each; a known,
unaddressed recall-depressing BED-granularity artifact; no runtime
comparison against any external caller; no environment lockfile.

## 16. Reproducibility

See `REPRODUCIBILITY.md`: exact versions, checksums, commands, hardware.

## 17. Discussion

The consistency of the cascade's near-null ΔF1 across genuinely independent
generalization axes (region, depth, sample) — with one honestly-reported
mixed result rather than a uniformly positive record — is the paper's
central empirical claim, and is deliberately reported without smoothing over
the one weaker cell. The genotype-layer result (Method C) is a cleaner,
larger, and more mechanistically understood effect than the cascade's own
ΔF1, and is arguably the stronger finding of the two. The low-mappability/
segdup weakness recurring across two independent samples (HG002, HG005) is
evidence that it is a structural property of the pileup-based evidence these
callers consume, not sample-specific noise — a natural target for future
work (e.g. local reassembly or haplotype-aware evidence in exactly those
regions, as DeepVariant/GATK already attempt via different mechanisms
[REFERENCE NEEDED]).

## 18. Conclusion

Under the tested scope — GIAB HG002/HG003/HG004/HG005, Illumina short-read,
GRCh38, SNP-only, the specific regions evaluated — a frozen confidence-router
cascade preserves Poisson-Binomial-level detection accuracy while
substantially reducing the fraction of loci requiring the expensive caller,
and a closed-form genotype posterior recovers standards-compliant genotype
accuracy without retraining. No broader claim (whole-genome, cross-platform,
clinical, or general-superiority) is supported by the evidence gathered.

## 19. References

`[REFERENCE NEEDED]` — this draft intentionally contains no fabricated
citations. Before submission, add: GIAB/NIST truth-set methodology (Zook et
al.), hap.py/vcfeval comparison methodology, DeepVariant, Clair3, and GATK
HaplotypeCaller primary papers, and any prior cascade/routing-classifier
literature the discussion section should be grounded against.
