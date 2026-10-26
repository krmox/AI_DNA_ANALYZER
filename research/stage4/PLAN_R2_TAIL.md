# Varion R2 — plan (draft, 2026-10-26; nothing implemented; R1.1 stays frozen)

## Owner decisions
* New version is named **R2**. R1.1 is not modified.
* Target: not worse than bcftools at every depth, and SNP F1 >= 0.99 at >= 30x (thresholds to be fixed in a preregistration BEFORE any run).
* Idea from the owner: a small, flexible **"final tail"** after the router, in C (ASM only if profiling demands it), costing <= ~2 s and little RAM.

## Evidence this rests on (Stage 3 diagnostics, HG004 chr22 30x s101; exploratory, post hoc)
* 71-91 % of FN/FP lie in GIAB alldifficultregions (21 % of truth SNP); outside them F1 ~ 0.996-0.998.
* 299/493 pure FP are shared with bcftools: depth 1.43x median, VAF ~0.28, alt MAPQ ~52 (collapsed paralogs).
* 194 Varion-only FP: VAF ~0.15, k ~4 (weak candidates accepted by R1). 138 FN: k>=3 rejected by R1.
* z-threshold is already at its oracle optimum (0.98452 vs 0.98453); QUAL is '.' in the VCF.
* Best-of-both ceiling on that run ~0.989: classifier work alone cannot reach 0.99 on chr22.

## Tail design (hypothesis, to be tested offline first)
Block 1 — difficulty detector (per candidate, cheap, no truth, no BED needed; optional researcher BED as extra input):
 local depth / contig-median depth, fraction of low-MAPQ reads, candidate clustering (nc25), VAF-vs-copy-number mismatch.
Block 2 — pattern layer: compare the candidate/local statistics with expected distributions at three scales
 (chromosome, researcher-defined region/gene, genome), e.g. expected variant density, het-VAF around 0.5 given local copy estimate,
 Ti/Tv. Output a bounded adjustment of the decision score. Researcher parameters: region list, scale weights, strength.
Cost model: candidates are ~44k of 40M loci per chr22, so the tail touches ~0.1 % of positions; 2 s and tens of MB are plausible in plain C.
ASM/AVX2 only for a measured hot loop; compiler -O3 is expected to suffice (state this honestly in the report).

## Risks to test, not assume
* Density/prior normalisation can suppress real variants in variant-dense regions (recall loss) -> report recall by stratum.
* Parameters chosen by the researcher must have frozen defaults set on dev data only.
* Training frame still has the truth-indel mask (candidates at/near truth indels never seen as negatives) -> retrain without mask.

## Data discipline
* Dev/validation (already inspected): HG004 chr20/21/22, HG002 chr20, HG003 chr8, HG002/HG003 chr22.
* Blind confirmation (do NOT inspect errors before R2 is frozen): chr18 of HG002/HG003/HG004, and HG005 (chr22, chr18).
* Preregistration first (endpoints, thresholds, ablations), then offline ablation (features / tail / genotype / training data), then C implementation, then a single confirmation run.

## Order of work and rough time (working estimate, not a promise)
1. Finish Stage 3 (running): ~5-8 h wall for remaining units + timing campaign + report.
2. Preregistration: ~0.5 day.
3. Offline ablations (Python on dumped features): 1-2 days.
4. C implementation + equivalence tests + timing guardrails: 1-2 days (ASM +1 day only if needed).
5. Confirmation run and report: ~1 day.
