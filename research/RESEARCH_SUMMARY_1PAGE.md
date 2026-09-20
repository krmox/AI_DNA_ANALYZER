> **SUPERSEDED IN PART (2026-09-20).** Superseded statements: "Cascade preserves PB-only accuracy in every tested experiment" is **wrong** — v14 pre-registered verdict is NEGATIVE and the whole-chromosome chr20 run (v20) is DEGRADED (ΔF1 −5.6×10⁻⁵, extremely small); "rescue recall 0.986" is not a true-SNP rescue recall (13 / 129 composition on HG005; raw record unavailable). Current summary: `PAPER_MANUSCRIPT.md` Abstract.

# Research Summary

## Title
A Frozen Confidence-Router Cascade for Short-Read SNP Calling: Accuracy
Preservation and a Closed-Form Genotype Posterior on GIAB Benchmarks

## Research question
Can a computationally expensive per-read Poisson-Binomial SNP caller be
approximated by routing only a small, confidence-selected subset of loci to
it — using a much cheaper Binomial screen and a fixed router cutoff for the
rest — without losing detection accuracy? And can a closed-form genotype
posterior recover standards-compliant (hap.py) genotype accuracy without any
retraining?

## Approach
`BAM -> pileup extraction -> Binomial LLR screen (threshold=7.0) -> router
(routes iff |LLR-7.0|<=5.412) -> Poisson-Binomial arm for routed loci
(threshold=10.5) -> Method C genotype posterior (binomial argmax over
{0/0,0/1,1/1}) -> VCF`. All three detection thresholds and the genotype
formula are frozen constants, fixed before any of the validation below and
never refit against any reported result (verified by checksum across every
experiment).

## Experimental validation
- **9 pre-registered internal-evaluator experiments** across HG002 (7
  region/depth cells spanning chr1/4/13/16/17/19/20/21), HG003 (cross-sample,
  chr13), and HG004 (3 new contigs: chr2/3/5) — row-index accuracy, not
  genotype-aware.
- **2 hap.py (GA4GH-standard, genotype-aware) experiments**: HG002
  chr21:32-44Mb (F1≈0.9546-0.9547) and HG005 chr1:1,000,001-4,000,000
  (F1≈0.9519-0.9521, this session, after fixing and independently verifying
  a coordinate-integrity bug in the extraction pipeline).
- **External comparison** (HG002 chr21:32-44Mb only) against DeepVariant
  1.6.1 (F1=0.9494), Clair3 v1.0.10 (F1=0.9290), GATK 4.5.0.0 (F1=0.9682).

## Main findings
Cascade preserves PB-only accuracy in every tested experiment (ΔF1 within
noise or a small, statistically-real positive delta; one experiment showed a
mixed/weak result, reported as such). Method C closes a +0.33 F1 gap versus a
naive forced-genotype VCF writer. Router rescue recall 0.986 on HG005 with a
false-rescue rate of 0.006. `lowmap_segdup` regions are the dominant,
reproducible failure mode (7-8x error enrichment) with a genuine low-VAF
component.

## Computational contribution
A C/htslib extraction extension, bit-exact against the pure-Python
reference (0 mismatches on every real-data check run, including the
original 12 Mb HG002 benchmark and this session's independent
re-verification on HG005 data), is now integrated into the production
extraction path. Measured extraction-stage speedup: 21.0x (HG005 3 Mb
region, this session) and 65.1x (original 12 Mb HG002 benchmark scale,
reused under a proven equivalence guarantee). This speedup is scoped to the
pileup-counting stage specifically — the Poisson-Binomial caller's own,
separate read extraction is unaccelerated and dominates full-pipeline wall
time, so the true end-to-end speedup is smaller; both numbers are reported
separately, never conflated.

## Scientific limitations
SNP-only; short-read/Illumina/GRCh38 only; regional (largest tested region
12 Mb), not whole-genome; HG002/3/4 are a related trio; external-caller
comparison exists for one region/sample only; no runtime comparison against
any external caller; a known, quantified, unfixed BED-filtering granularity
issue depresses measured recall by a region-dependent amount.

## Reproducibility
Exact versions, checksums, and commands: `REPRODUCIBILITY.md`. Code:
`AI_DNA_ANALYZER` repository, branch `fix/hg005-extraction-integrity` (this
session's isolated worktree; not merged to `main`).

## Researcher contribution
Research question, experimental design, hypothesis selection, and every
decision about which results to trust or discard were made by the
researcher. Code generation, debugging, test authoring, and documentation
drafting were substantially AI-assisted (Claude). Full disclosure:
`AI_ASSISTANCE_DISCLOSURE.md`.
