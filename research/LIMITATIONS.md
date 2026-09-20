> **SUPERSEDED IN PART (2026-09-20).** Superseded: "v18, v19 demonstrate the frozen cascade generalizes across trio members" overstates — HG002/3/4 are a related trio (weak independence). Chromosome-scale evidence now exists for one chromosome (HG002 chr20, same sample/depth). See `PAPER_MANUSCRIPT.md` §5.

# Scientific Limitations

Stated explicitly and not minimized.

## Scope of validation

- **SNP-only.** No indel, MNP, or structural-variant calling exists anywhere
  in this pipeline. Every indel row in `BENCHMARK_MASTER.csv` shows 0 calls
  by design, not by failure.
- **Short-read, Illumina-only.** No long-read (PacBio HiFi, ONT) data was
  ever used. No claim extends beyond Illumina short-read input.
- **GRCh38 only.** No T2T-CHM13 or other reference was evaluated.
- **GIAB high-confidence truth only.** Accuracy is measured only inside
  regions GIAB itself certifies as high-confidence; performance in
  GIAB-excluded regions (which include some of the hardest genomic sequence)
  is unknown by construction.
- **Regional, not whole-genome.** The largest single contiguous region ever
  evaluated is 12 Mb (`bench_v13`, HG002 chr21 pool) / 3 Mb per individual
  cell elsewhere. No chromosome-scale or genome-scale run has been completed
  for any sample. Whole-genome extrapolations in `PERFORMANCE_RESULTS.csv`
  are explicitly labeled `EXTRAPOLATION`, not measurement.
- **HG002/HG003/HG004 are a trio (father/mother/son) and are not
  independent samples** in the population-genetics sense — they share
  substantial genomic identity by descent. The HG003/HG004 "cross-sample"
  experiments (`v18`, `v19`) demonstrate the frozen cascade generalizes
  across trio members, which is a real and useful result, but is a weaker
  claim than generalization across an unrelated, ancestrally diverse cohort.
- **HG005 is the only non-trio validation sample**, and only one 3 Mb region
  of it has been evaluated with a corrected, verified pipeline. This is one
  additional data point, not a second independent cohort.
- **No population/ancestry diversity was tested.** All 4 samples (HG002-5)
  are GIAB reference materials; no claim about performance across broader
  human genomic diversity is supported.

## Model and pipeline limitations

- **Router rescue analysis is diagnostic only** and was never used to modify
  the router; whether the router's fixed cutoff is actually well-calibrated
  for samples/regions/depths outside those tested is not established one way
  or the other.
- **`lowmap_segdup` and low-VAF regions remain the dominant, reproducible
  weakness** across both HG002 and HG005 — this is the single most
  consistent negative finding in the entire project and should be treated as
  the primary caveat on any accuracy claim.
- **Method C's genotype model is a fixed 3-hypothesis binomial posterior**
  (0/0, 0/1, 1/1, eps=0.01) — it has no mechanism for detecting or modeling
  somatic mosaicism, CNV-driven allele-balance shifts, or ploidy other than
  diploid.
- **A window-level (not per-SNP) BED confidence filter measurably reduces
  recall** independent of detection quality (58% of HG005's FN have no
  candidate at all for this reason) — this is a known, quantified,
  unaddressed engineering choice (`PIPELINE_INTEGRITY_AUDIT.md` M-1), not a
  fundamental model limitation, but it currently depresses every recall
  number in this dossier by an unknown, region-dependent amount.

## Performance limitations

- **The C/htslib acceleration covers only the pileup-counting stage
  (`pileup_counts.load_counts`), not the full pipeline.**
  `read_level_pileup.load_reads` (feeds the Poisson-Binomial caller) is a
  separate extractor, unaccelerated, and is the larger cost in a run that
  performs both stages (HG005 3Mb: `load_reads` 696s vs `load_counts` 344s
  pre-acceleration) — so the true full-pipeline speedup for a run that
  recomputes PB calls from scratch is meaningfully smaller than the
  extraction-stage's own ~65x figure. Both numbers are reported separately
  in `PERFORMANCE_RESULTS.csv`; neither should be quoted for the other's
  scope. Separately, process-level chunking (a different, independent
  mechanism, ~3.1x at 8 workers, sub-linear) remains unintegrated.
- **All performance numbers are single-machine, single-hardware-configuration
  measurements** (one 4c/8t desktop CPU). No cloud, cluster, or alternative
  hardware benchmark exists.
- **No whole-genome runtime has ever been measured.** Any whole-genome time
  estimate in this dossier is a linear extrapolation from a 2-3 Mb
  measurement and should be treated with substantial skepticism (real-world
  scaling is very unlikely to be linear across genomic regions of widely
  varying complexity/coverage).

## Evaluation methodology limitations

- **Two internal evaluators with similar names but incompatible semantics
  have both been used historically** (`PIPELINE_INTEGRITY_AUDIT.md` M-2).
  This dossier labels every number's evaluator explicitly to prevent
  confusion, but any number lifted out of context from an underlying devlog
  file without that label risks being misinterpreted.
- **hap.py was only run for 2 experiments** (HG002 chr21:32-44M, HG005
  chr1:1,000,001-4,000,000) out of 9 total internal-evaluator experiments.
  The genotype-aware, standards-compliant accuracy figure is therefore only
  available for a small fraction of the tested scope.
- **No confidence interval or bootstrap was computed for the HG005 hap.py
  result** — it is a single point estimate on a single region/sample.
- **External caller comparison exists for exactly one region/sample/depth**
  (HG002 chr21:32-44M, 15x). No external-caller comparison exists for HG005,
  HG003, HG004, or any other region.
- **No `requirements.txt`/environment lockfile exists** in this repository
  (noted in `REPRODUCIBILITY.md`) — exact reproduction currently depends on
  the version table recorded by hand in that document, not on a pinned,
  automatically-installable environment.

## What this dossier does NOT claim

See `SCIENTIFIC_CLAIMS.md` for the explicit claim-by-claim classification.
In summary: no claim of whole-genome validation, clinical readiness,
universal or general superiority over any external caller, cross-platform
(long-read) validation, or population-scale generalization is made anywhere
in this research package.
