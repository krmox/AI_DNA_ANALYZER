# Architecture (R3.1)

Two independent callers share a repository and a front-end command. They have different evidence behind them.

```
                              varion  (Python, stdlib only)
                 ┌───────────────┬───────────────────────────────┐
                 │ varion snp    │ varion indel / depth / doctor │
                 ▼               ▼
        dnav2 (C++20, htslib)    nindels-call (C++20, htslib)
        SNP engine R2            indel engine R3.1
```

## SNP engine (R2, frozen)

`native_r2/` → `dnav2 run --r2`. Lineage: V1 (Python/C, fixed-error binomial screen → per-read Poisson-binomial
model on the few loci near the decision boundary) → V2 (C++ rewrite, byte-identical calls) → V2.x (+ cheap read-level
evidence and a 32-input logistic classifier) → R1.1 (BED-edge windows, truth-free frame) → R2 (self-calibration,
out-of-distribution guard, pattern tail, retrained model, `research/stage4/frozen_R2.model`). There is no neural
network in any version. R3.1 does not change this engine; `varion snp` only assembles its command line.

## Indel engine (R3.1)

`nindels/` → `nindels-call`. Pipeline per contig:

1. **Candidate discovery** from CIGAR I/D events (>= 2 supporting reads, alt fraction >= 0.03), left-normalised.
2. **Joint multi-allelic genotyping (A)** — `--genotyper multi`. Candidates whose repeat footprints overlap form one
   locus with alleles 1..K; every read is scored against every haplotype with a banded PairHMM; genotypes are
   unordered allele pairs, so a read from a third allele no longer counts as "reference". PASS = non-reference
   genotype with phred QUAL >= 20.
3. **Repeat error model (B)** — `--str-model`. Per (period, repeat-units) slippage rate ε estimated from the BAM being
   called (binomial-mixture EM over sampled repeat loci, **no truth labels**), turned into a position-specific
   gap-open probability inside each repeat in the PairHMM. B needs enough reads over repeats to estimate ε.
4. **Depth-adaptive switch (R3.1)** — `varion indel --mode auto`:
   - `nindels-call depth` samples 10-kb windows of each contig and reports the median depth and spread;
   - `varion/adaptive.py::choose_config` compares the contig median with `global_min_depth` from
     `varion/thresholds.json`; below it (or if no depth was ever found to favour B) it uses **A**, otherwise **A+B**;
   - inside an A+B contig, a **per-locus gate** (`--str-min-depth`, `local_min_depth`) applies B only at loci with
     at least that many reads, so low-coverage stretches of a high-coverage genome fall back to A;
   - `--mode A` / `--mode A+B` override the choice. The decision, its reason and the thresholds' provenance are written
     to the log, the JSON summary and the VCF header (`##varion_config_<contig>=`); B-scored loci carry the `STRB` flag.

Not used in R3.1 indel calling: the pattern memory, the learned logistic model and the reassessment gate described in
`nindels/SPEC_INTERPRETATION.md` (they need labelled training data and were not part of the evaluated configurations).

## Where the numbers come from

| Quantity | Produced by |
|---|---|
| depth, decision | `nindels-call depth`, `varion/adaptive.py` |
| switching thresholds | `nindels/scripts/r31/select_threshold.py` on HG004 chr20/chr21 only → `varion/thresholds.json` |
| benchmark tables | `nindels/scripts/r31/run_one.py` + `report_test.py` (hap.py 0.3.15, xcmp) |
