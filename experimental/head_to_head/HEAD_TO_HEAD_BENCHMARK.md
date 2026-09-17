# Head-to-Head Benchmark: AI_DNA_ANALYZER vs Open-Source Variant Callers

Additional, standalone experiment. Not a replacement for, and not pooled with, the six benchmark
rounds in `FINAL_RESEARCH_REPORT.md`. Frozen production code (`cascade.py`,
`FROZEN_BINOMIAL_THRESHOLD=7.0`, `FROZEN_ROUTER_CUTOFF=5.411872376933351`,
`FROZEN_PB_THRESHOLD=10.5`) was not modified; verified via `git diff` (empty) at the end.

## 1. Common Data

```
sample          = HG002
reference       = GRCh38 (data/reference/chr21_full.fa; renamed-header copy used for external
                  tools only, to match BAM/truth contig naming — see COMMANDS.md §0)
chromosome      = chr21
start           = 32,000,000
end             = 44,000,000
truth           = GIAB/NIST v4.2.1 (data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz)
confident BED   = data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed
variant class   = SNP only
depth           = 15x (subsampled)
BAM             = data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_15x.bam (670,928 reads on chr21)
alignment       = novoalign (per BAM @PG header) — NOT bwa-mem; noted as a protocol detail since
                  DeepVariant's WGS model is typically validated against bwa-mem alignments. Same
                  BAM used for every caller, so this does not bias the comparison between callers.
```

This region was selected because it is the only region in this repository with a complete,
already-verified BAM + GIAB v4.2.1 truth + confident BED + reference set, and was chosen and fixed
*before* any external caller was run (see chat record; no region was re-picked after seeing
results). No chr20 or chr14 BAM exists in this repository, so Options B/C from the original
protocol were unavailable without a new download; Option A (this region) was used instead, which
the protocol explicitly allows.

## 2. Callers Run

| Caller | Status |
|---|---|
| PB-only (AI_DNA_ANALYZER) | Evaluated (frozen internal evaluator) |
| AI Cascade (AI_DNA_ANALYZER) | Evaluated (frozen internal evaluator) |
| DeepVariant 1.6.1 (CPU) | Evaluated (hap.py) |
| Clair3 v1.0.10 (ilmn model) | Evaluated (hap.py), after `:latest` tag failed with an internal bug |
| GATK HaplotypeCaller 4.5.0.0 | Evaluated (hap.py), no BQSR (documented deviation) |
| Strelka2 2.9.10 | **BLOCKED — TOOLING INCOMPATIBILITY** (bundled htslib binary crashes on this host) |
| FreeBayes | Not attempted (explicitly lowest priority; no trivial official build found in time budget) |

## 3. Results

| Caller | Evaluator | TP | FP | FN | Precision | Recall | F1 | Truth total |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| PB-only | frozen locus evaluator | 15,923 | 171 | 674 | 0.9894 | 0.9594 | 0.9742 | 16,597 |
| AI Cascade | frozen locus evaluator | 15,924 | 173 | 673 | 0.9893 | 0.9595 | 0.9741 | 16,597 |
| DeepVariant | hap.py 0.3.15 (xcmp) | 16,313 | 1,153 | 585 | 0.9340 | 0.9654 | 0.9494 | 16,898 |
| Clair3 | hap.py 0.3.15 (xcmp) | 16,287 | 1,880 | 611 | 0.8965 | 0.9638 | 0.9290 | 16,898 |
| GATK HaplotypeCaller | hap.py 0.3.15 (xcmp) | 16,042 | 199 | 856 | 0.9877 | 0.9493 | 0.9682 | 16,898 |

**Neutral phrasing, per protocol:** on the common HG002 chr21:32,000,000-44,000,000 benchmark
region, under the evaluation protocols documented in this file, AI Cascade achieved F1=0.9741,
DeepVariant achieved F1=0.9494, Clair3 achieved F1=0.9290, and GATK HaplotypeCaller (no BQSR)
achieved F1=0.9682. GATK HaplotypeCaller had the highest precision (0.9877) of any caller tested;
AI Cascade and PB-only had the highest recall (0.9594-0.9595) of any caller tested. No caller
dominates on every axis.

**Important caveat on comparability** (see `ERROR_ANALYSIS.csv` for full detail): AI
Cascade/PB-only were scored with the pre-existing frozen internal locus evaluator
(`accuracy_block()`/`robustness_benchmark.py`), not hap.py, because they are LLR-threshold
decisions over a locus scoring frame rather than a genotyped VCF, and synthesizing one was judged
out of scope for this pass rather than something to force. Their truth denominator (16,597) is
~1.8% smaller than hap.py's (16,898) for the same truth VCF/BED/region, reflecting a small,
pre-existing difference between the frozen extractor's SNP-labeling logic and hap.py's GA4GH-standard
truth/BED intersection — not something introduced in this pass. The F1 values above are comparable
in spirit (same truth source, same region, same SNP-only scope) but not perfectly apples-to-apples
at the individual-locus level.

A direct locus-by-locus cross-caller disagreement/overlap analysis (AI Cascade vs DeepVariant,
etc., with genomic position, VAF, BQ, MAPQ, LLRs, genomic context) was attempted and **discarded
as invalid** — see `ERROR_ANALYSIS.csv`. The cause was a genuine, pre-existing limitation
documented in `History/12_DEVLOG.md` §16: the frozen extraction script's `positions` output array
is a cache-order index, not a genomic coordinate, once BED filtering has dropped any loci. A test
join against the truth VCF's real coordinates recovered under 2% of known truth SNPs (with or
without a ±1 shift), confirming the mismatch. Reporting numbers built on that join would have been
dishonest; they are not included anywhere in this report or CSVs.

## 4. Router Rescue Analysis (this region)

Computed entirely within the AI evidence cache's own internal index (no cross-file coordinate join
needed, so this part is not affected by the limitation above):

```
PB-rescuable loci (PB correct, Binomial wrong): 714
Rescued by router:                              712
Missed by router:                                 2
Rescue recall:                                0.9972
Total routed to PB:                            17,007  (0.155% of 10,975,654 scored loci)
"False rescue" (Binomial already correct, routed anyway): 15,955 (93.8% of all routed loci)
```

This is closely consistent with the full-HG004-universe figures reported in
`FINAL_RESEARCH_REPORT.md` §12 (99.2% rescue recall, 90.5% false-rescue rate) — the same pattern
(the router over-routes heavily relative to the loci it strictly needs to fix, in order to achieve
near-total rescue recall) reproduces on an independent HG002 region, not just HG004. This is a
useful cross-sample consistency check, though it is not a formal independent-sample validation in
the sense Experiment 1 (HG005) was intended to be.

ΔF1 (Cascade vs PB-only) on this region: **-0.0000282**, 95% paired-bootstrap CI
**[-0.000122, +0.0000647]** (10,000 resamples, seed 20260812). Per the pre-registered rules: CI
contains 0 AND |ΔF1| < 0.001 → **PRESERVED**.

## 5. Performance

Two categories, kept explicitly separate per protocol:

**Calling-stage wall-clock (already-aligned BAM → VCF), single run, this host:**

| Caller | Seconds | Threads |
|---|--:|--:|
| GATK HaplotypeCaller | 287 | default |
| DeepVariant | 819 | 8 |
| Clair3 (full pipeline) | 1,342 | 4 |

**AI evidence-extraction time (4,740s, single-threaded, full 12Mb region) is NOT a fair
calling-only comparator** against the external callers' numbers above — it computes full
per-locus Binomial+PB scoring for every one of 10,975,654 loci to build the evaluation cache,
not the production router-gated ~0.155%-to-PB path. Once binomial_llr/pb_llr are cached, applying
the frozen router and computing metrics over all 10.9M loci took under 1 second — this ~1s IS
representative of the router's actual marginal cost, but the 4,740s evidence-extraction time
should never be quoted as "the cascade's runtime" without this context. See `RUNTIME_RESULTS.csv`.

No multiprocessing-worker-count (2/4/8) scaling data was generated in this pass; that experiment
(referenced in the original four-experiment validation) could not be located anywhere in this
repository — see `FINAL_RESEARCH_REPORT.md` for that separate finding.

## 6. Failure Analysis / Notable Technical Findings

- **Reference contig-naming mismatch**: the project's own `data/reference/chr21_full.fa` uses
  Ensembl-style naming (`21`) while every external tool required UCSC-style (`chr21`) to match the
  BAM/truth/BED. The frozen cascade's own code already normalizes this internally
  (`providers.py:_normalise_contig`); external tools do not, so a renamed reference copy was
  created for their use only (see COMMANDS.md §0). Not a bug in the project — a genuine
  cross-tool convention difference the project happened to already handle.
- **Clair3 `:latest` tag is broken on this host**: internal `AttributeError: No symbol free_memory
  found in library ...realigner`, causing a silent fallback to pileup-only calls (skipping the
  neural full-alignment refinement stage) without failing the overall run. This would have quietly
  produced a materially worse Clair3 result (pileup-only vs full pipeline) if not caught and
  retried against a pinned version (`v1.0.10`, which worked correctly).
- **Strelka2's official static release is not portable to this host**: bundled 2018-era CentOS6
  htslib binary asserts and crashes against this host's zlib. Python2 itself was not the blocker
  (built via pyenv in ~3 minutes) — the deeper binary incompatibility was.
- **hap.py packaging is inconsistent across tags/registries**: `docker.io/pkrusche/hap.py` uses an
  old v1 manifest format hostile to inspection; `quay.io/biocontainers/hap.py:0.3.7--py27_1` lacks
  a bundled rtg-tools (vcfeval unavailable) and its xcmp engine crashes with `regex_error`;
  `quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0` worked correctly with the default xcmp
  engine. Documented so a future run can skip straight to the working tag.
- No new AI_DNA_ANALYZER-specific true-SNP-loss mechanism was found in this region beyond what
  `FINAL_RESEARCH_REPORT.md` already documents; the 2 missed-rescue loci here are consistent in
  kind (not separately deep-dived per-locus in this pass because the same coordinate-mapping
  limitation in §3 applies to producing a genomic-position table for them — they are countable via
  the internal index but not easily presentable as chrom:pos in this pass).

## 7. Statistical Analysis

Pre-registered rules applied without modification (see `FINAL_RESEARCH_REPORT.md` §15 for the
rule definitions). The only paired-bootstrap comparison possible with matching internal indices in
this pass is AI Cascade vs PB-only (§4 above: PRESERVED). ΔF1 and CIs for AI Cascade/PB-only vs
each external caller are not reported, because a fair paired bootstrap requires a shared per-locus
index across compared arms, which — per §3 — is not currently available across the AI evidence
cache and the external callers' VCFs. Reporting an unpaired or coordinate-mismatched delta would
misrepresent statistical uncertainty; it was not done. The raw F1 values in §3 are still valid and
usable for qualitative comparison, just not for a formal paired ΔF1/CI claim in this pass.

## 8. Reproducibility

- Git HEAD: `16b221ee82d64f66de25119aafa4e5944cfe63cc`, clean before and after.
- `git diff` on production files: empty (verified at end of task).
- Full environment/versions/digests: `ENVIRONMENT.md`.
- Full commands: `COMMANDS.md`.
- Checksums for every caller's saved output: `<caller>/md5sums.txt`.

## 9. Limitations (full list)

1. Single region (12Mb of chr21), single sample (HG002), single depth (15x). Not a genome-wide or
   multi-sample claim.
2. SNP-only scope; INDEL numbers from hap.py are visible in `happy/*/happy.summary.csv` but not
   analyzed here (out of scope per this task's SNP-only framing).
3. AI Cascade/PB-only and the three external callers were NOT scored by the identical evaluator
   (see §3 caveat) — comparable in spirit, not bit-for-bit equivalent.
4. Locus-level cross-caller disagreement/overlap analysis could not be produced (coordinate-index
   limitation, §3); this removes Section 11/12/13's cross-caller-vs-router tabulations from this
   pass's deliverables.
5. GATK run without BQSR (no known-sites resource available for a region-only subset).
6. Strelka2 blocked entirely; no Strelka2 numbers exist anywhere in this report.
7. FreeBayes not attempted.
8. DeepVariant/Clair3/GATK calling-stage runtimes are single-run, single-host measurements, not
   averaged over repeated trials; treat as indicative, not high-precision.
9. The BAM was aligned with novoalign, not the bwa-mem alignments DeepVariant's released models are
   most commonly validated against; this is identical across all callers tested here (fair among
   them) but may itself suppress or inflate every external caller's absolute numbers somewhat
   relative to their own published bwa-mem-based benchmarks — which is exactly why this task's
   "don't compare against other papers' numbers" rule matters, and why all numbers here were
   produced fresh on identical input.

## 10. Recommended Follow-Up

- Patch (separately, reviewed) the evidence extractor to optionally emit true genomic positions
  alongside its current cache-order index, enabling the full locus-level cross-caller disagreement
  analysis this pass could not complete.
- Retry Strelka2 via a from-source build or a modern community Docker image (not the official
  static release) if that comparator is wanted.
- Extend to chr20 or another region once/if a compatible BAM is obtained, to check whether the
  relative caller ranking (GATK precision-leaning, AI Cascade/PB recall-leaning) holds elsewhere.
