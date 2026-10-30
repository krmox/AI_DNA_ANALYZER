# nindels — implementation plan

Layout: `nindels/` (CMake, C++20, htslib only). One binary, `nindels-call`, with subcommands.
Python lives only under `nindels/scripts/` (download, vcfeval wrapper, stratified metrics, report).
Data and outputs go to `/mnt/hdd/nindels/{data,runs}`.

| # | Stage | Algorithms / data structures | Files | Tests (ctest, synthetic, no data) |
|---|---|---|---|---|
| 0 | Toolchain | htslib/samtools/bcftools 1.22 built in user space into `~/opt/nindels-tools`; rtg-tools 3.12.1 (bundled JRE) | `scripts/setup_tools.sh` | `--version` checks |
| 1 | Candidate discovery + normalization | chunked BAM scan (1 Mb + 1 kb padding) with a compact per-read record (pos, CIGAR ops, strand, MAPQ, clip, pointer to BQ); CIGAR I/D → left-shift against the reference (period roll); `unordered_map<AlleleKey,Acc>`; threads over chunks, deterministic merge by position | `io_bam.*`, `reference.*`, `candidates.*`, `normalize.*` | left-shift in homopolymers/repeats; reads with the same event in different CIGAR positions merge; chunk-boundary dedup |
| 2 | Evidence E_i | single sweep over the reads overlapping each candidate (sorted interval sweep); alt/ref/other classification from normalized CIGAR events; all features in spec-v2 §2 | `evidence.*` | hand-built reads → exact counts, SB, read-pos, clip |
| 3 | Pattern features P_i, complexity C_r | homopolymer run, minimal period, copy count, entropy; 2-bit packed context | `pattern.*`, `complexity.*` | AAAAAA → hp 6; CACACA → period 2, copies 3 |
| 4 | Similarity | Levenshtein (two-row DP, ≤ 81 bp, early cutoff not used, to stay exact), d_k, D, S | `similarity.*` | metric properties, empty-empty = 0, D ∈ [0,1] |
| 5 | Memory + index + statistics | bucket hash → vector of classes; leader clustering; medoid reservoir; binary memory file with header (sha256 of inputs, θ, s_*, α, β); top-K retrieval; P_context, N_eff_kish, N_eff_loci, S_max | `memory.*` | determinism (same input → identical bytes); retrieval == brute force on small sets |
| 6 | Labels | truth VCF (htslib) normalized and split; exact allele match inside the confident BED | `label.*` | match/mismatch fixtures |
| 7 | Count GL + diploid PairHMM GL + original Δℓ | 3-state PairHMM in log space (double, banded to ±(L+8) diagonals around the CIGAR placement, **ASSUMPTION**: the band is a speed choice, checked against unbanded on fixtures); MAPQ mixture; GL(0/0,0/1,1/1) | `pairhmm.*`, `genotype.*` | perfect read → strong preference; het mixture → GL01 max; banded == full on fixtures |
| 8 | Q_k | per-position aggregation from the HMM-best haplotype per read | `nucleotide.*` | synthetic conflict position gets low Q |
| 9 | Models | standardized logistic, IRLS + L2, missingness indicators + training masking; stage-1 P_base/P_evid; R_i, τ_review(A) | `model.*`, `reassess.*` | IRLS matches a closed-form 1-D case; A=0 vs 1 gives monotone reassessed counts |
| 10 | VCF writer + CLI | `candidates`, `label`, `build-memory`, `train`, `call` (with `--A`, `--no-memory`, `--neff kish|loci`, `--dl orig|diploid`, `--candidates-vcf` for shadow); INFO carries every feature, S_max, closest class ids, reassessed flag | `vcf_out.*`, `main.cpp` | VCF passes `bcftools view` |
| 11 | Glue + benchmark | download/subset scripts (pinned URLs, md5), downsampling with `samtools view -s` (fixed seed), bcftools baseline, vcfeval, stratification, `/usr/bin/time -v` ×3, `perf stat`/`perf record` | `scripts/*.sh`, `scripts/*.py` | — |
| 12 | Experiments | train → validation sweep (variants, A, ablations) → freeze → one test run → REPORT.md | `REPORT.md` | — |

Variants benchmarked against baseline-as-specified: `V-no-memory`, `V-dl-orig`, `V-neff-loci`,
`V-w-grid`, `V-lambda`, shadow-bcf, A sweep, plus a trivial threshold caller (n_alt ≥ 3 and f ≥ 0.2)
as a floor. If the learned model cannot beat that floor, the report will say so.

Order of work: stages 1–6 → first end-to-end with memory only (spec §17 recommendation) → measure
whether memory helps → stages 7–9 → full sweep → freeze → test.
