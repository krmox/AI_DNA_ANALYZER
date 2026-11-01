# nindels — Indel Pattern Engine (C++20)

An independent indel caller: candidates are discovered from the BAM, then each
candidate goes through the Pattern Engine (spec §13 `PROCESS_INDEL`).
Experimental. Benchmarks against bcftools on HG004 chr20 and HG001 chr19 are in `REPORT.md`, `EXPERIMENTS.md` and `results*`.

## Build

```bash
cmake -S . -B build -DHTSLIB_ROOT=$HOME/opt/nindels-tools   # omit HTSLIB_ROOT if htslib is in pkg-config
cmake --build build -j
./build/nindels-tests                                        # unit tests of every formula
```

Without htslib, the core library, the tests and the table subcommands
(`build-memory`, `fit`, `select-tau`) still build. `call` needs htslib.

## Where each formula lives

| Spec | Formula | Code |
|---|---|---|
| §2, §13.1–2 | normalization; P_i = (t, L, m, c, h, r) | `src/pattern.cpp` `left_normalize`, `extract_pattern` |
| §2.1, §3, §4 | E_i; f = n_alt/(n_alt+n_ref) (missing if 0/0); L_depth = ln(1+n_alt+n_ref); strand, read position, clipping, MAPQ, BQ, alignment | `src/evidence.cpp` `summarize_evidence` |
| §5.1 | d_type, d_length = min(1, \|ln(1+L_i)−ln(1+L_j)\|/s_L), d_sequence = Lev/max len; d_motif, d_hp, d_repeat | `src/similarity.cpp` `feature_distances` |
| §5.2–5.3 | D_ij = Σw d / Σw; S = 1 − D | `src/similarity.cpp` `pattern_distance`, `pattern_similarity` |
| §6, §14 | M_c = (μ_c, n⁺, n⁻, metadata); leader clustering; medoid; bucket index | `src/memory.cpp` `PatternMemory::build` |
| §7 | N⁺ = Σ S n⁺, N⁻ = Σ S n⁻, P_context = (N⁺+α)/(N⁺+N⁻+α+β) | `src/memory.cpp` `PatternMemory::query` |
| §7.1 | N_eff = W²/Σ S²(n⁺+n⁻) (Kish) and the distinct-locus version | `src/memory.cpp` `PatternMemory::query` |
| §8 | P(R\|H) PairHMM; ℓ(H) = Σ ln P(R\|H); Δℓ = ℓ(H_alt) − ℓ(H_ref) | `src/likelihood.cpp` `pairhmm_loglik`, `delta_ell_original` |
| §8 (correction) | diploid GL(g) = Σ ln[(2−g)/2 P(R\|H_ref) + g/2 P(R\|H_alt)], g ∈ {0/0, 0/1, 1/1} | `src/likelihood.cpp` `diploid_gl` |
| §9 | U = 1 − \|2P_base − 1\|; R = a_U U + a_C C + a_D D; τ_review(A) = τ_s − A(τ_s − τ_l) | `include/nindels/reassess.hpp` |
| §10 | C_r = v₁H + v₂R + v₃L + v₄G | `src/reassess.cpp` `region_complexity` |
| §11 | Q_k = w₁B + w₂M + w₃A + w₄S, evaluated against complete haplotypes | `src/reassess.cpp` `nucleotide_scores` |
| §12 | z = β₀ + Σβ X; P_final = 1/(1+e^−z); missingness indicators; PASS / LowProb / UNCERTAIN | `src/model.cpp`, `src/engine.cpp` |
| §13 | PROCESS_INDEL steps 1–17 | `src/engine.cpp` `process_indel` |

Every arbitrary constant is marked `ARBITRARY` in the headers and listed in
`SPEC_INTERPRETATION.md` §12.

## Run

```bash
# 1. candidate table with truth labels, every candidate reassessed (training data)
nindels-call call --bam train.bam --ref GRCh38.fa --region chr20 --truth truth.vcf.gz --bed conf.bed \
  --force-reassess --out-table train.tsv --threads 8
# 2. pattern memory from labelled training candidates (two-fold by 2-Mb block: --parity even|odd)
nindels-call build-memory --table train.tsv --out memory.tsv
# 3. table with memory features, then fit a variant (preset = feature subset)
nindels-call call ... --memory memory.tsv --force-reassess --out-table train_mem.tsv
nindels-call fit --table train_mem.tsv --preset full --out model.tsv
# 4. score validation, choose tau_call there
nindels-call call --bam val.bam ... --memory memory.tsv --model model.tsv --out-table val.tsv
nindels-call select-tau --table val.tsv
# 5. calls
nindels-call call --bam test.bam --ref GRCh38.fa --region chr19 --memory memory.tsv --model model.tsv \
  --A 0.5 --tau-call <from step 4> --out-vcf calls.vcf --threads 8
```

Variants are feature presets for `fit`:

| Preset | Features |
|---|---|
| `spec-x` | spec §12 X_i verbatim, including the original Δℓ |
| `full` | all evidence + memory + C_r + diploid HMM |
| `no-memory` | the same without the pattern memory |
| `dl-orig` | the original Δℓ in place of the diploid GL |
| `full-neff-loci` | distinct-locus N_eff in place of Kish |
| `evidence-only` | read evidence only |
| `stage1` | everything except the reassessment features (used as P_base) |

Shadow mode: `--candidates-vcf other.vcf.gz` evaluates another caller's indels
instead of discovering candidates. Its recall is not nindels' recall.
