# V2 programme — Routing, rescue and loss analysis (chr20, HG002, 300x)

Status: **analysis of frozen V1 outputs only.** No V1 code, threshold, cache, VCF, hap.py output or
manuscript was modified. Nothing here is a V2 result — this part needs no native engine and was done
while the root disk was full (V2 build blocked, see end).

Date: 2026-09-24. Scope: items 10–14 of the V2 brief.

## 0. Reproducibility

| Artefact | Path |
|---|---|
| Category definitions (written before per-locus classification) | `research/benchmark_chr20_300x/routing_analysis/CATEGORY_PREREGISTRATION.md` |
| Analysis script | `research/benchmark_chr20_300x/routing_analysis/routing_analysis.py` |
| Per-locus table (3,491 loci: 707 routed ∪ 892 rescuable ∪ 226 frame true-SNP losses ∪ 1,906 PB≠Binomial ∪ 1,593 cascade≠PB-only ∪ 1,000 random) | `routing_analysis/loci_of_interest.tsv` |
| All 1,321 hap.py SNP FN of the cascade, with "in locus table?" status | `routing_analysis/happy_snp_fn_cascade.tsv` |
| Population routing-dependence tables (all 56,242,697 frame loci) | `routing_analysis/routing_dependence.tsv` |
| Run summary + PB recompute check | `routing_analysis/stage1_summary.json` |

Inputs (read-only): `/mnt/archive/AI_DNA_ANALYZER_benchmark/ai_dna_analyzer_run/cache/chr20_300x_{counts,pb_llr,reads}.npz`,
V1 VCFs, hap.py outputs in `research/benchmark_chr20_300x/happy/`, NISTv4.2.1 truth + BED, GIAB
stratifications v3.1 (`data/strat_v31/{alldifficult,lowmap_segdup,tandemrepeats}.bed.gz`).
Runtime 5 min 16 s, peak RSS 17.8 GB.

Gate at script start (all passed, asserted in code): the script recomputes from the caches and
reproduces the frozen `pipeline_meta_chr20_300x.json` exactly — frame 56,242,697; routed 707;
rescuable 892; rescued 191.

## 1. Provenance finding: which PB variant produced the cached `pb_llr`

The streaming PB script used for the V1 300x run was **not saved** (the benchmark report describes it
but the file does not exist in the repo or archive). `poisson_binomial_llr` has two non-legacy modes
that differ at depth > 48 (`alt_index` recount vs `min(k, n_counted)`).

Recomputed on all 3,491 selected loci with the frozen functions from the read tensor:

| Variant | Loci | Bit-identical to cache | Mismatches | max |Δ| |
|---|---|---|---|---|
| `alt_index=candidate_alt(...)` | 3,491 | **yes** | 0 | 0.0 |
| `min(k, n_counted)` (no alt_index) | 3,491 | no | 2,400 | 453.1 |

→ The V1 300x `pb_llr` = frozen `poisson_binomial_llr(evidence, k, alt_index=alt_index)`. **This is
the reference V2 must reproduce.** (Evidence on 3,491 loci, not a proof over all 56.27M; the full
equivalence run of V2 will close that.)

## 2. Denominator note: "564 unnecessary routes"

No count of 564 exists in the frozen outputs. `unnecessary_pb_rate = 0.5644` is a rate over routed
loci; the frozen definition `routed & (Binomial correct)` gives **399** loci
(= 277 NEUTRAL_BOTH_CORRECT + 122 HARMFUL). The complement of rescues is 707 − 191 = **516**.
Both are reported below; 564 is not used.

## 3. The 707 routed loci (item 10)

Pre-registered, mutually exclusive outcomes:

| Outcome | n | Sub-type | hap.py query status (cascade) |
|---|---|---|---|
| RESCUE (B wrong, PB right) | **191** | 184 FP_AVOIDED, 7 FN_AVOIDED | 184 not emitted, 7 TP |
| HARMFUL (B right, PB wrong) = frozen "false rescue" | **122** | 122 FP_CREATED, 0 FN_CREATED | 122 FP |
| NEUTRAL_BOTH_CORRECT | **277** | — | 267 not emitted, 7 TP, 3 FP* |
| NEUTRAL_BOTH_WRONG | **117** | — | 110 FP, 7 not emitted |

\* 3 loci the detection layer scores correct but hap.py scores FP: label-level vs hap.py allele/genotype
matching differ; listed in the TSV.

Net effect of routing on this dataset: 191 fixed − 122 broken = **+69 detection-layer decisions**,
which is the whole cascade-vs-binomial difference (all other 56.24M loci use the Binomial call).

Pre-declared flags over the 707: DEPTH_CAPPED (depth > 48, PB saw a deterministic 48-read
subsample) **658 (93.1%)**; NEAR_PB_BOUNDARY (|PB−10.5| < 1) 51; NEAR_ROUTER_EDGE 75;
GENOTYPE_DIFF **0** (as expected structurally: Method C genotypes from counts, not LLR).

### Medians (IQR in TSV) by outcome

| Feature | RESCUE | HARMFUL | NEUTRAL_OK | NEUTRAL_WRONG |
|---|---|---|---|---|
| depth | 207 | 217.5 | 194 | 189 |
| VAF, full pileup (k/n) | 0.131 | 0.125 | 0.13 | 0.13 |
| VAF inside PB's 48-read subsample | **0.111** | **0.176** | 0.11 | 0.18 |
| ALT reads in subsample (k_tensor) | **4** | **6.5** | 3 | 6 |
| counted reads in subsample | 36 | 37 | 34 | 37 |
| mean BQ, ALT reads (subsample) | 18.3 | 18.7 | 17.9 | 18.7 |
| mean BQ, REF reads (subsample) | 33.6 | 34.1 | 33.5 | 34.4 |
| MAPQ (mean, full) | 70 | 70 | 70 | 70 |
| Binomial LLR | 9.25 | 4.33 | 4.25 | 9.70 |
| PB LLR | 2.75 | 18.06 | 1.97 | 17.73 |
| in GIAB alldifficult | 46.6% | 45.1% | 48.4% | 51.3% |
| in lowmap/segdup | 12.6% | 19.7% | 15.9% | 18.8% |

Descriptive observations (no causal claim):

1. **Routed loci are almost all low-VAF, low-BQ-ALT sites.** 648/707 (91.7%) have full-pileup VAF in
   [0.10, 0.15). At 300x with fixed ε = 0.01, Binomial LLR ≈ 7 ± 5.4 corresponds to that narrow VAF
   band. ALT-supporting reads have mean BQ ≈ 18 vs ≈ 34 for REF reads in every outcome group.
2. **BQ and MAPQ do not separate RESCUE from HARMFUL** (ALT-BQ gap > 10 in 78% vs 76%; MAPQ is 70,
   NovoAlign's maximum, in all groups).
3. **What separates them is the ALT count inside PB's 48-read subsample.** RESCUE/NEUTRAL_OK
   subsamples have lower VAF than the full pileup (0.11 vs 0.13); HARMFUL/NEUTRAL_WRONG subsamples have
   higher VAF (0.18 vs 0.13). PB's decision tracks the subsample; Binomial's tracks the full pileup.
   This is consistent with — but does not prove — that at 300x PB outcomes on these loci are driven
   substantially by which 48 reads the MAX_READS cap retains. Testing that requires PB on all reads,
   i.e. changing MAX_READS → **MODEL_CHANGE_REQUIRED**; not done.

## 4. The 191 rescues (item 11)

- 184 / 191 are **FP avoidance**: truth = ref, Binomial called (median LLR 9.25, VAF 0.131, depth
  207), PB rejected (median PB 2.75). 7 / 191 are **FN avoidance** (all truth 0/1, VAF 0.12–0.13,
  Binomial LLR 1.9–6.3, PB 11.2–38.3; positions in TSV and §3).
- Why Binomial "erred" (descriptive): at VAF ≈ 13% and depth ≈ 200, ~26 non-ref reads far exceed what
  ε/3 = 0.33% predicts, so the fixed-ε Binomial calls. Those reads have BQ ≈ 18 (p_err ≈ 1.6%);
  PB weighs them by their own error rate *and* only sees ~4 of them in its 48-read subsample.
- Summary: median depth 207; median VAF 0.131; median BQ (all counted) 31.7, ALT 18.3, REF 33.6;
  MAPQ 70 (IQR 70–70); alldifficult 46.6%, lowmap/segdup 12.6%, tandem repeats 19.9%.
- The same features do not distinguish rescues from harmful routes (§3); the rescue population is
  not characterised by a genomic context different from the other routed groups at this n.

## 5. True-SNP losses (item 12)

### 5a. Most cascade FN are not model decisions

hap.py cascade SNP FN = 1,321. Mapping every FN to the V1 locus table:

| Class | n | Share |
|---|---|---|
| Truth SNP position **absent from the extracted locus table** | **1,041** | 78.8% |
| – of which in a BED interval ≥ 64 bp (window-edge truncation) | 1,033 | |
| – of which in a BED interval < 64 bp | 8 | |
| In table, cascade did not emit a record | 242 | 18.3% |
| In table, cascade emitted a record hap.py scored FP (allele/GT mismatch) | 38 | 2.9% |

The locus table covers 56,266,816 of 56,916,239 chr20 BED bp (98.86%). The dropped 1.14% holds
1,041 truth SNPs that no V1 arm can call. This is the previously documented **M-1 window-vs-BED
artefact** (64-bp windows are only kept whole), not a Binomial/PB/router behaviour. Those FN are in
easy regions far more often (alldifficult 25.0%) than the in-table FN (93.2%).
Implication for V2: reproducing V1 exactly means reproducing this truncation. Fixing it would change
the evaluated locus set → **MODEL_CHANGE_REQUIRED** (extraction-scope change); to be offered as a
separate, labelled V2 variant, never mixed into the equivalence run.

### 5b. Frame-level loss groups (truth SNPs inside the locus table, n = 226 with a loss + 7 rescued)

| Group | n | median depth | median VAF | median BQ | median MAPQ | median Binom LLR | median PB LLR | alldifficult | routed |
|---|---|---|---|---|---|---|---|---|---|
| A LOST_BOTH | 198 | 108 | **0.008** | 36.3 | 57.6 | −42.4 | −14.9 | 98.0% | 3.5% |
| B LOST_CASCADE_ONLY | 19 | 146 | 0.100 | 36.6 | 54.5 | −21.3 | 18.6 | 100% | 0% |
| E LOST_PB_ONLY | 9 | 129 | 0.154 | 34.6 | 65.9 | 24.3 | 3.5 | 88.9% | 0% |
| C RESCUED_BY_PB | 7 | see §4 | 0.12–0.13 | | | 1.9–6.3 | 11.2–38.3 | | 100% |
| D FALSE_RESCUE_LOSS | **0** | — | — | | | | | | |

- Group D is empty: **no true SNP was lost because the router sent it to PB** on this dataset. All
  122 harmful routes are FP creations.
- Low-VAF (< 0.25, pre-declared): A 157/198, B 19/19, E 9/9. Group A is dominated by near-zero VAF in
  low-mappability/segdup sequence (97.5% lowmap_segdup, MAPQ median 57.6): evidence absent from the
  pileup, not mis-weighed. Truth GT: A 157 het / 41 hom-alt.
- The 19 B-group losses are not routed (Binomial LLR far below the band); only PB-only recovers them.
  B and E are n < 30 → **UNDERPOWERED** for any rate.

## 6. The non-rescue routes (item 13)

Of 707 routed: 516 non-rescue = 277 NEUTRAL_BOTH_CORRECT + 117 NEUTRAL_BOTH_WRONG + 122 HARMFUL.

| Pre-registered question | Answer |
|---|---|
| Binomial and PB gave the same decision | 394 (277 both correct, 117 both wrong) |
| PB changed the result without benefit | 122, all FP_CREATED |
| FP avoidance (outside rescues) | 0 — by definition every PB-changed decision is either RESCUE or HARMFUL |
| FP creation | 122 |
| Genotype-only difference | 0 (structural) |
| Near-threshold (|PB−10.5|<1), overlapping flag | 51 of 707 (TSV per outcome) |
| Other | none |

Frozen "unnecessary" (399) = 277 + 122.

## 7. Routing dependence (item 14)

Population: all 56,242,697 frame loci; full tables in `routing_dependence.tsv`.

- **VAF:** routing probability is 0 below VAF 0.10, **0.22 for VAF ∈ [0.10, 0.15)** (648 routed of
  2,946), 0.017–0.027 for [0.15, 0.35), ≈ 2e-4 for het/hom bands. Only the [0.10, 0.15) bin has
  ≥ 30 routed loci; every other VAF bin is **UNDERPOWERED** for outcome rates.
- **Depth:** P(route) falls monotonically from 1.0e-3 (< 50x) to 2.1e-6 (300–350x). Rescue share of
  routed per depth bin (bins with ≥ 30 routed): 10% (<50x), 19%, 21%, 31%, 31%, 26%, 32% (300–350x);
  ≥ 350x UNDERPOWERED.
- **BQ (mean over counted bases):** P(route) falls from 3.1e-2 (BQ 20–25) to 1.1e-6 (36–38).
  Bins < 20, 38–40, > 40 UNDERPOWERED.
- **MAPQ:** 99.7% of loci have mean MAPQ ≥ 61 (NovoAlign cap 70); every lower bin has < 30 routed →
  **UNDERPOWERED**. No MAPQ dependence can be stated from this BAM.
- **Stratification:** P(route) is 6.3× higher inside alldifficult (4.7e-5 vs 7.5e-6), 5.6× inside
  lowmap/segdup, 5.5× inside tandem repeats; outcome mixes inside vs outside are similar
  (rescue 26% vs 28%, harmful 16% vs 18% for alldifficult).

All of these are associations over one sample/one chromosome; no causal claim is made.

## 8. What this means for V2 (engineering, not science)

1. The cascade output depends on PB only at the 707 routed loci. V2 can compute PB for the cascade on
   routed loci only (identical math, identical result); full-genome PB is needed only for the
   PB-only control arm. This must be reported as such — not as a PB speed-up.
2. The V2 equivalence target is `poisson_binomial_llr(..., alt_index=...)` (§1).
3. The 3,491-locus set in `loci_of_interest.tsv` (includes all 707 routed, 191 rescues, all frame
   true-SNP losses, all PB/cascade disagreements, 1,000 seeded-random loci, seed 20260924) is the
   pre-full-run correctness gate for V2 (brief item 16).
4. MAX_READS = 48 and the 64-bp window scope are the two places where the data suggest scientific
   limitations; both are **MODEL_CHANGE_REQUIRED** and stay out of V2's equivalence build.

## Blocked / not done yet

- V2 native engine, profiling, benchmarks, hap.py re-run: **blocked** — root filesystem has 87 MB free;
  archive disk is mounted **read-only** (`mkdir` → "Read-only file system"), so it cannot host the
  build/worktree either.
- 1,000-random-locus and routed-locus BQ/MAPQ values above come from the V1 read tensor (≤ 48 reads);
  "full" columns come from the count matrix (all reads).
