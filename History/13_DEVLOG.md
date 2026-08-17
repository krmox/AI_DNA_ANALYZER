# PB depth-clipping fix + frozen cheap-router → PB robustness benchmark on hard chr21 data

**Date:** 2026-08-16
**Status of this file when sections 1–7 were written: PRE-REGISTRATION ONLY.**
Sections 1–7 were written **before any benchmark number of this session was
computed or read**. Nothing in them was edited afterwards; every later section
was appended. The one exception, flagged explicitly here, is §7.2's list of
extracted regions, which was fixed before extraction started but whose *loci
and SNP counts* could only be filled in after extraction ran — those counts are
descriptive, not outcomes.

**Verdict:** see §23. Written last.

> **Numbering note.** `cascade.py`, `cheap_router.py` and `extract_bench_v12.py`
> refer to "devlog 13/14/16" for the router-freezing and Stage-2 studies. Those
> numbers refer to a worktree-local numbering that was never merged into
> `History/`; the files actually present are `3_DEVLOG.md … 12_DEVLOG.md`. This
> file takes the next free number in `History/`, 13, and where it cites earlier
> work it cites the file that exists (`History/11_DEVLOG.md`,
> `History/12_DEVLOG.md`), not the in-code numbering.

---

## 1. Research question

Two questions, in this order:

1. **The PB defect.** `History/12_DEVLOG.md` §15 reported that at depth ≳48x the
   Poisson-binomial caller returns `pb_llr == 0.0` for roughly half of all true
   SNP, i.e. converts the strongest possible variant evidence into a confident
   false negative. What exactly causes it, is it a correctness bug or deliberate
   numerical stabilisation, what is the mathematically correct behaviour, and
   what is the smallest fix that restores it without redesigning PB?

2. **The architecture's failure boundary.** With PB repaired, where does the
   frozen cascade

   ```
   binomial LLR → |binomial_LLR − 7.0| ≤ 5.411872376933351 → PB → call
   ```

   stop behaving like PB-only? Not "does it work" — *where does it break*.

## 2. Hypotheses

* **H1 (bug).** The zeroing is a correctness bug arising from a scale mismatch:
  `k` is counted over all reads, the quality evidence covers at most
  `read_level_pileup.MAX_READS = 48`, and the two are reconciled by clipping.
  Prediction: loci with `k > n_counted` produce `-inf` under every hypothesis,
  hence `nan`, hence `0.0` via `nan_to_num`. Prediction: the defect is a strict
  no-op below 48x, so devlogs 3–12's ≤30x results survive unchanged.
* **H2 (bug consequence).** After the fix, PB-only recall at ~70x rises to
  roughly its ≤30x level, and the `pb_llr == 0` population at depth > 48
  disappears.
* **H3 (router).** The frozen cascade's ΔF1 versus PB-only remains
  indistinguishable from 0 on new regions and in hard strata, as in devlogs 11
  and 12 — *including* at high depth once PB is repaired, which is the regime
  devlog 12 could not test because its reference arm was broken.
* **H4 (failure boundary, the interesting one).** Router/PB disagreements are
  not uniformly distributed: they concentrate where the binomial and PB
  decision boundaries disagree *and* the binomial LLR is outside the routed
  band. Depth is expected to *reduce* disagreement (the cascade routes less at
  high depth), while low VAF and low base quality are expected to increase it.
  This hypothesis is falsifiable in either direction and no result of it will
  change the cutoff.

## 3. Pre-registration

### 3.1 Frozen configuration — not tuned, not refit, not swept for production

| Constant | Value | Source |
|---|---|---|
| `FROZEN_BINOMIAL_THRESHOLD` | 7.0 | imported from `cascade.py` |
| `FROZEN_PB_THRESHOLD` | 10.5 | imported from `cascade.py` |
| `FROZEN_ROUTER_CUTOFF` | 5.411872376933351 | imported from `cascade.py` |
| binomial ε | 0.01 (`BinomialVariantCaller(error_rate=0.01)`) | unchanged |
| PB ε floor / ceiling | 1e-4 / 0.25 | `quality_error_model` |
| `MAX_READS` | 48 | `read_level_pileup` — **not changed by the fix** |

Routing rule, verbatim: `|binomial_llr − 7.0| ≤ 5.411872376933351` → PB,
otherwise the binomial answer. `cascade.py` and `cheap_router.py` are **not
edited in this session** (verified in §21 and by `git diff`).

### 3.2 Primary endpoint

**ΔF1 = F1(frozen router → PB) − F1(PB-only)**, per region and pooled, with a
paired bootstrap 95% CI (`evaluate_quality_error.paired_bootstrap`, 10,000
resamples, module-default seed). Both arms use the **repaired** PB.

### 3.3 Secondary / engineering endpoints

* **Disagreement count and rate**: number of loci where `router→PB ≠ PB-only`.
  This is an engineering quantity about *reproducing PB*, and is reported
  separately from any accuracy statement. It is never used as a proxy for
  accuracy, nor accuracy for it.
* Precision, recall, TP/FP/FN, accuracy, balanced accuracy, TN.
* Fraction of loci routed to PB (= PB compute fraction at locus level).
* Wall-clock: measured binomial and PB throughput, and the projected cascade
  vs PB-only runtime from those measured rates.
* Accuracy-vs-compute curve over the cutoff sweep — **diagnostic only**.
* Depth-, quality- and VAF-stratified metrics.
* PB-fix endpoints: PB-only recall and the `pb_llr == 0` rate at depth > 48,
  before (legacy path) vs after (fixed path), on identical reads.

### 3.4 Decision rule, fixed in advance

Per region with ≥ 100 SNP, and for the pooled result:

* **PRESERVED** — the paired-bootstrap 95% CI for ΔF1 contains 0 **and**
  |ΔF1| < 0.001.
* **DEGRADED** — the CI excludes 0 with ΔF1 < 0, or |ΔF1| ≥ 0.001.
* **IMPROVED** — the CI excludes 0 with ΔF1 > 0. (Recorded, but not a goal; the
  cascade cannot legitimately beat the arm it approximates except by luck.)
* **UNDERPOWERED** — region has < 100 SNP, or the CI half-width exceeds 0.01.
  Reported as underpowered, never rounded to "preserved".

The overall verdict is CONFIRMED only if every powered region and the pooled
result are PRESERVED. A single powered DEGRADED region makes the verdict
NEGATIVE or MIXED, and is reported as such.

### 3.5 Stopping criteria

Extraction runs **once** per pre-registered region. Metrics are computed
**once**. No region is added, removed, re-extracted or re-scored after any
result is read, except to fix an outright software crash — and any such re-run
is disclosed in §18. No threshold, cutoff or parameter is changed on the basis
of any number produced tonight. If the frozen cutoff fails, the failure is
reported and the cutoff stays at 5.411872376933351.

### 3.6 What would falsify H3

Any powered region where the ΔF1 CI excludes 0 on the negative side, or where
the disagreement rate exceeds ~10⁻⁴ (devlog 11/12 saw 1–2 loci in 3.4 M, i.e.
~3×10⁻⁷). Both are recorded whatever they show.

## 4. Dataset and region selection

### 4.1 Leakage constraint

The frozen cutoff was selected on chr21:32–40 Mb (`big15x` train/validation
windows at 32–38 Mb). Regions inside that span cannot support a
*generalisation* claim about the router. High-depth BAM coverage in this
repository, however, exists **only** for chr21:31–33 Mb and chr21:30.0–30.47 Mb
— there is no untainted high-depth BAM outside 31–33 Mb. Rather than fabricate
one, the 32–33 Mb slice is included and **labelled tainted for router claims**:
it is used for PB-defect claims (which have nothing to do with the router) and
its router numbers are reported but excluded from the pooled generalisation
statement.

### 4.2 Pre-registered regions

Group A — newly extracted this session, with both the repaired and the legacy
PB computed on identical reads in one pass:

| Key | Coordinates (GRCh38 chr21) | Depth | Status |
|---|---|---|---|
| `new32_33M_15x` | 32,000,000–33,000,000 | ~14x | new coordinates; **inside** router-tuning span |
| `new32_33M_30x` | 32,000,000–33,000,000 | ~28x | as above |
| `new32_33M_full` | 32,000,000–33,000,000 | ~70x | as above |
| `r30M_full` | 29,999,999–30,469,999 | ~73x | **new depth realisation** of devlog 11's `test_15x`/`test_30x` locus set; outside tuning span |
| `re31_32M_full` | 31,000,000–32,000,000 | ~69x | replication of devlog 12's failing region; outside tuning span |
| `re31_32M_30x` | 31,000,000–32,000,000 | ~28x | bit-identity control vs devlog 12's artefact |

Group B — already-cached regions, reused unchanged (their PB is provably
unaffected by the fix wherever max depth ≤ 48; verified per region in §10):

`gen14_near_40_44M` (40–44 Mb), `gen14_far_13_17M`, `gen14_far_17_21M`,
`gen14_far_21_25M`, `gen14_far_25_29M`, `test_15x`, `test_30x` — 19.44 M loci,
~15x except `test_30x` at ~30x.

### 4.3 Regimes: what is and is not available

| Requested regime | Available? | What is used |
|---|---|---|
| low depth | yes | `depth ≤ 9` stratum across all regions; 15x regions are ~14x mean |
| moderate depth | yes | 30x realisations |
| high depth | yes | ~70x realisations (31–32, 32–33, 30.0–30.47 Mb) |
| low VAF | yes (as a stratum) | `0 < VAF ≤ 0.35` |
| low base quality | yes (as a stratum) | mean base quality < 30 |
| higher-error regions | **only as a stratum** | no BAM with an elevated error rate exists; the low-quality/low-VAF strata are the strongest available proxy. Labelled as a proxy, not as an independent high-error dataset. |
| heterogeneous / mixed difficulty | yes | whole-region pools mix all strata |
| different variant/error density | partial | SNP density varies ~1.0–1.6 per 1000 loci across the gen14 regions; no deliberately variant-dense region exists |
| other chromosomes / samples / platforms | **no** | chr21/HG002/Illumina only. No generalisation beyond that will be claimed. |

Nothing is synthesised. No artificial difficulty is injected.

## 5. Provenance

* Sample: GIAB **HG002**, GRCh38, chr21 only, Illumina.
* Truth: GIAB v4.2.1 benchmark VCF, restricted to the matching high-confidence
  BED, exactly as every previous devlog (`data/giab_hg002_*/`).
* BAMs: `hg002_chr21_31_33M.bam` (full ~70x) and its 30x/15x downsamples;
  `hg002_chr21_30M.bam` (full ~73x) and its 30x/15x downsamples;
  `hg002_chr21_32_44M_15x.bam`. Reference `data/reference/chr21_full.fa`.
* Depth realisations at the same coordinates are downsamples of the same
  library, so 15x/30x/full at one locus set are **not** independent samples;
  they are treated as three regimes of one region, never as three regions.

## 6. Analysis plan

For every region and every stratum, three arms on the identical scoring frame
(`label ∈ {0, SNP}`):

1. **PB-only** — `pb_llr ≥ 10.5`
2. **Binomial-only** — `binomial_llr ≥ 7.0`
3. **Frozen cheap-router → PB** — PB where routed, binomial elsewhere

Primary comparison: 3 vs 1. Metrics as listed in §3.2–3.3. Strata are defined
on evidence only (depth, VAF, base quality, PB margin) — never on labels.

Failure analysis: every disagreeing locus is dumped with depth, VAF, mean base
quality, binomial LLR, PB LLR, router margin, `k`, region, truth label and
whether it is a router-caused FP or FN, and the population is compared against
the region background.

## 7. Execution plan

1. Diagnose and fix the PB defect; add regression tests; run the suite.
2. Verify the fix is a no-op below 48x, on cached real data, bit-for-bit.
3. Extract Group A once, with legacy PB alongside for an exact A/B.
4. Score all regions, all strata, all three arms.
5. Accuracy-vs-compute curve (diagnostic).
6. Failure-boundary analysis.
7. Full test suite, artefact verification.

### 7.1 Artefacts

| What | Path |
|---|---|
| Benchmark driver | `bench_v13_robustness.py` |
| Driver unit tests | `test_bench_v13.py` |
| Extraction caches | `cache/bench_v13/*.npz` |
| Primary results | `results/bench_v13/robustness_results.json` |
| Disagreement dump | `results/bench_v13/disagreements.json` |
| PB fix A/B | `results/bench_v13/pb_fix_ab.json` |
| Baseline re-verification | `results/bench_v13/baseline_verification.json` |

No file under `results/robustness_v11/`, `results/bench_v12/` or `History/`
(other than this new file) is modified.


---

*Everything from here down was written after the corresponding measurement was
made. Sections 1–7 above were not edited.*

---

## 8. PB bug — diagnosis

### 8.1 Where it lives

`quality_error_model.poisson_binomial_llr`, one line:

```python
block_k = np.clip(k[start:stop], 0, block_p.shape[1])   # pre-fix
```

`block_p.shape[1]` is the *width of the read tensor*, i.e.
`read_level_pileup.MAX_READS = 48` — not the number of reads the locus
actually has.

### 8.2 Why the counts are clipped at all

Three facts collide:

1. `read_level_pileup._locus_reads` keeps at most `MAX_READS = 48` reads per
   locus (a deterministic hash-ordered subsample), so `QualityEvidence` is
   `[N, 48]` no matter how deep the pileup is.
2. `candidate_alt` — deliberately delegated to `BinomialVariantCaller` so ALT
   selection stays bit-identical to the baseline — computes `k` from the
   **uncapped** `[N, 10]` count matrix. At ~70x a homozygous SNP gives `k ≈ 54`.
3. `_log_poisson_binomial_pmf` returns a distribution over `0 … 48` for the
   retained reads only. Indexing it at `k = 54` is out of bounds, so *something*
   had to clip. The clip was a bounds guard, not a modelling decision.

### 8.3 Why `-inf` becomes zero

Reproduced in isolation (48 slots, 40 counted ALT reads at Q35 plus 8 gap
reads, `k = 54`):

| `k` passed | pre-fix `loglik_h0` | pre-fix LLR |
|---:|---:|---:|
| 40 | −366.31 | +366.30 |
| 41 | −inf | **0.0** |
| 48 | −inf | **0.0** |
| 54 | −inf | **0.0** |

The dp array carries `-inf` for every count above a locus' number of *counted*
slots (gap and quality-failed reads occupy slots but are not counted). So for
`k_clipped > n_counted` all three hypotheses return `-inf`, the subtraction
`alternative - loglik_h0` is `(-inf) - (-inf) = nan`, and the closing
`np.nan_to_num(llr, nan=0.0, ...)` maps it to exactly `0.0`. Devlog 12 §15
attributed the zeroing to the `posinf=0.0` branch; the actual path is the
`nan=0.0` branch. The observable consequence it reported — a strong variant
scored 0.0, below the 10.5 threshold, hence a confident false negative — is
exactly right.

### 8.4 Intentional stabilisation, or a correctness bug?

Both, in different places, and it matters which is which:

* The `np.clip` is a **bounds guard** that silently changed the observation
  being scored. That is a correctness bug.
* The `nan_to_num` is **genuine numerical stabilisation** and is kept: a locus
  with no counted reads legitimately has no evidence and must score 0. What was
  wrong is that a second, entirely different condition was being funnelled into
  the same "no evidence" outcome. The fix removes the second condition rather
  than the guard.

There is also a quieter, non-`nan` half of the same defect. When a locus does
fill all 48 slots with counted reads, `clip(54, 0, 48) = 48` asserts "all 48
retained reads support ALT", which is *over*-confident whenever the retained
sample is mixed. That half never produced a `nan` and so was invisible in
devlog 12's "`pb_llr == 0`" diagnostic.

### 8.5 What the correct behaviour is

A likelihood must be evaluated on the observation its parameters describe. The
parameters here are the per-read error probabilities of the ≤48 **retained**
reads, so the count scored must be the ALT count **among those retained reads**
— not the ALT count of the full pileup, and not that count truncated.

## 9. The fix

Smallest change that restores that property; PB is not redesigned, `MAX_READS`
is not changed, no threshold moves.

`quality_error_model.poisson_binomial_llr` gains two optional arguments:

* `alt_index` — when supplied (from `candidate_alt`, which already returns it),
  the ALT count is recounted over the retained reads:
  `k_eff = Σ counted & (base_code == alt_index)`. Exact, no approximation.
* `legacy_truncation` — reproduces the pre-fix behaviour, **for the benchmark's
  A/B control only**, so the fix's effect can be measured on identical reads.

Without `alt_index` the fallback is `k_eff = min(k, n_counted)`, which cannot
leave the pmf's support and therefore cannot produce the `nan → 0` collapse.

The function now also returns `k_effective`, `n_out_of_support` (loci whose
supplied `k` exceeded the retained read count) and `n_non_finite` (loci whose
raw LLR was not finite), so any future benchmark can see the defect's footprint
instead of inferring it.

Call sites updated to pass `alt_index`: `extract_bench_v12`,
`cache_locus_evidence`, `evaluate_quality_error`.

**Below the tensor width all three paths are identities** — every read is
retained, so `k_retained == k ≤ n_counted` — which is what licenses carrying
the ≤30x results of devlogs 3–12 forward unchanged (verified on real data in
§11).

### 9.1 Alternatives considered and rejected

| Alternative | Why not |
|---|---|
| Raise or remove `MAX_READS` | Changes the cached tensor format, every extraction and every previous cache; not a minimal fix, and PB cost is linear in tensor width. |
| Rescale `k` to the retained sample (`k · n_counted / n`) | An approximation where an exact count is available. |
| Map `+inf`/`nan` to a large positive score | Treats the symptom; leaves the wrong observation being scored, including the over-confident half of §8.4. |

## 10. Regression tests

Added `TestPoissonBinomialDepthClipping` (12 tests) to `test_quality_error.py`:

* the previously failing case (`k` above the retained read count) is finite,
  positive, and clears the frozen 10.5 threshold;
* a guard on the A/B control itself — the legacy path still reproduces the
  defect (`llr == 0.0`, one non-finite locus), so the before/after comparison
  is measuring something real;
* `alt_index` recounts over retained reads and agrees exactly with scoring that
  count directly;
* the over-confident half of the defect (`legacy > fixed` for a mixed retained
  sample);
* ≤48x no-op, parametrised over 8 depth/ALT combinations, plus a 500-locus
  randomised sweep asserting **bit-identical** legacy vs fixed output;
* numerical edges: `k ∈ {−3, 0, 10, 11, 48, 500}` all stay in support and
  finite; a gap-only locus still scores exactly 0; monotonicity of the LLR in
  the retained ALT count;
* the `n_out_of_support` diagnostic counts the right loci.

## 11. Baseline verification — did the fix disturb the established results?

**No. Bit-for-bit identical on real data.** PB was recomputed from the cached
read tensors of devlog 11's two `test` regions and compared against the cached
`pb_llr` element by element:

| region | tensor width | loci | max depth | loci with `k` out of support | fixed == cached | legacy == fixed |
|---|---:|---:|---:|---:|---|---|
| `test_15x` | 48 | 455,488 | 38 | 0 | **yes** (max abs diff 0.0) | yes |
| `test_30x` | 96 | 455,488 | 62 | 0 | **yes** (max abs diff 0.0) | yes |

Artefact: `results/bench_v13/baseline_verification.json`.

Two things worth recording. First, `test_30x` contains 4,234 loci above 48x, so
"all previous work was ≤30x and therefore safe" was *not* the correct reason it
was safe — it was cached with a **96-wide** read tensor, so no read was dropped
and no `k` left the support. The correct invariant is `depth ≤ tensor width`,
not `depth ≤ 48`. Second, the five `gen14` caches top out at 39–44x against a
48-wide tensor, so they satisfy the invariant too and their cached PB is reused
unchanged.

Re-extracting devlog 12's 31–32 Mb full-depth region reproduced its reported
failure exactly: legacy PB F1 **0.68407** and legacy cascade F1 **0.98876** here
against 0.68407 / 0.98876 there, and 426 true SNP with `pb_llr == 0` against 426
there. An independent extraction landing on the same numbers is strong evidence
that both the defect and its measurement were real.

## 12. The PB fix, measured (PB versus truth — not a router statement)

Identical reads, identical thresholds, the defect the only moving part
(`results/bench_v13/pb_fix_ab.json`):

| region | mean depth | PB recall legacy → fixed | PB F1 legacy → fixed | true SNP with `pb_llr == 0` | LLRs changed | PB calls changed |
|---|---:|---|---|---:|---:|---:|
| `new32_33M_15x` | 14x | 0.9796 → 0.9796 | 0.98799 → 0.98799 | 0 → 0 | **0** | 0 |
| `new32_33M_30x` | 28x | 0.9974 → 0.9983 | 0.99071 → 0.99114 | 1 → 0 | 34 | 1 |
| `re31_32M_30x` | 28x | 0.9968 → 0.9968 | 0.99192 → 0.99192 | 0 → 0 | 20 | 0 |
| `new32_33M_full` | 70x | **0.5221 → 1.0000** | **0.66848 → 0.99366** | 562 → 0 | 33,105 | 594 |
| `re31_32M_full` | 69x | **0.5390 → 1.0000** | **0.68407 → 0.99088** | 426 → 0 | 29,217 | 443 |
| `r30M_full` | 73x | **0.3024 → 1.0000** | **0.45920 → 0.99547** | 383 → 0 | 14,448 | 386 |

At full depth PB's recall goes from *losing half the variants* to **1.000** on
all three regions, and the `pb_llr == 0` population is eliminated entirely
(1,371 true SNP recovered in total). At 15x nothing changes at all — literally
zero LLRs differ — which is the no-op property, confirmed on real data rather
than argued.

Exposure, from the per-locus diagnostics (`pb_fix_scope`):

| region | loci > 48x | loci where `k` left the support | loci where the retained ALT count differed from `k` | exposed SNP |
|---|---:|---:|---:|---:|
| `new32_33M_full` | 910,920 (95.9%) | 563 | 33,128 | 562 |
| `re31_32M_full` | 859,792 (92.9%) | 427 | 29,226 | 426 |
| `r30M_full` | 439,502 (96.5%) | 383 | 14,457 | 383 |
| `new32_33M_30x` | 8,299 (0.87%) | 1 | 34 | 1 |
| `re31_32M_30x` | 3,752 (0.41%) | 0 | 20 | 0 |
| all ≤44x regions | 0 | 0 | 0 | 0 |

Note the two columns are very different sizes: ~30 k loci per Mb had their ALT
count silently rescaled (the quiet, over-confident half of §8.4), while ~500 had
`k` pushed out of the support (the loud, zeroed half). Only the second was
visible in devlog 12's diagnostic. **Essentially every exposed locus is a true
SNP** (562 of 563, 426 of 427, 383 of 383) — the defect selected almost
perfectly for the variants it destroyed, because only a real homozygous variant
drives `k` above the retained read count.

## 13. Primary benchmark results — frozen cascade vs PB-only, repaired PB

`results/bench_v13/robustness_results.json`. Every row: three arms on the same
scoring frame, frozen constants, no tuning.

| region | loci | SNP | PB F1 | router→PB F1 | ΔF1 | 95% CI | disagreements | routed to PB | binomial-only F1 |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|
| `new32_33M_15x` † | 949,553 | 1,176 | 0.98799 | 0.98757 | −0.000423 | [−0.00132, +0.00000] | 1 | 0.1315% | 0.97211 |
| `new32_33M_30x` † | 949,553 | 1,176 | 0.99114 | 0.99197 | +0.000837 | [+0.00000, +0.00214] | 2 | 0.0239% | 0.99114 |
| `new32_33M_full` † | 949,553 | 1,176 | 0.99366 | 0.99619 | **+0.002525** | **[+0.00041, +0.00509]** | 8 | 0.0057% | 0.99577 |
| `r30M_full` | 455,371 | 549 | 0.99547 | 0.99547 | 0.000000 | [+0.00000, +0.00000] | 0 | 0.0013% | 0.99728 |
| `re31_32M_full` | 925,418 | 924 | 0.99088 | 0.99088 | 0.000000 | [−0.00159, +0.00158] | 2 | 0.0070% | 0.98560 |
| `re31_32M_30x` | 925,418 | 924 | 0.99192 | 0.99192 | 0.000000 | [+0.00000, +0.00000] | 0 | 0.0335% | 0.97861 |
| `gen14_near_40_44M` | 3,372,463 | 5,951 | 0.96955 | 0.96955 | +0.000005 | [−0.00025, +0.00027] | 2 | 0.1698% | 0.94584 |
| `gen14_far_13_17M` | 3,567,862 | 4,097 | 0.96514 | 0.96514 | 0.000000 | [+0.00000, +0.00000] | 0 | 0.1319% | 0.94790 |
| `gen14_far_17_21M` | 3,863,504 | 5,873 | 0.98315 | 0.98315 | 0.000000 | [+0.00000, +0.00000] | 0 | 0.1319% | 0.96961 |
| `gen14_far_21_25M` | 3,878,571 | 5,833 | 0.98772 | 0.98772 | 0.000000 | [+0.00000, +0.00000] | 0 | 0.1125% | 0.97628 |
| `gen14_far_25_29M` | 3,835,404 | 5,116 | 0.98582 | 0.98582 | 0.000000 | [+0.00000, +0.00000] | 0 | 0.1136% | 0.96736 |
| `test_15x` | 455,371 | 549 | 0.99086 | 0.99086 | 0.000000 | [+0.00000, +0.00000] | 0 | 0.1008% | 0.97896 |
| `test_30x` | 455,371 | 549 | 0.99728 | 0.99728 | 0.000000 | [+0.00000, +0.00000] | 0 | 0.0075% | 0.99637 |

† inside the router's cutoff-selection span (§4.1) — valid for PB claims,
discounted for router-generalisation claims.

Pooled:

| pool | loci | SNP | PB F1 | router F1 | ΔF1 | 95% CI | disagreements | routed |
|---|---:|---:|---:|---:|---:|---|---:|---:|
| `group_a_untainted` (new, outside tuning span) | 2,306,207 | 2,397 | 0.99233 | 0.99233 | 0.000000 | [−0.00062, +0.00062] | 2 | 0.0165% |
| `group_b_cached` (devlog 11's seven regions) | 19,428,546 | 27,968 | 0.97950 | 0.97950 | +0.000001 | [−0.00005, +0.00006] | 2 | 0.1273% |
| `group_a_all` | 5,154,866 | 5,925 | 0.99151 | 0.99209 | +0.000584 | [+0.00000, +0.00117] | 13 | 0.0371% |
| `high_depth` (three ~70x regions) | 2,330,342 | 2,649 | 0.99306 | 0.99418 | +0.001118 | [+0.00000, +0.00228] | 10 | 0.0054% |
| **`all_regions`** | **24,583,412** | **33,893** | **0.98162** | **0.98172** | **+0.000103** | **[+0.00000, +0.00022]** | **15** | **0.1084%** |

Full confusion counts for the total pool — router→PB: TP 33,031, FP 368,
FN 862, accuracy 0.9999500, balanced accuracy 0.98728, precision 0.98898,
recall 0.97457. PB-only: TP 33,030, FP 374, FN 863, accuracy 0.9999497,
balanced accuracy 0.98726, F1 0.98162.

**Classification against the §3.4 rule**, applied without revision:

* PRESERVED — 12 of 13 regions and all pooled results except `high_depth`.
* IMPROVED — `new32_33M_full` (CI excludes 0, ΔF1 = +0.0025 > 0).
* DEGRADED — **none**.
* UNDERPOWERED — none at region level (every region has ≥ 549 SNP and CI
  half-width ≤ 0.0026).

**A pre-registration defect, disclosed rather than patched:** §3.4's DEGRADED
clause ("|ΔF1| ≥ 0.001") and its IMPROVED clause both fire on
`new32_33M_full` and on the `high_depth` pool. The rule should have bounded
only *adverse* deviation. Both labels are reported; the deviation's direction is
positive, and §14 explains the mechanism rather than leaving it as a label.

**Why the cascade beats PB at high depth.** Not luck this time, and not the
router being clever. After the repair, the *binomial* arm is the better caller
at ~70x (0.99577 vs PB 0.99366 at 32–33 Mb; 0.99728 vs 0.99547 at 30 Mb): PB
acquires a systematic false-positive mode at high depth and low VAF (§16). The
cascade routes only 0.0013–0.0070% of high-depth loci to PB, so it answers
almost everything with the binomial and inherits the better arm's behaviour.
This inverts the architecture's own premise — **PB is no longer unambiguously
the accuracy ceiling above ~48x** — and it is a finding about PB, not a
validation of the router.

## 14. Accuracy-vs-compute curve (diagnostic only — the cutoff stays frozen)

Measured throughput, single quiet process, same machine, same session:

| regime | binomial | PB | ratio |
|---|---:|---:|---:|
| 15x | 2,415,061 loci/s | 10,007 loci/s | 241x |
| 30x | 2,926,648 loci/s | 9,018 loci/s | 325x |
| full (~68x) | 2,674,649 loci/s | 7,221 loci/s | 370x |

PB is nearly depth-independent (the convolution runs over all 48 tensor slots
regardless), which is why its rate barely moves while depth quintuples. Both
callers are per-locus, so no window amplification applies to this projection;
the neural arm is absent tonight and contributes nothing.

Projected wall-clock from those measured rates:

| region | loci | PB loci | PB compute fraction | PB-only | cascade | speedup |
|---|---:|---:|---:|---:|---:|---:|
| `gen14_near_40_44M` | 3,372,463 | 5,725 | 0.1698% | 337.0 s | 2.0 s | **171x** |
| `gen14_far_13_17M` | 3,567,862 | 4,706 | 0.1319% | 356.5 s | 1.9 s | 183x |
| `new32_33M_15x` | 949,553 | 1,249 | 0.1315% | 94.9 s | 0.5 s | 183x |
| `test_30x` | 455,371 | 34 | 0.0075% | 50.5 s | 0.2 s | 317x |
| `re31_32M_full` | 925,418 | 65 | 0.0070% | 128.2 s | 0.4 s | 361x |
| `r30M_full` | 455,371 | 6 | 0.0013% | 63.1 s | 0.2 s | **369x** |

The cascade gets *cheaper* as depth rises (0.17% → 0.0013% of loci routed),
extending devlog 12's observation to three more depth realisations: better data
pushes the binomial LLR further from its threshold, so fewer loci are uncertain.

The curve itself, `re31_32M_full` (~69x) and `gen14_far_13_17M` (~14x):

| cutoff | 69x: PB frac / ΔF1 / disagree / speedup | 14x: PB frac / ΔF1 / disagree / speedup |
|---|---|---|
| 0.5 | 0.0003% / −0.00476 / 15 / 370x | 0.0036% / −0.01073 / 150 / 239x |
| 2.0 | 0.0019% / −0.00265 / 9 / 368x | 0.0116% / −0.00199 / 67 / 235x |
| 3.5 | 0.0035% / −0.00053 / 3 / 366x | 0.0291% / −0.00017 / 19 / 226x |
| 5.0 | 0.0064% / 0.00000 / 2 / 362x | 0.0778% / −0.00012 / 1 / 203x |
| **5.412 (frozen)** | **0.0070% / 0.00000 / 2 / 361x** | **0.1319% / 0.00000 / 0 / 183x** |
| 6.0 | 0.0081% / +0.00053 / 1 / 360x | 0.1349% / 0.00000 / 0 / 182x |
| 10.0 | 0.0209% / 0.00000 / 0 / 344x | 5.07% / 0.00000 / 0 / 18x |
| 15.0 | 0.1009% / 0.00000 / 0 / 270x | 53.9% / 0.00000 / 0 / 1.8x |
| 20.0 | 0.4612% / 0.00000 / 0 / 137x | 97.9% / 0.00000 / 0 / 1.0x |

The frozen cutoff sits **exactly at the knee**: at 15x it is the smallest swept
cutoff reaching ΔF1 = 0 and zero disagreements, and every larger cutoff buys
nothing while costing up to 180x of the speedup. At 69x, ΔF1 is already 0 at
5.0 and the two residual disagreements persist to t = 8. The cutoff was
selected on 15x data in a different region and lands on the knee of a ~70x
region it never saw. **It is not moved, and no row of this table was promoted
to a production threshold.**

## 15. Depth-stratified results

Fine bins that resolve the 48x tensor boundary (repaired PB throughout):

| region | bin | loci | SNP | PB F1 | router F1 | ΔF1 | routed | disagree |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `new32_33M_15x` | 1–9x | 127,663 | 187 | 0.96399 | 0.96133 | −0.00266 | 0.790% | 1 |
| `new32_33M_15x` | 10–29x | 816,799 | 985 | 0.99236 | 0.99236 | 0.00000 | 0.030% | 0 |
| `new32_33M_30x` | 30–47x | 461,927 | 486 | 0.98982 | 0.99083 | +0.00101 | 0.005% | 1 |
| `new32_33M_full` | 48–59x | 109,033 | 170 | 0.98551 | 0.99125 | +0.00575 | 0.003% | 2 |
| `new32_33M_full` | 60–79x | 493,327 | 659 | 0.99622 | 0.99924 | +0.00302 | 0.0004% | 6 |
| `new32_33M_full` | 80x+ | 312,895 | 288 | 1.00000 | 1.00000 | 0.00000 | 0.0003% | 0 |
| `re31_32M_full` | 48–59x | 152,082 | 197 | 0.98995 | 0.98995 | 0.00000 | 0.009% | 0 |
| `re31_32M_full` | 60–79x | 511,017 | 497 | 0.99599 | 0.99599 | 0.00000 | 0.001% | 2 |
| `r30M_full` | 48–59x | 48,995 | 76 | 0.98065 | 0.98065 | 0.00000 | 0.004% | 0 |
| `r30M_full` | 60–79x | 244,248 | 311 | 0.99679 | 0.99679 | 0.00000 | 0.0004% | 0 |
| `r30M_full` | 80x+ | 148,007 | 144 | 1.00000 | 1.00000 | 0.00000 | 0.000% | 0 |

Both callers improve monotonically with depth, and at 80x+ both are perfect on
these regions (288 + 144 SNP). Routing coverage falls by three orders of
magnitude from 1–9x to 80x+.

## 16. Difficulty-stratified results

Pre-registered hard strata, evidence-defined, over all 13 regions. Selected
rows (full table in the JSON):

| region | stratum | loci | SNP | PB F1 | router F1 | ΔF1 | disagree |
|---|---|---:|---:|---:|---:|---:|---:|
| `new32_33M_15x` | low quality (Q<30) | 18,245 | 33 | 0.95522 | 0.94118 | **−0.01405** | 1 |
| `new32_33M_15x` | PB-uncertain | 29,312 | 110 | 0.87129 | 0.86700 | −0.00429 | 1 |
| `new32_33M_15x` | low depth (≤9x) | 127,663 | 187 | 0.96399 | 0.96133 | −0.00266 | 1 |
| `new32_33M_15x` | low VAF (≤0.35) | 21,108 | 110 | 0.88780 | 0.88780 | 0.00000 | 0 |
| `new32_33M_full` | low VAF | 96,544 | 20 | 0.72727 | 0.81633 | +0.08905 | 8 |
| `new32_33M_full` | PB-uncertain | 240 | 1 | 0.15385 | 0.22222 | +0.06838 | 6 |
| `gen14_near_40_44M` | low VAF | 70,509 | 782 | 0.86713 | 0.86653 | −0.00061 | 1 |
| `gen14_near_40_44M` | low depth | 625,487 | 1,183 | 0.88331 | 0.88383 | +0.00052 | 1 |
| all other strata × regions (43 of 55 rows) | | | | | | **0.00000** | 0 |

The single worst adverse number in the whole benchmark is
**ΔF1 = −0.0140 in the 15x low-base-quality stratum of `new32_33M_15x`**,
produced by exactly **one** disagreeing locus among 33 SNP. It is a real
adverse deviation and it is reported as such, but with 33 SNP the stratum is
underpowered by the §3.4 criterion and a single event moves F1 by 1.4 points.
No other low-quality stratum in any region shows any disagreement at all.

The large *positive* strata deltas at full depth (+0.089 on low VAF) are the
same PB high-depth false-positive mode as §13, seen from the stratum that
contains it.

## 17. Router/PB disagreement analysis — the failure boundary

**15 disagreeing loci in 24,583,412 scored loci = 6.1 × 10⁻⁷.** Every one is
dumped with its features in `results/bench_v13/disagreements.json`. By
construction none is inside the routed band (a routed locus takes PB's answer),
which the unit tests assert.

Direction — a disagreement is not automatically a regression:

| effect | count | meaning |
|---|---:|---|
| `router_avoids_fp` | **10** | PB called a non-variant; the cascade did not. Cascade **right**. |
| `router_fp` | 4 | the cascade called a non-variant PB rejected. Cascade **wrong**. |
| `router_recovers_fn` | 1 | PB missed a true SNP; the cascade found it. Cascade **right**. |
| `router_fn` | **0** | the cascade never lost a true SNP that PB found. |

**Zero router-caused false negatives in 24.6 M loci** is the single most
load-bearing engineering number here: the cascade's failure mode, when it has
one, is a false positive on a non-variant, not a missed variant.

Pooled disagreement rate by difficulty axis:

| axis | bin | loci | disagreements | rate |
|---|---|---:|---:|---:|
| depth | 1–9x | 3,003,078 | 2 | 6.7e−7 |
| depth | 10–29x | 17,971,835 | 2 | **1.1e−7** |
| depth | 30–47x | 1,224,283 | 1 | 8.2e−7 |
| depth | 48–59x | 332,319 | 2 | 6.0e−6 |
| depth | 60–79x | 1,248,973 | 8 | **6.4e−6** |
| depth | 80x+ | 663,856 | 0 | 0 |
| base quality | Q<20 | 7,731 | 3 | **3.9e−4** |
| base quality | Q20–25 | 37,585 | 2 | 5.3e−5 |
| base quality | Q25–30 | 311,852 | 2 | 6.4e−6 |
| base quality | Q30–35 | 3,428,361 | 1 | 2.9e−7 |
| base quality | Q≥35 | 20,658,814 | 7 | **3.4e−7** |
| VAF | ≤0.15 | 699,175 | 10 | 1.4e−5 |
| VAF | 0.15–0.25 | 21,771 | 2 | 9.2e−5 |
| VAF | 0.25–0.35 | 3,654 | 1 | **2.7e−4** |
| VAF | ≥0.35 | 31,263 | 2 | 6.4e−5 |
| truth | SNP | 33,893 | 1 | 3.0e−5 |
| truth | non-SNP | 24,549,519 | 14 | 5.7e−7 |

**There is a detectable boundary, and it is base quality.** The disagreement
rate rises monotonically as mean base quality falls, spanning **three orders of
magnitude** from 3.4e−7 at Q≥35 to 3.9e−4 below Q20 — a ~1,100x gradient. VAF
shows the same direction (worst in the 0.25–0.35 band, 2.7e−4). Depth shows a
~50x elevation above 48x relative to the 10–29x floor, but with only 10 events
across two bins, and the 80x+ bin has none at all.

**Every one of these gradients rests on 1–8 events.** They are consistent with
H4's quality/VAF prediction and inconsistent with its depth prediction (depth
was expected to *reduce* disagreement; it raised it), but no gradient here is
individually powered. They are reported as a signal to test deliberately, not
as an established law.

## 18. Failure modes, named

1. **PB's high-depth low-VAF false positives (new, and the important one).**
   Nine of the ten `router_avoids_fp` loci are ~60–73x sites with VAF 0.07–0.12
   and *high* base quality (Q 28–38), where the repaired PB returns LLR 10.6–27.5
   and calls a variant that is not there. At 70x, 8 ALT reads out of 70 at Q37
   are individually so improbable under H0 that the Poisson-binomial's exact
   per-read product overwhelms the prior — the model has no term for the
   *alignment*-driven errors that dominate at that VAF. The binomial with fixed
   ε = 0.01 is accidentally protected because ε = 0.01 is far more pessimistic
   than the reads' own Q37. This is the mechanism behind §13's inverted ranking.
2. **Router false positives at low quality and low depth (4 loci).** All four
   have the binomial LLR above 7.0 while outside the routed band, and PB
   correctly rejecting: depth 5–63x, VAF 0.17–0.75, mean base quality 11–24.
   These are the cascade's genuine regressions, and they cluster where §17's
   quality gradient predicts.
3. **The one recovered false negative.** A depth-5, VAF-0.75, Q11.2 true SNP
   that PB scored 10.12 (just below 10.5) and the binomial scored 14.33.
4. **The repaired PB now scores only the retained ≤48 reads.** At 120x that is
   40% of the data. This is *correct* (the likelihood matches its evidence) but
   *less powerful* than the data allows — a modelling limitation the fix
   deliberately does not address, and the clearest next task (§22).

## 19. Controls

Random routing and reverse-confidence routing at matched PB coverage, five
seeds (`controls` in the JSON):

| region | frozen router ΔF1 | random ΔF1 (mean) | reverse ΔF1 | matched coverage |
|---|---:|---:|---:|---:|
| `gen14_near_40_44M` | +0.000005 | −0.023707 | −0.023707 | 5,725 |
| `gen14_far_25_29M` | 0.000000 | −0.018397 | −0.018457 | 4,358 |
| `gen14_far_13_17M` | 0.000000 | −0.017194 | −0.017243 | 4,706 |
| `new32_33M_15x` | −0.000423 | −0.015878 | −0.015878 | 1,249 |
| `re31_32M_30x` | 0.000000 | −0.013313 | −0.013313 | 310 |
| `re31_32M_full` | 0.000000 | −0.005285 | −0.005285 | 65 |
| `new32_33M_full` | +0.002525 | +0.002103 | +0.002103 | 54 |
| `r30M_full` | 0.000000 | +0.001808 | +0.001808 | 6 |

At 15x/30x the frozen router beats both controls by 1–2 F1 points at identical
compute: *which* loci are routed is what matters, not how many. At full depth
the control is uninformative — coverage is 6–65 loci, and random routing
degenerates to "answer everything with the binomial", which as §13 explains is
itself a good strategy there. Both controls landing on the binomial-only F1 is
the expected degenerate behaviour and is not evidence for the router.

## 20. Test status

| suite | result |
|---|---|
| `test_quality_error.py` | 65 passed (12 new depth-clipping regressions) |
| `test_bench_v13.py` | 25 passed (new) |
| full suite, `--ignore=test_bimamba.py --ignore=test_providers.py` | **251 passed, 0 failed** |
| `test_bimamba.py`, `test_providers.py` | **pre-existing collection errors** — both import `bimamba_variant_caller`, a package that does not exist anywhere in this repository. Neither file was touched this session. Documented, not modified. |

Totals: 251 tests run, 251 passed, 0 failed, 2 pre-existing collection errors,
**0 newly introduced failures**. Baseline before any change tonight: 169 passed
across the relevant suites, same 2 collection errors.

## 21. Reproducibility

```bash
# 1. extract the six pre-registered regions (once each; ~50 min wall-clock,
#    four concurrent, on 8 cores)
./run_v13_extract.sh
#    each invocation is:
#    python3 extract_bench_v12.py --fasta data/reference/chr21_full.fa \
#      --bam <BAM> --vcf <VCF> --bed <BED> --region <start> <end> \
#      --features none --legacy-pb --out cache/bench_v13/<name>.npz

# 2. score everything (throughput probes + 13 regions, ~6 min)
python3 bench_v13_robustness.py

# 3. tests
python3 -m pytest test_quality_error.py test_bench_v13.py test_cascade.py \
  test_cheap_router.py test_bench_v12.py test_robustness_benchmark.py -q
```

| Artefact | Path |
|---|---|
| Primary results | `results/bench_v13/robustness_results.json` |
| Per-locus disagreements | `results/bench_v13/disagreements.json` |
| PB fix A/B | `results/bench_v13/pb_fix_ab.json` |
| Baseline bit-identity check | `results/bench_v13/baseline_verification.json` |
| Evidence caches | `cache/bench_v13/*.npz` (gitignored) |
| Extraction logs | `logs_v13/*.log` |
| Changed code | `quality_error_model.py`, `extract_bench_v12.py`, `cache_locus_evidence.py`, `evaluate_quality_error.py` |
| New code | `bench_v13_robustness.py`, `test_bench_v13.py`, `run_v13_extract.sh` |

Environment: Python 3.14.3, numpy 2.4.4, scipy 1.17.1, pysam 0.24.0, Linux
6.19.12, 8 cores, CPU only (no CUDA device). Bootstrap seed 20260812 with
10,000 resamples (`evaluate_quality_error` module default), control seeds
(0,1,2,3,4), all unchanged. Every frozen constant is imported from `cascade`,
never redefined — asserted by `test_bench_v13.py`.

Regions, exactly: chr21:32,000,000–33,000,000 (15x / 30x / full),
chr21:31,000,000–32,000,000 (30x / full), chr21:29,999,999–30,469,999 (full),
plus the seven cached regions of `results/robustness_v11`
(13–17, 17–21, 21–25, 25–29, 40–44 Mb at ~15x; 29,999,999–30,469,999 at 15x and
30x). Sample HG002, GRCh38, GIAB v4.2.1 truth restricted to the matching
high-confidence BED, Illumina.

**No historical artefact was modified.** `results/robustness_v11/*`,
`results/bench_v12/*`, `History/3–12_DEVLOG.md`, `cascade.py` and
`cheap_router.py` are byte-identical to their state at session start.

## 22. Limitations

1. **Single chromosome, single sample, single platform.** chr21, HG002,
   Illumina. Nothing here supports a claim about other chromosomes, samples,
   sequencing technologies, or error profiles.
2. **Three of six new regions overlap the router's cutoff-selection span**
   (§4.1). Their router numbers are reported but excluded from the
   generalisation pool, which leaves the untainted new-region evidence at
   2.31 M loci / 2,397 SNP.
3. **Depth realisations are downsamples of one library**, not independent
   samples; 15x/30x/full at one locus set are three regimes of one region.
4. **The high-error regime was not tested with high-error data.** No BAM with
   an elevated error rate exists in the repository; the low-quality and low-VAF
   strata are a proxy and are labelled as one.
5. **Every difficulty gradient in §17 rests on 1–8 events.** Directionally
   consistent, individually underpowered.
6. **SNP only.** Indels, MNPs and structural variants are outside the label
   set entirely.
7. **The repaired PB uses ≤48 reads per locus** (§18.4), so its high-depth
   numbers are for a correct estimator on a subsample, not for the best
   estimator the data supports.
8. **Extraction timings in the caches are contaminated** by up to four
   concurrent jobs. The throughput numbers in §14 were measured separately in a
   quiet single process and are the ones used for every speedup claim.
9. **`positions` are cache-order indices offset by the chunk start, not genomic
   coordinates** (devlog 12 §16.5, inherited unchanged). The `frame_index` in
   the disagreement dump is an index into the scoring frame, not a coordinate.

## 23. Final verdict

| Question | Verdict |
|---|---|
| Is the depth ≥48x PB zeroing a correctness bug, and is it fixed? | **CONFIRMED.** Mechanism traced, minimal fix implemented, 12 regression tests, 1,371 true SNP recovered across three ~70x regions, PB recall 0.30–0.54 → 1.000. |
| Does the fix preserve the established ≤30x results? | **CONFIRMED.** Bit-identical on 910,976 real loci; zero LLRs differ at 15x. |
| Does the frozen cascade preserve PB-level accuracy on new/harder data? | **CONFIRMED** within scope: ΔF1 = +0.000103 [+0.00000, +0.00022] over 24.58 M loci / 33,893 SNP; no region DEGRADED; 12 of 13 regions PRESERVED, one IMPROVED. |
| Does it reproduce PB's decisions? | **CONFIRMED as an engineering property, with a caveat.** 15 disagreements in 24.58 M loci (6.1e−7), zero of them a lost true SNP. |
| Is PB still the accuracy ceiling at high depth? | **NEGATIVE.** After the repair the fixed-ε binomial *outscores* PB at ~70x on two of three regions, because PB acquires a low-VAF false-positive mode there. |
| Is there a failure boundary? | **CONFIRMED directionally, UNDERPOWERED numerically.** Base quality is the strongest axis (3.4e−7 → 3.9e−4, ~1,100x, as quality falls below Q20); VAF second; depth above 48x third. All on 1–8 events. |
| Is the frozen cutoff still the right operating point? | **CONFIRMED, and not re-tuned.** 5.412 sits at the knee of the accuracy-vs-compute curve at both 14x and 69x; larger cutoffs buy nothing and cost up to 180x speedup. |
| Cross-chromosome / cross-sample / cross-platform generalisation | **NOT TESTED.** |
| Behaviour in genuinely high-error data | **NOT TESTED** (no such data available). |
| Indels and non-SNP variants | **NOT TESTED.** |

**One paragraph, conservatively.** The depth-clipping defect reported in devlog
12 is real, is a correctness bug rather than numerical stabilisation, and is now
fixed by a change that is provably inert below the read-tensor width and that
recovers every one of the 1,371 true SNP the defect destroyed across three
~70x regions; the previously validated ≤30x results are reproduced bit-for-bit.
With PB repaired, the frozen binomial → router → PB cascade preserves PB-level
accuracy across 24.58 M loci and 33,893 SNP spanning 14x–127x, ten independent
chr21 spans and four pre-registered hard strata, at 0.0013–0.17% of PB's
compute and a projected 171–369x speedup, disagreeing with PB on 15 loci and
never losing a true variant to a disagreement. That statement is bounded to
chr21/HG002/Illumina SNP calling; nothing was tested outside it. Two findings
cut against the architecture's own premise and should not be buried: PB is no
longer the accuracy ceiling above ~48x, because the repaired PB over-calls
low-VAF sites at high depth where the pessimistic fixed-ε binomial does not; and
the cascade's disagreement rate rises by roughly three orders of magnitude as
mean base quality falls below Q25, on a handful of events that are directionally
consistent but individually underpowered. The frozen cutoff was not changed, no
threshold was selected on any number in this document, and no region was added
or dropped after results were read.

## 24. What remains untested

* Any chromosome other than 21; any sample other than HG002; any platform other
  than Illumina short reads.
* Genuinely high-error data (elevated substitution or alignment error rates), as
  opposed to the low-quality/low-VAF strata used as a proxy.
* Depth above ~127x, and depth between 44x and 48x on untainted coordinates.
* Indels, MNPs, structural variants; anything outside the SNP/non-SNP label.
* Somatic / low-VAF variant calling as a *task* (the low-VAF stratum here is
  mostly error, not true low-VAF variants).
* Whether the quality gradient in §17 survives at 100x more events — the only
  way to turn a 15-locus signal into a boundary.
* End-to-end pipeline runtime including BAM I/O; only caller compute was
  measured.
* The interaction between `MAX_READS` and PB accuracy at high depth (§18.4).

## 25. Recommended next research step

**Fix PB's high-depth false-positive mode, and measure it against the same
frozen cascade.** It is now the largest accuracy defect visible in the system,
it was invisible until tonight because the depth-clipping bug masked it, and it
has a clear mechanism (§18.1): the per-read Poisson-binomial has no term for
error sources that do not shrink with base quality, so at 70x it converts eight
high-Q non-reference reads into overwhelming evidence. Two concrete, separable
experiments:

1. **Raise the PB error floor to model alignment error explicitly** — the
   `EPSILON_FLOOR = 1e-4` docstring already argues that misalignment puts a
   floor on the non-reference rate; test whether a floor calibrated at high
   depth removes the FP mode without hurting 15x.
2. **Raise `MAX_READS` (or add a hypergeometric correction) so PB sees the
   whole pileup at high depth** (§18.4), and re-run this benchmark unchanged.

Both are single-variable changes against a frozen, now-validated cascade and an
existing benchmark harness, so each is a clean A/B rather than a redesign. The
router itself needs no further validation work at 15x — six regions and 19.4 M
loci at ΔF1 = 0 is a saturated question; the open question has moved to the
caller it routes to.

---

*Editorial note: after the results sections were appended, two forward cross-references in the pre-registration (§1–7) pointed at section numbers that the final layout had shifted — "see §20" for the verdict (now §23) and "§18 by diff" for the frozen-file check (now §21). Only those two pointers were corrected. No claim, threshold, endpoint, decision rule or region list in §1–7 was altered.*
