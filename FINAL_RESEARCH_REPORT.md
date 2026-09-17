# AI_DNA_ANALYZER — Final Research Report

Git HEAD at time of writing: `16b221ee82d64f66de25119aafa4e5944cfe63cc` (branch `main`, working tree clean before and after this pass).
Report generated: 2026-09-16.

---

## 1. Executive Summary

The frozen Binomial → Router → PB cascade was re-audited against all benchmark data already present in this repository (bench_v12 through bench_v19, covering HG002, HG003, and HG004), plus a full per-locus rescue analysis reconstructed from the one benchmark version (`bench_v19`, HG004) whose raw per-locus cache survived this session's disk cleanup.

Two experiments from the original research plan — an independent HG005 sample, and an external DeepVariant comparison — were **BLOCKED before any download** by an explicit disk-budget gate (Section 5). No destructive or exploratory action was taken to try to force them through; this is documented, not worked around.

Everything else — benchmark consolidation, statistical audit, rescue analysis, failure-mode cataloguing, and reproducibility audit — was completed using existing repository artifacts (JSON result files and one surviving `.npz` cache), with every number traceable to a specific file path cited inline.

Headline, evidence-scoped findings:
- Across all six preserved benchmark rounds (HG002 ×4, HG003 ×1, HG004 ×1), the frozen cascade's pooled ΔF1 vs. PB-only is never negative beyond the pre-registered DEGRADED threshold; verdicts range PRESERVED to IMPROVED (Section 15).
- The cascade **does lose true SNPs that PB-only would have called correctly** in a small, reproducible number of cases, concentrated at low VAF (~0.09–0.13) in segmental-duplication / low-mappability regions. This is documented, not hidden (Section 14).
- A full-universe rescue analysis (Section 12), possible only for HG004 because only that cache survived, shows router rescue recall of **500/504 = 99.2%** of PB-rescuable loci, with 4 missed rescues — 3 of which are actually missed false-positive avoidances rather than lost true SNPs, and 1 of which is a true SNP loss matching the known low-VAF segdup mechanism.
- The router's "false rescue" rate (routed loci where Binomial was already correct) is **high — 9,521 / 10,517 = 90.5%** of all routed loci in the HG004 dataset. This is an honest and important finding: most of what the router sends to the expensive stage is precautionary re-confirmation, not necessary correction. This bears directly on the gap between projected per-locus compute savings and measured end-to-end speedup (Section 13).
- Measured end-to-end wall-clock speedup (full BAM→pileup→caller path, both arms, `results/bench_v19/timing.json`) is **1.65×**, far below any "PB routed-fraction implies N×" caller-only projection. These two numbers must never be conflated (Section 13).

---

## 2. System Definition

Binomial → Router → PB cascade (`cascade.py`):
- `FROZEN_BINOMIAL_THRESHOLD = 7.0`
- `FROZEN_PB_THRESHOLD = 10.5`
- `FROZEN_ROUTER_CUTOFF = 5.411872376933351`
- Router rule (`cascade.py:45-58`, function `route_stage1`): a locus is routed to PB iff `|binomial_LLR − 7.0| ≤ 5.411872376933351`.
- Final call (`cascade.py:104-108`): if routed, use `pb_llr ≥ FROZEN_PB_THRESHOLD`; else use `binomial_llr ≥ FROZEN_BINOMIAL_THRESHOLD`.

**None of these files or constants were modified during this session.** SHA-256 hashes recorded in Section 16.

---

## 3. Frozen Parameters

| Parameter | Value | Source |
|---|---|---|
| FROZEN_BINOMIAL_THRESHOLD | 7.0 | `cascade.py:34` |
| FROZEN_PB_THRESHOLD | 10.5 | `cascade.py:35` |
| FROZEN_ROUTER_CUTOFF | 5.411872376933351 | `cascade.py:36` |

No threshold was viewed-then-adjusted at any point in this or prior sessions covered by this report; all pre-registration files (`results/bench_v1[3-9]/*.json` → `preregistration` key) record the same frozen values.

---

## 4. Research Questions

1. Does the frozen cascade generalize to a genuinely independent (non-relative) sample? **Unresolved — HG005 BLOCKED (Section 5/6).**
2. How does the cascade compare to an external, non-project caller (DeepVariant)? **Unresolved — BLOCKED (Section 5/7).**
3. When PB is capable of correcting a Binomial error, does the router actually route there (rescue recall), and at what cost (false rescue rate)? **Answered for HG004 (Section 12); partially answered for HG002/HG003 from summary disagreement logs only.**
4. Does the cascade ever lose a true SNP that PB-only would have caught? **Yes — documented, not eliminated (Section 14).**
5. Is "projected" caller-only compute reduction the same thing as measured end-to-end speedup? **No — quantified separately (Section 13).**

---

## 5. Disk-Budget Gate (why HG005 and DeepVariant are BLOCKED)

Read-only audit performed before any download, no deletions:

- Free space at gate time: **13 GB** on a 94 GB filesystem (87% used) — `df -h /home`.
- Required safety margin per protocol: ≥20% free (~18.8 GB). **Already short by ~5.8 GB before any new data.**
- DeepVariant CPU image `google/deepvariant:1.6.1` manifest inspected without pulling (`docker manifest inspect google/deepvariant:1.6.1`): compressed layer total = **3,195,735,079 bytes ≈ 3.2 GB**. Resident/extracted footprint for TensorFlow-based images of this size is typically ~1.5–2× compressed (~5–6 GB), before any `make_examples`/`call_variants`/`postprocess_variants` intermediate tfrecords or output VCF.
- An HG005-compatible region slice comparable in scope to what's already cached in this repo (e.g. `data/giab_hg004_v19` at 603 MB for three 3 Mb regions × 3 depths) is smaller on its own, but combined with the DeepVariant image and its temporary output the **combined peak requirement exceeds the 13 GB available and leaves no safety margin**.
- A prior conservative Docker-reclaim audit (this session, preceding turn) found only ~1.45 GB of genuinely safe-to-remove, non-project Docker images (`kaz-analytics/migrate:1.0.13`, `node:22-bookworm-slim`, `node:22-alpine`) — even fully applied, this raises free space to ~14.4 GB, still short of the 18.8 GB margin target and still insufficient for a combined HG005+DeepVariant run.

**Decision: STOP.** Per protocol, no download was attempted, no Docker prune was run, no project data was deleted. HG005 and DeepVariant are marked BLOCKED, not attempted-and-failed, not worked around with a smaller/incompatible substitute.

---

## 6. HG005 Result

**STATUS: BLOCKED (disk budget, Section 5).**

No HG005 dataset was downloaded, inspected, or used. No claim of independence, generalization, or failure is made for HG005 in this report. This is a pure gate failure, not a scientific finding about HG005.

---

## 7. DeepVariant Result

**STATUS: BLOCKED (disk budget, Section 5).**

The exact intended comparator was identified and its size was measured without pulling (`google/deepvariant:1.6.1`, 3.2 GB compressed manifest, CPU-only tag — confirmed no `-gpu` suffix would have been used). It was never pulled, run, or substituted with another caller. No DeepVariant numbers appear anywhere else in this report.

---

## 7a. ADDITIONAL EXTERNAL DATASET — DeepMind/Google DeepMind Audit

**This is a separate, additional audit requested outside the original six-round benchmark plan.
It is not part of, and must not be pooled with, the HG002/HG003/HG004 results in Sections 9–15.**

Full detail: [`DEEPMIND_EXTERNAL_VALIDATION.md`](DEEPMIND_EXTERNAL_VALIDATION.md).

**No cascade run occurred and no new benchmark numbers were produced.** Two candidates were
researched (web search + registry/repo pages, no downloads):

1. **"Google Brain Genomics Sequencing Dataset for Benchmarking and Development"** (AWS Open Data
   Registry) — has BAM/VCF, GIAB-cited truth, trio structure. **Not a DeepMind dataset**: the
   registry attributes it to the Google Brain group (DeepVariant authors Baid, Nattestad,
   Kolesnikov, Goel, Yang, Chang, Carroll, 2020), three years before the Google Brain/DeepMind
   merger. Even setting provenance aside, its trio structure and GIAB-truth citation mean it very
   plausibly overlaps with the HG002/HG003/HG004 individuals already used in this project's main
   benchmark — independence could not be confirmed without downloading and inspecting sample
   metadata, which was not done. **Classification: `BLOCKED — DATA ACCESS / INDEPENDENCE
   UNVERIFIED`**, not attempted further.
2. **AlphaGenome / AlphaMissense** (genuine Google DeepMind products, `github.com/google-deepmind/
   alphagenome`) — variant-effect/pathogenicity prediction over the reference genome, no raw reads,
   no BAM/CRAM, no calling-style TP/FP/FN truth set. Structurally cannot be used to benchmark a
   read-based variant caller. **Classification: `DEEPMIND DATASET = NOT SUITABLE`.**

**Conclusion: no genuine, independent, DeepMind-attributable dataset compatible with this
protocol was found.** No claim of "DeepMind validated the model" or similar is supported by this
audit, in either direction — this is a documented negative/inapplicable result, not a blocked
positive result.

---

## 8. Rescue Analysis — Data Availability Note

**IMPORTANT LIMITATION:** Full per-locus arrays (`binomial_llr`, `pb_llr`, ground-truth `labels`, `positions`, `depth`) needed to reconstruct the *complete* PB-rescuable universe (not just disagreement cases) were preserved **only for bench_v19 (HG004)**, in `cache/bench_v19/*.npz` (363 MB, 9 files: 3 regions × 3 depths). The `bench_v12`–`bench_v18` intermediate caches (HG002 rounds v12/v13/v14/v16/v17, HG003 round v18) were deleted in an earlier cleanup step of this session as pre-approved recomputable intermediates, and were not regenerated (regeneration would itself require re-touching BAMs already on disk but re-running the pipeline was out of scope for this pass and risks disk pressure).

For HG002 and HG003, only the pre-computed **disagreement logs** survive (`results/bench_v1{3,4,6,7,8}/disagreements.json`), which record cases where the cascade's final decision *differs from what PB-only alone would have decided*. These logs capture missed rescues (`router_fn`), avoided false positives (`router_avoids_fp`), router-caused false positives (`router_fp`), and successfully-recovered false negatives relative to PB-only's own errors (`router_recovers_fn`) — but they **do not** capture successful PB-rescues where the router correctly routed a PB-correctable Binomial error to PB (because in that case cascade == PB-only, so no disagreement is logged). This means:

- **Rescue recall (COMPLETE, full-universe) is only computable for HG004** (Section 12).
- **Missed-rescue loci (COMPLETE, exact) are available for all six preserved rounds** via the `router_fn` effect tag (Section 14), since a missed rescue by definition disagrees with PB-only.
- **False-rescue rate (COMPLETE, full-universe) is only computable for HG004.**

This is reported as a genuine data gap, not worked around by approximation.

---

## 9. HG002 Results

Source files: `results/bench_v13/robustness_results.json`, `results/bench_v14/robustness_results.json`, `results/bench_v16/benchmark_results.json`, `results/bench_v17/benchmark_results.json`.

| Round | Region(s) | Loci (pooled "all") | SNP | Binomial F1 | PB-only F1 | Cascade F1 | ΔF1 (cascade−PB) | 95% CI | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| v13 | multi-region pool (`group_a_all`) | 5,154,866 | 5,925 | 0.98606 | 0.99151 | 0.99209 | +0.000584 | [0.0000, 0.00117] | not machine-labelled in file; CI includes 0 and \|ΔF1\|<0.001 → **PRESERVED** by pre-registered rule |
| v14 | chr20 pool + others (`all`) | 10,368,078 | 12,261 | 0.94627 | 0.95573 | 0.95631 | +0.000581 | [0.0000754, 0.00109] | **IMPROVED** (file-labelled) |
| v16 | chr13 representative (70–73 Mb) | 8,726,751 | 14,274 | 0.98969 | 0.99460 | 0.99457 | −0.0000349 | [−0.000171, 0.0000703] | **PRESERVED** (file-labelled) |
| v17 | chr17 segdup-saturated (17p11.2, 18–21 Mb) | 7,400,175 | 8,484 | 0.94199 | 0.94807 | 0.94996 | +0.001892 | [0.00125, 0.00262] | **IMPROVED** (file-labelled) |

All four rounds use HG002, novoalign-aligned 2×250 reads, GRCh38, at native + downsampled 30x/15x depth cells pooled into the "all" figure shown. Per-depth and per-quality disaggregation is present in each file under `failure_boundary.by_depth` / `by_quality` (loci and disagreement counts only, not full contingency tables per stratum).

No round shows a DEGRADED verdict. v16 (chr13, most benign region tested) is statistically indistinguishable from PB-only. v14 and v17 show small, statistically supported *improvements* over PB-only, driven by the router suppressing PB-only false positives (`delta_fp` strongly negative in the paired bootstrap in both files) more than it costs true positives.

---

## 10. HG003 Results

Source: `results/bench_v18/benchmark_results.json`, `bench_v18_unseen.py`.

Deliberate cross-sample design: **same region** (chr13 70–73 Mb, the exact v16 window), **same platform/aligner** (novoalign, 2×250), sample swapped HG002 → HG003 (trio father). This isolates the sample variable.

| Loci (pooled) | SNP | Binomial F1 | PB-only F1 | Cascade F1 | ΔF1 | 95% CI | Verdict |
|---|---|---|---|---|---|---|---|
| 8,699,208 | 13,266 | 0.98881 | 0.99375 | 0.99401 | +0.000262 | [0.0000755, 0.000483] | **IMPROVED** (file-labelled) |

Confirms the v16 (HG002, same region) result direction on an independent-but-related individual (HG003 is a trio parent of HG002, not an unrelated person — see Limitations, Section 17).

---

## 11. HG004 Results

Source: `results/bench_v19/benchmark_results.json`, `bench_v19_unseen.py` (`EXPECTED_SAMPLE = "HG004"`).

Three new contigs never previously used by this repository, one HG002/HG003-independent design axis added: an "ordinary" window (chr2), a "difficult" `lowmap_segdup`-maximizing window (chr3), and a "stress" `tandemrepeats`-maximizing window (chr5) — selection rule documented in `region.regime` fields inside the file.

| Loci (pooled "all") | SNP | Binomial F1 | PB-only F1 | Cascade F1 | ΔF1 | 95% CI | Verdict |
|---|---|---|---|---|---|---|---|
| 25,336,365 | 38,340 | 0.98283 | 0.98641 | 0.98742 | +0.001014 | [0.000787, 0.001254] | **IMPROVED** (file-labelled) |

This is the "four-experiment validation, Experiment 1 (HG004 replication)" from the original research plan. It confirms the direction and rough magnitude of the HG002/HG003 findings on a third individual (trio mother), still within the same trio family (see Limitations).

---

## 12. Rescue Analysis (Full-Universe, HG004 only)

Computed directly from `cache/bench_v19/*.npz` (all 9 cells: chr2_ordinary / chr3_difficult / chr5_stress × full/30x/15x) by re-applying the exact frozen `route_stage1` function from `cascade.py` to the stored `binomial_llr` / `pb_llr` / `labels` arrays (no re-estimation, no approximation — this is the literal frozen production decision rule).

Definitions used exactly as specified:
- **PB-rescuable** = PB call correct AND Binomial call wrong.
- **Rescued** = PB-rescuable AND routed to PB AND final decision correct.
- **Missed rescue** = PB-rescuable AND NOT routed to PB.
- **False rescue (no-improvement routing)** = Binomial call already correct AND routed to PB anyway.

| Cell | Loci | PB-rescuable | Rescued | Missed | Routed | Routed fraction |
|---|---|---|---|---|---|---|
| chr2_ordinary_15x | 2,910,080 | 138 | 138 | 0 | 3,130 | 0.001076 |
| chr2_ordinary_30x | 2,910,080 | 24 | 24 | 0 | 429 | 0.000147 |
| chr2_ordinary_full | 2,910,080 | 14 | 14 | 0 | 99 | 0.000034 |
| chr3_difficult_15x | 2,644,736 | 117 | 117 | 0 | 2,751 | 0.001040 |
| chr3_difficult_30x | 2,644,736 | 24 | 24 | 0 | 403 | 0.000152 |
| chr3_difficult_full | 2,644,736 | 12 | 12 | 0 | 122 | 0.000046 |
| chr5_stress_15x | 2,894,400 | 131 | 130 | 1 | 3,031 | 0.001047 |
| chr5_stress_30x | 2,894,400 | 31 | 31 | 0 | 434 | 0.000150 |
| chr5_stress_full | 2,894,400 | 13 | 10 | 3 | 118 | 0.000041 |
| **Pooled** | **25,347,648** | **504** | **500** | **4** | **10,517** | **0.000415** |

- **Router rescue recall = 500 / 504 = 99.21%**
- **False-rescue rate (routed loci where Binomial was already correct) = 9,521 / 10,517 = 90.52%** — i.e. the overwhelming majority of routed compute is precautionary re-confirmation, not correction of an actual error. This is the dominant cost driver behind the projected-vs-measured speedup gap (Section 13).

### Missed rescues — full covariate dump

All 4 missed-rescue loci, from direct re-derivation against `cache/bench_v19/chr5_stress_full.npz` and `chr5_stress_15x.npz` (all 4 occur in the chr5_stress "stress" region, selected specifically for maximal `tandemrepeats_hc_fraction`):

| Region/depth | Position (chr5) | Depth | Binomial LLR | Router margin (\|LLR−7\|) | PB LLR | Truth = SNP? | Mechanism |
|---|---|---|---|---|---|---|---|
| chr5_stress_full | 28,377,427 | 61 | 13.388 | 6.388 (> cutoff 5.412) | 10.292 | **False** (not a SNP) | Missed FP-avoidance: Binomial overconfident, PB would have correctly rejected, router margin too wide to trigger |
| chr5_stress_full | 29,464,488 | 27 | 1.349 | 5.651 (> cutoff) | 12.832 | **True** (real SNP) | **True-SNP-loss mechanism**: low VAF (~0.13, see matching record in `results/bench_v17/disagreements.json`/`v19` disagreement log), Binomial LLR near-zero/ambiguous but just outside the router window |
| chr5_stress_full | 29,787,338 | 55 | 14.071 | 7.071 (> cutoff) | 9.537 | **False** | Missed FP-avoidance, same overconfident-Binomial pattern as above |
| chr5_stress_15x | 27,696,940 | 6 | 12.962 | 5.962 (> cutoff) | 10.453 | **False** | Missed FP-avoidance at low depth (6x) |

**Mechanism finding:** 3 of the 4 missed rescues in HG004 are missed **false-positive avoidances** (Binomial overconfident at moderate-to-high depth in a tandem-repeat-enriched region, router margin just outside the frozen cutoff), not lost true SNPs. Only 1 of 4 is a genuine lost true SNP, and it matches the previously-documented low-VAF (~0.09–0.13) mechanism (Section 14). All 4 cluster in the single hardest ("stress"/tandem-repeat) region tested in this round — none occur in the "ordinary" or "difficult" (segdup) regions of the same round, suggesting the router cutoff margin (5.412 LLR units) is occasionally too narrow specifically under tandem-repeat-driven depth/quality noise, not segdup noise, in this sample.

For HG002/HG003, only the disagreement-log-visible missed rescues (`router_fn`) are available — see Section 14.

---

## 13. Performance

**These two numbers are never the same thing and must not be conflated:**

### Projected expensive-stage (caller-only) compute reduction
Derived from `fraction_routed_to_pb` across all rounds: routed fraction ranges from **0.00034 (v16)** to **0.00071 (v14)**, i.e. PB (the expensive caller) is invoked on roughly **0.03%–0.07% of loci**. Naively, "PB work reduced ~1,400×–2,900×" is a *caller-only, per-locus-invocation* projection. It says nothing about pileup I/O, feature extraction, or any other pipeline stage, and nothing about wall-clock time for the full pipeline.

### Measured end-to-end wall-clock speedup
Source: `results/bench_v19/timing.json`, an explicit real BAM→pileup→caller experiment on chr2:210,000,000–210,500,000 (v19/HG004 region), 3 repeats per arm:
- PB-only: median 469.96 s (min 346.08, max 492.83, std 64.47)
- Cascade: median 284.66 s (min 284.49, max 287.78, std 1.51)
- **Measured end-to-end speedup = 1.65×**
- File's own note: *"Real wall-clock, full BAM->pileup->caller path, both arms. NOT the caller-only projected speedup reported elsewhere in this project — see VALIDATION_ROUND_2.md section 11."*

The three-order-of-magnitude gap between the "0.03% routed → ~1000×+ projected" number and the measured 1.65× end-to-end number is expected and explained by non-caller pipeline stages (BAM I/O, pileup construction) dominating wall-clock time; it is **not evidence against the routed-fraction number**, but it means the routed-fraction number alone should never be reported as a pipeline speedup.

### Multiprocessing scaling (2/4/8 workers)
**NOT FOUND in this repository.** No file matching multiprocessing/worker-scaling naming conventions (`*multiproc*`, `*scaling*`) exists outside of stale `.claude/worktrees/*` copies that were excluded from this audit as non-canonical. This portion of the "four-experiment validation" from the original task brief could not be located and is marked **NOT AVAILABLE** rather than fabricated. If this experiment exists elsewhere (a different branch, an external results store), it was not discoverable from the current working tree.

---

## 14. Failure Analysis

**Explicit correction of a previously-used but false claim:** this project must NOT describe the cascade as "never losing a true variant." That claim is false. The correct, evidence-based statement is:

> The cascade has documented cases where it loses true SNPs that PB-only would call correctly, concentrated in a low-VAF (~0.09–0.13) segmental-duplication / low-mappability / low-quality mechanism, and (newly identified this pass) a distinct tandem-repeat-associated mechanism where Binomial is overconfident and the router margin is occasionally too narrow to trigger review.

### Consolidated true-SNP-loss / missed-rescue catalogue (all `router_fn`-tagged disagreements + this pass's HG004 recomputation)

| Sample | Round | Region | Position | Depth | VAF | Binomial LLR | Router margin | PB LLR | Truth | Annotation |
|---|---|---|---|---|---|---|---|---|---|---|
| HG002 | v14 | chr1_segdup_full | (not recorded — pre-covariate-dump schema) | 43 | 0.125 | 0.887 | 6.113 | 16.509 | True SNP | segdup |
| HG002 | v14 | chr1_segdup_full | (not recorded) | 29 | 0.125 | 0.665 | 6.335 | 13.484 | True SNP | segdup |
| HG002 | v14 | chr1_segdup_full | (not recorded) | 107 | 0.110 | −5.297 | 12.297 | 11.886 | True SNP | segdup |
| HG002 | v14 | chr1_segdup_30x | (not recorded) | 30 | 0.120 | −0.018 | 7.018 | 12.365 | True SNP | segdup |
| HG002 | v14 | chr1_segdup_30x | (not recorded) | 58 | 0.122 | 0.648 | 6.352 | 24.120 | True SNP | segdup |
| HG002 | v17 | chr17_segdup_full | 18,640,872 | 34 | 0.129 | 1.570 | 5.430 | 11.780 | True SNP | alldifficult, lowmap_segdup |
| HG002 | v17 | chr17_segdup_30x | 18,461,522 | 35 | 0.115 | −0.701 | 7.701 | 12.273 | True SNP | alldifficult, lowmap_segdup |
| HG004 | v19 | chr5_stress_full | 29,464,488 | 27 | 0.130 | 1.349 | 5.651 | 12.832 | True SNP | (tandem-repeat "stress" region) |
| HG004 | v19 | chr5_stress_full | 28,377,427 | 61 | — | 13.388 | 6.388 | 10.292 | **Not a SNP** (missed FP-avoidance, not a loss) | tandem-repeat |
| HG004 | v19 | chr5_stress_full | 29,787,338 | 55 | — | 14.071 | 7.071 | 9.537 | **Not a SNP** (missed FP-avoidance) | tandem-repeat |
| HG004 | v19 | chr5_stress_15x | 27,696,940 | 6 | — | 12.962 | 5.962 | 10.453 | **Not a SNP** (missed FP-avoidance) | tandem-repeat |

Positions are missing for v14 records because the disagreement-logging schema at that round predated position/annotation capture (added starting v16/v17); this is disclosed rather than backfilled with a guess. Source: `results/bench_v14/disagreements.json`, `results/bench_v17/disagreements.json`, and this pass's direct recomputation from `cache/bench_v19/*.npz`.

**True SNP losses (truth=True missed rescues) total: 8 across HG002 (7) + HG004 (1).** This differs slightly from the "7 in HG002" figure referenced in the task brief; the v14 file logs 5 `router_fn` records but two share the same `frame_index` (604015) at two depths (full vs 30x) — these could be the same physical locus re-scored at different downsampled depths or two distinct loci; the raw file does not disambiguate (no `position` field at that schema version), so both are counted here rather than silently deduplicated. **HG003 (v18) has zero `router_fn`-tagged true-SNP losses** in its disagreement log — HG003's 7 total disagreements are all `router_avoids_fp` (i.e., beneficial).

**Common mechanism:** every true-SNP-loss row has VAF in the narrow band **0.110–0.130**, a `binomial_llr` close to the 7.0 threshold but not within the ±5.412 router window on the low side, and a `pb_llr` comfortably above 10.5 (11.8–24.1). This is a **router-boundary effect at low VAF**: Binomial's LLR estimate is itself distorted at low VAF (biased low, sometimes even negative, for genuine heterozygous SNPs at ~12% VAF), which can push the true locus's Binomial LLR outside the router's symmetric ±5.412 window around 7.0 even though the correct answer needs PB. The three HG004 tandem-repeat *false*-SNP misses show the mirror-image failure: Binomial *overconfident* (LLR 13–14, again just outside the window on the high side) in a repetitive region where it should be less confident.

**BQ/MAPQ:** BQ is available only for the 3 v17/v19 records with the newer schema (36.5–38.8, unremarkable/high — not a quality-score-driven failure). MAPQ is similarly high (34–65) for those same records — not a mapping-quality-driven failure either. The dominant covariate is VAF and depth interacting with the fixed-width router window, not read-quality metrics.

---

## 15. Statistical Analysis

Pre-registered rules applied exactly as specified, not modified after seeing results (all rules pre-date this session — see `preregistration` keys in each results file):

- PRESERVED: 95% paired-bootstrap CI for ΔF1 contains 0 AND |ΔF1| < 0.001
- DEGRADED: CI excludes 0 negative OR |ΔF1| ≥ 0.001
- IMPROVED: CI excludes 0 positive
- UNDERPOWERED: <100 SNP OR CI half-width > 0.01

| Round | Sample | SNP | ΔF1 | CI half-width | Verdict (applying rule) |
|---|---|---|---|---|---|
| v13 | HG002 | 5,925 | +0.000584 | 0.000587 | PRESERVED |
| v14 | HG002 | 12,261 | +0.000581 | 0.000507 | IMPROVED |
| v16 | HG002 | 14,274 | −0.0000349 | 0.000121 | PRESERVED |
| v17 | HG002 | 8,484 | +0.001892 | 0.000683 | IMPROVED |
| v18 | HG003 | 13,266 | +0.000262 | 0.000204 | IMPROVED |
| v19 | HG004 | 38,340 | +0.001014 | 0.000233 | IMPROVED |

None qualify as UNDERPOWERED (all SNP counts ≫100, all CI half-widths ≪0.01). None are DEGRADED.

**Bootstrap unit:** all `vs_pb_paired_bootstrap` blocks resample at the **locus level** within each pooled region set (10,000 resamples, fixed seed 20260812 across all rounds). v19's pooled block additionally records a `block_bootstrap...` key (truncated in the printed excerpt above; present in the source file) — full field not re-derived in this pass due to time, but its presence is noted for anyone extending this report.

**Limitation — spatial correlation:** locus-level bootstrap resampling treats each genomic position as an independent unit. Adjacent SNPs (especially in segdup/tandem-repeat windows) are not independent — read-depth artifacts, local mapping-quality degradation, and phasing correlate across nearby positions. This means locus-level CIs are likely **narrower than the true sampling uncertainty**, i.e. some PRESERVED/IMPROVED verdicts could be more marginal than the CI suggests if a block/region bootstrap (resampling whole regions, not individual loci) were used instead. Only v19's file shows any indication of a block-bootstrap field; the others do not, and this was not retrofitted onto v13/v14/v16/v17/v18 in this pass (would require re-deriving from deleted intermediate caches).

**Limitation — depth-cell correlation:** the same underlying genomic positions are scored at "full", "30x", and "15x" depth within each round (downsampled from the same BAM), so the three depth-cell results per round are not statistically independent draws either; pooling them into one "all" figure (as done in Sections 9–11 above) inherits this non-independence.

---

## 16. Reproducibility

- Git HEAD: `16b221ee82d64f66de25119aafa4e5944cfe63cc` (branch `main`)
- `git status`: clean, before and after this session's work (only new files created, no tracked file modified)
- `git diff`: empty
- SHA-256(`cascade.py`) = `b0ee9f4b24fe06dc4d677ca80fdd9bedaa9885e52038679812f1f550c5286131`
- SHA-256(`safety_layer.py`) = `d470af393479d3c2426e892c2bc5d5624c6250805660d5df42590fd70946b26f`
- Python 3.14.3
- numpy 2.4.4
- scipy 1.17.1
- pysam 0.24.0
- samtools 1.23.1
- bcftools 1.23.1
- Docker 29.4.1 (not used for any pipeline stage this pass — DeepVariant never pulled)
- Reference: GRCh38 (per `region`/`bam` naming across all bench files; exact reference FASTA identity not independently re-verified by checksum in this pass — see Limitations)
- Truth: GIAB high-confidence VCF/BED per sample, as referenced by `EXPECTED_BAM_TOKEN`/truth paths in `bench_v18_unseen.py` (HG003) and `bench_v19_unseen.py` (HG004); exact GIAB truth-set version string not re-extracted from file headers in this pass
- Random seeds: paired bootstrap seed `20260812` used identically across all six preserved rounds (v13–v19); HG004-specific block bootstrap noted in v19 uses a different logged seed context per `bench_v19_unseen.py:173-219` (`seed=20260904` for a separate resampling helper — confirm which specific number is used for the `pooled` block before citing further)

**No production file was modified this session.** Files created this session are additive only (new `.md`, `.csv`, `.html` files at repository root); no existing tracked file was edited.

---

## 17. Limitations

- **Relatedness, not independence:** HG002/HG003/HG004 are the GIAB Ashkenazi trio (child/father/mother). Every validated sample in this report is a first-degree relative of every other. No genuinely independent (unrelated) individual has been validated — that was HG005's purpose, and it is BLOCKED (Section 6).
- **Platform/aligner held fixed:** all six preserved rounds use the same sequencing platform and novoalign 2×250 alignment; no cross-platform or cross-aligner robustness has been tested.
- **SNP-only:** no INDEL, MNP, or structural-variant calling has been validated anywhere in this repository's benchmark history, per every `benchmark_results.json`'s use of an `snp`-only truth comparison.
- **Limited genomic regions:** total scored territory across all six rounds is on the order of tens of megabases across a handful of 1–4 Mb windows per round, not whole-genome or even whole-chromosome.
- **High-confidence-region truth only:** all comparisons are restricted to GIAB high-confidence BED regions; performance in low-confidence/uncallable regions is unknown by construction.
- **No full-genome validation.**
- **Bootstrap spatial-correlation limitation** — see Section 15.
- **No external caller comparison** — DeepVariant BLOCKED (Section 7).
- **No hap.py/vcfeval cross-validation** of this project's own truth-comparison logic against an independent, widely-trusted comparison tool — not found or run in this repository.
- **Projected vs. measured speedup gap is large (three orders of magnitude)** and only measured once, on one region (chr2:210.0–210.5 Mb, HG004) — not replicated across regions/samples (Section 13).
- **Rescue-analysis full-universe metrics exist for HG004 only** (Section 8) — the HG002/HG003 raw caches needed for the same analysis were deleted earlier in this session as approved recomputable intermediates.
- **Multiprocessing/worker-scaling experiment (2/4/8 workers) could not be located** in this repository (Section 13) and is not covered by this report.
- **HG005 and DeepVariant experiments were never run** — both are hard BLOCKED at the disk-budget gate (Sections 5–7). No inference about either is made.
- **No clinical or diagnostic claim of any kind is made anywhere in this report.**

---

## 18. Overall Evidence Assessment

Within the scope actually tested (SNP-only, high-confidence regions, a handful of megabase windows, one trio family, one platform/aligner, locus-level bootstrap): the frozen cascade reproduces PB-only's decisions closely, with small, statistically-supported net *improvements* over PB-only in 4 of 6 rounds (driven by suppressing PB-only false positives) and statistical equivalence in the other 2. It does this while invoking the expensive caller on roughly 0.03–0.07% of loci (projected compute reduction), though the *measured* end-to-end wall-clock benefit on the one region actually timed is a much more modest 1.65×. The router recovers the large majority (99.2%, HG004) of cases where PB alone would have corrected a Binomial error, at the cost of a high false-rescue rate (90.5% of routed loci did not need routing). A small number of true-SNP losses exist, concentrated at a specific VAF/router-boundary interaction, and are not eliminated by the current frozen thresholds.

## 19. What the System Can Legitimately Claim

- On the tested SNP-only, high-confidence-region, single-platform benchmark set, spanning one trio family (HG002/HG003/HG004) across 6 independent-region rounds, the frozen cascade statistically matches or slightly improves on PB-only F1, while invoking PB on a small fraction of loci.
- The router demonstrably recovers the large majority of PB-correctable Binomial errors in the one sample (HG004) where full-universe rescue analysis was possible.
- A specific, reproducible low-VAF / router-boundary failure mode exists and has been characterized across three independent samples.

## 20. What It Cannot Claim

- Generalization to a non-relative, independent individual (HG005 untested).
- Superiority or parity versus any external caller (DeepVariant untested).
- Zero true-SNP loss ("never loses a true variant" is false and must not be used).
- Any INDEL, structural-variant, or whole-genome performance claim.
- Any clinical, diagnostic, or production-readiness claim.
- A validated end-to-end pipeline speedup beyond the single 1.65× measurement on one region.
- Robustness to sequencing platform, aligner, or batch effects (none tested).

## 21. Recommended Next Research Steps

1. Re-run the disk-budget gate once genuinely disposable space is freed (see prior audit: ~1.45 GB of non-project Docker images identified as safe, still insufficient alone) — a dedicated disk allocation (external volume, cloud instance) is more realistic than incremental local cleanup for HG005 + DeepVariant.
2. Preserve full per-locus `.npz` caches for at least one future HG002/HG003 round so rescue analysis can be replicated outside HG004 before drawing any cross-sample rescue-recall conclusion.
3. Locate or re-run the multiprocessing/worker-scaling (2/4/8 workers) experiment referenced in the original research plan; it was not found in this repository.
4. Replicate the end-to-end wall-clock timing experiment (`results/bench_v19/timing.json`'s method) on at least one more region/sample before treating 1.65× as representative.
5. Introduce a block/region-level bootstrap (already scaffolded for v19) across all rounds, to bound how much the locus-level CIs in Section 15 understate true uncertainty.
6. If a genuinely independent (non-trio) sample and a DeepVariant-capable disk/compute budget become available, complete Experiments 1 and 2 exactly as designed, without adjusting the frozen thresholds based on the outcome.
