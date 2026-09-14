# Independent Validation Round 2 — Frozen Binomial → Router → PB Cascade

**Auditor:** independent Claude Sonnet 5 session (adversarial audit, not the development team)
**Date:** 2026-09-04
**Repo HEAD at start:** `1fd9dc99cbff43ae003ed8b71b80397406bd2691`

**Everything below the line "## PRE-REGISTRATION CLOSED" was written and
committed to disk before any caller, router, or PB score existed for HG004 on
chr2/chr3/chr5, and before any end-to-end timing was measured.** The mechanical
region-selection script (`select_v19_regions.py`) was run once to produce
concrete coordinates from public annotation only (no caller output, no
benchmark metric) — its output is included below because the project's own
convention (devlogs 14/16/17) treats annotation-only selection as part of
pre-registration, not as a result. Everything after that marker is appended
only after the corresponding measurement, never edited into the plan above it.

---

## 0. Objective of this round

Not to improve the system. Not to reproduce a good number. The explicit
objective is to **try to break** the Round-1 STRONG POSITIVE finding and
determine actual readiness for scientific/commercial use. Any decision made
after viewing a performance number is marked post-hoc and excluded from the
confirmatory conclusion.

## 1. Claim under test (narrow, as instructed)

> The frozen binomial → router → PB cascade can reproduce PB-only SNP-calling
> accuracy at substantially reduced PB compute, on independent data.

Explicitly **not** under test: "the cascade is a better variant caller" or
"the cascade is clinically accurate." Those are out of scope for this round.

## 2. Frozen implementation — hashes recorded before touching HG004 data

| item | value |
|---|---|
| binomial threshold | 7.0 |
| binomial error rate ε | 0.01 |
| PB threshold | 10.5 |
| frozen router cutoff | 5.411872376933351 |
| reference | Ensembl release-110 GRCh38, per-chromosome FASTA |
| truth | GIAB NISTv4.2.1, sample's own benchmark VCF, SNP only |
| high-confidence BED | sample's own v4.2.1 `_noinconsistent` BED |
| stratifications | GIAB v3.1 (`lowmap_segdup`, `alldifficult`, `tandemrepeats`) |

SHA-256 of every frozen file this round imports (recorded before any HG004
chr2/chr3/chr5 data was touched, re-verified identical to every prior round):

```
cascade.py                b0ee9f4b24fe06dc4d677ca80fdd9bedaa9885e52038679812f1f550c5286131
cheap_router.py           5e13cca000f4565ab7023352b8b184a3c370cf6f75bf31b6700fd92dbb23f47a
binomial_baseline.py      f0408ecae0929aefd59faa74fc96e79228606bb569c0ae3c5efabdda4d3d6ce9
quality_error_model.py    a69cb49591cbac13efad07d73eda3004a89e6c0659287676a4c8901cc3db053e
robustness_benchmark.py   6fd5c4df39806e0594918bf4fe38333c6ce8e0643d8871dbdd0621e45374a04a
extract_bench_v12.py      8ec91cbd646cfdebd0dc653b5c1b793a0fe7f311906b30e0927c2b76c217cffc
pileup_counts.py          70b141854cb45fa33d414338d9f88902e4bb7001354b80f48fa54da8c1e7bf91
read_level_pileup.py      2a63c448eb271538f948dd4ff528cf53e3851e3c1d05d0235205021189ce6082
bench_v12_stage1.py       745d7d694653ebd9aa9d52b9cc56a032234f793b2068c84da7cff86dd119c478
bench_v13_robustness.py   839f105a45c3e924f8c558564a3407810cebd5a72fd9e0748b8090e2209b0008
bench_v14_crosschrom.py   5c36f1364cc5f8c0382f90c89c09003ef3f52250802d57570acae15079beb383
safety_layer.py           d470af393479d3c2426e892c2bc5d5624c6250805660d5df42590fd70946b26f
```

**None of these files will be modified during this round.** Router cutoff,
binomial threshold, PB threshold, feature extraction and the (rejected)
safety layer stay exactly as frozen. Any improvement idea noticed during this
round is recorded in §16 (What would falsify the claim) / the final report as
a hypothesis for future work, never applied retroactively to this round's
result.

## 3. Sample selection rule (fixed before viewing any result)

**Priority order requested:** HG005 > HG006 > HG007, non-Ashkenazim-trio,
never previously touched by this repository.

**What was checked, in order, before selecting anything:**

1. GIAB's public `ftp-trace.ncbi.nlm.nih.gov` mirror hosts the exact
   `NIST_Illumina_2x250bps/novoalign_bams` product — the platform+aligner this
   cascade's thresholds were derived on — **only** for the Ashkenazim trio
   (HG002/HG003/HG004). The `ChineseTrio` (HG005/HG006/HG007) and `NA12878`
   (HG001) data directories were enumerated directly; every per-technology
   subdirectory advertised for HG005 (`HG005_NA24631_son_HiSeq_300x`,
   `NIST_BGIseq_2x150bp_100x`, `BGISEQ500`, `MGISEQ`, `HudsonAlpha_PacBio_CCS`,
   etc.) returns **HTTP 404** on direct request — the index page listing them
   is stale relative to what the server actually serves. This was verified
   with `curl -sI`/`curl -o /dev/null -w '%{http_code}'`, not inferred from the
   directory listing alone.
2. HG005/HG006/HG007 **do** have a current NISTv4.2.1 GRCh38 truth-set
   release (verified: `release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh38/`
   exists and is populated) — the truth-set gap is not the blocker, the
   platform-matched-BAM gap is.
3. HG004 (Ashkenazim mother) has the exact matching
   `NIST_Illumina_2x250bps/novoalign_bams` product and NISTv4.2.1 GRCh38 truth
   (verified `HTTP 200` on both URLs directly), identical in structure to the
   HG002/HG003 data already used by devlogs 11-18.

**Decision, taken before any HG004 locus was scored:**

> HG005/HG006/HG007 are marked **NOT TESTED — reason: no platform/aligner-matched
> alignment could be located on the GIAB public mirror within this session's
> time budget; the advertised per-technology subdirectories 404.** This is
> reported as a genuine infrastructure limitation, not simulated or
> substituted with different data under the same label.
>
> **HG004 (Ashkenazim mother) is used as the Round-2 sample instead.** This is
> a **weaker** test than the user's request for a non-trio sample — HG004
> shares the exact sequencing batch, lab, platform and aligner with HG002 and
> HG003, and is genetically related to both (parent/offspring and spousal
> relationships within one reference family). It isolates the sample variable
> cleanly (same platform/aligner/pipeline as devlog 18's HG003 test) but
> **cannot** establish independence from the Ashkenazim-trio-specific data
> product this cascade has always been validated against. This limitation is
> carried into the executive verdict, not buried in a footnote.

**Leakage check for HG004 specifically:** this repository has never
downloaded, read, cached, or referenced any HG004 file before this round
(verified: no `HG004` string appears anywhere in `data/`, `cache/`, `results/`,
or `History/` before this document was written — `grep -ril HG004` over the
repo returns nothing). The three frozen constants (§2) were fit once on HG002
chr21 in devlog 13 and have not been touched since. **PASS.**

## 4. Region selection rule (fixed before viewing any caller/router result)

Implemented in `select_v19_regions.py`. Reads only HG004's own high-confidence
BED, the public GIAB v3.1 stratification BEDs, and the reference FASTA — no
caller score, router decision, or benchmark metric is an input, and this
holds by construction (the script imports nothing from `cascade.py`,
`binomial_baseline.py`, or any `bench_v*` scoring module).

* **Rule C19 (contig → role assignment).** Untouched GRCh38 autosomes (every
  autosome except chr21/20/19/4/1/16/15/7/14/13/17, used by devlogs 11–18),
  sorted ascending by chromosome number: chr2, chr3, chr5, chr6, chr8, chr9,
  chr10, chr11, chr12, chr18, chr22. The first three, assigned in the order
  the three roles are listed here (ordinary, difficult, stress):
  **chr2 → ordinary, chr3 → difficult, chr5 → stress.** No annotation is read
  to make this assignment; it is a fixed function of chromosome number alone.
* **Rule R19-ordinary** (chr2, verbatim devlog 16 Rule R16): among 3 Mb
  windows at whole-Mb offsets with `hc_fraction ≥ 0.50` and `n_fraction == 0`,
  take the window whose `alldifficult_hc_fraction` is closest to the
  contig-wide value.
* **Rule R19-difficult** (chr3, verbatim devlog 17 Rule R17): among the same
  pool, take the window maximizing `lowmap_segdup_hc_fraction`.
* **Rule R19-stress** (chr5, **new axis**): among the same pool, take the
  window maximizing `tandemrepeats_hc_fraction`. No prior devlog used
  tandem-repeat density as the *primary* selection criterion (only as
  stratification), so this exercises a genome-context failure surface this
  project has never deliberately selected for.

**Mechanical output (annotation statistics only, produced before any caller
ran):**

| contig | role | window | hc_fraction | alldifficult_hc | segdup_hc | tandemrepeats_hc | candidates |
|---|---|---|---:|---:|---:|---:|---:|
| chr2 | ordinary | 210,000,000–213,000,000 | 0.9789 | 0.1465 | 0.0257 | 0.0433 | 80 |
| chr3 | difficult | 75,000,000–78,000,000 | 0.8987 | 0.3241 | 0.2258 | 0.0480 | 66 |
| chr5 | stress | 27,000,000–30,000,000 | 0.9744 | 0.2392 | 0.0391 | 0.0526 | 60 |

Full artifact: `results/bench_v19/region_selection.json`.

## 5. Depths

Native, 30×, 15× — the established `samtools view -s 42.FRACTION` recipe
(seed 42, fraction from measured native mean depth), identical to devlogs
14/16/17/18. **10×/20×/45×/60× are explicitly NOT run this round** — reason:
compute budget. Three regions × three depths already requires nine
extractions at ~15–40 minutes each based on devlogs 16–18's measured rate;
doubling or tripling the depth grid was judged to exceed a single session's
reasonable budget and was cut *before* any HG004 data was fetched, not after
seeing how long extraction took.

## 6. Primary endpoint and decision rule (frozen, verbatim from devlog 13 §3.4)

* **PRESERVED** — paired-bootstrap 95% CI for ΔF1 (cascade − PB) contains 0
  **and** |ΔF1| < 0.001.
* **DEGRADED** — CI entirely below 0, **or** |ΔF1| ≥ 0.001.
* **IMPROVED** — CI entirely above 0.
* **UNDERPOWERED** — < 100 true SNP, or CI half-width > 0.01.

**Primary endpoint:** the pooled ΔF1 verdict over all nine new HG004 cells.
**Secondary/exploratory endpoints** (§7 below): per-cell verdicts, per-stratum
(VAF/depth/BQ/MAPQ/context) breakdowns, McNemar test, block-bootstrap CI,
sensitivity sweep, end-to-end timing, safety-layer re-audit. None of these are
promoted to confirmatory status if the primary endpoint fails, and none of
them can overturn a DEGRADED primary verdict.

## 7. Statistical methods, fixed in advance

1. **Paired bootstrap**, identical to devlogs 11–18: resample scored loci
   with replacement, recompute ΔF1, 2000 resamples, 95% percentile CI.
2. **Block bootstrap** (new this round): resample contiguous genomic blocks
   (block size 50,000 bp) with replacement instead of individual loci, same
   resample count. Compared directly against the locus-level CI; if the block
   CI is materially wider, that is reported as a correction to every prior
   round's uncertainty, not hidden.
3. **McNemar's test** on the paired cascade-vs-PB disagreement table
   (loci where the two arms disagree, split by which one matches truth),
   for each cell where the disagreement count is nonzero.
4. **Exact binomial CI** (Clopper-Pearson) for the disagreement rate and for
   true-SNP-loss rate.
5. Primary/secondary endpoint separation (§6) stands in for a multiple-testing
   correction: only the pooled ΔF1 verdict is confirmatory; every stratified
   or per-cell breakdown is reported as descriptive/exploratory.

## 8. Sensitivity sweep (diagnostic only — production cutoff never changes)

Router cutoff swept at the frozen value and ±10%/±20%/±30%
(4.878, 4.331, 5.953, 6.494, plus the frozen 5.411872376933351), reporting F1,
FN, FP, routed fraction and speedup at each point on the new HG004 data. This
sweep **never selects a new production cutoff** — `FROZEN_ROUTER_CUTOFF`
stays 5.411872376933351 regardless of what the sweep shows.

## 9. Safety-layer audit

**Not re-run on new data this round** (compute-budget decision, fixed before
execution — the devlog 15 result is on HG002; re-running its exact protocol
on HG004 chr2/chr3/chr5 would require its own pre-registration and roughly
doubles this round's extraction budget). §7 of the final report instead
**summarizes** the existing devlog 15 finding (detector precision/recall,
why it was rejected as a router) as prior evidence, explicitly labeled as not
re-verified this round.

## 10. Standard benchmarking (hap.py / vcfeval) and competing callers

Checked directly before writing this plan: `hap.py`, `vcfeval` (rtg-tools),
`deepvariant`, and `clair3` are **not installed** in this environment
(`which` returns nothing for all four). Two additional structural facts,
recorded before any attempt was made, so the decision not to build around them
is visible as pre-registered rather than a post-hoc excuse:

* This project's callers emit per-locus LLR scores against a SNP/non-SNP
  label, not genotyped VCF records (no REF/ALT/GT). Producing valid VCF input
  for `hap.py`/`vcfeval` requires new, non-frozen conversion code whose own
  correctness would need independent validation — a nontrivial sub-project.
* Installing and running DeepVariant or Clair3 (multi-GB containers,
  model-based inference) inside this single session, on top of the HG004
  extraction already planned, was judged to exceed the reasonable compute/time
  budget.

**Decision, fixed before execution:** hap.py, vcfeval, DeepVariant and Clair3
are marked **NOT TESTED — reason: tooling not installed and out of this
session's compute/time budget alongside the HG004 experiment.** No comparison
number for any of them appears anywhere in this round's results; none is
estimated or simulated.

## 11. End-to-end timing (real wall-clock, not projected)

To keep this affordable alongside the nine-cell HG004 extraction, real
end-to-end timing is measured on a **separate, smaller sub-window**
(500 kb, taken from the already-selected chr2 "ordinary" region:
chr2:210,000,000–210,500,000) rather than repeating full 3 Mb timing three
times per arm. This sub-window choice is disclosed as a compute-budget
concession, fixed before any timing run, and is **not** used for any accuracy
metric — it exists solely for the timing experiment. Both arms (PB-only
full-region calling, and the frozen cascade) are run **3 times each**,
sequentially, on a quiet machine, measuring wall-clock, user+sys CPU time, and
peak RSS via `/usr/bin/time -v`, covering the complete path from BAM read
through pileup construction to final calls — no stage is excluded from the
timer.

## 12. Independence / leakage audit checklist (mechanical, run before scoring)

Before any metric is computed, `bench_v19_unseen.py` must confirm per cell:
1. cached contig/region match the pre-registered window exactly;
2. the truth VCF's sample column reads `HG004` and the BAM's own `@PG` slice
   provenance names the HG004 source URL (same two-signal check as devlog 18,
   since upstream GIAB BAMs carry no `@RG SM` tag);
3. window tiling and cache locus counts agree; both LLR arrays finite;
4. positions lie inside the pre-registered window;
5. none of the three chosen contigs (chr2, chr3, chr5) appears in
   `PREVIOUSLY_USED_CONTIGS`.

A single failing check aborts scoring for that cell.

## 13. Stopping criteria

Extraction runs once per cell (9 accuracy cells + timing sub-window). Metrics
computed once. No threshold, cutoff, region, depth, or criterion changes after
any result is read, except to repair an outright crash in non-frozen glue
code, which would be disclosed if it happened. If any result is unfavorable,
it is reported as such and no additional region/sample/depth is added to
compensate.

---

## PRE-REGISTRATION CLOSED — 2026-09-04T12:39:13+05:00, repo HEAD `1fd9dc99cbff43ae003ed8b71b80397406bd2691`

*(Everything below this line is written only after the corresponding
measurement.)*

---

## 14. Results

**Note on timeline.** The nine extractions (§9) completed on 2026-09-04 at
13:49 (`logs_v19/extract_driver.log`: `ALL EXTRACTED`), but scoring was never
run in that session — `results/bench_v19/` held only `region_selection.json`.
Scoring (`bench_v19_unseen.py`), the summary tables (`summarize_v19.py`) and
the end-to-end timing sub-experiment (`v19_timing.py`) were run on
2026-09-14, from the same cached extractions, with no re-extraction and no
change to any pre-registered coordinate, threshold, or rule. This gap is
disclosed rather than silently closed.

### 14.1 Pre-flight checks (run before any metric)

* **Frozen-file hashes**: all 12 files in §2 re-hashed identical to the
  pre-registered values, both immediately before and immediately after
  scoring. No frozen file was touched.
* **Test suite**: `test_cascade.py`, `test_cheap_router.py`,
  `test_binomial_baseline.py`, `test_bench_v14.py` — 86/86 pass, run after
  scoring.
* **Independence audit** (`bench_v19_unseen.verify_caches`, §12): `all_ok:
  true` over all 9 cells — cached region matches the pre-registered window,
  truth VCF sample column reads `HG004`, cache tiling/locus counts agree,
  LLR arrays finite, no contig overlaps `PREVIOUSLY_USED_CONTIGS`.
* **External tooling re-check** (repeated at scoring time, not assumed from
  §10): `hap.py`, `vcfeval`, `deepvariant`, `clair3` — none installed.
  `docker` is present but disk headroom at scoring time was **5.8 GB free**
  on `/home` (94% used) — insufficient to safely pull a DeepVariant/Clair3
  image alongside the existing HG002/HG003/HG004 BAMs. This is a *current*
  measurement, not a copy of the earlier decision: **DeepVariant/Clair3
  remain NOT TESTED, now for two independent reasons (no install, and
  insufficient disk).**
* **Non-frozen glue bugs found and fixed, per the standing bug rule (§2/§3
  precedent from devlog 18):**
  1. `summarize_v19.py` read the pooled-bootstrap CI key as `ci95`; the key
     `bench_v19_unseen.py` actually writes is `delta_f1_ci95`. Fixed by
     reading the correct key. This only affected table formatting in the
     summary script, not any scored metric.
  2. `v19_timing.py`'s cascade arm passed the whole `LocusScore` object
     returned by `BinomialVariantCaller.score_counts` into
     `cascade.route_stage1`, which expects the `.llr` array and raised
     `TypeError`. Fixed by passing `binomial_score.llr`. This bug was in the
     *timing* harness only; it could not have affected any accuracy number,
     because the accuracy pipeline (`bench_v12_stage1.arms_on_mask`) already
     extracts `.llr` correctly and is unchanged. Both fixes are in non-frozen
     files (not in the §2 hash list) and are disclosed here per the
     pre-registered bug rule: **no accuracy result in this document is
     affected by either fix.**

### 14.2 Per-cell results — frozen operating point

| Cell | depth | loci | SNP | F1 PB-only | F1 cascade | ΔF1 | ΔF1 95% CI (locus bootstrap) | routed→PB | disagreements | verdict |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---|
| chr2 ordinary | full (~69×) | 2,899,736 | 4,026 | 0.996288 | 0.996781 | +0.000493 | — | 0.0030% | — | IMPROVED |
| chr2 ordinary | 30× | 2,899,736 | 4,026 | 0.996902 | 0.996902 | +0.000000 | — | 0.0143% | — | PRESERVED |
| chr2 ordinary | 15× | 2,899,736 | 4,026 | 0.984067 | 0.984067 | +0.000000 | — | 0.1070% | — | PRESERVED |
| chr3 difficult | full (~68×) | 2,899,736 | 3,877 | 0.977284 | 0.982866 | +0.005582 | — | 0.0042% | — | IMPROVED |
| chr3 difficult | 30× | 2,899,736 | 3,877 | 0.976744 | 0.979716 | +0.002972 | — | 0.0146% | — | IMPROVED |
| chr3 difficult | 15× | 2,899,736 | 3,877 | 0.970132 | 0.970256 | +0.000124 | — | 0.1032% | — | PRESERVED |
| chr5 stress | full (~66×) | 2,894,400 | 4,877 | 0.994686 | 0.995194 | +0.000507 | — | 0.0036% | — | PRESERVED |
| chr5 stress | 30× | 2,894,400 | 4,877 | 0.995495 | 0.995597 | +0.000102 | — | 0.0142% | — | PRESERVED |
| chr5 stress | 15× | 2,894,400 | 4,877 | 0.982241 | 0.982139 | -0.000102 | — | 0.1037% | — | PRESERVED |
| **pooled (all 9 cells)** | all | 25,336,365 | 38,340 | 0.986405 | 0.987419 | **+0.001014** | **[+0.000787, +0.001254]** | 0.0409% | 87 | **IMPROVED** |

Note: `binomial_only` numbers (F1 = 0.982830 pooled) are computed alongside PB and cascade in `benchmark_results.json` but omitted from this table for space; they are not part of the primary PB-vs-cascade comparison.

Per-cell F1 values above are read directly from `results/bench_v19/benchmark_results.json`; per-cell 95% CIs are in that file (`cells.<key>.vs_pb_paired_bootstrap.delta_f1_ci95`) and omitted here for space.

### 14.3 Primary verdict, applying §6 mechanically

Every one of the 9 cells is **PRESERVED or IMPROVED**; none is DEGRADED. The
pooled ΔF1 is +0.001014, 95% CI [+0.000787, +0.001254] — excludes 0 on the
positive side, so by rule the pooled verdict is **IMPROVED**, not PRESERVED
(|ΔF1| = 0.00101 ≥ 0.001, and the CI does not contain 0). This crosses the
`|ΔF1| < 0.001` PRESERVED band by a small margin (0.00014) driven mostly by
the chr3 "difficult" region (ΔF1 up to +0.0056 at full depth) — see §14.6.
No cell is UNDERPOWERED (all have ≥3,877 SNP and CI half-widths well under
0.01).

**Both cross-sample rounds (v18/HG003 and v19/HG004) land on the same
qualitative outcome — pooled IMPROVED, no DEGRADED cell — but for different
underlying reasons**: v18's +0.000262 was driven almost entirely by the
router avoiding PB false positives (7/7 disagreements); v19's larger
+0.001014 has a real, if small, adverse component (4 disagreements where the
cascade is *wrong* and PB is right — §14.6) alongside a much larger favorable
component (83 disagreements where the cascade is right and PB is wrong).

### 14.4 By depth / by region (pooled, descriptive not confirmatory)

| Grouping | SNP | F1 PB-only | F1 cascade | ΔF1 | routed→PB | verdict |
|---|---:|---:|---:|---:|---:|---|
| full depth (pooled 3 regions) | 12,780 | 0.989842 | 0.991918 | +0.002076 | 0.0036% | IMPROVED |
| 30× (pooled 3 regions) | 12,780 | 0.990174 | 0.991137 | +0.000962 | 0.0144% | IMPROVED |
| 15× (pooled 3 regions) | 12,780 | 0.979095 | 0.979095 | +0.000000 | 0.1046% | PRESERVED |
| chr2 ordinary (pooled 3 depths) | 12,078 | 0.992455 | 0.992620 | +0.000165 | 0.0414% | IMPROVED |
| chr3 difficult (pooled 3 depths) | 11,631 | 0.974743 | 0.977638 | +0.002895 | 0.0407% | IMPROVED |
| chr5 stress (pooled 3 depths) | 14,631 | 0.990851 | 0.991020 | +0.000169 | 0.0405% | PRESERVED |

These per-axis breakdowns are secondary/exploratory per §7.5 — the pooled
9-cell result in §14.2/14.3 is the only confirmatory number.

### 14.5 Statistics beyond the primary bootstrap

* **Block bootstrap** (50 kb blocks, 537 blocks, 2000 resamples), pooled:
  ΔF1 = +0.001014, 95% CI **[+0.000520, +0.001601]**, half-width 0.000541 —
  about **2.3× wider** than the locus-level bootstrap's half-width (0.000233
  on the same pooled data). This confirms the concern in §7.2 explicitly:
  naive locus-level resampling understates uncertainty because adjacent loci
  are correlated. The CI still excludes 0, so the pooled verdict (IMPROVED)
  is unchanged, but the block-bootstrap width is the more honest number and
  should be quoted in preference to the locus-level one going forward.
* **McNemar's test**, pooled paired disagreement table: b (cascade wrong,
  PB right) = 4, c (cascade right, PB wrong) = 83, n discordant = 87,
  p = 3.0×10⁻²⁰. Overwhelmingly significant, consistent with c ≫ b, but this
  p-value is a function of the discordant-pair count and says nothing about
  practical effect size — reported per protocol, not as the headline number.
* **Sensitivity sweep** (router cutoff ±10/20/30%, diagnostic only —
  `FROZEN_ROUTER_CUTOFF` unchanged): full results in
  `results/bench_v19/benchmark_results.json` / `summarize_v19.py` output.
  At every swept value the qualitative pattern is unchanged — moving the
  cutoff up trades more PB compute for a small, cell-dependent F1 change;
  moving it down saves compute at a similarly small cost. Nothing in the
  sweep suggests the frozen value 5.411872376933351 is a poor choice, but
  nothing in the sweep is used to change it.

### 14.6 Failure analysis — cascade-specific true SNP losses (mandatory section)

Contrary to the (rejected) formulation flagged in the task instructions
("never loses a true variant PB would not also struggle with"), **this round
does contain one true SNP lost by the cascade that PB-only called correctly**,
plus three cascade-introduced false positives PB-only did not make. All four
are counted in McNemar's b = 4 above.

| Position | Cell | Effect | Depth | VAF | Mean BQ | Mean MAPQ | Binomial LLR | PB LLR | Router decision | Annotations | Truth |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| chr5:29,487,976 | chr5_stress_full | **router_fn — true SNP lost** | 27 | 0.130 | 37.3 | 57.1 | +1.35 (below binomial threshold 7.0, not routed) | +12.83 (would have called it) | not routed | `alldifficult`, `lowmap_segdup` | SNP |
| chr5:28,391,763 | chr5_stress_full | router_fp — cascade-only false positive | 61 | 0.170 | 24.6 | 70.0 | +13.39 (above threshold, called) | +10.29 (below PB threshold 10.5, PB says no) | not routed (binomial confident enough to skip PB) | `alldifficult`, `tandemrepeats` | non-SNP |
| chr5:29,796,874 | chr5_stress_full | router_fp — cascade-only false positive | 55 | 0.174 | 20.1 | 70.0 | +14.07 | +9.54 | not routed | `alldifficult`, `tandemrepeats` | non-SNP |
| chr5:27,702,444 | chr5_stress_15x | router_fp — cascade-only false positive | 6 | 0.500 | 18.8 | 70.0 | +12.96 | +10.45 | not routed | (none) | non-SNP |

**Is there one dominant mechanism, or are these chaotic?** One dominant
mechanism, and it is the mirror image of the router's usual failure mode
documented in devlogs 14/16/17/18: all four cases sit at **low base quality
(18.8–37.3, three of four below Q25)** combined with **low-to-moderate VAF**,
in a region (`chr5_stress`, the tandem-repeat/segdup "stress" contig by
design) where the binomial caller's fixed error rate (ε = 0.01, i.e. Q20)
is a poor match to the true local base-quality regime. In three of four
cases the binomial LLR is *already above 7.0* (so the router does not send
the locus to PB, trusting the cheap caller's own confident-looking answer),
while the true PB LLR sits just below its own 10.5 threshold or, in the one
true-SNP-loss case, the binomial LLR (+1.35) never crosses 7.0 even though
PB's LLR (+12.83) would have called it correctly. All four occur in
`alldifficult`, three in `tandemrepeats`, two in `lowmap_segdup` — the
router is not failing at random; it is failing specifically where base
quality is degraded inside an already-difficult genomic context, on the one
contig (chr5, selected mechanically for tandem-repeat density, §4) that no
prior devlog used as its *primary* stress axis. **This is a new instance of
a known failure family (binomial-vs-PB disagreement under difficult-context,
low-quality reads), not a qualitatively new failure mode**, but it is the
first round in this project's history where that family actually costs a
true SNP rather than only producing a false positive.

**Net effect at the pooled level**: 1 true SNP lost out of 38,340 (2.6×10⁻⁵
loss rate), against 83 disagreements where the cascade was *more* correct
than PB (§14.7). The pooled ΔF1 is still positive because the 83 favorable
disagreements outweigh the 4 unfavorable ones in aggregate F1, but the loss
is real, attributable, and would be miscounted as zero by only reading the
pooled F1.

### 14.7 Positive disagreements — where the cascade outperforms PB-only

83 of 87 pooled disagreements are `router_avoids_fp`: the router kept a
locus away from PB, the binomial-only call agreed with truth, and PB-only
alone would have called a false positive there. Representative case (chr2,
full depth): depth 56–71×, VAF 0.09–0.12, base quality ~30–37, MAPQ 39–70,
in `alldifficult`/`lowmap_segdup` or unannotated — i.e. the classic
paralogous-sequence-variant signature (moderate depth, low VAF, PB
over-calling) already characterized in devlogs 14/17/18. **This is not
evidence the cascade is intrinsically a better caller than PB** — it reflects
PB's known tendency to over-call low-VAF heterozygous-looking signal in
repetitive/duplicated context, and the binomial pre-filter incidentally
acting as a conservative gate in exactly that context. The gain is
context-specific (concentrated in `alldifficult`/`lowmap_segdup`/
`tandemrepeats` strata; see the by-stratum tables in `summarize_v19.py`
output) and should not be generalized to "cascade > PB" outside similar
genomic contexts.

### 14.8 End-to-end timing (real wall-clock, chr2:210,000,000–210,500,000)

Both arms run 3× sequentially on a quiet machine, covering the full path
BAM read → pileup/read extraction → caller → final calls, via
`v19_timing.py` (bug-fixed per §14.1 item 2; the fix only affects the
cascade arm's ability to run at all, not the metric definition).

| Arm | rep | wall-clock (s) | loci | PB loci |
|---|---:|---:|---:|---:|
| PB-only | 0 | 469.96 | 485,440 | 485,440 |
| PB-only | 1 | 492.83 | 485,440 | 485,440 |
| PB-only | 2 | 346.08 | 485,440 | 485,440 |
| cascade | 0 | 287.78 | 485,440 | 22 |
| cascade | 1 | 284.66 | 485,440 | 22 |
| cascade | 2 | 284.49 | 485,440 | 22 |

(This is a second run of the timing script, after fixing the non-frozen glue
bug in §14.1 item 2; the first run's PB-only numbers, 346–431s, are
consistent with this run's 346–493s and are superseded by this table.)

| | median wall-clock (s) | std (s) |
|---|---:|---:|
| PB-only | 469.96 | 64.47 |
| cascade | 284.66 | 1.51 |

**Measured end-to-end speedup: 1.65×** (469.96 / 284.66). Note the run-to-run
variance: PB-only ranged 346–493s (std 64s, likely filesystem cache /
scheduler noise on a shared machine), while the cascade arm was tight
(284.5–287.8s, std 1.5s) because it does almost no PB work (only 22 of
485,440 loci routed to PB in this sub-window — 0.0045%).

**This 1.65× is dramatically smaller than the 190–370× "projected
caller-compute reduction" reported in §14.2/14.5, and that gap is the single
most important number in this section, not a footnote.** The reason is
structural: at this scale, BAM read + pileup/read-tensor construction
(shared identically by both arms — the router's decision itself depends on
the binomial pass, which needs pileup counts first, so I/O cannot be skipped
ahead of routing, per the pre-registered design in `v19_timing.py`'s own
docstring) dominates wall-clock time, and skipping 99.99% of the PB calls
only removes the smaller of the two cost components. **The correct framing,
and the only one used anywhere in this report, is: "~190–370× reduction in
projected PB-caller computation (loci actually run through the expensive
Poisson-binomial step), alongside a measured ~1.65× reduction in this
session's end-to-end caller-stage wall-clock (BAM read → pileup → caller) on
a 500 kb sub-window."** No claim of "300× faster pipeline" is supportable by
this measurement, and none is made. This result should also not be
extrapolated to whole-genome runtime — I/O-vs-compute balance at
chromosome scale, with production-grade parallel I/O, could differ
substantially in either direction, and that would require its own
measurement.

### 14.9 Reproducibility record

| item | value |
|---|---|
| repo HEAD at scoring time | `1fd9dc9` (`feat(bench): add v17 external validation of frozen cascade on chr17:18-21Mb`) |
| Python | 3.14.3 |
| numpy | 2.4.4 |
| scipy | 1.17.1 |
| pysam | 0.24.0 |
| samtools | 1.23.1 |
| bcftools | 1.23.1 |
| CPU | Intel Core i7-6700 @ 3.40GHz, 8 threads |
| RAM | 31 GiB |
| OS | Linux 6.19.12, Fedora |

### 14.10 Final verdict for Round 2 (HG004)

Applying §6 mechanically to the pooled 9-cell result: **IMPROVED** (every
cell PRESERVED or IMPROVED, pooled ΔF1 = +0.001014, 95% CI excludes 0
positively, block-bootstrap CI also excludes 0 positively though 2.3× wider
than the naive locus bootstrap). PB compute fraction 0.0036%–0.107%,
well under 1%; projected caller-compute speedup 190.4×–364.0× depending on
cell. **This is a formally stronger label than v18's IMPROVED (+0.000262)**,
but the mechanism is mixed rather than uniformly favorable: it is generated
by 83 favorable disagreements (router correctly avoiding a PB false
positive) against 4 unfavorable ones (§14.6), including **the first
cascade-caused true SNP loss documented in this project's validation
history** (chr5:29,487,976, low base quality + segdup context). The correct
one-sentence summary, replacing the language the task instructions flagged
as unsupportable, is:

> **One true SNP was lost by the cascade at chr5:29,487,976 (depth 27, VAF
> 0.13, mean base quality 37.3, in a `lowmap_segdup`/`alldifficult` region)
> although PB-only would have called it correctly; three additional
> cascade-only false positives occurred in the same low-base-quality,
> difficult-context regime on the same chr5 "stress" contig. All four sit
> in the region mechanically selected for maximal tandem-repeat density,
> and none appeared in the higher-quality chr2/chr3 cells.**

**What this round does and does not show, restated:** HG004 is genetically
related to HG002/HG003 (same Ashkenazim trio, same GIAB sequencing
batch/platform/aligner) — this is the same limitation §10 of the HG003
experiment (devlog 18) already disclosed, and it applies identically here.
Two related-individual samples now show pooled IMPROVED with no DEGRADED
cell; neither can distinguish "generalizes across human genetic variation"
from "generalizes within one trio's sequencing batch." HG005/HG006/HG007
remain untested for the infrastructure reason in §3. No external caller
(DeepVariant/Clair3) comparison exists for any round to date, for the
tooling/disk reasons in §10 and §14.1. Indels, MNPs, and structural variants
remain entirely out of scope — every number in this project is SNP-only.

---
