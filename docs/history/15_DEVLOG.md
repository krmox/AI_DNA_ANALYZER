# Can the frozen cascade's one failure mode be detected cheaply? — an ε-sensitivity safety layer

**Date:** 2026-08-17

**Status of sections:** §1–§8 are the pre-registration. They were written and
saved **before any caller, router, safety rule or metric was run on any chr16 /
chr15 / chr7 / chr14 locus**, and before the reference FASTAs for those contigs
had finished downloading. Everything from §9 onwards was written after the
corresponding measurement. §7.3 discloses, in advance, exactly which parts of
the design were informed by looking at the *validation* data — because they
were, and hiding that would make the pre-registration worthless.

Devlogs 3–14 are not edited by this experiment. `cascade.py`, `cheap_router.py`,
`binomial_baseline.py` and `quality_error_model.py` are not edited by this
experiment; §16 verifies that mechanically.

---

## 1. Research question

Devlog 14 validated the frozen cascade

```
binomial LLR (ε = 0.01) → route to PB if |binomial_LLR − 7.0| ≤ 5.411872376933351 → call
```

across four new chromosomes and 10.37 M loci, and found exactly one failure
mode. Quoting devlog 14 §17.3 verbatim, because this whole experiment is a
response to that one sentence:

> The frozen cascade departs from PB exactly at loci where the fixed-ε binomial
> is **confidently wrong** — and the router, which measures *uncertainty* rather
> than *error*, cannot see them by construction.

Concretely: 5 true SNP that PB-only finds were lost by the cascade, all inside
`chr1_segdup`, all at VAF 0.11–0.13 with high base quality. Widening the router
cutoff does not reach them (devlog 14 §13) — they have *large* router margins
precisely because the binomial is confident.

The question here is **not** whether the router can be improved. The router is
frozen. The question is:

> Can a stage placed *after* the frozen router detect "the binomial is
> confidently wrong" from statistics the pipeline has already computed, cheaply
> enough that the cascade's ~170–400× compute advantage survives?

A repair that costs a large fraction of PB compute is not a repair; it is a
slower way of running PB. That is stated now, before any number exists, and it
is written into the acceptance criteria in §7.5 as criterion **A3**.

## 2. Hypotheses

* **H1 (the signal exists and is cheap).** The failure regime is ε
  misspecification, not an irreducible property of the counts. The binomial's
  confidence at the failing loci is an artefact of the single constant
  ε = 0.01; re-evaluating the *same* likelihood ratio at bracketing error rates
  flips the call at those loci and at almost nowhere else. If so, the detector
  costs two extra `binom.logpmf` evaluations over `(k, n)` the caller has
  already produced — no pileup re-read, no read tensor, no per-read quality, and
  no PB call to decide whether PB should run.
* **H2 (the budget is small).** Because the flip condition is a *sensitivity*
  condition rather than a *support* condition, the escalated set is far smaller
  than any raw support threshold (`k ≥ t`) achieving the same repair, and small
  relative to the frozen router's own routed fraction (7.0×10⁻⁴ on the devlog 14
  data).
* **H3 (it transfers).** The regime is a property of the *statistic*
  (fixed-ε binomial vs per-read quality model) and of *segmental duplication*,
  not of chr1:120–121 Mb. So a rule frozen on the devlog 14 cells should repair
  the same failure on segdup windows of chromosomes it has never seen.
* **H4 (the repair is not free, and the cost is precision).** Devlog 14 §17.1
  found the cascade's departures from PB are *net favourable*: 26 avoided PB
  false positives against 5 lost true SNP, which is why its pooled ΔF1 vs PB is
  **positive** (+0.000581). Any layer that escalates the low-VAF/high-quality
  corner will hand some of those avoided false positives back to PB. Prediction:
  the safety layer raises recall, lowers precision, and lands ΔF1 between
  PB-only and the frozen cascade rather than above both.

H4 is a prediction that the intervention is *worse on the primary metric of
devlogs 11–14* while being better on the failure mode that motivated it. It is
recorded now so that it cannot be presented later as an unexpected nuance.

## 3. Existing baseline — frozen, verified, unmodified

Nothing below is re-derived, re-fit or re-selected. Verified byte-identical to
`HEAD` (commit `3837b61`) before any work in this session; §16 repeats the check
afterwards.

| Constant | Value | Source |
|---|---|---|
| `FROZEN_BINOMIAL_THRESHOLD` | 7.0 | imported from `cascade.py` |
| `FROZEN_PB_THRESHOLD` | 10.5 | imported from `cascade.py` |
| `FROZEN_ROUTER_CUTOFF` | 5.411872376933351 | imported from `cascade.py` |
| binomial ε | 0.01, `error_floor` 1e-3 | `BinomialVariantCaller` defaults |
| PB ε floor / ceiling | 1e-4 / 0.25 | `quality_error_model` |
| `MAX_READS` | 48 | `read_level_pileup` |
| PB implementation | the repaired one (devlog 13 §9) | `quality_error_model.poisson_binomial_llr` |

Baseline integrity, checked at the start of this session and recorded here:

```
git status --porcelain           -> (empty; clean tree at 3837b61)
cascade.py            sha256 b0ee9f4b24fe06dc…  IDENTICAL-TO-HEAD
cheap_router.py       sha256 5e13cca000f4565a…  IDENTICAL-TO-HEAD
binomial_baseline.py  sha256 f0408ecae0929aef…  IDENTICAL-TO-HEAD
quality_error_model.py sha256 a69cb49591cbac13… IDENTICAL-TO-HEAD
robustness_benchmark.py sha256 6fd5c4df39806e05… IDENTICAL-TO-HEAD
pytest (8 relevant suites)  -> 190 passed
```

The routing rule is used verbatim through `cascade.route_stage1`; the safety
layer takes that function's output mask as an *input* and returns a mask
**disjoint** from it, so no locus the router declined can be un-declined and no
locus it routed can be un-routed. This is asserted at runtime inside
`safety_layer.safety_escalate` and tested in `test_safety_layer.py`.

## 4. Motivation from the discovered failure regime

Devlog 14 §17 characterised 41 router/PB disagreements over 10.37 M loci:

| effect | count | what it means for the cascade |
|---|---:|---|
| `router_avoids_fp` | 26 | favourable — PB false positive dodged |
| `router_fp` | 8 | unfavourable — cascade false positive |
| `router_fn` | **5** | **unfavourable — true SNP lost, the failure being repaired** |
| `router_recovers_fn` | 2 | favourable |

and located both unfavourable populations on one axis: the fixed-ε binomial is
confidently wrong, over-calling below ~Q25 and under-calling above ~Q30, most
strongly at VAF ≲ 0.25 and high depth, concentrated in segmental duplication.

§9 of this devlog re-derives that characterisation independently from the cached
per-locus evidence rather than taking devlog 14's summary on trust, because §2's
whole design rests on the mechanism being ε misspecification and not something
else that merely correlates with it.

## 5. Data provenance

Two disjoint sets. Nothing is synthesised; no artificial difficulty is injected;
no read is modified.

### 5.1 Validation set — the devlog 14 cells (already used, already read)

The 12 cached cells of devlog 14: `chr20_neutral`, `chr19_gcrich`,
`chr4_atrich`, `chr1_segdup`, each at native (~70x) / 30x / 15x. 10,368,078
scored loci, 12,261 SNP. Provenance is devlog 14 §5.5 and is unchanged. These
cells contain the 5 known failures and were fully inspected by devlog 14, so
they are **discovery data** and cannot support any generalisation claim. They
are used for one purpose: choosing the rule.

### 5.2 Test set — new contigs, untouched until §11

* Sample: GIAB **HG002** (NA24385 son), GRCh38 — held constant deliberately.
* Reads: `HG002.GRCh38.2x250.bam`, NIST Illumina 2×250 novoalign, sliced
  remotely per region with `samtools view -b <URL> <contig>:<start+1>-<stop>`.
  Same URL and recipe as devlogs 12–14.
* Truth: GIAB **v4.2.1** benchmark VCF, region-sliced with `bcftools view -r`.
* High-confidence BED:
  `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`.
* Reference: Ensembl release-110 GRCh38 per-chromosome FASTA (contigs `16`,
  `15`, `7`, `14`) — same source and release as every existing
  `data/reference/chr*_full.fa`.
* Difficulty annotation: GIAB genome-stratifications **v3.1**, GRCh38 —
  `lowmap_segdup`, `alldifficult`, `tandemrepeats`, the same three files
  devlog 14 used.
* Retrieved 2026-08-17. Exact commands: `fetch_v15_reference.sh`,
  `fetch_v15_data.sh`.

Cross-sample and cross-platform transfer remain untested, as in devlog 14, and
no claim will be made about them.

### 5.3 Test region selection — mechanical and results-blind

Implemented in `select_v15_regions.py`, which reads only the GIAB
high-confidence BED, the GIAB v3.1 stratification BEDs and the reference FASTA.
It cannot see a caller score, a router decision, a safety-layer decision or a
metric. Its contig rule was run to completion **before any FASTA for a new
contig existed**; its window rules run before any caller.

* **Rule C (which contigs).** Among the GRCh38 autosomes excluding every contig
  any previous devlog used (`chr21`, `chr20`, `chr19`, `chr4`, `chr1`), rank by
  `|HC ∩ lowmap_segdup| / |HC|` over the whole contig. Top three are the hard
  test contigs; the bottom one is the ordinary-sequence control contig. Ties
  broken by the smaller contig number.
* **Rule H (hard window)** — verbatim from devlog 14 §5.2, applied to each hard
  contig: among 1 Mb windows at whole-Mb offsets with `hc_fraction ≥ 0.50` and
  `n_fraction = 0`, take the one maximising the segdup fraction of its
  high-confidence bases.
* **Rule O (neutral window)** — verbatim from devlog 14 §5.2, applied to the
  control contig: among windows with `hc_fraction ≥ 0.90`,
  `difficult_fraction ≤ 0.15`, `n_fraction = 0`, take the one closest to the
  contig midpoint.

Rule H is applied to three contigs rather than one because the failure being
repaired occurs only in segmental duplication: a test set that cannot contain
the failure mode cannot test a repair for it. The single Rule O window exists so
the layer's **false-escalation cost** is measured on ordinary sequence, which is
where almost all of the genome lives. This is a deliberate enrichment of the
test set for the failure regime and it is a limitation as well as a design
choice — it is recorded as such in §17.

Rule C output, computed before the pre-registration was finalised and before any
new FASTA existed (`results/bench_v15/contig_selection.json`):

| contig | segdup fraction of HC bases | role |
|---|---:|---|
| chr16 | 0.0910 | hard test window (Rule H) |
| chr15 | 0.0891 | hard test window (Rule H) |
| chr7 | 0.0872 | hard test window (Rule H) |
| … 14 other contigs … | 0.0861 → 0.0334 | not used |
| chr14 | 0.0334 | ordinary control window (Rule O) |

**Disclosed rule revision:** the first version of `CANDIDATE_CONTIGS` listed
`chr2…chr18, chr22` and so failed to exclude `chr4`, which devlog 14 used. It
was corrected to exclude every previously used contig. `chr4` ranked 12th of 18
under the buggy version and was selected neither as a hard contig nor as the
control, so the corrected and uncorrected rules return **the same four contigs**;
the fix changes no coordinate. Recorded because pre-registration hygiene is
worth more than a tidy narrative.

Each region is extracted at three depths — full (native), 30x, 15x — by
`samtools view -s` with a fixed seed, exactly as devlog 14 §5.4. Depth
realisations of one region are downsamples of one library and are treated as
three regimes of one region, never as three independent regions. **12 test
cells: 4 regions × 3 depths.**

## 6. Independence audit — planned checks

Run mechanically in §11 and recorded in
`results/bench_v15/independence_audit.json`:

1. **Chromosome.** Every test cell's contig ∈ {chr16, chr15, chr7, chr14}; the
   intersection with {chr21, chr20, chr19, chr4, chr1} — every contig used by
   devlogs 3–14 — is empty. Coordinate overlap with any previous training,
   validation, tuning or test region is therefore impossible, not merely
   unlikely.
2. **Coordinate.** Each cell's cached `region` equals its pre-registered span,
   every reconstructed locus coordinate lies inside it, and the locus count
   matches an independent re-tiling of the windows.
3. **Sample.** HG002 in both sets, deliberately. This is a *dependence* and is
   declared: the experiment tests transfer across chromosome and structure, not
   across sample.
4. **Dataset provenance.** Same BAM URL, same truth VCF release, same HC BED,
   same reference release, same stratification release as validation — so the
   only thing that differs between validation and test is *which sequence*.
5. **Selection blindness.** `select_v15_regions.py` imports nothing from any
   caller, router or benchmark module; asserted by an import-graph check.
6. **Rule freezing.** The frozen rule's parameters are written into
   `safety_layer.FROZEN_SAFETY_RULE` and its selection record into
   `results/bench_v15/rule_selection.json` **before** any test cell is
   extracted; the file's mtime versus the test caches' mtimes is recorded.

## 7. Pre-registered protocol

### 7.1 Candidate feature set

Only statistics the pipeline already computes, or functions of them that cost
`O(1)` vectorised arithmetic per locus:

`k` (alternate support), `n` (usable observations), `VAF = k/n`, `depth`,
`quality_sum` and the mean base quality derived from it, `mapping_sum` and the
mean MAPQ derived from it, `gap_reads`, `insertion_reads`, and the fixed-ε
`binomial_llr` itself.

Explicitly **excluded**: `pb_llr` and anything derived from the read tensor or
per-read qualities (using PB to decide whether to run PB is circular and is
banned by the objective); neural models of any kind; any external annotation
(segdup/difficult BEDs are used for *stratified reporting* only, never as a rule
input, because a production caller cannot assume a stratification BED and
because "escalate all of segdup" is a trivially expensive non-answer).

### 7.2 Candidate rule families

1. **Family S — ε-sensitivity bracket.** Escalate when the binomial's call flips
   between two fixed error rates: `(llr(ε_lo) ≥ 7.0) ≠ (llr(ε_hi) ≥ 7.0)`.
   Grid: `(ε_lo, ε_hi) ∈ {(1e-3,1e-2), (1e-3,3e-2), (1e-3,5e-2), (1e-4,1e-2),
   (1e-4,3e-2), (1e-4,5e-2), (1e-4,1e-1), (3e-3,3e-2), (1e-2,5e-2)}`.
2. **Family Sq — quality-derived-ε straddle.** Escalate when the call differs
   between fixed ε = 0.01 and the existing quality-derived ε of
   `BinomialVariantCaller(quality_derived_epsilon=True)`, floor ∈ {1e-4, 1e-3,
   3e-3}. This uses an already-implemented v2 caller, not a new one.
3. **Family K — raw support threshold.** Escalate when `k ≥ t`, `t ∈ {2,3,4,5,6}`.
   The dumb baseline the clever rules must beat on budget.
4. **Family KV — support × VAF band.** `k ≥ t` and `VAF ∈ [lo, hi]`,
   `t ∈ {3,4}`, band ∈ {[0.05,0.35], [0.05,0.25], [0.08,0.30]}.
5. **Family SK — conjunction.** Family S ∧ `k ≥ k_min`, `k_min ∈ {1,2,3,4}`.

All are threshold rules and boolean combinations over ≤ 3 features, evaluated
with integer/float comparisons on arrays the caller already holds. No family is
trained; each is a small enumerated grid.

### 7.3 What was informed by the validation data — disclosed in advance

The families and grids in §7.2 were **not** chosen blind. Before writing this
section, the validation cells were inspected: the 41 disagreements were
re-dumped with their evidence, and the five families above were evaluated on
them. That is what validation data is *for*, and the resulting numbers are
reported in §9–§10 as validation results, not as evidence of anything.

The consequence, stated plainly: **the validation numbers in §9–§10 are
optimistically biased and no claim in §18 may rest on them.** Everything load-
bearing rests on §11–§15, computed once on cells that did not exist when this
section was written.

### 7.4 Rule selection — decided now, executed in §10

Pooled over the 12 validation cells, on the standard scoring frame
(`label ∈ {0, SNP}`), for every candidate in §7.2 compute: escalated fraction of
all loci; total PB fraction (router ∪ escalation); `lost_true_snp`, the number
of true SNP called by PB-only but not by the arm; and ΔF1 vs PB-only.

**Feasible set** = candidates satisfying all of

* (a) `lost_true_snp = 0` (the frozen cascade's value on validation is 5);
* (b) `ΔF1 vs PB-only ≥ 0` — the layer may not make the cascade worse than the
  arm it approximates;
* (c) escalated fraction ≤ 5×10⁻⁴, so total PB compute stays under ~2× the
  frozen router's 7.0×10⁻⁴.

**FROZEN-FULL** = the feasible candidate minimising escalated fraction. Ties
broken by: narrower ε bracket (`ε_hi/ε_lo`), then larger `k_min`, then earlier
in the grid order of §7.2.

**FROZEN-LEAN** = the same procedure with (a) relaxed to `lost_true_snp ≤ 1`.
Reported alongside FROZEN-FULL as the minimum-budget point of the
accuracy-vs-compute curve. If the two coincide, only one is reported.

If the feasible set for FROZEN-LEAN is empty, **no rule is frozen**, the test set
is not touched, and the verdict is NEGATIVE or ARCHITECTURAL FAILURE per §7.6.

Both rules are written into `safety_layer.FROZEN_SAFETY_RULE` /
`results/bench_v15/rule_selection.json` before a single test cell is extracted,
and are not changed afterwards for any reason.

### 7.5 Comparison arms and acceptance criteria

Five arms on the identical scoring frame, per test cell and pooled:

* **A. Binomial-only** — `binomial_llr ≥ 7.0`.
* **B. PB-only** — `pb_llr ≥ 10.5`. The reference arm.
* **C. Frozen router → PB** — the production candidate. **Primary baseline.**
* **D. Frozen router → safety layer (FROZEN-FULL) → PB.**
* **D-lean. Frozen router → safety layer (FROZEN-LEAN) → PB.**

Acceptance criteria, evaluated on the **pooled test set**, fixed now:

* **A1 — the failure is repaired.** `lost_true_snp(D) < lost_true_snp(C)` and
  `lost_true_snp(D) ≤ ½·lost_true_snp(C)`. If `lost_true_snp(C) = 0` on the test
  set, A1 is **UNINFORMATIVE** — the test could not observe the failure mode —
  and the verdict is NULL/INCONCLUSIVE regardless of the other criteria.
* **A2 — accuracy is not sacrificed.** The paired-bootstrap 95% CI for
  ΔF1(D vs PB-only) has lower bound > −0.001 (non-inferiority to PB at the
  equivalence bound inherited from devlogs 13–14).
* **A3 — the compute advantage survives.** Pooled total PB fraction of D ≤ 2× that
  of C, **and** projected end-to-end speedup of D ≥ 0.8× that of C.
* **A4 — the benefit is structure, not budget.** D captures strictly more of C's
  lost true SNP than random routing and than reverse-confidence routing at
  matched additional PB coverage (5 seeds each).

ΔF1(D vs C) is a **reported secondary endpoint, not an acceptance criterion**;
§2 H4 predicts it is negative and it is reported whatever it is.

### 7.6 Verdict mapping, fixed now

* **A. STRONG POSITIVE** — A1, A2, A3, A4 all hold, and the ΔF1(D vs PB) CI
  contains 0 or lies above it.
* **B. WEAK POSITIVE** — A1, A2, A3 hold but A4 is underpowered or not
  demonstrated; or A1 holds with fewer than 3 observed C-failures on test.
* **C. NULL / INCONCLUSIVE** — A1 is UNINFORMATIVE (no C-failures on test), or
  the observed event counts cannot separate the arms.
* **D. NEGATIVE** — A1 fails with C-failures present, or A2 fails, while A3 holds.
* **E. ARCHITECTURAL FAILURE** — A3 fails: the regime can be detected only by
  escalating enough loci to destroy the compute advantage.

The verdict is read off these rules, not chosen.

### 7.7 Stopping criterion

Extraction runs **once** per test cell. Metrics are computed **once**. No region
is added, removed, re-extracted or re-scored after any test result is read,
except to fix an outright software crash — any such re-run is disclosed in §19.
No threshold, cutoff, bracket or `k_min` is changed on the basis of any test
number. If the frozen rule fails, the failure is reported and the rule stays as
frozen. The experiment stops when §§9–17 are complete; no further optimisation
cycle is run regardless of outcome.

### 7.8 Negative-result criterion

A negative result is a result. If the frozen rule does not reduce lost true SNP
on the test set, or reduces them only by escalating more than 2× the router's PB
budget, that is reported as NEGATIVE or ARCHITECTURAL FAILURE in §18 with the
mechanism analysed in §16, and the production recommendation stays "frozen
cascade, unchanged, with a documented failure mode".

## 8. Required metrics and controls

Per arm, per test cell and pooled: F1, precision, recall, TP, FP, FN, TN,
accuracy, balanced accuracy, ΔF1 vs PB-only with paired-bootstrap 95% CI
(`evaluate_quality_error.paired_bootstrap`, 10,000 resamples, module-default
seed — same estimator and seed as devlogs 11–14), fraction of loci sent to PB,
fraction sent by safety escalation, `lost_true_snp` vs PB-only, router/PB
disagreement rate, safety-layer capture rate (fraction of C's PB-disagreements
the layer escalates), false-escalation rate (escalated loci at which PB and the
binomial agree, so the escalation changed nothing), PB compute, total compute,
projected wall-clock and speedup vs PB-only.

Curves and stratifications: accuracy-vs-compute over the full §7.2 grid;
PB-compute-vs-coverage; depth-, VAF- and base-quality-stratified metrics; and
`lowmap_segdup` / `alldifficult` / `tandemrepeats` stratification from the
external GIAB v3.1 BEDs.

Controls, per cell, at matched *additional* PB coverage (the size of D's
escalated set), 5 seeds: **random routing**, **reverse-confidence routing**
(loci furthest from the binomial threshold), plus the existing frozen-router,
PB-only and binomial-only arms. The safety rule must show its benefit is not
merely the consequence of spending more PB.

---

*(End of pre-registration. Everything below was written after the corresponding
measurement.)*

---

## 9. Failure-regime characterisation — re-derived, not assumed

Recomputed from the cached per-locus evidence of the 12 validation cells by
`bench_v15_safety.py --characterise`
(`results/bench_v15/validation_characterisation.json`). Devlog 14 §17's summary
is *not* taken on trust, because this whole experiment's design rests on the
mechanism being ε misspecification rather than something that merely correlates
with it.

10,368,078 loci, 12,261 SNP, 7,298 routed (7.04×10⁻⁴), **41** cascade/PB
disagreements (rate 3.95×10⁻⁶). The four populations reproduce devlog 14 §17.1
exactly: 26 `router_avoids_fp`, 8 `router_fp`, **5 `router_fn`**, 2
`router_recovers_fn`.

### 9.1 The four populations, by every available variable

Values are min / median / max. `router_margin` is `|binomial_llr − 7.0|`, the
only quantity the frozen router can see. `ε-shift` is
`llr(ε=1e-4) − llr(ε=5e-2)`, the sensitivity statistic the safety layer keys on.

| variable | `router_fn` (5) | `router_avoids_fp` (26) | `router_fp` (8) | background, unrouted (10.36 M) |
|---|---|---|---|---|
| binomial LLR | −5.3 / 0.65 / 0.89 | −28.3 / −3.9 / 1.3 | 12.5 / 14.2 / 27.0 | −51.9 / −16.4 / … |
| PB LLR | 11.9 / 13.5 / 24.1 | 11.2 / 15.8 / 26.7 | −5.5 / 6.1 / 10.3 | — |
| VAF | 0.110 / 0.122 / 0.125 | 0.068 / 0.107 / 0.125 | 0.154 / 0.194 / 0.500 | 0 / 0 / 1.0 |
| mean base quality | 37.1 / 38.6 / 38.8 | 26.8 / 37.7 / 38.7 | **9.6 / 22.0 / 24.1** | — |
| mean MAPQ | 47.2 / 57.3 / 66.8 | 39.7 / 59.7 / 70.0 | 68.4 / 70.0 / 70.0 | — |
| depth | 29 / 43 / 107 | 44 / 64 / 157 | 13 / 58 / 88 | 0 / 24 / 260 |
| k | 3 / 4 / 10 | 4 / 5 / 12 | 3 / 8.5 / 11 | 0 / 0 / 93 |
| **router margin** | **6.1 / 6.4 / 12.3** | 5.7 / 10.9 / 35.3 | 5.5 / 7.2 / 20.0 | 5.4 / 23.4 / 522.8 |
| **ε-shift** | **17.7 / 23.6 / 58.5** | 22.9 / 29.3 / 70.6 | 18.6 / 51.1 / 66.7 | **−12.9 / −1.1 / 581.1** |

Cell distribution: all 5 `router_fn` in `chr1_segdup` (3 full, 2 at 30x); 24 of
26 `router_avoids_fp` likewise; the 8 `router_fp` spread over `chr20_neutral`
(3), `chr19_gcrich` (4) and `chr4_atrich` (1).

### 9.2 Which variables actually distinguish failures from ordinary loci

1. **Not the router's own statistic.** Every one of the 13 unfavourable
   disagreements (`router_fn` + `router_fp`) lies **outside** the frozen band —
   0 inside, 13 outside — with margins 5.5 to 20.0 against a cutoff of 5.412.
   This is the direct measurement behind devlog 14 §17.3's claim, and it is the
   reason widening the cutoff cannot work: the failures are on the *far* side of
   the band, so reaching the worst of them requires a cutoff near 20, which
   routes a large multiple of the loci.
2. **VAF and base quality separate the two failure directions, not failures from
   background.** `router_fn` is low-VAF/high-Q; `router_fp` is
   moderate-VAF/**low-Q** (median Q22, i.e. a true error rate ~2.5× the
   assumed 0.01). Both are ε misspecification, in opposite directions. But
   millions of ordinary loci also sit at low VAF, so neither variable alone is
   a usable rule.
3. **ε-sensitivity does separate them.** Every failure population has a strongly
   positive ε-shift (min 17.7) while the unrouted background has median −1.1.
   That is the signal, and it is exactly the quantity a fixed-ε caller cannot
   report about itself.
4. **Mapping quality does not separate them.** `router_fn` MAPQ 47–67,
   `router_fp` 68–70: overlapping, and no threshold on MAPQ divides them from
   ordinary segdup sequence.

### 9.3 Is the regime broad or localised?

**Localised in space, narrow in statistics, but not a single-locus accident.**
It is confined to segmental duplication for the *false-negative* direction (5/5)
but *not* for the false-positive direction (8/8 outside `chr1_segdup`, including
on the neutral chr20 control). So it is not "a chr1 artefact"; it is a property
of the fixed-ε caller that segmental duplication makes frequent.

### 9.4 Enrichment of the ε-flip signal

| candidate signal (unrouted loci) | loci flagged | fraction | captures `router_fn` | captures `router_fp` | enrichment for unfavourable disagreements |
|---|---:|---:|---:|---:|---:|
| flip between ε 1e-3 and 5e-2 | 392 | 3.78×10⁻⁵ | 4 / 5 | 6 / 8 | **20,300×** |
| flip between ε 1e-4 and 5e-2 | 2,703 | 2.61×10⁻⁴ | **5 / 5** | 6 / 8 | 3,250× |

H1 survives its first test: the signal exists, it is a function of `(k, n)`
alone, and it concentrates the failure population by three to four orders of
magnitude.

## 10. Validation — candidate grid and the frozen decision rule

**These numbers are discovery-set numbers and are optimistically biased (§7.3).
Nothing in the final verdict rests on them.**

All 50 candidates of the §7.2 grid, pooled over the 12 validation cells
(`results/bench_v15/rule_selection.json`). Reference points: PB-only F1
0.955726; frozen cascade F1 0.956306 (ΔF1 vs PB **+0.000581**), losing **5** true
SNP that PB-only finds; routed fraction 7.039×10⁻⁴.

Feasibility (§7.4): (a) `lost_true_snp = 0`, (b) ΔF1 vs PB ≥ 0, (c) escalated
fraction ≤ 5×10⁻⁴. **9 of 50** candidates are feasible for FROZEN-FULL, **19**
for FROZEN-LEAN.

Selected rows, ordered by budget:

| candidate | escalated | escalated fraction | PB budget × frozen | ΔF1 vs PB | lost true SNP | feasible |
|---|---:|---:|---:|---:|---:|---|
| Sq(floor=3e-3) | 16 | 1.5×10⁻⁶ | 1.00 | +0.000505 | 4 | no |
| S[3e-3, 3e-2] | 32 | 3.1×10⁻⁶ | 1.00 | +0.000624 | 4 | no |
| **Sq(floor=1e-3)** | **292** | **2.82×10⁻⁵** | **1.04** | **+0.000122** | **1** | **LEAN** |
| S[1e-3, 1e-2] | 313 | 3.02×10⁻⁵ | 1.04 | +0.000122 | 1 | lean |
| S[1e-3, 5e-2] | 392 | 3.78×10⁻⁵ | 1.05 | +0.000358 | 1 | lean |
| **SK[1e-4, 5e-2] & k ≥ 3** | **1,277** | **1.232×10⁻⁴** | **1.17** | **+0.000086** | **0** | **FULL** |
| KV(k ≥ 3, VAF∈[0.08,0.30]) | 1,636 | 1.58×10⁻⁴ | 1.22 | +0.000082 | 0 | yes |
| S[1e-4, 5e-2] (no k gate) | 2,703 | 2.61×10⁻⁴ | 1.37 | +0.000086 | 0 | yes |
| KV(k ≥ 3, VAF∈[0.05,0.35]) | 4,076 | 3.93×10⁻⁴ | 1.56 | +0.000004 | 0 | yes |
| K(k ≥ 3) — the dumb baseline | 15,625 | 1.51×10⁻³ | **3.14** | +0.000000 | 0 | no (budget) |
| K(k ≥ 2) | 37,494 | 3.62×10⁻³ | **6.14** | +0.000000 | 0 | no (budget) |

Three things worth stating from this table:

* **H2 holds on validation.** The dumb support threshold `k ≥ 3` also reaches
  zero lost SNP, but needs **12.2× more escalations** than the ε-bracket rule to
  do it (1.51×10⁻³ vs 1.232×10⁻⁴) and triples the PB budget. Sensitivity beats
  support, which is the whole argument for the layer being cheap.
* **The `k ≥ 3` conjunction is free.** `SK[1e-4,5e-2]&k≥3` and the ungated
  `S[1e-4,5e-2]` have *identical* TP/FP/FN, but the gate removes 53% of the
  escalations. Single-read flips never change an answer.
* **Every rule that repairs the failure costs precision** — the ΔF1 vs PB of the
  feasible set (+0.000004 to +0.000086) is *below* the frozen cascade's
  +0.000581. Exactly as H4 predicted, the repair hands back some of the PB false
  positives the cascade was dodging.

### 10.1 The frozen decision rule

Selected mechanically by §7.4, written to
`results/bench_v15/rule_selection.json` and
`safety_layer.FROZEN_SAFETY_RULE` **before any test cell was extracted**:

```
FROZEN-FULL   escalate  <=>  (llr(k, n; ε=1e-4) >= 7.0) != (llr(k, n; ε=5e-2) >= 7.0)
                             AND k >= 3
                             AND the frozen router did not already route the locus

FROZEN-LEAN   escalate  <=>  (llr(k, n; ε=0.01) >= 7.0) != (llr(k, n; ε=10^(-meanQ/10), floored at 1e-3) >= 7.0)
                             AND the frozen router did not already route the locus
```

`llr` is the project's own binomial likelihood ratio, asserted equal to
`BinomialVariantCaller.score_counts` to 1e-9 by `test_safety_layer.py`. The
router is untouched: both rules take its mask as an input and return a set
disjoint from it.

Freeze timestamp and the test caches' mtimes are recorded in
`results/bench_v15/independence_audit.json` (§11).

## 11. Test-set extraction and independence audit

12 cells extracted once each by `run_v15_extract.sh` (same extractor, same
flags, same resume guard as devlog 14), 15:55–16:23 on 2026-08-17.
`results/bench_v15/independence_audit.json`: **all 12 cells OK**.

| cell | loci | SNP | native mean depth | achieved 30x | achieved 15x |
|---|---:|---:|---:|---:|---:|
| `chr16_segdup` | 815,398 | 617 | 62.19 | 30.04 | 15.06 |
| `chr15_segdup` | 724,684 | 414 | 55.33 | 30.06 | 14.99 |
| `chr7_segdup` | 648,437 | 1,325 | 62.63 | 30.10 | 15.08 |
| `chr14_neutral` | 977,716 | 1,265 | 73.90 | 30.08 | 14.98 |

**9,498,705 scored loci, 10,863 SNP** across all 12 cells.

Audit results, per §6:

1. **Chromosome.** Every cell sits on chr16 / chr15 / chr7 / chr14. Intersection
   with {chr21, chr20, chr19, chr4, chr1} — every contig devlogs 3–14 used — is
   **empty**. Coordinate overlap with any earlier region is impossible.
2. **Coordinate.** Each cell's stored `region` equals its pre-registered span;
   every reconstructed locus coordinate lies inside it; the locus count matches
   an independent re-tiling of the windows; all cells are 64-locus aligned; no
   non-finite binomial or PB score anywhere.
3. **Sample.** HG002 in both splits — a declared *dependence*, not an
   independence claim. This experiment tests transfer across chromosome and
   structure only.
4. **Provenance.** Same BAM URL, truth VCF release, HC BED, reference release
   and stratification release as the validation split, so the only difference
   between the splits is which sequence.
5. **Selection blindness.** `select_v15_regions.py` imports only `gzip`, `json`,
   `numpy` and `pysam` — no caller, router, safety or benchmark module.
6. **Rule freezing.** `rule_selection.json` records `frozen_at
   2026-08-17T15:52:06+05:00`; every test cache was written between 16:04 and
   16:22, i.e. **after** the rule was frozen, and the audit asserts this per cell
   (`extracted_after_rule_was_frozen: true` × 12). One caveat is recorded in §16:
   a `git stash -u`/`pop` at 16:01 reset the *mtimes* of the results JSONs, so
   the mtime ordering holds but the in-file `frozen_at` field is the
   authoritative timestamp.

## 12. Test results — the frozen safety layer on untouched data

Full tables in `results/bench_v15/summary.md`; raw numbers in
`results/bench_v15/test_results.json`.

### 12.1 All five arms, pooled test set (9,498,705 loci, 10,863 SNP)

| arm | F1 | precision | recall | TP | FP | FN | accuracy | balanced acc. | ΔF1 vs PB | lost true SNP | PB fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **A** binomial-only | 0.923862 | 0.915184 | 0.932707 | 10,132 | 939 | 731 | 0.99982419 | 0.966304 | −0.001217 | 67 | 0 |
| **B** PB-only | 0.925079 | 0.911683 | 0.938875 | 10,199 | 988 | 664 | 0.99982608 | 0.969385 | 0 | 0 | 1.0 |
| **C** frozen router → PB | **0.927766** | 0.917004 | 0.938783 | 10,198 | **923** | 665 | 0.99983282 | 0.969343 | **+0.002686** | **1** | 5.809×10⁻⁴ |
| **D** safety FROZEN-FULL | 0.925079 | 0.911683 | 0.938875 | 10,199 | 988 | 664 | 0.99982608 | 0.969385 | **+0.000000** | **0** | 6.879×10⁻⁴ |
| **D-lean** safety FROZEN-LEAN | 0.925955 | 0.913472 | 0.938783 | 10,198 | 966 | 665 | 0.99982829 | 0.969341 | +0.000875 | 1 | 6.044×10⁻⁴ |

Paired-bootstrap 95% CIs on ΔF1 vs PB-only: **C** [+0.001968, +0.003414] (excludes
0 — the frozen cascade genuinely beats the arm it approximates); **D**
[−0.000209, +0.000209]; **D-lean** [+0.000396, +0.001388].

The single most important row-pair: **arm D reproduces PB-only exactly** —
identical TP, FP, FN on all 9.5 M loci — while arm C beats PB by +0.002686.

### 12.2 Per-cell

| cell | depth | SNP | F1 B | F1 C | F1 D | ΔF1(D−PB) | ΔF1(D−C) | lost C | lost D | routed | escalated | PB × C |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `chr16_segdup` | full | 617 | 0.894955 | 0.896266 | 0.895695 | +0.000741 | −0.000570 | **1** | 0 | 0.0283% | 278 | 2.20 |
| `chr16_segdup` | 30x | 617 | 0.892410 | 0.893901 | 0.892410 | 0 | −0.001491 | 0 | 0 | 0.0536% | 100 | 1.23 |
| `chr16_segdup` | 15x | 617 | 0.891156 | 0.891915 | 0.891156 | 0 | −0.000758 | 0 | 0 | 0.1533% | 1 | 1.00 |
| `chr15_segdup` | full | 414 | 0.779859 | 0.796651 | 0.781690 | +0.001831 | **−0.014961** | 0 | 0 | 0.0302% | 137 | 1.63 |
| `chr15_segdup` | 30x | 414 | 0.783751 | 0.791315 | 0.783751 | 0 | −0.007563 | 0 | 0 | 0.0386% | 84 | 1.30 |
| `chr15_segdup` | 15x | 414 | 0.779577 | 0.781523 | 0.779577 | 0 | −0.001947 | 0 | 0 | 0.1025% | 8 | 1.01 |
| `chr7_segdup` | full | 1,325 | 0.917730 | 0.920014 | 0.916755 | −0.000975 | −0.003259 | 0 | 0 | 0.0210% | 130 | 1.96 |
| `chr7_segdup` | 30x | 1,325 | 0.918534 | 0.924454 | 0.918534 | 0 | −0.005920 | 0 | 0 | 0.0335% | 86 | 1.40 |
| `chr7_segdup` | 15x | 1,325 | 0.919942 | 0.922291 | 0.919942 | 0 | −0.002349 | 0 | 0 | 0.1261% | 28 | 1.03 |
| `chr14_neutral` | full | 1,265 | 0.998028 | 0.997634 | 0.998028 | 0 | **+0.000394** | 0 | 0 | 0.0020% | 119 | 6.95 |
| `chr14_neutral` | 30x | 1,265 | 0.996848 | 0.996848 | 0.996848 | 0 | 0 | 0 | 0 | 0.0133% | 44 | 1.34 |
| `chr14_neutral` | 15x | 1,265 | 0.988889 | 0.988889 | 0.988889 | 0 | 0 | 0 | 0 | 0.1061% | 1 | 1.00 |

`chr14_neutral_full` is the one cell where the safety layer *helps* relative to
the cascade (+0.000394): there the binomial was confidently wrong in the
over-calling direction, and escalating fixed it. Every segdup cell moves the
other way.

Note the `PB × C` column on the two lowest-routing cells: `chr14_neutral_full`
routes only 0.0020% of loci, so 119 escalations multiply its PB budget by 6.95×
— but 6.95 × 0.0020% is still 0.0139%. Ratios are misleading at these
magnitudes, which is why §7.5's A3 is evaluated on the pooled budget and on wall
clock rather than per cell.

### 12.3 The failure mode on independent data

**The frozen cascade lost exactly 1 true SNP in 9,498,705 loci** (in
`chr16_segdup` at native depth), against 5 in 10,368,078 validation loci. The
safety layer recovered it: arm D loses **0**. But one event is one event, and
§7.6 caps the verdict accordingly.

Cascade/PB disagreements on test: **80** (rate 8.42×10⁻⁶, vs 3.95×10⁻⁶ on
validation), of which 74 lie outside the router band. Their composition is what
decides everything below.

## 13. What the safety layer actually did

Pooled, FROZEN-FULL (`results/bench_v15/test_results.json` → `pools.all.safety`):

| quantity | value |
|---|---:|
| escalated loci | 1,016 |
| escalated fraction of all loci | 1.070×10⁻⁴ |
| total PB fraction (router ∪ escalation) | 6.879×10⁻⁴ |
| PB budget vs frozen cascade | **1.18×** |
| unrouted cascade/PB disagreements | 80 |
| of those, escalated (**capture rate**) | 74 (**92.5%**) |
| cascade-lost true SNP captured | 1 / 1 |
| false escalations (PB and binomial agreed anyway) | 942 (**92.7%**) |
| calls changed vs frozen cascade | 74 |
| … of which **fixed** a wrong call | **5** |
| … of which **broke** a right call | **69** |

**The detector works and the intervention does not.** The layer finds 92.5% of
the loci where the cascade departs from PB, at 1.18× the PB budget and ~0
compute — H1 and H2 confirmed on independent data. But of the 74 calls it
changed, 69 were calls the cascade had *right* and PB had *wrong*. Its net
effect is to convert the cascade into PB-only, which is exactly what the pooled
F1 shows.

This is H4, and it is stronger on test than on validation: the favourable /
unfavourable ratio among the cascade's departures from PB was 26:5 on validation
and **69:5** here.

FROZEN-LEAN is worse in the way its budget predicts: 223 escalations, 53.8%
capture, **0** of the cascade's lost true SNP captured, 0 calls fixed, 43 broken.
It pays a small accuracy cost for no repair at all.

## 14. Compute — the layer is genuinely almost free

Measured on this experiment's own data, single quiet process per cell:

| stage | throughput | relative |
|---|---:|---:|
| binomial caller | 2,971,220 loci/s | 1× |
| **safety layer (production, `k`-gated)** | **122,978,607 loci/s** | **0.02×** the binomial caller |
| safety layer (naive dense form — *not used*) | 2,058,718 loci/s | 1.44× the binomial caller |
| PB (repaired) | 6,041 loci/s | 492× the binomial caller per locus |

The gated form evaluates the two `logpmf` passes only at loci with `k ≥ k_min`
(a locus with `k = 0` cannot produce a positive call at any ε, so it cannot
flip), which is ~0.1% of loci. The mask is bit-identical to the dense
computation — asserted by `test_safety_layer.py::test_sparse_escalation_matches_dense`
— so this is an evaluation-order optimisation, not a different rule. **The naive
implementation would have cost more than the entire binomial caller; the gated
one costs 2% of it.**

Projected pooled wall clock:

| accounting | PB-only | frozen cascade | safety FULL | safety LEAN |
|---|---:|---:|---:|---:|
| **caller-only** (devlogs 11–14 convention) | 1,045.9 s (1×) | 4.1 s (**255.9×**) | 4.3 s (**243.0×**, 0.949× of C) | 4.2 s (248.8×, 0.972× of C) |
| **end-to-end** (incl. pileup construction) | 5,257.1 s (1×) | 1,316.1 s (**4.0×**) | 1,316.7 s (4.0×, **1.000×** of C) | 1,316.3 s (4.0×, 1.000× of C) |

Two things must be said about this table rather than left implicit:

1. **The layer costs essentially nothing.** 0.2 s over 9.5 M loci caller-only;
   0.6 s end-to-end. A3 passes under *both* accountings, which is why the
   ambiguity in §7.5's wording ("end-to-end") did not have to be adjudicated:
   0.949× and 1.000× are both ≥ 0.8.
2. **The project's headline speedup is a caller-only number.** End-to-end, with
   the measured pileup construction both arms pay, the frozen cascade's advantage
   over PB-only is **4.0×**, not 256×. That is not a finding about the safety
   layer — it applies identically to devlogs 11–14 — but it is measured here for
   the first time and it belongs in the record. The extraction timings feeding it
   come from four concurrent extractions and are upper bounds, which makes 4.0× a
   conservative floor.

## 15. Controls — the effect is structure, not budget

Matched *extra* PB coverage: the frozen router's routed set is held fixed and
only the additional 1,016 loci differ. 5 seeds for random.

| arm | F1 | ΔF1 vs PB | TP | FP | FN | lost true SNP |
|---|---:|---:|---:|---:|---:|---:|
| frozen cascade (no extra PB) | 0.927766 | +0.002686 | 10,198 | 923 | 665 | 1 |
| **safety layer** | 0.925079 | +0.000000 | 10,199 | 988 | 664 | **0** |
| reverse-confidence escalation | 0.927766 | +0.002686 | 10,198 | 923 | 665 | 1 |
| random escalation, seeds 0–4 | 0.927766 | +0.002686 | 10,198 | 923 | 665 | 1 |

Random escalation std across 5 seeds: **0.000000**. Both controls changed
**not one call** at the same budget. Spending 1,016 loci of extra PB compute at
random, or on the loci the binomial is *most* confident about, is exactly
equivalent to not spending it. The safety layer at the identical budget changed
74 calls and recovered the lost variant.

So A4 passes decisively, and the interpretation is unambiguous: everything the
safety layer does — the repair *and* the damage — comes from its structure, not
from the additional compute.

## 16. Engineering, QC and integrity

* **Test suite.** 190 tests in the 8 directly relevant suites passed before any
  work began; **394 pass after** (`pytest -q`, excluding `test_bimamba.py` and
  `test_providers.py`, which fail to *collect* at `HEAD` too — they import a
  `bimamba_variant_caller` package that does not exist in this flat layout. That
  was verified against a clean `HEAD` checkout, so it is pre-existing and
  unrelated).
* **New tests.** `test_safety_layer.py` (34) and `test_bench_v15.py` (80),
  covering: the frozen constants are imported and never assigned; the frozen rule
  matches the on-disk selection record; the layer's LLR equals
  `BinomialVariantCaller` to 1e-9 at five ε values; `quality_epsilon` equals the
  v2 caller's `_epsilon`; escalation is disjoint from routing for **every** rule
  in the grid; the router's mask is never mutated; raising `k_min` only shrinks
  and widening the bracket only grows the escalated set; the sparse and dense
  evaluations are bit-identical; `k = 0` can never escalate at any ε; every A1–A4
  criterion and every §7.6 verdict branch; and signature/source introspection
  proving no truth label, coordinate, PB score or difficulty annotation can reach
  a routing decision.
* **No existing behaviour changed.** `git status --porcelain` shows **zero
  modified tracked files** — the session is purely additive.
  `cascade.py`, `cheap_router.py`, `binomial_baseline.py`,
  `quality_error_model.py`, `robustness_benchmark.py`, `bench_v12_stage1.py`,
  `bench_v13_robustness.py`, `bench_v14_crosschrom.py` and
  `extract_bench_v12.py` are all byte-identical to `HEAD` (`3837b61`).
* **No previous artefact or devlog touched.** `results/bench_v14/*` and
  `History/14_DEVLOG.md` retain their 2026-08-17 10:29–10:41 mtimes.
* **JSON validity.** All six `results/bench_v15/*.json` parse.
* **No orphaned processes.** All extraction, fetch and scoring processes exited;
  the two shell wait-loops were terminated explicitly and nothing remains.
* **Disclosed incident.** At 16:01 a `git stash -u` / `git stash pop` was run to
  verify that the two collection errors above pre-exist at `HEAD`. This was a
  careless thing to do while an extraction was running. Consequences, checked:
  the BAMs and `.npz` caches are `.gitignore`d (`data/`, `*.npz`) so they were
  never stashed; the running extractions were unaffected and all 12 completed;
  the empty `cache/bench_v15/` directory was removed by the stash and recreated
  (the extractor `mkdir -p`s it anyway); and the `results/bench_v15/*.json`
  mtimes were reset to 16:01:23. The freeze-ordering claim in §11 therefore rests
  on the in-file `frozen_at` field (15:52:06) as well as on the reset mtimes
  (16:01:23), both of which precede every test cache (16:04–16:22). No result
  was altered. Recorded because an undisclosed timestamp perturbation in a
  pre-registered experiment would be indistinguishable from tampering.

## 17. Stratified analysis

`results/bench_v15/stratified.json`, produced by `analyse_v15_strata.py` (post-hoc
reporting on the frozen arm; it cannot change a decision). Pooled over all 12
test cells. "changed" = calls arm D makes differently from arm C.

### 17.1 By VAF

| VAF bin | loci | SNP | escalated | F1 B | F1 C | F1 D | lost C | lost D | changed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0–0.15 | 9,483,072 | 577 | **894** | 0.016416 | 0.015198 | 0.016484 | **1** | 0 | **70** |
| 0.15–0.25 | 3,600 | 113 | 113 | 0.290909 | 0.287582 | 0.289474 | 0 | 0 | 4 |
| 0.25–0.35 | 1,753 | 504 | 9 | 0.788177 | 0.788177 | 0.788177 | 0 | 0 | 0 |
| 0.35–0.45 | 1,572 | 1,456 | 0 | 0.964558 | 0.964558 | 0.964558 | 0 | 0 | 0 |
| 0.45–0.55 | 2,581 | 2,299 | 0 | 0.986225 | 0.986225 | 0.986225 | 0 | 0 | 0 |
| 0.55–0.75 | 1,570 | 1,542 | 0 | 0.992278 | 0.992278 | 0.992278 | 0 | 0 | 0 |
| ≥0.75 | 4,557 | 4,372 | 0 | 0.994153 | 0.994153 | 0.994153 | 0 | 0 | 0 |

**Above VAF 0.35 the layer never fires.** 88% of escalations and 95% of changed
calls are at VAF < 0.15 — the regime devlog 14 predicted, replicated exactly.

### 17.2 By mean base quality

| meanQ bin | loci | SNP | escalated | F1 B | F1 C | F1 D | lost C | lost D | changed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0–20 | 278,780 | 311 | 5 | 0.000000 | 0.000000 | 0.000000 | 0 | 0 | 0 |
| 20–25 | 13,473 | 22 | 145 | 0.476190 | 0.465116 | 0.476190 | 0 | 0 | 2 |
| 25–30 | 130,403 | 162 | 348 | 0.827795 | 0.827795 | 0.827795 | 0 | 0 | 0 |
| 30–35 | 1,500,126 | 1,661 | 237 | 0.947653 | 0.949081 | 0.947653 | 0 | 0 | 5 |
| 35–60 | 7,575,923 | 8,707 | 281 | 0.941458 | 0.944655 | 0.941458 | **1** | 0 | **67** |

The escalation *rate* is far higher at Q20–30 (1.1% and 0.27% of those bins) than
at Q35–60 (0.004%), which is the ε-mismatch signal behaving exactly as designed —
the assumed ε = 0.01 is worst calibrated there. But the *consequential*
escalations are in the Q35–60 bin, which is where the low-VAF segdup corner
lives.

### 17.3 By depth

| depth bin | loci | SNP | escalated | F1 B | F1 C | F1 D | lost C | lost D | changed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1–5 | 157,935 | 194 | 0 | 0.496403 | 0.496403 | 0.496403 | 0 | 0 | 0 |
| 5–10 | 365,354 | 484 | 0 | 0.881907 | 0.881907 | 0.881907 | 0 | 0 | 0 |
| 10–15 | 924,610 | 1,182 | 0 | 0.948696 | 0.948696 | 0.948696 | 0 | 0 | 0 |
| 15–20 | 1,148,464 | 1,405 | 8 | 0.965420 | 0.965420 | 0.965420 | 0 | 0 | 0 |
| 20–30 | 1,719,218 | 1,995 | 85 | 0.959552 | 0.960488 | 0.959552 | 0 | 0 | 4 |
| ≥30 | 4,908,670 | 5,300 | **923** | 0.936420 | 0.941440 | 0.936420 | **1** | 0 | **70** |

**Below 15x the layer is inert** (8 escalations in 2.6 M loci, 0 changed calls).
The regime is a high-depth phenomenon: only at ≥30x is a handful of low-VAF
non-reference reads enough evidence to move either caller.

By depth realisation (pooled): full-depth 664 escalations / 86% capture / 5 fixed
/ 31 broken; 30x 314 / 100% / 0 fixed / 28 broken; 15x 38 / 100% / 0 fixed / 10
broken.

### 17.4 By external difficulty annotation (reporting only — never a rule input)

| stratum | loci | SNP | escalated | lost C | lost D |
|---|---:|---:|---:|---:|---:|
| in `lowmap_segdup` | 4,725,945 | 5,055 | 525 | 1 | 0 |
| not in `lowmap_segdup` | 4,772,760 | 5,808 | 491 | 0 | 0 |
| in `alldifficult` | 5,317,956 | 5,721 | 704 | 1 | 0 |
| not in `alldifficult` | 4,180,749 | 5,142 | 312 | 0 | 0 |
| in `tandemrepeats` | 538,581 | 663 | 258 | 258 → n/a | 0 |
| not in `tandemrepeats` | 8,960,124 | 10,200 | 758 | 1 | 0 |

Escalations are enriched 2.3× per locus inside `alldifficult` and 8.0× inside
`tandemrepeats`, but the layer fires in ordinary sequence too (491 escalations
outside segdup) — it keys on the statistic, not on the annotation, which is what
§7.1 required of it.

## 18. Failure analysis — why a working detector is not a working repair

The layer met every engineering goal it was set:

* it detects the regime — **92.5%** of the cascade's departures from PB, and the
  one true SNP the cascade lost;
* it is cheap — **1.18×** the PB budget, **0.02×** the binomial caller's compute,
  no measurable change in end-to-end wall clock;
* its effect is structural, not budgetary — matched-budget random and
  reverse-confidence escalation changed **zero** calls.

And it made the pipeline **worse**: ΔF1(D − C) = **−0.002686** pooled, against a
0.001 equivalence bound. The reason is a single fact this experiment measured and
the previous one only hinted at:

> In segmental duplication, **PB is not a better caller than the binomial — it is
> a worse one**, and the frozen cascade's "errors" relative to PB are mostly it
> being right.

The evidence:

1. Arm C beats arm B by +0.002686 [+0.001968, +0.003414] on the pooled test set —
   statistically significant and 2.7× the equivalence bound. The cascade beats the
   reference it approximates.
2. Of the 80 cascade/PB disagreements, the cascade is right **69** times and wrong
   **5** (6 are indeterminate under the frame). The favourable ratio was 26:5 on
   validation and is 69:5 here.
3. The layer captured 74 of the 80 and, by escalating them, adopted PB's answer at
   all of them: 5 fixes, 69 breakages, net −64 correct calls.
4. Precision is where it lands: 0.917004 → 0.911683, exactly PB-only's. Recall is
   essentially unchanged (0.938783 → 0.938875, i.e. the one recovered SNP).

So the failure is **not** insufficient signal — the signal is excellent — and
**not** insufficient data at the detection level. It is architectural, but not in
the way §7.6's category E means it (compute): the layer's *action* is wrong. It
routes to an oracle that is not an oracle in the regime it routes in. Escalating
"the binomial may be confidently wrong" to PB assumes PB is right there; in
segdup, PB's own failure mode — trusting per-read base quality with no term for
mapping error, so a few paralogous high-quality reads look like overwhelming
evidence (devlog 13 §18.1) — fires on exactly the same loci.

Two secondary limitations of the evidence, stated plainly:

* **Power for A1 is very weak.** The test set contained **one** cascade failure.
  The rule recovered it, and the controls did not, which is real but is one event.
  §7.6 caps the verdict at WEAK POSITIVE for precisely this reason, and that cap
  is doing real work here rather than being a formality.
* **The damage estimate is much better powered than the benefit estimate.** 69
  broken calls versus 5 fixed is not a marginal count, and it replicates the
  validation direction. The conclusion "this layer costs more accuracy than it
  buys" is the well-supported half of this experiment.

FROZEN-LEAN adds nothing to the analysis: at 53.8% capture it caught none of the
cascade's lost SNP, fixed 0 calls and broke 43. It is strictly dominated.

## 19. Limitations

1. **One cascade failure in the test set.** The primary endpoint A1 rests on a
   single event. Everything about the *repair* is therefore weakly powered;
   everything about the *cost* is not.
2. **The test set is deliberately enriched** for segmental duplication (3 of 4
   regions, 6.57 M of 9.50 M loci). This was necessary — an unenriched test set
   would have contained no failures at all — but it means the pooled numbers are
   not a genome-wide expectation. The `chr14_neutral` arm is the only
   ordinary-sequence estimate here, and there the layer is neutral-to-slightly-
   positive.
3. **One sample, one platform.** HG002, Illumina 2×250, novoalign, GRCh38.
   Nothing here supports a claim about other samples, aligners or technologies.
4. **Four 1 Mb windows, 10,863 SNP.** Wide in structure, narrow in quantity.
5. **Depth realisations are downsamples of one library**, not independent runs.
6. **GIAB truth is least reliable exactly where this experiment lives.** In
   segdup, F1 ≈ 0.78–0.92 for *every* arm; some of what is scored as a false
   positive is truth-set error. The claim "the cascade is right and PB is wrong 69
   times" inherits that uncertainty, though it would take a very unusual truth-set
   bias to reverse a 69:5 ratio.
7. **SNP only.** Indels, MNPs and SVs are outside the label set, and they matter
   more in segdup than elsewhere.
8. **PB above 48x still sees at most `MAX_READS = 48` reads** (devlog 13 §18.4).
9. **The rule family is one family.** ε-sensitivity, support and VAF-band rules
   were searched; many cheap statistics were not (strand balance, per-read MAPQ
   distribution, local haplotype context). A different cheap signal might select a
   *different* subset of the disagreements — specifically, the 5 where PB is right
   rather than the 69 where it is wrong. This experiment shows that
   ε-sensitivity does not separate those two populations; it does not show that
   nothing can.
10. **The 170–400× figure is caller-only.** Measured end-to-end with pileup
    construction, the frozen cascade's advantage over PB-only on this data is 4.0×.
    This applies to devlogs 11–14 equally and is not a criticism of the safety
    layer.

## 20. Final verdict

Applying §7.6 mechanically to the pooled test set (`results/bench_v15/summary.md`,
computed by `bench_v15_safety.acceptance` / `verdict_from`, not by hand):

**FROZEN-FULL** — `SK[1e-4, 5e-2] & k ≥ 3`

| criterion | result | evidence |
|---|---|---|
| A1 failure repaired | **PASS** | 1 → 0 lost true SNP |
| A2 non-inferior to PB | **PASS** | ΔF1 95% CI [−0.000209, +0.000209], bound −0.001 |
| A3 compute preserved | **PASS** | PB budget 1.18× (≤2); end-to-end 1.000× of cascade, caller-only 0.949× (both ≥0.8) |
| A4 beats matched-budget controls | **PASS** | random best 1 lost, reverse 1 lost, safety **0** lost |

All four pass; the test set contained only 1 cascade failure (< 3), so §7.6 caps
the verdict:

> ## **B. WEAK POSITIVE**

**FROZEN-LEAN** — `Sq(floor=1e-3)`: A1 **FAIL** (1 → 1), A4 **FAIL**, A2/A3 pass →
> ## **D. NEGATIVE**

### 20.1 What that verdict does and does not say

The central question was:

> "Can we cheaply detect the specific regime where the frozen binomial→router
> architecture is confidently wrong, and selectively send those loci to PB
> without sacrificing the ~170–400× compute advantage?"

Answering it in three parts, because the parts have different answers:

* **Can we detect it cheaply? — Yes, convincingly.** An ε-sensitivity bracket over
  `(k, n)` captures 92.5% of the cascade's departures from PB and 100% of its lost
  true variants, at 1.07×10⁻⁴ of loci escalated, 3,250× enrichment, and 2% of the
  binomial caller's compute. This is the solid, well-powered result of the
  experiment, and it replicated from validation to independent chromosomes.
* **Without sacrificing the compute advantage? — Yes.** 1.18× the PB budget,
  0.949× the caller-only speedup, 1.000× end-to-end. Category E does not apply.
* **Should we send those loci to PB? — No.** Doing so converts the cascade into
  PB-only and costs −0.002686 F1 against a production baseline that beats PB by a
  statistically significant margin. 69 correct calls broken to fix 5.

So the honest one-line summary is not "the safety layer works" and not "the
safety layer fails" but: **the detector works, the escalation target is wrong.**

### 20.2 Production recommendation

**Do not adopt the safety layer.** Keep the production candidate exactly as it is:

```
binomial (ε = 0.01) → |binomial_LLR − 7.0| ≤ 5.411872376933351 → PB → call
```

with its documented failure mode: ~1 true SNP lost per 10⁷ loci, concentrated at
VAF < 0.15, base quality ≥ Q35, depth ≥ 30x, inside segmental duplication. On
this evidence that failure mode is **cheaper to accept than to repair by
escalation to PB**.

The frozen rule is retained in `safety_layer.py` as a *detector*, not as a
router: it is the only cheap way found so far to enumerate the loci where the two
callers disagree, which is exactly the candidate pool a future third opinion
would need.

## 21. Recommended next experiment

This experiment's result reframes the problem, so the obvious next step changed:

The bottleneck is no longer *finding* the disputed loci — that is solved, at
1.07×10⁻⁴ of the genome and negligible cost. The bottleneck is that **neither
existing caller is trustworthy on them**: the binomial wins 69 to 5 against PB
there, but both are far from the truth set (F1 ≈ 0.78–0.92 in these regions).

The right next question is therefore: *given the ~10⁻⁴ locus pool this detector
produces, can a third opinion beat both callers on it?* That pool is small enough
that an expensive model — including the neural arm this project has been building
toward — is affordable on it, which is the first time that has been true. The
key design constraint from this experiment is that such a model must be evaluated
against **truth**, and against **both** the binomial and PB, never against PB as a
reference, because PB is the weaker of the two callers exactly where the pool
lives.

A necessary preliminary, cheaper than any of that: establish whether the 69:5
ratio is real or a truth-set artefact, by checking those loci against an
independent truth source or a mapping-quality-aware caller. If it is an artefact,
this devlog's verdict flips, and the layer becomes worth adopting.

## 22. Reproducibility

Everything below was run on 2026-08-17 from a clean `HEAD` = `3837b61`.

```bash
# 1. contig selection (annotation only; no FASTA needed)
python select_v15_regions.py --contigs-only     # -> results/bench_v15/contig_selection.json

# 2. reference FASTAs for the chosen contigs (Ensembl release-110 GRCh38)
./fetch_v15_reference.sh                        # -> data/reference/chr{16,15,7,14}_full.fa

# 3. window selection (Rules H and O)
python select_v15_regions.py                    # -> results/bench_v15/region_selection.json

# 4. failure-regime characterisation, validation split only
python bench_v15_safety.py --characterise       # -> results/bench_v15/validation_characterisation.json

# 5. FREEZE THE RULE, validation split only, before any test data exists
python bench_v15_safety.py --select             # -> results/bench_v15/rule_selection.json

# 6. test data
./fetch_v15_data.sh                             # -> data/giab_hg002_v15/*_full.bam, *.vcf.gz, *_highconf.bed
./prepare_v15_depths.sh                         # -> *_30x.bam, *_15x.bam  (samtools view -s, seed 42)
./run_v15_extract.sh                            # -> cache/bench_v15/*.npz  (12 cells, once each)

# 7. score the frozen rule on the untouched test set
python bench_v15_safety.py --test               # -> results/bench_v15/{test_results,independence_audit}.json
python summarize_v15.py                         # -> results/bench_v15/summary.md
python analyse_v15_strata.py                    # -> results/bench_v15/stratified.json

# 8. tests
python -m pytest -q --ignore=test_bimamba.py --ignore=test_providers.py   # 394 passed
```

**Environment.** Python 3.14, numpy, scipy, pysam; samtools/bcftools for slicing;
Linux 6.19.12, 8 cores, 31 GB RAM. All randomness is seeded: bootstrap seed and
resample count are `evaluate_quality_error`'s module defaults (10,000 resamples),
control seeds are 0–4, `samtools view -s` seed is 42, the synthetic-fixture seeds
in the tests are fixed literals.

**Artefacts.**

| file | contents |
|---|---|
| `results/bench_v15/contig_selection.json` | Rule C ranking of all 18 candidate autosomes |
| `results/bench_v15/region_selection.json` | Rules H/O window statistics and the 4 chosen windows |
| `results/bench_v15/validation_characterisation.json` | §9 failure-regime re-derivation |
| `results/bench_v15/rule_selection.json` | all 50 candidates, the §7.4 criterion, the frozen rules, `frozen_at` |
| `results/bench_v15/independence_audit.json` | §11 per-cell audit, freeze-vs-extraction ordering |
| `results/bench_v15/test_results.json` | every arm, cell, pool, control, curve, stratum, timing |
| `results/bench_v15/stratified.json` | §17 VAF / quality / depth strata for the safety arm |
| `results/bench_v15/summary.md` | the rendered tables quoted in §12–§16 |
| `cache/bench_v15/*.npz` | the 12 extracted test cells (12 × ~12 MB) |

**Code added** (nothing existing modified): `safety_layer.py`,
`bench_v15_safety.py`, `select_v15_regions.py`, `summarize_v15.py`,
`analyse_v15_strata.py`, `fetch_v15_reference.sh`, `fetch_v15_data.sh`,
`prepare_v15_depths.sh`, `run_v15_extract.sh`, `test_safety_layer.py`,
`test_bench_v15.py`.
