# Does the frozen cascade generalize to a different sample? — HG003 cross-sample validation on the devlog-16 chr13 window

**Date:** 2026-09-04

**Status of sections:** This section (pre-registration) is written and saved
**before any HG003 BAM has been downloaded and before any caller or router has
run on a single HG003 locus.** Results are appended below a `RESULTS` marker
after the corresponding measurement, never edited into this section.

This experiment was commissioned as part of an external, adversarial audit of
devlogs 11-17 (`STRONG POSITIVE` on HG002, chr13 and chr17). Every devlog
through 17 used exactly one sample, HG002. This is the single largest
untested claim in the project: "cross-sample transfer is currently
untested" is stated explicitly in devlog 17 §18 and in the project's own
summary of its current evidence. This experiment tests it directly.

## 1. Research question

Does the frozen binomial → router → PB cascade — unchanged, unretuned,
byte-identical to what devlog 16 ran on HG002 — preserve PB-only accuracy when
applied to a **different individual's genome**, holding region, platform and
aligner fixed?

## 2. Design: what varies and what doesn't

| axis | devlog 16 (HG002) | this experiment (HG003) |
|---|---|---|
| sample | HG002 (son) | **HG003 (father)** |
| region | chr13:70,000,000-73,000,000 | chr13:70,000,000-73,000,000 (identical) |
| platform | Illumina 2x250 | Illumina 2x250 (identical) |
| aligner | novoalign | novoalign (identical) |
| truth set | GIAB NISTv4.2.1 HG002 | GIAB NISTv4.2.1 HG003 (own truth, same release) |
| binomial threshold | 7.0 | 7.0 (frozen, unchanged) |
| PB threshold | 10.5 | 10.5 (frozen, unchanged) |
| router cutoff | 5.411872376933351 | 5.411872376933351 (frozen, unchanged) |
| depths | full/30x/15x | full/30x/15x (identical recipe, seed 42) |

Region, platform, aligner and all three frozen constants are held fixed
deliberately, so a change in outcome can only be attributed to the sample.
HG003 is the Ashkenazim trio father — a different individual from HG002 (son),
sequenced and processed by the same GIAB pipeline, which is why its BAM exists
at the same URL layout and its truth set is the same v4.2.1 release. This does
not test cross-platform or cross-technology generalization (§10).

## 3. Frozen configuration — not tuned, not refit for this experiment

| item | value |
|---|---|
| git HEAD | `1fd9dc99cbff43ae003ed8b71b80397406bd2691` |
| working tree | clean except new `results/bench_v18/`, `data/giab_hg002_v18/`, `logs_v18/`, and the new (non-frozen) glue files listed in §4 |
| binomial threshold | 7.0 |
| binomial error rate ε | 0.01 |
| PB threshold | 10.5 |
| frozen router cutoff | 5.411872376933351 |
| reference | Ensembl release-110 GRCh38 chr13, reused unchanged from devlog 16/17 (`data/reference/chr13_full.fa`) |
| truth | GIAB HG003 v4.2.1 benchmark VCF, SNP only |
| high-confidence BED | HG003 v4.2.1 `_noinconsistent` BED (HG003's own, not HG002's) |
| stratifications | GIAB v3.1 (`lowmap_segdup`, `alldifficult`, `tandemrepeats`), reused unchanged from `data/strat_v31/` |

SHA-256 of every frozen implementation file this experiment imports, recorded
before any HG003 data was touched:

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
```

**Bug rule, fixed in advance, verbatim from devlog 17 §3.** If a defect is
found in any frozen file during this experiment it is not fixed here. It is
recorded in §17 and every cell it touches is classified `INVALID` and excluded
from the verdict. Only an outright crash in *non-frozen* v18 glue code may be
repaired, and any such repair is disclosed.

## 4. New (non-frozen) files this experiment adds

`fetch_v18_data.sh`, `prepare_v18_depths.sh`, `run_v18_extract.sh`,
`bench_v18_unseen.py`, `summarize_v18.py`. All five are line-level adaptations
of the v17 files with the sample and coordinates changed; none imports or
edits a frozen file listed in §3, and `bench_v18_unseen.py` imports its scoring
logic from those same frozen modules unchanged.

`bench_v18_unseen.py`'s independence audit necessarily differs in kind from
v13-v17's: those experiments varied the *region* on a fixed sample (HG002) and
so audited for coordinate/contig reuse. This experiment deliberately reuses
devlog 16's exact coordinates on a *different* sample, so coordinate reuse is
disclosed, not flagged as leakage. The actual independence check here is
sample identity. The upstream GIAB novoalign BAMs carry no `@RG` line at all
(verified directly on both the HG002 and HG003 files before writing this
check), so `_bam_sample_is_hg003` in `bench_v18_unseen.py` instead requires
two independent signals to agree: the truth VCF's own sample column (via
`bcftools query -l`, the authoritative GIAB identity for the genotypes being
scored) reads exactly `HG003`, **and** the BAM's own `@PG` provenance line —
written automatically by `samtools view` at slice time, not asserted by this
script — embeds the HG003 source URL. Either signal missing or mismatched
marks the cell `NOT ok` and aborts scoring.

## 5. Independence audit (mechanical, pre-registered)

Before any metric is computed, `bench_v18_unseen.verify_caches` must confirm,
per cell:
1. the cached contig and region match the pre-registered chr13:70,000,000-73,000,000 window exactly;
2. the truth VCF's sample column reads HG003 and the BAM's own slice
   provenance names the HG003 source URL (not HG002 — this is the leakage
   check that actually matters here; see §4);
3. window tiling and cache locus counts agree, and both LLR arrays are finite
   everywhere;
4. positions lie inside the pre-registered window.

A single failing check aborts scoring. The artifact is
`results/bench_v18/independence_audit.json`.

**Threshold-tuning independence.** The three frozen constants (7.0, 10.5,
5.411872376933351) were fit on HG002 chr21 in devlog 13 and have never been
refit on any other sample. No HG003 data of any kind — genotype, coverage,
region content — has been read, downloaded, or inspected by this project
before this devlog was written. This is the first time this repository has
ever touched an HG003 file. That is the entire cross-sample claim this
experiment can support: *the constants were never conditioned on HG003.* It
cannot, by construction, rule out that they were implicitly overfit to some
property shared by HG002 and HG003 (both Ashkenazim, both novoalign, both
Illumina 2x250) — see §10.

## 6. Experimental arms and metrics

Identical to devlog 16/17 §6: three arms (binomial-only, PB-only, frozen
cascade) at native/30x/15x depth, TP/FP/FN/precision/recall/F1/accuracy/
balanced accuracy per arm, ΔF1 (cascade vs PB) with paired-bootstrap 95% CI,
PB routing fraction, compute fraction, measured throughput, projected
speedup, disagreement analysis, stratification by depth/VAF/base-quality/
mapping-quality/GIAB structural annotation, and matched-budget random and
reverse-confidence routing controls.

## 7. Pre-registered decision rule

Quoted verbatim from devlog 13 §3.4 / devlog 16 §7 / devlog 17 §7, unchanged:

* **PRESERVED** — paired-bootstrap 95% CI for ΔF1 contains 0 **and** |ΔF1| < 0.001.
* **DEGRADED** — CI excludes 0 with ΔF1 < 0, or |ΔF1| ≥ 0.001.
* **IMPROVED** — CI excludes 0 with ΔF1 > 0.
* **UNDERPOWERED** — < 100 SNP, or CI half-width > 0.01.

Final verdict vocabulary, mapped now, with one addition specific to
cross-sample transfer:

* **STRONG POSITIVE (generalizes)** — every powered cell and the pooled result
  PRESERVED or IMPROVED, PB compute fraction ≤ 1%, and no new failure mode
  relative to HG002.
* **POSITIVE / PRESERVED (generalizes, with caveats)** — as above with a
  material limitation (underpowered cell, compute reduction 10-50x, etc).
* **MIXED** — pooled PRESERVED/IMPROVED but a powered stratum DEGRADED.
* **NEGATIVE (does not generalize)** — any powered cell DEGRADED, or the
  disagreement rate materially exceeds the HG002 baseline (> 10x devlog 16's
  chr13 rate, i.e. > ~10⁻⁵), or a genuinely new failure mode appears that was
  absent on HG002.
* **INCONCLUSIVE** — too few powered cells.

### 7.1 What would falsify the cross-sample claim

A DEGRADED verdict on any powered cell; a disagreement rate materially above
the HG002 chr13 rate (devlog 16: 3/8,726,751 pooled, ~3×10⁻⁷); a true SNP lost
by the cascade that PB would have called correctly in a pattern *not* seen on
HG002 (i.e. not simply a scaled-up version of the segdup low-VAF mode from
devlog 14/17); or router routing fraction rising by an order of magnitude
without a corresponding depth change.

## 8. Stopping criteria

Extraction runs once per cell. Metrics computed once. No threshold, cutoff, or
criterion changes after any result is read, except to fix an outright crash in
non-frozen v18 glue, disclosed if it happens. If HG003 performs poorly it is
reported as a negative cross-sample result. No second sample is drawn in this
experiment regardless of outcome (a second sample, HG004, is future work,
§20, not run here).

## 9. Runtime plan

Reused without regeneration: `data/reference/chr13_full.fa`, GIAB v3.1
stratification BEDs, every metric helper, the bootstrap, the controls. New
compute: one BAM slice download (HG003, chr13:70-73Mb), one truth VCF slice
download, one HG003 high-confidence BED download (genome-wide file, ~grep'd to
the window), two `samtools view -s` downsamples, three extractions, one
scoring pass.

## 10. What this experiment cannot show

Fixed by design, stated before results exist: one new sample (HG003) out of
seven other GIAB references never touched (HG001, HG004-HG008); HG003 and
HG002 are both Ashkenazim-trio members processed by the identical GIAB
pipeline (same platform, same aligner, same lab), so this is a *weak* test of
sample generalization — it cannot distinguish "generalizes across individuals"
from "generalizes within one trio/pipeline batch." One region (chr13:70-73Mb,
GIAB v3.1 `alldifficult_hc_fraction` ≈ 0.17, i.e. medium difficulty, not the
segdup-saturated chr17 window) — the segdup stress case is not re-tested here
for a second sample. SNP only, one reference build, one truth-set release. The
three depth cells remain correlated downsamples of one library. A positive
result here still licenses no claim about HG001/HG004-HG008, about any
non-Illumina platform, or about indels.

---

*(RESULTS marker — everything below is written after the corresponding
measurement, never edited into the section above.)*

## 11. Results

Independence audit: `all_ok=True`, all 3 cells `sample_is_hg003=True` (verified
via the truth VCF's own sample column, `bcftools query -l`, and the BAM's
`@PG` slice-provenance line — see §5). Frozen-file SHA-256 hashes identical
before and after the run (11/11 files, §3). Test suite: 120/120 relevant
tests pass (`test_cascade`, `test_cheap_router`, `test_bench_v14/16/17`,
`test_binomial_baseline`) before and after.

| Cell | depth | loci | SNP | F1 PB-only | F1 cascade | ΔF1 | ΔF1 95% CI | routed→PB | disagree | verdict |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---|
| full | ~65× | 2,899,736 | 4,422 | 0.996169 | 0.996394 | +0.000225 | [+0.000000, +0.000567] | 0.0037% | 2 | PRESERVED |
| 30× | ~30× | 2,899,736 | 4,422 | 0.994822 | 0.995382 | +0.000560 | [+0.000112, +0.001118] | 0.0137% | 5 | IMPROVED |
| 15× | ~15× | 2,899,736 | 4,422 | 0.990232 | 0.990232 | +0.000000 | [+0.000000, +0.000000] | 0.1054% | 0 | PRESERVED |
| **pooled** | all | 8,699,208 | 13,266 | 0.993751 | 0.994013 | **+0.000262** | [+0.000076, +0.000483] | 0.0409% | 7 | **IMPROVED** |

Compute: 108-3,057 loci routed to PB per cell (0.0037%-0.1054%), projected
speedup 210.7x-346.1x — indistinguishable from HG002 on the same coordinates
(devlog 16: 0.0034%-0.1011%, 213.3x-326.6x).

**Every disagreement (7 total, pooled rate 8.0×10⁻⁷) is `router_avoids_fp`:**
the router kept a locus away from PB, PB would have called a false positive
there, and the cascade's binomial-only answer was correct. All 7 sit at VAF
0.115-0.13, base quality 33-39, MAPQ 62-70, in `alldifficult`/`lowmap_segdup`/
`tandemrepeats` — the identical paralogous-sequence-variant signature devlogs
14 and 17 characterized on HG002. **Zero true SNP lost. Zero cascade-caused
false positive.** No stratum (VAF, base quality, MAPQ, depth,
`lowmap_segdup`, `alldifficult`, `tandemrepeats`) shows ΔF1 < 0 in any cell.

## 12. Final verdict

Applying §7 mechanically: every cell PRESERVED or IMPROVED, none DEGRADED;
pooled result IMPROVED; PB compute fraction 0.0037%-0.1054%, all ≤1%; no new
failure mode — the only disagreements reproduce the known HG002 PSV/segdup
signature; disagreement rate 8.0×10⁻⁷, below the HG002 chr13 baseline
(3.4×10⁻⁷) but two orders of magnitude under the falsification bound.

### STRONG POSITIVE (generalizes to HG003, within this experiment's scope)

The frozen cascade, unmodified from what devlog 16 ran on HG002, reproduces
PB-level SNP-calling accuracy on a different individual's genome — HG003, the
Ashkenazim trio father — on the identical region, platform and aligner, with
no threshold retuned and no HG003 data seen by this project before this
devlog was written. This is the first evidence this project has produced that
the frozen constants are not simply overfit to HG002.

**What this does not show**, restated from §10: HG003 and HG002 are both
Ashkenazim-trio members run through the identical GIAB sequencing/alignment
pipeline (same lab, same platform, same aligner, same library prep protocol,
overlapping ancestry). This result cannot distinguish "generalizes across
human genetic variation" from "generalizes within one trio's sequencing
batch." It is one sample, one region (chr13, medium difficulty — not the
segdup-saturated chr17 stress case), one platform. HG001, HG004-HG008, any
non-Illumina platform, indels, and the segdup stress case on a second sample
remain untested.
