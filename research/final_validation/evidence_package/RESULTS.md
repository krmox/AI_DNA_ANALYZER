# Expanded benchmark report — frozen Varion R2 on HG001 and HG007 (chr17, chr19)

**Read first.** HG001 and HG007 were **previously evaluated** by this project (10/20/30× campaigns, 2026-10-28/29). This is an **expanded validation / replication, NOT an independent blind-genome test**; no unseen genome was available (`PHASE2_BLOCKER.md`). Preregistration: `expanded_benchmark_preregistration.md` (sha256 + timestamp in `expanded_prereg.sha256`; Addendum A documents one deviation). Everything below is SNP-only, GIAB v4.2.1 high-confidence regions, GRCh38, Illumina HiSeq 2×148, two contigs per sample.

## 1. What ran
| Part | Provenance | Status |
|---|---|---|
| 10×, 20×, 30× (seeds 101–103), R2 / R1.1 / bcftools 1.21 | EXISTING (re-analysed with new CIs; HG007 CIs are new) | complete |
| **15×** (seeds 101–103) and **50×** (seed 101), same 3 callers, 4 cells | **NEW**: 48 runs, all exit 0, hap.py OK, 0 metric inconsistencies | complete |
| DeepVariant 1.10.0, 30×, seed 101 | **NEW**: 2 of 4 cells (HG007 chr19, HG001 chr19) | **chr17 cells UNASSESSED: stopped by user decision (too slow)** |
Frozen R2 identity verified before and after (commit 3cdfe5a, 21/21 sources, binary `6c675174…`, model `649215d5…`; `expanded/tests/identity_{before,after}.txt`). HG001 natives re-fetched; identity by BED-mean depth (Addendum A): chr19 267.296853257729 vs 267.2968532577294; chr17 279.737449352576 vs 279.7374493525764.

## 2. Primary endpoint — SNP F1 at 30× (95% block-bootstrap CI, 2-Mb blocks, 4000 resamples)
| Comparison | HG001 | HG007 | Pooled (4 cells) |
|---|---|---|---|
| R2 − R1.1 | +0.00109 [+0.00043, +0.00210] | +0.00160 [+0.00072, +0.00290] | **+0.00134 [+0.00073, +0.00211]** |
| R2 − bcftools | +0.00229 [+0.00114, +0.00423] | +0.00268 [+0.00109, +0.00514] | **+0.00249 [+0.00141, +0.00398]** |
**P1 (≥ +0.003 over R1.1): FAIL** for both samples and pooled. R2's gain over R1.1 is real in the sense that the CI excludes 0 (P1b PASS) but is about 45% of the preregistered target. Δ are absolute F1 differences (0.00134 = 0.13 F1 points).

## 3. Secondary endpoints
| Endpoint | Result |
|---|---|
| **P2** non-inferior to bcftools (lower bound ≥ −0.0005), 10/15/20/30/50× | **PASS at all five depths**, per sample and pooled; lower bound > 0 everywhere (superior). Pooled Δ: 10× +0.0351, 15× +0.0134, 20× +0.0062, 30× +0.0025, 50× +0.0018 [+0.0006, +0.0034]. 50× has one seed (block variation only). The advantage shrinks monotonically with depth. |
| **P3** F1 ≥ 0.990 at 30× and 50× | R2 reaches it in all 8 cells. bcftools also reaches it on chr19 (both samples, 30×/50×) and HG001 chr17 50×, but not on HG007 chr17 (30× and 50×) or HG001 chr17 30× → "reached by R2 only" there (not a data-ceiling claim). |
| **P4** no regression vs R1.1 at 10×/15× | **PASS**; pooled Δ 10× +0.0009 [+0.0006, +0.0014], 15× +0.0016 [+0.0011, +0.0023]. |
| **Resource guardrails** R2 vs R1.1 (NEW rows) | CPU and RSS within ±2.1% everywhere. Wall: 7 of 8 NEW conditions PASS; **HG001 chr19 50× FAIL** (+18.2%, +3.4 s; CPU −0.7%; single seed, 10 s runs). Treated as a FAIL by rule; with equal CPU it looks like timing noise but this was not repeated. |
Per-cell means: `ENDPOINTS_EXPANDED.md`; per-run TP/FP/FN/P/R/F1: `expanded_results.csv`; CIs: `BOOTSTRAP_EXPANDED.tsv`.

## 4. Speed and memory (this machine: i7-6700, 8 threads; NEW runs, median)
| Caller | Wall | CPU | Peak RSS |
|---|---|---|---|
| R2 | 10.5 s | 74.5 s | 216 MB |
| R1.1 | 10.5 s | 74.0 s | 215 MB |
| bcftools 1.21 (32 shards) | 59 s | 431 s | 117 MB |
| DeepVariant 1.10.0 (2 runs, CPU) | 1,688–1,879 s | 10,659–10,885 s | 8.9–9.2 GB |
R2 ≈ 5.6× faster than bcftools (wall), ≈ 5.8× less CPU, ≈ 1.8× its memory; ≈ 160–180× faster than DeepVariant on wall time and ≈ 145× less CPU, ≈ 42× less memory. Hardware- and pipeline-specific; DeepVariant here is CPU-only without tuning of its sharding.

## 5. DeepVariant (separate table, descriptive, same hap.py/truth/BED/reference/BAM/seed 101)
| Cell (30×, chr19) | R2 F1 | R1.1 | bcftools | **DeepVariant** | DV − R2 |
|---|---|---|---|---|---|
| HG001 | 0.99315 [0.99069, 0.99525] | 0.99256 | 0.99198 | **0.99559** [0.99317, 0.99779] | +0.00244 |
| HG007 | 0.99379 [0.99239, 0.99507] | 0.99284 | 0.99281 | **0.99652** [0.99498, 0.99786] | +0.00273 |
**DeepVariant was more accurate than R2** in both cells (higher precision AND recall; 35–44% fewer errors (FP+FN: 921→594 on HG001, 843→473 on HG007)). The individual intervals overlap (single seed, not paired-tested here), so the size of the gap is uncertain, but the direction is the same in both cells. **No superiority claim over DeepVariant is possible**; the supported statement is that R2 is far cheaper but not as accurate as DeepVariant on these two cells. chr17 cells: UNASSESSED. DeepVariant also calls indels (not scored here).

## 6. Answers
- **Ran:** 48 new three-caller runs + 2 DeepVariant runs; all completed; 0 harness inconsistencies.
- **Passed:** P2 (all depths), P3, P4, P1b (improvement over R1.1 is real); resource guardrails (except one wall-time cell).
- **Failed:** P1 (+0.003 over R1.1); one wall-time guardrail cell.
- **Unassessed:** DeepVariant on the two chr17 cells; DeepVariant at 10×/20×; native depth for new data; INDEL performance (Varion does not call INDELs); any statement about unseen genomes.
- **Unseen genome still needed:** **yes.** None of the 7 germline GIAB samples with a GRCh38 v4.2.1 truth set is untouched; a new GIAB germline benchmark (e.g., HG008) or an external dataset would be required.

## 7. Defensible statements
1. On HG001 and HG007 chr17/chr19, SNP-only, R2 is non-inferior and in fact superior to bcftools 1.21 at 10–50×, by a margin that falls from +0.035 F1 (10×) to +0.002 (30–50×).
2. R2 improves slightly but consistently on R1.1 (+0.001 to +0.002 F1 at 10–50×); the preregistered +0.003 was not reached.
3. R2 reaches F1 ≥ 0.990 at ≥ 30× in all tested cells at ~10 s per contig on 8 threads, ~216 MB.
4. DeepVariant 1.10.0 was ~0.0025 F1 higher on two cells at ~150× the compute.
**Not supported:** "revolution", genome-wide performance, generalization to other samples/platforms, superiority over DeepVariant, INDEL performance.
