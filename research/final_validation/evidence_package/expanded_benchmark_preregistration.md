# Expanded benchmark — preregistration

Written 2026-10-29 BEFORE any new run, download or new analysis for this benchmark. The existing audit (`research_history.md`, `claims_audit.md`, `PHASE2_BLOCKER.md`) is unchanged and remains in force.

## 0. What this is and is not
- **NOT an independent blind-genome validation.** HG001 and HG007 were **previously evaluated** by this project (campaigns of 2026-10-28/29: 10/20/30× on chr17 and chr19). They were never used in R1.1/R2 development, but their 10/20/30× results are already known to the analysts. This is an **expanded validation / replication**. No genuinely unseen genome is used. `PHASE2_BLOCKER.md` still stands.
- Honest status of the endpoints: the 10/20/30× data and point estimates already exist (seen). Only the **15× and 50× conditions, the HG007 confidence intervals, the pooled analysis, the runtime summary and the DeepVariant arm are new**. Endpoint rules and margins below are fixed before the new data are produced; for 10/20/30× they are retrospective and are labelled that way in the report (`EXISTING`, not `NEW`).
- Frozen R2, R1.1 and their models are not changed. No tuning on any result.

## 1. Data (existing, reused)
| Item | HG001 / NA12878 | HG007 / NA24695 |
|---|---|---|
| Native BAM | GIAB HiSeq 2×148 novoalign 300×, GRCh38 analysis set, streamed per contig from the URL in `/mnt/hdd/AI_DNA_ANALYZER_hg001/logs/source_url.txt`. Natives were deleted after the first campaign; **re-fetched** (the only download), size/sha256 must equal `independent_genome/logs/fetch_sha256.log`, else stop. | existing natives on HDD (`/mnt/hdd/AI_DNA_ANALYZER_hg007`), no download |
| Truth / BED | GIAB v4.2.1 GRCh38, `HG001_…benchmark.{vcf.gz,bed}` (local sha256 in `SHA256SUMS.local.txt`) | `HG007_…benchmark.{vcf.gz,bed}` |
| Reference | GRCh38 no-alt + hs38d1 analysis set `.fna` (same file for all) | same |
| Regions | chr17 and chr19 (confident BED, whole contig) | chr17 and chr19 |
Cells = 2 samples × 2 contigs = 4. Contigs fixed now (the existing ones); no contig is added or dropped after results.

## 2. Conditions
| Depth | Seeds | Status |
|---|---|---|
| 10×, 20×, 30× | 101, 102, 103 | **EXISTING** (reused unchanged) |
| 15× | 101, 102, 103 | **NEW** |
| 50× | 101 (one seed, as in the earlier preregistration) | **NEW** |
Downsampling: `samtools view -b -s SEED.FFFF`, FFFF = round(10000 × target / BED-mean native depth); BED-mean depth from `samtools bedcov` of the native BAM (same recipe as existing runs). Native-depth runs are out of scope for new data (HG007 native results exist and are reported as EXISTING).

## 3. Callers
| Caller | Version / command | Role |
|---|---|---|
| Varion R2 | frozen commit `3cdfe5a`, `dnav2` sha256 `6c675174…`, `frozen_R2.model` sha256 `649215d5…`, `--mode cascade --r2 --threads 8` | subject |
| Varion R1.1 | `dnav2` sha256 `e83d87fd…`, `R1.model` `63816707…`, cascade, 8 threads | directly comparable |
| bcftools 1.21 | `mpileup -q20 -Q13 -d10000 -a FORMAT/AD \| call -mv \| view -v snps -i 'QUAL>=30'`, 32 shards, 8 parallel (`research/stage3/scripts/bcftools_sharded.sh`) | directly comparable (same version as the existing HG001/HG007 runs) |
| DeepVariant 1.10.0 | conditional, see §7 | separate comparison, only if feasible |
Hashes of R2 (21 sources, binary, model) are verified before and after execution; any mismatch invalidates the run set.

## 4. Evaluation (identical for all callers and conditions)
hap.py 0.3.15 (xcmp engine, conda env `hp`), `hap.py TRUTH QUERY -f CONFIDENT_BED -r REF -l CONTIG --threads 4`, SNP rows with `Filter=PASS`. TP = TRUTH.TP, FN = TRUTH.FN, FP = QUERY.FP; precision/recall/F1 as reported by hap.py and re-derived from TP/FP/FN (a mismatch > 1e-6 is a harness defect). **INDEL results are reported separately; Varion is SNP-only so no INDEL claim is made.** Scores produced by different evaluators are never put in one table.

## 5. Endpoints and decision rules (fixed now)
Unit of analysis: SNP F1 averaged over seeds within a cell, then over cells (contigs). Δ = difference in absolute F1 (not relative).

**Confidence interval.** Percentile 95% paired block bootstrap, 2-Mb blocks (`floor((pos−1)/2e6)`), 4000 resamples, RNG seed 20261027 (as in `research/stage4/PREREG_R2.md`). Resampling is per contig; one draw of blocks is shared across that contig's seeds and callers; the statistic is recomputed from summed TP/FP/FN per condition then averaged. Sample-level interval: HG001 (chr17+chr19) and HG007 (chr17+chr19) separately; pooled interval: all 4 cells. Pooled = two samples, **not a population interval**; it quantifies spatial variation inside four contigs only. Seeds are re-subsamples of one BAM and are not replicates of the genome.

**Primary endpoint.** SNP F1 at 30×: Δ(R2 − R1.1) and Δ(R2 − bcftools), per sample and pooled.
- **P1 (gain over R1.1):** PASS if pooled Δ ≥ +0.003 and CI lower bound > 0 (target carried over from PREREG_R2). Reported per sample too. *At writing, the existing point estimates (+0.0023 / +0.0011) already fall below +0.003; P1 is expected FAIL and will be reported as such.*
- **P1b (any improvement over R1.1):** PASS if CI lower bound > 0 (descriptive; not the preregistered target).

**Secondary endpoints.**
- **P2 (non-inferiority to bcftools):** at each depth d ∈ {10, 15, 20, 30, 50}, PASS if the lower 95% bound of Δ(R2 − bcftools) ≥ **−0.0005** (margin carried over from PREREG_R2; not new). UNASSESSED for a depth if any of the 4 cells lacks all three callers' OK results at that depth. 50× has one seed: its CI reflects only block variation. Superiority over bcftools at a depth: lower bound > 0 (stated separately from non-inferiority).
- **P3 (absolute level):** F1(R2) ≥ 0.990 at 30× and at 50×, per cell; claimed only if bcftools also reaches it on that cell (ceiling check), otherwise reported as "reached by R2 only".
- **P4 (low-depth non-regression vs R1.1):** Δ(R2 − R1.1) ≥ −0.002 at 10× and 15× (point estimate), CI reported.
- Per-depth precision, recall, F1, TP/FP/FN per seed and caller are tabulated (no pooling of TP/FP/FN across contigs for the headline).

**Resource guardrails (R2 vs R1.1, same BAM):** wall ≤ +10% or ≤ +2 s (whichever larger), CPU ≤ +10%, peak RSS ≤ +15%, evaluated on the median over seeds per condition from `wait4` rusage. Wall time is valid only for runs executed when no other benchmark job is running (all downloads finish before any timed run; run sequentially; `loadavg` before each run recorded). Runs affected by an overlap are flagged and not used for wall ratios. Single-machine, single-hardware results; no claim beyond this machine (i7-6700, 8 threads).

**Result states.** PASS / FAIL (rule evaluated and violated) / UNASSESSED (data missing or run failed) / EXISTING (data already seen before this preregistration). Nothing is passed by implication.

## 6. Failed runs, reruns, exclusions
Every failed run is logged in `FAILED_RUNS.md` with the cause and kept. A harness defect (not a caller defect) may be fixed and the affected conditions rerun in full; the failed rows are archived, not deleted. A caller failure is not rerun to obtain a better number (one deterministic retry only for infrastructure faults such as a remount, documented). No seed, contig or depth is removed after seeing results. Missing calls and filtering follow the scorer: monoallelic and multiallelic sites are as handled by hap.py xcmp; no extra normalisation beyond hap.py's.

## 7. DeepVariant arm (conditional, separate table)
Feasibility is audited first (`deepvariant_feasibility.md`). It runs only if (a) the official image `google/deepvariant:1.10.0` (the version used in `benchmark_chr20_300x`, digest recorded at pull) can be unpacked and executed without root/Docker/Podman, and (b) a smoke test on a small slice succeeds. Scope if feasible: **30×, seed 101, on the 4 cells** (primary DV comparison), `run_deepvariant --model_type=WGS --regions CONTIG --num_shards=8`, same BAMs as the other callers, same hap.py/truth/BED/reference. DV is given the whole contig; Varion's call domain is the confident BED; hap.py scores only inside the BED for all, so this difference is documented but unavoidable. Optional 10× and 20× (seed 101) only if time allows; not required for any pass. Per-run abort limit: 4 h wall or > 24 GB RSS → recorded as FAILED, not retried with other settings.
DV endpoints are descriptive: SNP F1 / P / R with block-bootstrap CIs (single seed). **No non-inferiority or superiority claim vs DeepVariant is preregistered**; none will be made. Where DV cannot be run, its comparison is reported UNASSESSED with the precise blocker.

## 8. Deliverables
`expanded_benchmark_preregistration.md` (this), `expanded_benchmark_report.md`, `expanded_results.csv`, `expanded_runtime_memory.csv`, `deepvariant_feasibility.md`, `limitations.md`, `reproduction_guide.md` in `research/final_validation/evidence_package/`; scripts, logs, raw `RESULTS.jsonl` in `research/final_validation/expanded/`. Nothing is committed.

## Addendum A (2026-10-29, written after the HG001 re-fetch and BEFORE any new caller run; changes no endpoint)
The §1 rule "re-fetched BAM sha256 must equal the first-campaign value" is **not satisfiable literally**: `samtools view -o` writes its own output path into the `@PG` header line, and the re-fetch used a different path (NVMe instead of HDD), so the whole-file sha256 differs (chr19 10,293,591,311 B vs 10,293,591,289 B; chr17 15,173,222,007 vs 15,173,221,985 — both +22 B, the path-length difference). Reads come from the same remote GIAB BAM and the same contig. Replacement identity check (fixed now): (i) header differs only in that `@PG samtools` CL path; (ii) BED-mean depth of the re-fetched native BAM equals the value recorded in the first campaign (`meta/HG001_*_native_cov.json`) to ≥ 12 significant digits. chr19: 267.296853257729 vs 267.2968532577294 → satisfied. chr17: verified before use (result in the report). Downsampling flags are derived from that recorded depth, so the same reads are selected; the downsampled BAM sha256 values will differ for the same header reason and are not used as an identity test. Stop rule: if chr17's depth does not match, HG001 chr17 is excluded and reported UNASSESSED.
