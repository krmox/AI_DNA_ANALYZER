# chr20 300x Post-v1.0 Benchmark: AI_DNA_ANALYZER vs DeepVariant vs GATK HaplotypeCaller

Date: 2026-09-23 / 2026-09-24
Status: **same-sample, different-chromosome replication / benchmark — NOT an independent validation**
(HG002 was used elsewhere in this project's threshold-fitting history; chr20 itself was not the
fitting region, but is the same sample. See Limitations, §15.)

Frozen v1.0 thresholds, production algorithm, manuscript, and figures were **not modified** by
this benchmark. This is a separate, post-v1.0 engineering/benchmark exercise.

---

## 1. Dataset

- **Sample:** HG002 / NA24385 (GIAB Ashkenazim Trio, son)
- **BAM:** `HG002.GRCh38.300x_chr20.bam` — official NCBI/GIAB storage, NovoAlign V3.02.07,
  Illumina HiSeq2500 v1 Rapid, 2×148bp PCR-free paired-end, ~300x nominal genome-wide coverage.
  Full download/validation record: `/mnt/archive/AI_DNA_ANALYZER_benchmark/HG002_GRCh38_chr20/dataset_manifest.txt`.
  - Size: 12,100,052,820 bytes; local SHA256 `2ce70e1a5f08c90cbc7fd3d336deab9ab8f0709ab3f00bd7cc78bdbed8482965`
    (no official checksum published for this chr20-sliced file — see that manifest's Open Concerns).
  - 124,888,117 total reads; 123,038,598 mapped to chr20; index-validated, random-access-validated.
  - **No `@RG` line** in the header — required GATK-specific preprocessing, see §6 and
    `CALLER_SPECIFIC_PREPROCESSING.md`. AI_DNA_ANALYZER and DeepVariant used the file unmodified.

## 2. Reference

`GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna` — official NCBI
`seqs_for_alignment_pipelines.ucsc_ids/` release.
- md5 `a056c57649f3c9964c68aead3849bbf8`, verified against official `md5checksums.txt`.
- 2580 contigs — exact match to the BAM header's `@SQ` count.
- chr20 sequence byte-identical (sequence-only md5 `b18e6c531b0bd70e949a7fc20859cb01`) to this
  project's pre-existing `data/reference/chr20_full.fa`.
- Chosen over the alt-containing `full_plus_hs38d1_analysis_set.fna.gz` because the BAM's own
  `@PG` line names `GRCh38_full_plus_hs38d1_analysis_set_minus_alts` — "no_alt" = "minus_alts".

## 3. Truth set

GIAB HG002 GRCh38 **NISTv4.2.1**:
- `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` (+`.tbi`)
- Source: `ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/`
- v5.0q exists as a newer draft/preview release; v4.2.1 chosen as the last fully-finalized NIST
  release (user-confirmed choice, see conversation).
- chr20 truth SNP bases (this project's own labelling, provider-classified): reported per-run
  below; hap.py's own `TRUTH.TOTAL` differs slightly due to hap.py's independent normalization
  (71,333 SNP loci in the confident regions per hap.py).

## 4. Confident regions

`HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`, chr20 subset extracted and verified
**byte-identical** to this project's pre-existing `data/giab_hg002_chr20/chr20_highconf.bed`
(10,192 intervals).

## 5. Software versions

| Component | Version | Verified via |
|---|---|---|
| AI_DNA_ANALYZER | commit `9c0933e4edc9bbdb2df64878e58e5d6101ce9d79`, worktree `.claude/worktrees/opt-chr20` (branch `fix/hg005-extraction-integrity`) — **not on `main`** | `git log -1` |
| DeepVariant | 1.10.0 | `run_deepvariant --version` → "DeepVariant version 1.10.0" |
| GATK | 4.6.2.0 | `gatk --version` → HTSJDK 4.2.0, Picard 3.4.0 |
| hap.py | 0.3.15 (`quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0`) | `hap.py --version` |
| Clair3 | not run — see §6 | — |
| gcc/g++ | 15.3.1 20260722 (Red Hat) | |
| CMake | 3.31.11 | |
| Python | 3.14.3 | |
| samtools/bcftools/htslib | 1.23.1 | |

Full environment record: `ENVIRONMENT.md` in this directory.

## 6. Exact commands (representative; full commands in `logs/*.log`)

**AI_DNA_ANALYZER** (native C backend `pileup_native.c`/`reads_native.c`, built via `native/build.sh`):
1. `cache_counts.py`-equivalent: `load_counts(fasta, bam, vcf, bed, (0, 64444167), contig="chr20", return_positions=True)`
2. Read-level tensor: `load_reads(fasta, bam, vcf, bed, (0, 64444167), contig="chr20", max_reads=48)`
3. PB LLR (streaming, see §13): custom script reading `reads.npz`'s `reads.npy` zip member in
   1M-row chunks via `zipfile`, calling `extract_quality_evidence` / `candidate_alt` /
   `poisson_binomial_llr` verbatim, unmodified, from `quality_error_model.py`.
4. `binomial_llr = cheap_router.binomial_llr(counts)` (cheap, closed-form, unmodified).
5. Frozen cascade + Method C genotyping: `experimental/chr20_validation/run_chr20_300x_pipeline.py`
   — every statement in `extract_records`/`optimized_method_c`/`write_vcf_c`/`HEADER_TMPL` copied
   verbatim from the frozen `experimental/chr20_validation/run_chr20_pipeline.py`; only I/O paths differ.

**GATK:**
```
gatk AddOrReplaceReadGroups -I HG002.GRCh38.300x_chr20.bam \
  -O gatk_input_with_rg/HG002.GRCh38.300x_chr20.withRG.bam \
  -RGID HG002_chr20_300x -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM HG002 --CREATE_INDEX true
gatk --java-options "-Xmx20g" HaplotypeCaller \
  -R GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
  -I HG002.GRCh38.300x_chr20.withRG.bam -O gatk_chr20_300x.vcf.gz -L chr20 --native-pair-hmm-threads 6
```

**DeepVariant:**
```
/opt/deepvariant/bin/run_deepvariant --model_type=WGS \
  --ref=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
  --reads=HG002.GRCh38.300x_chr20.bam --output_vcf=deepvariant_chr20_300x.vcf.gz \
  --regions=chr20 --num_shards=8 --sample_name=HG002
```

**hap.py** (default xcmp engine — see §13 for why not `--engine=vcfeval`):
```
hap.py truth.vcf.gz query.vcf.gz -f chr20_highconf.bed -r <reference> -o <prefix> -l chr20 --threads 8
```

## 7. Hardware

Intel Core i7-6700 @ 3.40GHz, 8 logical cores, 31 GiB RAM, Linux 6.19.12-200.fc43.x86_64
(Fedora 43). No GPU used by DeepVariant ("Could not find cuda drivers").

## 8. AI_DNA_ANALYZER results

- **Routing:** 56,242,697 loci scored (SNP/ref frame); only **707 loci (0.001%)** fell within the
  frozen router band (`|binomial_LLR − 7.0| ≤ 5.411872376933351`) and were routed to PB. At 300x
  depth the fixed-epsilon binomial LLR is almost always far from the decision boundary, so the
  cascade result is **almost identical to binomial-only** — this is a genuine, expected property
  of the frozen router at much higher depth than its calibration regime, not a bug.
- **Cascade VCF:** 72,661 SNP calls (789 forced 0/0→0/1 by Method C's floor).
- **PB-only VCF:** 72,872 SNP calls (1,403 forced 0/0→0/1).
- **Detection-layer raw accuracy** (pre-hap.py, this project's own metric): binomial-only 0.999949,
  PB-only 0.999947, cascade 0.999951.
- **Rescue composition:** 892 PB-rescuable loci; 191 rescued by the router (recall 0.214); false
  rescue rate 0.173; 56.4% of routed loci were already correct under binomial alone.
- AI_DNA_ANALYZER in this pipeline calls **SNPs only** — no indel calling. Comparisons below are
  SNP-only for this caller; DeepVariant/GATK also called indels, reported separately.

## 9. Competitor results (raw caller output, pre-hap.py)

- **DeepVariant:** 13,855+ variant records written per shard progression; final chr20 VCF non-empty, valid.
- **GATK:** 137,027 total variant records (SNPs + indels combined), depths 270–410x observed,
  consistent with 300x nominal coverage.

## 10. hap.py results (against NISTv4.2.1 truth + confident BED, default xcmp engine)

| Run | Type | TRUTH.TOTAL | TP | FN | QUERY.TOTAL | FP | Recall | Precision | F1 |
|---|---|---|---|---|---|---|---|---|---|
| **AI_DNA_ANALYZER cascade** | SNP | 71,333 | 70,012 | 1,321 | 72,661 | 2,649 | 0.981481 | 0.963543 | **0.972429** |
| AI_DNA_ANALYZER cascade | INDEL | 11,256 | 0 | 11,256 | 0 | 0 | 0.0 | 0.0 | n/a (not called) |
| AI_DNA_ANALYZER pb-only | SNP | 71,333 | 70,019 | 1,314 | 72,872 | 2,853 | 0.981579 | 0.960849 | 0.971104 |
| **DeepVariant 1.10.0** | SNP | 71,333 | 71,036 | 297 | 90,363 | 43 | 0.995836 | 0.999395 | **0.997613** |
| DeepVariant 1.10.0 | INDEL | 11,256 | 11,220 | 36 | 21,527 | 16 | 0.996802 | 0.998628 | 0.997714 |
| **GATK 4.6.2.0** | SNP | 71,333 | 71,056 | 277 | 114,696 | 911 | 0.996117 | 0.987347 | **0.991712** |
| GATK 4.6.2.0 | INDEL | 11,256 | 11,224 | 32 | 23,485 | 44 | 0.997157 | 0.996242 | 0.996699 |

Full `.summary.csv`, `.extended.csv`, `.metrics.json.gz`, `.roc.*.csv.gz`, `.vcf.gz` for every run
are in `happy/{ai_dna_analyzer,ai_pb_only,deepvariant,gatk}/`.

**QUERY.TOTAL** differs substantially by caller because DeepVariant/GATK's QUERY.TOTAL includes a
large `QUERY.UNK` component (multi-allelic / representation-normalized records hap.py couldn't
uniquely classify: 19,234 for DeepVariant, 42,700 for GATK SNP rows) — this is a hap.py
normalization artifact of each caller's own VCF representation style, not a call-count in the
same sense as AI_DNA_ANALYZER's un-normalized QUERY.TOTAL (Frac_NA=0.0, no UNK records at all,
since AI_DNA_ANALYZER emits only simple biallelic SNP records).

## 11. Runtime comparison (wall-clock, this hardware, this run only)

| Stage | Caller | Time |
|---|---|---|
| End-to-end | GATK HaplotypeCaller | ~3h 43m (18:23 → 22:06) |
| End-to-end | DeepVariant | 1h 27m 53s |
| Extraction (counts, native C) | AI_DNA_ANALYZER | 13.2 min (790.5s), 56.27M loci |
| Extraction (read tensor, native C) | AI_DNA_ANALYZER | 109.2 min (6553.8s), 21.6 GB raw |
| PB LLR (Poisson-binomial, streaming) | AI_DNA_ANALYZER | 185.6 min (11137.8s), 56,266,816 loci |
| Binomial LLR (closed-form) | AI_DNA_ANALYZER | negligible (<1s, vectorized over cached counts) |
| Router + genotyping + VCF write | AI_DNA_ANALYZER | 96.0s (both cascade and pb-only VCFs) |
| **AI_DNA_ANALYZER total** | | **~5h 10m** (dominated by PB LLR + read-tensor extraction) |
| hap.py evaluation (×4 runs) | hap.py | 1:08–3:44 each (not part of any caller's runtime) |
| GATK-specific: AddOrReplaceReadGroups | preprocessing | 20.12 min (separately logged, not counted in GATK's own runtime above) |

**AI_DNA_ANALYZER's PB LLR stage is, by a wide margin, the dominant cost** (185.6 of ~310 total
minutes, ~60%) — consistent with the frozen code's own documentation calling PB "the expensive
part." This is measured, not a projected/theoretical speedup claim.

## 12. Memory comparison

| Caller | Peak RSS (measured) |
|---|---|
| AI_DNA_ANALYZER, read-tensor extraction (bulk, in-RAM) | ~21.6 GB (matches raw array size) |
| AI_DNA_ANALYZER, PB LLR (bulk-load attempt) | **OOM-killed twice** (see §13) |
| AI_DNA_ANALYZER, PB LLR (streaming, final working version) | ~8.0–8.2 GB, stable |
| GATK HaplotypeCaller | `-Xmx20g`/`-Xmx24g` heap cap set; observed container RSS ~1GB at one point mid-run (not a true peak measurement — not continuously profiled) |
| DeepVariant | not continuously profiled; ran successfully within 31GB system RAM alongside other load |

No peak-RSS profiler (e.g. `/usr/bin/time -v` or cgroup memory.peak) was attached to the GATK/DeepVariant
docker runs; only system-level `free -h` snapshots and the AI_DNA_ANALYZER Python process's own RSS
were tracked. This is a real gap — flagged here rather than presenting an invented number.

## 13. Error/failure notes (nothing hidden)

1. **GATK: no `@RG` in BAM.** First HaplotypeCaller attempt crashed
   (`IllegalStateException: the sample list cannot be null or empty`). Fixed by running
   `AddOrReplaceReadGroups` on a **separate copy** of the BAM (original untouched); documented
   fully in `CALLER_SPECIFIC_PREPROCESSING.md`.
2. **Docker/containerd storage on NTFS.** `/mnt/archive` is NTFS (`ntfs3`), which cannot host
   overlayfs (layer extraction failed on Linux-only path/xattr semantics — e.g. `App::Cpan.3`
   filenames, `dpkg/info/*.list`). Fixed with a 60GB ext4 loopback image mounted from a file on
   the NTFS drive (`container_storage.img` → `container_storage_ext4/`), Docker/containerd
   repointed there. No system-root or existing-user-data changes.
3. **Root filesystem exhaustion.** Root disk (`/`) was independently near-full
   (`/home/mark` 88GB including pre-existing podman storage, `/var/lib/flatpak` 9GB) — unrelated
   to this benchmark, not touched.
4. **PB LLR: two OOM kills.** Loading the full 21.6GB read tensor via `np.load(...)['reads']`
   (which fully decompresses the npz member into one array) plus the 56.27M×44-channel feature
   tensor this project's `cache_locus_evidence.py` also builds by default, exceeded available
   RAM+swap twice (once via the OS-level OOM killer, once via Claude Code's own memory-pressure
   background-task reaper). Root-caused (not guessed) via `journalctl` OOM-kill entries and by
   noting the process died exactly during npz decompression. **Fix:** (a) skip the 44-channel
   feature tensor entirely (not needed — only `pb_llr` feeds the frozen cascade, the tensor is
   for an unused neural stage), and (b) stream the `reads.npy` zip member directly via Python's
   `zipfile` + `numpy.lib.format` header parsing, in 1M-row chunks, never materializing more than
   one chunk (~8GB peak, vs ~30GB+ for the bulk approach). This is an I/O engineering fix; the
   per-locus math (`extract_quality_evidence`, `candidate_alt`, `poisson_binomial_llr`) is
   byte-for-byte the same frozen function calls either way.
5. **hap.py `--engine=vcfeval` unavailable.** The `quay.io/biocontainers/hap.py:0.3.15` image
   does not bundle `rtg` (RTG Tools); `vcfeval` mode failed with `rtg: command not found`
   (return code 127). All four hap.py runs used hap.py's **default xcmp engine** instead — a
   tool-configuration fact, not a scientific substitution; documented here for transparency.
6. **`cache_counts.py`'s CLI has no `--contig` flag** and defaults to `chr21`, silently producing
   0 loci against our chr20 data on the first attempt. Caught immediately (0-loci output is an
   obvious signal), fixed by calling `load_counts(..., contig="chr20")` directly instead of the
   CLI wrapper.
7. `HaplotypeCaller`'s first successful attempt was interrupted mid-run once by Claude Code's
   background-task memory-pressure reaper (unrelated to GATK itself — the *docker container*
   survived independently of the killed orchestrating shell wrapper, and completed normally once
   noticed).

No result in this report was fabricated or silently patched over; every failure above is either
fixed via a documented, verbatim-preserving engineering change, or left as an open gap (§12).

## 14. Raw artifact locations

- BAM/reference/truth/BED: `/mnt/archive/AI_DNA_ANALYZER_benchmark/{HG002_GRCh38_chr20,reference,truth_HG002_GRCh38_v4.2.1}/`
- AI_DNA_ANALYZER caches + VCFs: `/mnt/archive/AI_DNA_ANALYZER_benchmark/ai_dna_analyzer_run/`
  (`chr20_300x_counts.npz`, `chr20_300x_reads.npz`, `chr20_300x_pb_llr.npz`,
  `ai_cascade_chr20_300x.vcf.gz`, `ai_pb_only_chr20_300x.vcf.gz`, `pipeline_meta_chr20_300x.json`)
- GATK: `/mnt/archive/AI_DNA_ANALYZER_benchmark/gatk_run/`,
  `/mnt/archive/AI_DNA_ANALYZER_benchmark/HG002_GRCh38_chr20/gatk_input_with_rg/`
- DeepVariant: `/mnt/archive/AI_DNA_ANALYZER_benchmark/deepvariant_run/`
- Container storage engineering: `/mnt/archive/AI_DNA_ANALYZER_benchmark/container_storage.img`,
  `container_storage_ext4/`
- hap.py outputs: `research/benchmark_chr20_300x/happy/{ai_dna_analyzer,ai_pb_only,deepvariant,gatk}/`
- Logs: `research/benchmark_chr20_300x/logs/*.log`
- Pipeline script (verbatim-logic adapter): `.claude/worktrees/opt-chr20/experimental/chr20_validation/run_chr20_300x_pipeline.py`
- Environment record: `research/benchmark_chr20_300x/ENVIRONMENT.md`
- Caller-specific preprocessing record: `research/benchmark_chr20_300x/CALLER_SPECIFIC_PREPROCESSING.md`

## 15. Limitations

- **Not independent validation.** chr20 was not the region used for threshold fitting
  (`threshold_fit_span: chr21:32-40 Mb` per the project's own prior chr20 validation script), but
  HG002 as a *sample* has been used elsewhere in this project's history. This run is honestly
  labeled **same-sample, different-chromosome replication / benchmark**.
- **AI_DNA_ANALYZER is SNP-only** in this pipeline — no indel model is exercised here. All
  INDEL rows for it in §10 are structural zeros (nothing was attempted), not a failure.
- **`MAX_READS=48` read-level cap.** The frozen read-level tensor keeps at most 48 reads/locus
  (a documented, deterministic, allele-unbiased subsample), calibrated in the code's own comment
  for "max depth 43 (train) and 38 (test)". At 300x, this means the PB stage sees roughly
  **48/~300 ≈ 16%** of available reads per locus — a real scope limitation of the frozen model at
  this depth, not something this benchmark introduced or could fix without touching frozen code
  (which was explicitly out of scope).
- **PB-routed fraction is tiny at 300x (0.001%)**, so the cascade result is statistically close
  to binomial-only on this data; this benchmark cannot speak to PB/router behavior in the regime
  it was actually designed for (lower depth, wider router engagement) — see `bench_v20_chr20.py`
  and related devlogs for that regime's own prior results.
- **No official checksum** for the chr20-sliced BAM (only the parent 560GB file); provenance
  rests on official-server-path + structural validation (see the original `dataset_manifest.txt`).
- **hap.py used the default xcmp engine, not vcfeval**, because `rtg` wasn't available in the
  chosen container image. xcmp and vcfeval do not always produce identical results on complex
  regions; this affects all three callers equally (compared under the same engine), but is not
  the field's most commonly cited hap.py configuration.
- **Peak-memory figures for GATK and DeepVariant are not rigorously measured** (no cgroup/`time -v`
  profiling attached); only AI_DNA_ANALYZER's own Python-process RSS was tracked precisely.
- Clair3 was **not run**: per the original plan's rule, it is not included in the primary
  numerical comparison automatically because it is designed primarily for long-read calling, and
  no verified Illumina+GRCh38 Clair3 model/configuration was confirmed available in this session's
  container image search. Marking as `NOT_APPLICABLE_FOR_THIS_INPUT` rather than forcing a
  possibly-miscalibrated run.

## 16. Final factual comparison

No subjective "winner" is declared. Measured SNP F1 on this chr20/300x/NISTv4.2.1 replication run,
default-xcmp-engine hap.py, same BAM/reference/truth/confident-BED for all three:

| Caller | SNP F1 | SNP Recall | SNP Precision | Indel F1 | Total runtime | Notes |
|---|---|---|---|---|---|---|
| AI_DNA_ANALYZER (cascade) | 0.972429 | 0.981481 | 0.963543 | — (not called) | ~5h10m (PB LLR-dominated) | Frozen thresholds unmodified; PB engaged on 0.001% of loci at this depth |
| DeepVariant 1.10.0 | 0.997613 | 0.995836 | 0.999395 | 0.997714 | 1h27m53s | CPU-only (no GPU detected) |
| GATK 4.6.2.0 | 0.991712 | 0.996117 | 0.987347 | 0.996699 | 3h43m (+20min RG preprocessing) | Required @RG addition to run at all |

These are the measured values under the conditions and caveats documented above (§13, §15).
Differences in QUERY.TOTAL/FP driven by hap.py's UNK-record normalization (§10) mean precision
numbers are not perfectly apples-to-apples without deeper per-caller VCF-representation
normalization, which was not performed beyond what hap.py does internally.
