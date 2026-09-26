# Environment manifest — chr20 300x benchmark (post-v1.0)

Date: 2026-09-23

## AI_DNA_ANALYZER version under test

- **Location:** `.claude/worktrees/opt-chr20` git worktree (branch `fix/hg005-extraction-integrity`)
  — NOT `main`. The frozen cascade/native-backend code that scores real BAM data
  currently lives only here (per project memory: native backend was never merged
  to `main`).
- **Commit:** `9c0933e4edc9bbdb2df64878e58e5d6101ce9d79` (2026-09-17 17:59:21 +0500,
  "FINAL test and validation")
- **Pipeline entrypoint:** `experimental/chr20_validation/run_chr20_pipeline.py`,
  fed by a native-backed count/LLR cache built via `cache_counts.py` /
  `cache_locus_evidence.py` (native C extension `native/pileup_native.c`,
  `native/reads_native.c`).
- **Frozen thresholds (verified read from `cascade.py`, never refit here):**
  - `FROZEN_BINOMIAL_THRESHOLD = 7.0`
  - `FROZEN_PB_THRESHOLD = 10.5`
  - `FROZEN_ROUTER_CUTOFF = 5.411872376933351`

## Toolchain

- gcc/g++: 15.3.1 20260722 (Red Hat 15.3.1-1)
- CMake: 3.31.11
- Python: 3.14.3
- samtools: 1.23.1
- bcftools: 1.23.1
- htslib: 1.23.1 (system libhts.so.3, via pysam's bundled headers)

## Competitor callers (Docker images, exact digests)

- **DeepVariant 1.10.0** — `google/deepvariant:1.10.0`
  digest `sha256:` (see `docker images --digests`), verified via
  `run_deepvariant --version` -> "DeepVariant version 1.10.0"
- **GATK 4.6.2.0** — `broadinstitute/gatk:4.6.2.0`
  digest `sha256:71b17ee42d149e8ec112603f5305c873ab60d93949ef8bb62a4fff85427f56fb`,
  verified via `gatk --version` -> "The Genome Analysis Toolkit (GATK) v4.6.2.0",
  HTSJDK 4.2.0, Picard 3.4.0
- **hap.py 0.3.15** — `quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0`
  digest `sha256:d63b963a6cb01b4830393b22369e7b91d298e4156dde353739e74e4cfa4f96d0`,
  verified via `hap.py --version` -> "Hap.py"
- **Clair3** — not pulled yet; per plan, only attempted if a short-read mode
  is genuinely applicable, else recorded as `NOT_APPLICABLE_FOR_THIS_INPUT`.

## Hardware / OS

- CPU: Intel(R) Core(TM) i7-6700 @ 3.40GHz, 8 logical cores (`nproc`)
- RAM: 31 GiB total
- Kernel: Linux 6.19.12-200.fc43.x86_64 (Fedora Linux 43)

## Container storage note (engineering detail, not part of the scientific result)

`/mnt/archive` (the benchmark HDD) is formatted **NTFS**, which cannot host
Docker/containerd's overlayfs snapshotter (layer extraction fails on Linux-only
path/xattr semantics). Fix: a 60 GB ext4 image file
(`/mnt/archive/AI_DNA_ANALYZER_benchmark/container_storage.img`) mounted via
loop device at `container_storage_ext4/`, with Docker's `data-root` and
containerd's `root` repointed there. All container images/layers/temp files
live inside this ext4 loopback filesystem, itself a file on the HDD — nothing
was written to the system root disk for this benchmark, and no existing data
(BAM/reference/truth, `/home`, `/var/lib/flatpak`) was modified or deleted.
This loop mount does not persist across reboot unless added to `/etc/fstab`
(intentionally left out — out of scope for this benchmark).

## Reference & truth set (already documented in dataset_manifest.txt, repeated here for the full benchmark record)

- Reference: `GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna`
  (md5 `a056c57649f3c9964c68aead3849bbf8`, matches official NCBI
  `md5checksums.txt`; 2580 contigs, matches BAM header `@SQ` count exactly;
  chr20 sequence byte-identical (sequence-only md5
  `b18e6c531b0bd70e949a7fc20859cb01`) to the project's pre-existing
  `data/reference/chr20_full.fa`)
- Truth: GIAB HG002 GRCh38 **NISTv4.2.1**
  (`HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` + `.tbi` +
  `_noinconsistent.bed` confident regions), downloaded from
  `ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/`.
  chr20 subset of the confident-regions BED is byte-identical to the
  project's pre-existing `data/giab_hg002_chr20/chr20_highconf.bed`.
- BAM: `HG002.GRCh38.300x_chr20.bam` (see `dataset_manifest.txt` in the
  sibling `HG002_GRCh38_chr20/` directory on the HDD for full validation).
  Note: this differs from the BAM previously used for other chr20 work in
  this project (`data/giab_hg002_chr20/chr20_full.bam`, 2x250bp NovoAlign,
  NIST_Illumina_2x250bps release) — per instructions, only the freshly
  validated 300x 2x148bp novoalign BAM is used as input to any caller in
  this benchmark.
