# Commands — Head-to-Head Benchmark

All work was done in an isolated git worktree (`git worktree add /tmp/head_to_head_work HEAD`,
detached at `16b221ee82d64f66de25119aafa4e5944cfe63cc`) with `data/` symlinked in (not copied) to
avoid duplicating ~11GB of gitignored BAM/reference data. Nothing was committed. Final
deliverables were copied into `experimental/head_to_head/` in the main working tree afterward
(still uncommitted).

## 0. Reference contig-name adapter (required preprocessing, not a cascade change)

The project's `data/reference/chr21_full.fa` uses Ensembl-style contig naming (`>21`), while the
GIAB HG002 BAM/truth/BED all use UCSC-style (`chr21`). `providers.py:_normalise_contig` already
handles this transparently for the frozen cascade's own internal use, but external tools
(DeepVariant/Clair3/GATK) require an exact name match. A renamed COPY (same sequence, header
changed) was created for external-tool use only; the production reference file was never touched.

```bash
sed 's/^>21 />chr21 /' data/reference/chr21_full.fa > shared_ref/chr21_chrname.fa
samtools faidx shared_ref/chr21_chrname.fa
gatk CreateSequenceDictionary -R shared_ref/chr21_chrname.fa
```

## 1. DeepVariant

```bash
docker pull google/deepvariant:1.6.1
docker run --rm \
  -v <shared_ref>:/ref -v <data/giab_hg002_chr21_12Mb>:/bam -v <out>:/out \
  google/deepvariant:1.6.1 \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WGS \
  --ref=/ref/chr21_chrname.fa \
  --reads=/bam/hg002_chr21_32_44M_15x.bam \
  --regions="chr21:32000000-44000000" \
  --output_vcf=/out/deepvariant.vcf.gz \
  --sample_name=HG002 --num_shards=8 --logging_dir=/out/logs
docker rmi google/deepvariant:1.6.1   # cleanup per sequential-disk protocol
```

## 2. Clair3

First attempt (`hkubal/clair3:latest`) hit `AttributeError: No symbol free_memory found in library
/opt/bin/preprocess/realign/realigner` during the full-alignment realignment step and silently
fell back to pileup-only output (documented, then discarded). Retried with a pinned version:

```bash
docker pull hkubal/clair3:v1.0.10
docker run --rm \
  -v <shared_ref>:/ref -v <data/giab_hg002_chr21_12Mb>:/bam -v <out>:/out \
  hkubal/clair3:v1.0.10 \
  /opt/bin/run_clair3.sh \
  --bam_fn=/bam/hg002_chr21_32_44M_15x.bam \
  --ref_fn=/ref/chr21_chrname.fa \
  --threads=4 --platform=ilmn --model_path=/opt/models/ilmn \
  --ctg_name=chr21 --sample_name=HG002 --output=/out
docker rmi hkubal/clair3:v1.0.10   # cleanup
```
Model selected before running: `ilmn` (bundled Illumina short-read model), the only match for our
platform among {hifi, hifi_revio, hifi_sequel2, ilmn, ont, ont_guppy5, r1041_*}.

## 3. GATK HaplotypeCaller

```bash
curl -sL -o gatk.zip https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
unzip gatk.zip   # static release, no Docker; Java 25 (OpenJDK) already present, worked without issue
gatk/gatk AddOrReplaceReadGroups -I hg002_chr21_32_44M_15x.bam -O hg002_rg.bam \
  -RGID HG002 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM HG002   # original BAM had no @RG at all
samtools index hg002_rg.bam
gatk/gatk HaplotypeCaller -R shared_ref/chr21_chrname.fa -I hg002_rg.bam \
  -L chr21:32000000-44000000 -O gatk_hc.vcf.gz
```
No BQSR applied — no known-sites resource (dbSNP/Mills indels) was on hand for a chr21-only
subset, and downloading one was judged not justified for a single-region experiment. Documented
protocol deviation from GATK best practices, not something applied silently.

## 4. Strelka2 — BLOCKED

```bash
curl -sL -o strelka2.tar.bz2 https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
tar xjf strelka2.tar.bz2
pyenv install 2.7.18   # bundled workflow scripts require python2, absent on this host
strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
  --bam <bam> --referenceFasta <ref> --region chr21:32000000-44000000 --runDir <out>
# -> CONFIGURATION ERROR: htsfile -h <bam> fails:
#    "bgzf.c:347: bgzf_hopen: Assertion `compressBound(0xff00) < 0x10000' failed."
```
Bundled `libexec/htsfile` (2018 CentOS6 static build) is binary-incompatible with this host's
zlib. No system `htsfile`/rtg-tools alternative was available to substitute. Not patched/rebuilt
from source (out of scope). Marked `BLOCKED — TOOLING INCOMPATIBILITY`.

## 5. AI Cascade / PB-only evidence extraction (frozen code, unmodified)

```bash
python3 extract_bench_v12.py \
  --fasta data/reference/chr21_full.fa \
  --bam data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_15x.bam \
  --vcf data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz \
  --bed data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed \
  --region 32000000 44000000 --features none \
  --out cache_h2h/hg002_32_44M_15x.npz
```
This is the same script/protocol used for the pre-registered bench_v12 study
(`History/12_DEVLOG.md` §17), just pointed at the full 32-44Mb span instead of a 1Mb sub-slice.
Frozen constants (`FROZEN_BINOMIAL_THRESHOLD`, `FROZEN_PB_THRESHOLD`, `FROZEN_ROUTER_CUTOFF`) are
imported from `cascade.py`, never redefined. Router/accuracy logic reused from
`robustness_benchmark.py` (`route_mask`, `accuracy_block`) exactly as `bench_v12_stage1.py` does —
see `analyze_cascade.py` in this directory for the exact reuse.

## 6. hap.py (common evaluator for external callers)

```bash
docker pull quay.io/biocontainers/hap.py:0.3.7--py27_1     # first attempt
# -> vcfeval engine: rtg-tools binary not bundled in this build (path guess failed, and
#    `find / -iname '*rtg*'` inside the image found nothing)
# -> default xcmp engine: quantify step failed with `regex_error` (old Boost-regex-era binary
#    bug, reproducible even with LC_ALL=C/LANG=C)
docker rmi quay.io/biocontainers/hap.py:0.3.7--py27_1
docker pull quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0   # newer tag, same xcmp engine
hap.py <truth.vcf.gz> <query.vcf.gz> -f <highconf.bed> -r shared_ref/chr21_chrname.fa -o <out>/happy
docker rmi quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0    # cleanup
```
Run once per external caller (DeepVariant, Clair3, GATK), same truth/BED/reference/region each
time. This is the GA4GH-standard evaluator requested as first priority; it worked once the correct
tag was found, so the "write our own adapter" fallback was not needed for the external-caller
comparisons. AI-Cascade/PB-only were NOT run through hap.py (see `ERROR_ANALYSIS.csv`,
`CAVEAT_ON_CROSS_METHOD_COMPARISON`) — they were evaluated with the pre-existing frozen internal
locus evaluator instead, since they are LLR-threshold decisions over a scoring frame, not a
genotyped VCF, and synthesizing one was judged out of scope.

## 7. AI Cascade rescue/accuracy analysis

See `analyze_cascade.py` (in this directory) — reuses `route_mask`/`accuracy_block` from
`robustness_benchmark.py` and `paired_bootstrap` from `evaluate_quality_error.py`, all frozen,
unmodified, imported not reimplemented.

## Sequential disk protocol log

| Step | Free space before | Action | Free space after |
|---|---|---|---|
| baseline | 13GB | — | 13GB |
| DeepVariant pull | 13GB | `docker pull google/deepvariant:1.6.1` (3.2GB compressed) | 5.8GB |
| DeepVariant run+save+cleanup | 5.8GB | ran, saved VCF, `docker rmi` | 13GB |
| Clair3 pull (:latest) | 13GB | `docker pull hkubal/clair3:latest` | 9.3GB |
| Clair3 :latest run (failed) + cleanup | 9.3GB | ran (bug), saved partial output, `docker rmi` | 12GB |
| Clair3 pull (v1.0.10) | 12GB | `docker pull hkubal/clair3:v1.0.10` | 9.8GB |
| Clair3 v1.0.10 run+save+cleanup | 9.8GB | ran (success), saved VCF, `docker rmi` | 13GB |
| Strelka2 static binary | 13GB | 31MB tarball (tmpfs `/tmp`, not `/home`) | 13GB (unaffected) |
| GATK static zip | 13GB | 741MB zip (tmpfs `/tmp`, not `/home`) | 13GB (unaffected) |
| hap.py pull (0.3.7, then 0.3.15) | 12GB | ~0.4GB each, sequential, both removed after use | 12GB |

Note: Strelka2/GATK/hap.py tool downloads and the cascade extraction's working files were placed
under `/tmp` (tmpfs, RAM-backed, 16GB cap, separate from the `/home` 94GB disk this whole protocol
is about), so they did not compete with the `/home` disk budget at all — only the three Docker
image pulls (DeepVariant, Clair3 ×2, hap.py ×2) touched `/home`, and each was removed immediately
after its run per the sequential protocol above.
