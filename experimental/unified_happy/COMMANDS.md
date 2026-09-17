# Commands — Unified hap.py Benchmark

All commands run directly in the main checkout (not a worktree), per task
instructions, since `experimental/`/`data/` are untracked and not
worktree-materializable. No tracked file was modified.

## 1. Sanity check (Step 1)

```bash
python3 experimental/unified_happy/sanity_check.py
```
See `SANITY_CHECK.md` for full output and verdict (PASS).

## 2. Recover true positions across the full region (Step 2)

```bash
python3 experimental/unified_happy/recover_full_positions.py
# -> experimental/unified_happy/recovered_full.npz
# 24 chunks x 500kb, chr21:32,000,000-44,000,000, ~1242s wall clock
```

## 3. Build VCFs from frozen calls + recovered positions (Step 2)

```bash
python3 experimental/unified_happy/build_vcfs.py
# consistency check: recomputed frame-filtered labels vs
# experimental/head_to_head/ai_cascade_per_locus.npz's stored labels
#   -> label agreement: 1.000000 (0 mismatches out of 10975654)
#   -> snp array exact match: True
# -> ai_pb_only.vcf (16094 records), ai_cascade.vcf (16097 records)
```

## 4. Normalize + hap.py (Steps 3-4)

```bash
bash experimental/unified_happy/run_happy.sh
```
Internally:
```bash
bcftools sort -Oz -o ai_pb_only.vcf.gz ai_pb_only.vcf
bcftools index -t ai_pb_only.vcf.gz
bcftools sort -Oz -o ai_cascade.vcf.gz ai_cascade.vcf
bcftools index -t ai_cascade.vcf.gz
# (No bcftools norm applied -- matches experimental/head_to_head/COMMANDS.md
# section 6, which also applied no norm step to DeepVariant/Clair3/GATK VCFs
# before hap.py; `tabix` binary is absent on this host, `bcftools index -t`
# produces an equivalent .tbi index.)

docker pull quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0
# Digest: sha256:d63b963a6cb01b4830393b22369e7b91d298e4156dde353739e74e4cfa4f96d0
# (confirmed identical to experimental/head_to_head/COMMANDS.md section 6)

docker run --rm \
  -v <data/giab_hg002_chr21_12Mb>:/truth \
  -v <experimental/unified_happy>:/query \
  -v <experimental/head_to_head/shared_ref>:/ref \
  -v <experimental/unified_happy/happy/ai_pb_only>:/out \
  quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0 \
  hap.py /truth/hg002_chr21_32_44M.vcf.gz /query/ai_pb_only.vcf.gz \
  -f /truth/hg002_chr21_32_44M_highconf.bed -r /ref/chr21_chrname.fa \
  -o /out/happy
# (repeated identically for ai_cascade.vcf.gz)

docker rmi quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0   # cleanup
```
Mount layout confirmed by inspecting `experimental/head_to_head/happy/deepvariant/happy.runinfo.json`'s
`final_args` (`vcf1=/truth/...`, `vcf2=/query/deepvariant/...`, `ref=/ref/...`,
`reports_prefix=/out/happy`) so the invocation structure is identical to the
one already used for the three external callers.

## Sequential disk protocol log (this pass)

| Step | Free space before | Action | Free space after |
|---|---|---|---|
| baseline | 12GB | — | 12GB |
| hap.py pull (0.3.15) | 12GB | `docker pull` (~0.4GB) | 12GB (rounds to same GB) |
| both hap.py runs | 12GB | ran, saved outputs, checksummed | 12GB |
| `docker rmi` cleanup | 12GB | removed image | 12GB |

No DeepVariant/Clair3 images were pulled (external caller results reused
unchanged from `experimental/head_to_head/`).

## 5. Runtime, unified table, disagreements (Steps 5-7)

See `runtime.csv`, `UNIFIED_HAPPY_BENCHMARK.csv`/`.md`, `DISAGREEMENTS.csv`.
