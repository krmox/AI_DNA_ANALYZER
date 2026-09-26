# Caller-specific input preprocessing (per benchmark rule: log separately, never silently)

## GATK HaplotypeCaller: added @RG read-group line

**Problem:** `HG002.GRCh38.300x_chr20.bam` (the official, validated, untouched input —
see `dataset_manifest.txt`) has **no `@RG` line** in its header. This was already
flagged as an open concern during the original download validation. GATK
HaplotypeCaller refuses to run without one:

```
java.lang.IllegalStateException: the sample list cannot be null or empty
	at ...GenotypingEngine.<init>...
	at ...HaplotypeCallerEngine.initializeActiveRegionEvaluationGenotyperEngine...
```//
(full trace in `logs/gatk_run.log`, first attempt, 2026-09-23 18:01)

**Fix applied:** `gatk AddOrReplaceReadGroups` on a **separate copy** of the BAM,
written to `HG002_GRCh38_chr20/gatk_input_with_rg/HG002.GRCh38.300x_chr20.withRG.bam`.
The original `HG002.GRCh38.300x_chr20.bam` used for AI_DNA_ANALYZER and
DeepVariant is never modified.

Exact command:
```
gatk AddOrReplaceReadGroups \
  -I HG002.GRCh38.300x_chr20.bam \
  -O gatk_input_with_rg/HG002.GRCh38.300x_chr20.withRG.bam \
  -RGID HG002_chr20_300x -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM HG002 \
  --CREATE_INDEX true
```

`-RGSM HG002` assigns the sample name from the official GIAB source directory
(dataset provenance), since it cannot be read from within the original BAM
(no @RG there to begin with — see the open concern already logged in
`dataset_manifest.txt` at download time). `-RGPL ILLUMINA` matches the
documented sequencing technology (HiSeq2500). `-RGID`/`-RGLB`/`-RGPU` are
arbitrary single-sample placeholders — GATK requires them to be present and
non-empty, but their specific values do not affect single-sample germline
calling.

**Scope of this change:** read-group metadata addition only. No read is
added, removed, remapped, or resorted; no base is altered. This is add-only
metadata required purely by the caller's implementation, not a change to
the evidence.
