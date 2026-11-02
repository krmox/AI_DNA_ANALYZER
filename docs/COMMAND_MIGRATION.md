# Command migration guide (to R3.1)

R3.1 adds one front-end command, `varion`. **No old command was removed or changed**: the old binaries and
scripts keep their arguments and outputs. Where you used to assemble several steps by hand, `varion` now does them.
The only change to old binaries is additive (new subcommands/flags in `nindels-call`, one extra `STRB` VCF flag).

Install the front end with `pip install -e .` (no third-party dependencies) or run `python3 -m varion`.

## Indels

| Before (indel engine, branch `nindels-test`) | Now | Notes |
|---|---|---|
| `nindels-call estimate-str --bam B --ref R --region chr20 --out S.tsv` then `nindels-call call --bam B --ref R --region chr20 --force-reassess --genotyper multi --str-model S.tsv --out-vcf o.vcf` | `varion indel --bam B --ref R --out-dir out/ --region chr20 --mode A+B` | `--mode auto` (default) picks A or A+B from sequencing depth and logs the choice |
| `nindels-call call ... --genotyper multi` (A) | `varion indel ... --mode A` | |
| `nindels-call call ...` (default, "baseline", biallelic) | unchanged: `nindels-call call ...` | kept only as the frozen reference behaviour; not recommended (7x more false positives at 67x on HG004 chr20) |
| (none) | `varion depth --bam B --region chr20` | sampled-window depth profile (JSON) |
| (none) | `nindels-call depth`, `nindels-call contigs` | new engine subcommands used by `varion` |
| (none) | `nindels-call call ... --str-min-depth N` (needs `--str-model`) | per-locus gate: B only where the locus has >= N reads |
| (none) | `nindels-call call ... --summary-json J` | machine-readable counts |
| scripts `nindels/scripts/{run_regression,run_validation}.sh` | unchanged; superseded for new work by `nindels/scripts/r31/` | |

Outputs of `varion indel` (all in `--out-dir`): `<sample>.indels.vcf.gz` (+ `.tbi` when `bgzip`/`tabix` exist),
`<sample>.indels.summary.json` (schema `varion.indel.summary/1`), `<sample>.indels.log`. The VCF header records the
configuration chosen per contig (`##varion_config_<contig>=...`); loci that used the repeat model carry the `STRB` flag.

## SNPs

| Before | Now | Notes |
|---|---|---|
| `native_r2/build/dnav2 run --bam B --ref R --bed BED --out-dir D --contig chr20 --contig-len N --sample S --tag T --threads 8 --mode cascade --model research/stage4/frozen_R2.model --r2 --ood auto` | `varion snp --bam B --ref R --bed BED --out-dir D --contig chr20 --sample S` | the wrapper looks up `--contig-len` from the `.fai` and uses the frozen R2 model and the R2 flags above |
| `native_v2x/build/dnav2 run ...` (V2 / V2.x) | unchanged | frozen; reproduces V1 byte-for-byte without `--model` |
| `bash verify.sh [--data]` | unchanged | |
| V1 Python scripts (`bench_v*.py`, `extract_bench_v12.py`, ...) | unchanged | historical, see `docs/HISTORY_V1_V2_README.md` |

`varion snp` needs a regions BED (the engine scores only positions inside it) and the engine built with
`make -C native_r2 engine`. SNP and indel calling are separate runs and separate VCFs.

## Deprecated

Nothing is removed. Soft-deprecated (still works, no longer recommended for new analyses):

- calling indels with the default `nindels-call call` (baseline biallelic genotyper): use `varion indel`;
- hand-assembling `estimate-str` + `call --str-model` pipelines: use `varion indel --mode A+B`.
