# Changelog

Research software; no stability guarantees. Versions before 3.1.0 had no changelog; the entries below for them are
reconstructed from git tags and are deliberately brief.

## 3.1.0 — tag `3.1.0` @ `f374301`, 2026-11-03 (Zenodo 10.5281/zenodo.23856530)

Built from public `main` @ `8297691` (which is `3.0.0` @ `910685b` plus documentation/layout commits) by re-applying the
R3.1 implementation (originally developed on local branch `r3.1`, commit `31e14df`) and two corrective patches.
The `r3.1` branch itself was not published (it differs from `main` in 519 files, mostly unrelated).

### Included

- **`varion` front end** (`varion/`, stdlib-only Python): `varion indel`, `varion depth`, `varion snp`, `varion doctor`.
  Old commands are unchanged (`docs/COMMAND_MIGRATION.md`).
- **`nindels` R3.1 engine changes:** per-locus repeat-model gate (`--str-min-depth`), `depth` and `contigs`
  subcommands, `--summary-json`, `STRB` VCF flag.
- **Depth-dependent configuration choice** (`varion/adaptive.py`, `varion/thresholds.json`): contig depth < 67.37x -> A,
  >= 67.37x -> A+B, no per-locus gate. `--mode A` / `--mode A+B` override it.
- **R3.1 protocol and selection record:** `docs/R3.1_PROTOCOL.md` (pre-registered, unmodified),
  `results/r31/selection.json` (unmodified), experiment scripts `nindels/scripts/r31/`,
  `docs/R3.1_THRESHOLD_SELECTION.md` (new, see below), `docs/ARCHITECTURE.md`, `docs/COMMAND_MIGRATION.md`.
- **Two `nindels` fixes** written after a code audit (not present on the original `r3.1` branch):
  1. multi-allelic VCF records are now emitted in minimal form (previously 746 of 27,106 records on HG001 chr19 30x
     were non-minimal; `bcftools norm` rewrote them; calls are unchanged after normalisation: 30,160 split records,
     0 differences);
  2. an invalid `--region` (reversed, start 0, beyond contig length) is an error with non-zero exit status instead of an
     empty VCF with exit 0 (regression test `nindels/tests/test_cli_region.sh`).
- Tests: `nindels-tests` 160 (measured: 153 on `main`, 156 on `r3.1` @ `31e14df`, 160 here; the +4 comes with fix 1); `tests_varion` 22.
- Metadata: README, `CITATION.cff` (version 3.1.0, concept DOI), `docs/INDEX.md`, this file.

### Not included

- The R3.1 **held-out test results**: `docs/R3.1_FINAL_TEST_REPORT.md` and `results/r31/test_report*`,
  `dev_report*` are not part of this release. No statement about R3.1 on held-out data is made here.
- Raw data, per-run VCFs and hap.py outputs of any experiment (archived separately as a Zenodo dataset, supplementary to this release: DOI 10.5281/zenodo.23938931, https://zenodo.org/records/23938931).
- Any change to the R2 SNP engine, to the thresholds, to the selection output or to the protocol.
- The Zenodo record of version 3.0.0 (DOI 10.5281/zenodo.23855773) does not contain R3.1 or `varion`; the 3.1.0
  release is archived separately as DOI 10.5281/zenodo.23856530.

### Known issues and honesty notes

- **Protocol deviation.** The frozen thresholds do not follow protocol step 6 (no local-depth bin qualified, so the
  protocol prescribes A at every depth; the freeze script checked only the global threshold). The frozen A+B rule is
  not protocol-compliant and is weakly supported. Details, data and limits: `docs/R3.1_THRESHOLD_SELECTION.md`.
  `--mode A` is the protocol-compliant behaviour; `--mode auto` (default) is not.
- A claimed beforehand check that PASS-only scoring reproduces earlier numbers (protocol amendment 1) has no record
  in the repository.
- `nindels/scripts/r31/datasets.json` and several scripts hard-code machine-specific paths.
- `verify.sh` does not run `tests_varion` or the CLI region test; run them as shown in `README.md`.
- Outputs were checked with htslib 1.21; earlier published results used 1.23.x (identical on the cell tested here).
- Reproduction check (HG001 chr19 30x, configuration A): the original `r3.1` build reproduces the recorded VCF body
  exactly (27,106 records, 0 differing lines); the 3.1.0 build differs in the 746 multi-allelic records only by
  representation (fix 1), with identical calls after normalisation.

## 3.0.0 — tag `3.0.0` @ `910685b`, 2026-11-01 (Zenodo 10.5281/zenodo.23855773)

Portability work (shared htslib resolution, grouped `verify.sh`, conda/Docker recipes), R1.1/R2 SNP engine, `nindels`
Baseline/A/B/A+B. Does not contain R3.1 or `varion`.

## 2.0.0 — tag `v2.0.0`, 2026-09-26

V2 native engine, V2.x accuracy layer.

## 1.0.0 — tag `v1.0.0`, 2026-09-22

V1 publication package.
