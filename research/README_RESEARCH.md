> **SUPERSEDED IN PART (2026-09-20).** Index note: `PAPER_MANUSCRIPT.md` is the current manuscript; `FINAL_FORENSIC_AUDIT.md` and `FINAL_FREEZE_REPORT.md` describe the 2026-09-20 audit and freeze; other files may contain superseded statements (see `LIMITATIONS_AND_OPEN_QUESTIONS.md` section D).

# AI_DNA_ANALYZER — Research Dossier

Independent research project investigating whether a computationally cheap
Binomial screening stage, gated by a confidence router, can approximate a
more expensive Poisson-Binomial SNP caller's accuracy on GIAB short-read
data — and, downstream, whether a closed-form genotype posterior (Method C)
recovers standards-compliant (hap.py) genotype accuracy without retraining
or tuning.

## How to read this dossier

Start with `RESEARCH_SUMMARY_1PAGE.md` for the shortest honest overview.
Then:

- `METHODS.md` — architecture, evaluation methodology, statistical protocol.
- `RESULTS.md` — consolidated findings, cross-referencing the CSVs below.
- `SCIENTIFIC_CLAIMS.md` — **read this before quoting any number.** Every
  claim this project could plausibly make is classified
  SUPPORTED/PARTIALLY SUPPORTED/NOT SUPPORTED/NOT YET TESTED.
- `LIMITATIONS.md` — explicit, not minimized.
- `PIPELINE_INTEGRITY_AUDIT.md` — a read-only engineering audit, including
  one CRITICAL finding (fixed this session, the HG005 coordinate bug) and a
  self-correction: an earlier pass of this dossier wrongly marked the "~65x
  C/htslib speedup" as unsubstantiated (a real implementation existed in a
  worktree that was not checked); it has since been independently
  re-verified and integrated into the production extraction path.
- `REPRODUCIBILITY.md` — exact versions, checksums, commands.
- `PEER_REVIEW_SELF_AUDIT.md` — an adversarial self-review against 18
  standard reviewer questions.

Data: `BENCHMARK_MASTER.csv` (every experiment, one row per arm),
`ERROR_ANALYSIS.csv`, `RESCUE_ANALYSIS.csv`, `PERFORMANCE_RESULTS.csv`,
`EXPERIMENT_MANIFEST.json` (machine-readable index).

Narrative/outreach material: `PAPER_DRAFT.md`, `OUTREACH_PACKAGE.md`,
`AI_ASSISTANCE_DISCLOSURE.md`.

## One-sentence summary

A frozen Binomial->Router->Poisson-Binomial SNP-calling cascade with a
closed-form genotype posterior layer preserves PB-only accuracy (ΔF1 within
noise across every tested experiment) while routing under 0.1% of loci to
the expensive arm, validated on GIAB HG002/HG003/HG004 (internal evaluator,
9 pre-registered region/depth/sample cells) and HG002/HG005 (hap.py,
genotype-aware, 2 regions) — with a dominant, reproduced weakness in
low-mappability/segmental-duplication regions and no whole-genome or
cross-platform validation yet performed.

## Status

`RESEARCH PACKAGE STATUS: PARTIALLY COMPLETE` — see the final section of
`PEER_REVIEW_SELF_AUDIT.md` for the exact breakdown of what is and is not
ready.
