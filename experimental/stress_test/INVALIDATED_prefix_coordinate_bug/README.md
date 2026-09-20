# INVALIDATED — do not cite as results

Everything in this directory was produced by the HG005 chr1:1–25 Mb run made **before** the coordinate fix
(`positions = chunk_start + arange(n_rows)`, see manuscript §3.9 C1). Coordinates written to the VCFs were
wrong; hap.py therefore scored F1 = 0.0085 (TRUTH.TOTAL 28,435). That number is an artefact, **not a
performance measurement**.

Moved here on 2026-09-20 (originally `experimental/stress_test/{happy,ai_*_C.vcf*}` and
`HG005_MAX_FEASIBLE_STRESS_TEST.*`) so it cannot be mistaken for the post-fix result. Nothing was edited or
deleted; contents are byte-identical. The valid post-fix HG005 results (hap.py F1 0.9521 / 0.9519 on
chr1:1–4 Mb) are recorded in `research/` (manuscript Table 3, `CLAIM_EVIDENCE_MAP.csv`); their raw post-fix
artefacts (cache, region-scoped hap.py output) are NOT present in this repository (lost with a wiped
temporary worktree) — see manuscript "Data and code availability" and `research/FINAL_FREEZE_REPORT.md`.
