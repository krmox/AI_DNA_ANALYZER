# Phase 2 gate — blocked (2026-10-29)

Requirement: a GIAB genome never used for development, tuning, debugging, feature selection, or previous evaluation, with GRCh38 reads and a germline benchmark.

| Candidate | Why it fails |
|---|---|
| HG002, HG003, HG004 | development/tuning/R2 TRAIN+VALID |
| HG005 | debugged on (HG005 extraction bug); used in stress test |
| HG006 | R2 blind-new run on chr16 |
| HG007, HG001 | evaluated in this session's campaigns (independent of *development*, but not unseen) |
| HG008 | `data_somatic/HG008/NIST/` lists only tumor bulk/clones; no germline benchmark found on the GIAB FTP |
| HG009 | tumor/normal, no germline v4.2.1 GRCh38 truth |

The GIAB FTP `release/` tree has only AshkenazimTrio, ChineseTrio and NA12878_HG001 with germline benchmarks. No eligible sample exists on that site.

Decision (user, 2026-10-29): Phase 1 only, then stop. No preregistration file was written and no new benchmark was run.
Smallest unblocking steps: (a) relax the rule to "sample-seen, region-unseen" (HG005/HG006 on untouched chromosomes), (b) treat HG001/HG007 as a replication with a fresh preregistration, or (c) wait for a GIAB germline benchmark of HG008.
