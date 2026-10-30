# Failed or skipped runs

No failed runs recorded in RESULTS.jsonl.

Skipped by design (see REPORT.md): native 300x depth (impractical: ~3x the 100x cost for every caller and hap.py; the protocol targets 10–30x); 15x and 50x depths (PREREG P2/P4 depths) not in the requested matrix.

## Driver-side problems during the campaign (no caller or model failure)
All were bugs in `run_campaign.py` (my driver), found on the first condition(s), never in the frozen caller; affected rows were archived (not deleted) and the conditions re-run from scratch with the fixed driver:
* `FAILED_RUNS_driver_bug_archived.jsonl` — 3 rows (chr19 10x s101): `HAPPY_FAILED: 'QUERY.TP'` (hap.py summary has no QUERY.TP column; QTP is now QUERY.TOTAL−FP−UNK).
* `FAILED_RUNS_driver_bug_archived2.jsonl` — 3 rows: the first fix accidentally commented out the TP/FN/FP assignments; all rows discarded and the campaign restarted before any other condition was kept.
* Machine restart during chr17 (after the 10x s102 condition): the half-written downsampled BAM `HG001_chr17_10x_s103.bam` was deleted and the condition re-run (`logs/matrix_chr17_resume.log`). No result from the interrupted condition was kept.
* An initial fetch attempt used a wrong index filename (script typo) and produced an empty 24 KB BAM; it was deleted and the fetch re-run (no data used).
Note: timing of chr19 conditions overlapped with the chr17 download (disk/network contention); chr17 conditions ran with no concurrent download.
