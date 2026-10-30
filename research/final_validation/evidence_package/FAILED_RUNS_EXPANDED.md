# Failed / interrupted / skipped runs — expanded benchmark

| Item | What happened | Cause | Counted as |
|---|---|---|---|
| DeepVariant HG007 chr17 30× s101 | Started 2026-10-29 ~18:0x, killed at 18:29 | **Stopped by user decision** (too slow). Not a caller or harness failure. Partial output (call.log only) kept; partial downsampled BAM and DV temp files deleted. | UNASSESSED |
| DeepVariant HG001 chr17 30× s101 | Never started | same user decision | UNASSESSED |
| DeepVariant 10×/20× (optional) | Not run | optional in preregistration; not needed for any pass | not run |
| Whole-file sha256 of re-fetched HG001 BAMs | differs from first campaign | `@PG` header path (+22 B); identity verified by BED-mean depth instead (Addendum A) | documented deviation |
| Caller or harness failures among the 48 new three-caller runs | none (all rc 0, all hap.py OK) | – | – |
| DeepVariant HG007 chr19 and HG001 chr19 | completed OK | – | reported |
