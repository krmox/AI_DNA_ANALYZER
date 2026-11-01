## [2026-08-10] — Option 3 (gamma reduction): negative result on both seeds. Part B: recommend per-checkpoint threshold selection.

### Hypothesis:
FOCAL_GAMMA=3.0 rewards increasing confidence on already-correctly
classified tokens by construction, and is the most directly implicated
cause of the Normal-class confidence tail drifting upward across
epochs (diagnosed 2026-08-08, reproduced on clean data 2026-08-10).
Lowering it should slow or remove that drift.

### What was done:
- Added CLI overrides to run_train.py (--gamma, --seed, --max-epochs,
  --tag), all defaulting to prior hardcoded behavior — a bare
  `python run_train.py` invocation is unchanged.
- --tag additionally enables per-epoch checkpoint snapshots
  (checkpoints/epoch{N}_{tag}.pt) alongside the existing best+last
  files, fixing last session's blind spot where only 2 of 6 epochs
  survived to disk. Logging-only change, does not alter training.
- Gamma set to 1.0 (down from 3.0) — meaningfully below even the
  standard focal-loss default of 2.0, while staying a real focal loss
  (gamma>0) rather than degenerating to plain weighted CE at gamma=0,
  which would have confounded "less focusing" with "no focusing."
- Two seeds (20260810, 20260811), 6 epochs each, single variable
  (FOCAL_GAMMA only) — same discipline as the 2026-08-07 SNP-weight
  experiment. No other hyperparameter touched.

### Results (Facts):

Both seeds: AUC/AP = 1.0000 at every epoch — 18/18 checkpoints checked
across this session and the prior one never show ranking degradation.
SNP support = 381, consistent with the prior session.

**Seed A — F1@0.5 / Normal-tail q9999 / max, by epoch:**
| Epoch | F1@0.5 | q9999 | max |
|---|---|---|---|
| 1 | 0.9943 | 0.216 | 0.229 |
| 2 | 0.8816 | 0.082 | 0.095 |
| 3 | 0.6761 | 0.204 | 0.325 |
| 4 (peak) | 0.4438 | 0.835 | 0.918 |
| 5 | 0.5567 | 0.540 | 0.744 |
| 6 | 0.6167 | 0.536 | 0.696 |

**Seed B — same columns:**
| Epoch | F1@0.5 | q9999 | max |
|---|---|---|---|
| 1 | 0.9943 | 0.028 | 0.030 |
| 2 | 0.6684 | 0.105 | 0.115 |
| 3 | 0.5372 | 0.246 | 0.291 |
| 4 | 0.6268 | 0.511 | 0.705 |
| 5 (peak) | 0.3923 | 0.853 | 0.956 |
| 6 | 0.5625 | 0.706 | 0.936 |

Reference, prior session's gamma=3.0 (single unseeded run): epoch 1
max=0.051; epoch 6 q9999=0.747, max=0.821.

### Issues found:
1. **Gamma reduction does not slow the tail drift — negative on both
   seeds.** Both gamma=1.0 seeds reach worse-or-comparable tail
   extremes than the gamma=3.0 baseline, and reach them earlier
   (epoch 4-5 vs epoch 6). F1@0.5 is also noisier (non-monotonic)
   rather than calmer. Whatever drives Normal's confidence tail
   upward survives a 3x cut to the focal-loss focusing exponent —
   this rules out focusing strength as the sole or primary cause.
2. **Caveat on the comparison:** the gamma=3.0 reference point is
   single-seed/unseeded, not a matched 2-seed control at that gamma.
   A same-seed gamma=3.0 vs gamma=1.0 comparison would be a cleaner
   attribution check — not run this session (would be a second
   variable, out of scope for a single-variable session).
3. **AUC/AP have zero discriminating power at this task scale.**
   18/18 checkpoints checked so far score exactly 1.0000 regardless
   of gamma, seed, or epoch — useful as a leak/ranking-collapse
   diagnostic, but cannot select between checkpoints.

### Decision:
Option 3 (gamma-based recalibration) is **ruled out**, not fixed —
a genuine negative result, not an inconclusive one. FOCAL_GAMMA is
not the lever controlling tail drift.

**Part B (checkpoint-selection criterion) — recommended: (ii)
per-checkpoint optimal threshold.**
- (i) fixed 0.5: disproven three times now (prior session + both
  seeds this session) — measures threshold drift, not model quality.
- (iii) AUC/AP-based: disqualified by its own success — can't
  discriminate when every checkpoint scores 1.0000. Keep reporting
  it as a diagnostic (it would catch a real ranking collapse), but
  it can't be the selector.
- (ii) per-checkpoint optimal threshold: the only option with a
  verified positive result — the prior session's Option 2 sweep
  showed epoch 1 and epoch 6 reach identical F1 (0.9943) at their
  own best thresholds (0.50 vs 0.83). `threshold_selection.py`
  already implements this cheaply.

Recommendation only — **not implemented this session**. Wiring this
into checkpoint-selection logic is a distinct single-variable change
deserving its own isolated verification run.

### Open for next session:
- Implement and isolated-verify Part B's per-checkpoint
  optimal-threshold selection logic.
- Matched-seed gamma=3.0 vs gamma=1.0 comparison, if tighter
  attribution is wanted before trying a different lever.
- Next candidates for the actual drift-control lever: CLASS_WEIGHT_ALPHA
  / CLASS_WEIGHT_CAP, weight decay — each a separate single-variable
  session.
- Standing question: AUC/AP = 1.0000 on every one of 18 checkpoints
  regardless of configuration suggests this task, at this data scale
  (381 SNP support), may be too easily separable to meaningfully
  stress-test threshold/calibration logic. Worth considering whether
  future validation needs a larger or harder slice before trusting
  checkpoint-selection conclusions drawn here.
- All prior-week checkpoints (best_model_buggy_data_2026-08-07_
  phantom_era.pt, best_model.pt, last_checkpoint.pt, epoch1_halved_
  snp_weight_0.9246.pt, ...seed2_0.9524.pt) verified untouched
  (md5 + mtime) after this session's run.
