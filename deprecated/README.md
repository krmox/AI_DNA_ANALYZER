# deprecated/ — abandoned neural-model line (kept as negative-result evidence)

The project began as a neural sequence model; on real reads every neural stage was null or negative
against the statistical baselines. Nothing in V1/V2/V2.x/R2/nindels uses these files.

- `neural/` — training, evaluation, sanity-check and summary scripts, read-level/LLR model variants,
  and their tests (`test_bimamba.py`, `test_read_level.py`). Run from the repository root, e.g.
  `python3 deprecated/neural/train_read_level.py …`; `_bootstrap.py` / `conftest.py` add the root to `sys.path`.
  Pre-existing failures (unchanged): `test_bimamba.py` (imports the old package name
  `bimamba_variant_caller`), one test in `test_read_level.py`.
- `checkpoints/` — model checkpoints (`readlevel`, `residual_binomial`, `rawpileup_llr`, `rawpileup_15x`);
  moved from `checkpoints_<name>/`. Nothing deleted.

Checkpoint directories that are still referenced by code kept at the root remain there:
`checkpoints/` (default output of `run_train.py`), `checkpoints_rawpileup/` (`train_raw_pileup.py`) and
`checkpoints_pb_residual/` (read by `bench_v12_stage23.py`, whose recorded output embeds that path).
Model files that other modules still import (`model.py`, `model_raw_pileup.py`, `model_pb_residual.py`,
`model_residual_binomial.py`, `run_train.py`, `train.py`, `train_raw_pileup.py`) also remain at the root.
