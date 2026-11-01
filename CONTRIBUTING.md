# Contributing

Varion is a frozen research record plus a small amount of active code (`nindels/`, `native_r2/`).
Scientific results, thresholds and anything under `research/`, `results/`, `experimental/` are
evidence: do not edit them to "tidy up" — add new files instead.

1. **Set up** — follow *Installation* in `README.md` (conda env or venv + system htslib 1.23.x).
   `bash scripts/hts_env.sh report` shows which htslib is used.
2. **Verify before and after your change** — `bash verify.sh full` must end with
   `COMPLETE VERIFICATION PASSED`. Extra Python tests: `python3 -m pytest tests benchmarks --continue-on-collection-errors` (`tests/test_providers.py` fails to collect, as before the reorganisation).
3. **Run scripts from the repository root.** Data and result paths are relative to it.
4. **Moving files** — V1 modules at the root are imported by name by frozen scripts; do not move them.
   If you move anything else, update every reference and add a row to `docs/PATH_MAP.md`.
5. **New experiments** go in a new directory with a README stating the exact command, data identifiers and
   commit; keep development, validation and held-out data apart.
6. **Never change** frozen constants, benchmark outputs or test assertions to make a result pass.
7. Layout and status of each directory: `docs/INDEX.md`.
