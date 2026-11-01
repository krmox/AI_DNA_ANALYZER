# tests/ — additional Python unit tests

The **verification suite** that `verify.sh` runs (`test_cascade`, `test_cheap_router`, `test_bench_v12`,
`test_robustness_benchmark`, the three native-equivalence tests and `test_hg005_extraction_integrity`)
stays at the repository root next to the modules it tests. This directory holds the remaining
unit tests of V1 components: `binomial_baseline`, `quality_error_model`, `locus_evidence`, `providers`,
data loading, mock BAM/VCF generators and `model_residual_binomial`.

Run from the repository root: `python3 -m pytest tests`. `conftest.py` puts the root on `sys.path`.
Tests for the benchmark rounds live in `benchmarks/`, tests for abandoned neural code in `deprecated/neural/`.
Known pre-existing failures (unchanged by the reorganisation): `test_providers.py` (collection error).
