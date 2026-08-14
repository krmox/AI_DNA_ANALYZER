"""A pre-router that decides where expensive calling is worth paying for.

The problem this exists to solve
--------------------------------
``History/12_DEVLOG.md`` established that routing 3% of loci to the Mamba
reproduces the full neural caller, but saves only 0.44% of end-to-end runtime,
because the router read the Poisson-binomial LLR and the PB prior *is* the
expensive component:

    Poisson-binomial prior     8,025 loci/s   (CPU, 48-step exact convolution)
    fixed-epsilon binomial 1,156,000 loci/s   (CPU, closed form)
    PB-residual Mamba      1,140,000 loci/s   (GPU)

A router that consumes the PB LLR cannot avoid computing the PB LLR. Anything
built here must therefore decide *before* PB runs, from quantities the pipeline
already has.

The cheap-information boundary, stated precisely
------------------------------------------------
Every feature in this module is a function of the ``[N, 10]`` pileup count
matrix (``pileup_counts.COUNT_COLUMNS``) and nothing else. That matrix is
produced by the single pass over the BAM that *all* arms already pay for, so
these features are free in the marginal sense that matters: no second pass, no
read tensor, no per-read arithmetic.

What that boundary excludes is as important as what it includes. Strand
balance, read-position summaries and the reference-versus-alternate quality
split -- the features the 44-channel representation added and experiment 8
proved decisive -- all require the read tensor and cost 73,270 loci/s to build.
They are 9x cheaper than PB but not free, and they are deliberately outside the
counts-only router. Only ``quality_sum`` and ``mapping_sum`` survive, as
depth-pooled means; these are exactly the pooled summaries ``locus_evidence``
criticised as too coarse for *calling*, which is not an objection to using them
for *triage*.

One optional feature crosses a second boundary and is kept separable: the
fixed-epsilon binomial LLR. It is a closed-form function of the same counts at
1.16 M loci/s, i.e. 144x cheaper than PB, so a router that reads it is still a
genuine pre-router -- but it is a *statistic*, not raw evidence, so routers with
and without it are reported separately rather than merged.

Training target: distillation, not supervision
----------------------------------------------
The router is trained to predict which loci the PB-margin oracle would route,
computed on the **training split only**. It never sees a truth label, a VCF, or
a variant call -- not even on train. Two consequences:

* the leakage audit is trivial: the router cannot leak labels it was never
  shown;
* the router's ceiling is explicitly the oracle it distils, which is the right
  framing -- the question is how much of the oracle's routing quality survives
  losing access to the oracle's input.

PB is computed on the training split once, offline, to build that target. That
is a one-off training cost and appears in no inference path.

Leakage boundary at inference
-----------------------------
``score`` takes a count matrix and returns an uncertainty. No PB LLR, no Mamba
output, no label, no position, no test-derived quantile. Enforced by
``test_cheap_router``.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np

from binomial_baseline import BinomialVariantCaller
from quality_error_model import candidate_alt

#: Counts-only feature layout. Order is fixed and asserted by the unit tests.
CHEAP_FEATURE_NAMES: tuple[str, ...] = (
    "k", "n", "ref_count", "other_count", "log_depth",
    "vaf", "ref_fraction", "other_fraction", "gap_fraction", "insertion_fraction",
    "mean_base_quality", "mean_mapping_quality", "k_small", "has_alt",
)

#: Normalisation constants, matching ``locus_evidence`` so the two
#: representations are on comparable scales.
DEPTH_SCALE = 4.0
PHRED_SCALE = 40.0
MAPPING_QUALITY_SCALE = 60.0
K_SCALE = 10.0

#: Fraction of the training split the distillation target marks as "the oracle
#: would route this". Chosen a priori as a round number comfortably above the
#: operating points of interest, so the router learns the shape of the
#: uncertain region rather than one particular coverage.
ORACLE_TARGET_FRACTION = 0.05


def cheap_features(counts: np.ndarray) -> np.ndarray:
    """Build the counts-only feature matrix.

    Every column is O(1) per locus over the already-materialized count matrix:
    no read is touched, no convolution is run.

    Args:
        counts: ``[N, 10]`` count matrix following
            ``pileup_counts.COUNT_COLUMNS``.

    Returns:
        ``[N, len(CHEAP_FEATURE_NAMES)]`` float32 features.
    """
    counts = np.asarray(counts, dtype=np.float64)
    base_counts = counts[:, 0:4]
    reference_index = counts[:, 9].astype(int)
    depth = counts[:, 6]
    gap = counts[:, 4]
    insertion = counts[:, 5]

    _, k, n = candidate_alt(base_counts, reference_index)
    k = k.astype(np.float64)
    n = n.astype(np.float64)
    safe_n = np.maximum(n, 1.0)
    safe_depth = np.maximum(depth, 1.0)

    # The reference-supporting count is the reference base's own column, except
    # at reference N where no column is the reference.
    reference_count = np.where(
        reference_index > 0,
        np.take_along_axis(base_counts, np.maximum(reference_index - 1, 0)[:, None],
                           axis=1).reshape(-1),
        0.0)
    other_count = np.maximum(n - k - reference_count, 0.0)

    columns = {
        "k": k / K_SCALE,
        "n": np.log1p(n) / DEPTH_SCALE,
        "ref_count": np.log1p(reference_count) / DEPTH_SCALE,
        "other_count": np.log1p(other_count) / DEPTH_SCALE,
        "log_depth": np.log1p(depth) / DEPTH_SCALE,
        "vaf": np.where(n > 0, k / safe_n, 0.0),
        "ref_fraction": np.where(n > 0, reference_count / safe_n, 0.0),
        "other_fraction": np.where(n > 0, other_count / safe_n, 0.0),
        "gap_fraction": np.where(depth > 0, gap / safe_depth, 0.0),
        "insertion_fraction": np.where(depth > 0, insertion / safe_depth, 0.0),
        # Depth-pooled quality summaries: free, because the pileup pass already
        # accumulated the sums.
        "mean_base_quality": np.where(depth > 0, counts[:, 7] / safe_depth, 0.0) / PHRED_SCALE,
        "mean_mapping_quality": np.where(depth > 0, counts[:, 8] / safe_depth, 0.0)
        / MAPPING_QUALITY_SCALE,
        "k_small": np.minimum(k, K_SCALE) / K_SCALE,
        "has_alt": (k > 0).astype(np.float64),
    }

    features = np.empty((counts.shape[0], len(CHEAP_FEATURE_NAMES)), dtype=np.float32)
    for position, name in enumerate(CHEAP_FEATURE_NAMES):
        features[:, position] = columns[name]
    if not np.all(np.isfinite(features)):
        raise ValueError("cheap features contain a non-finite value")
    return features


def binomial_llr(counts: np.ndarray,
                 caller: BinomialVariantCaller | None = None) -> np.ndarray:
    """The fixed-epsilon binomial LLR: the cheap approximation to PB.

    Args:
        counts: ``[N, 10]`` count matrix.
        caller: Configured caller; the project's v1 default if omitted.

    Returns:
        ``[N]`` log-likelihood ratios.
    """
    caller = caller or BinomialVariantCaller(error_rate=0.01)
    return caller.score_counts(counts[:, 0:4], counts[:, 9].astype(int)).llr


# ------------------------------------------------------------- analytical ----

def binomial_margin_uncertainty(counts: np.ndarray, binomial_threshold: float,
                                llr: np.ndarray | None = None) -> np.ndarray:
    """Trivial rule A: distance from the *cheap* caller's decision boundary.

    The pre-PB analogue of the winning post-PB router. If this works, no
    learned model is needed and the honest recommendation is to use it.

    Args:
        counts: ``[N, 10]`` count matrix.
        binomial_threshold: The binomial caller's validation-frozen threshold.
        llr: Precomputed binomial LLR, to avoid recomputing it.

    Returns:
        ``[N]`` uncertainty, higher meaning closer to the cheap boundary.
    """
    values = binomial_llr(counts) if llr is None else np.asarray(llr)
    return -np.abs(values - float(binomial_threshold))


def alt_count_uncertainty(counts: np.ndarray) -> np.ndarray:
    """Trivial rule B: loci with a small non-zero alternate count.

    The crudest possible triage, and the one a bioinformatician would write by
    hand: k = 0 is decided, k >= half the depth is decided, and everything
    between is worth a second look. Expressed as ``-|vaf - 0.5| `` restricted to
    loci carrying at least one alternate read, so it ranks rather than merely
    partitions.

    Args:
        counts: ``[N, 10]`` count matrix.

    Returns:
        ``[N]`` uncertainty, higher meaning more ambiguous alternate evidence.
    """
    counts = np.asarray(counts, dtype=np.float64)
    _, k, n = candidate_alt(counts[:, 0:4], counts[:, 9].astype(int))
    vaf = np.where(n > 0, k / np.maximum(n, 1), 0.0)
    return np.where(k > 0, 1.0 - np.abs(vaf - 0.5) * 2.0, -1.0)


def depth_uncertainty(counts: np.ndarray) -> np.ndarray:
    """Trivial rule C: shallow loci first. A control, not a proposal."""
    return -np.asarray(counts, dtype=np.float64)[:, 6]


# ----------------------------------------------------------------- learned ----

@dataclass
class LogisticRouter:
    """A logistic regression over the counts-only features.

    Small by construction: one weight per feature plus an intercept, so its
    inference is a single matrix-vector product and its cost is dominated by
    the feature build rather than the model.

    Attributes:
        weights: ``[F]`` coefficients.
        intercept: Scalar bias.
        feature_names: Column names the weights correspond to.
        use_binomial: Whether a binomial-LLR column was appended.
        mean: Feature means used for standardization.
        scale: Feature standard deviations used for standardization.
    """

    weights: np.ndarray
    intercept: float
    feature_names: tuple[str, ...]
    use_binomial: bool
    mean: np.ndarray
    scale: np.ndarray

    @property
    def parameter_count(self) -> int:
        """Number of learned parameters, standardization constants included."""
        return int(self.weights.size + 1 + self.mean.size + self.scale.size)

    def design_matrix(self, counts: np.ndarray, binomial_threshold: float | None = None,
                      llr: np.ndarray | None = None) -> np.ndarray:
        """Assemble the standardized design matrix for a count matrix.

        Args:
            counts: ``[N, 10]`` count matrix.
            binomial_threshold: Required when the router uses the binomial LLR.
            llr: Precomputed binomial LLR, optional.

        Returns:
            ``[N, F]`` standardized features.
        """
        matrix = cheap_features(counts).astype(np.float64)
        if self.use_binomial:
            values = binomial_llr(counts) if llr is None else np.asarray(llr)
            margin = -np.abs(values - float(binomial_threshold))
            matrix = np.column_stack([matrix, np.clip(values, -50.0, 50.0) / 10.0,
                                      np.clip(margin, -50.0, 0.0) / 10.0])
        return (matrix - self.mean) / self.scale

    def score(self, counts: np.ndarray, binomial_threshold: float | None = None,
              llr: np.ndarray | None = None) -> np.ndarray:
        """Router uncertainty: the model's logit, higher meaning "route me".

        Args:
            counts: ``[N, 10]`` count matrix.
            binomial_threshold: Required when the router uses the binomial LLR.
            llr: Precomputed binomial LLR, optional.

        Returns:
            ``[N]`` uncertainty scores.
        """
        return self.design_matrix(counts, binomial_threshold, llr) @ self.weights \
            + self.intercept


def oracle_target(pb_llr: np.ndarray, pb_threshold: float,
                  fraction: float = ORACLE_TARGET_FRACTION) -> np.ndarray:
    """The distillation target: would the PB-margin oracle route this locus?

    Args:
        pb_llr: ``[N]`` Poisson-binomial LLRs for the **training** split.
        pb_threshold: PB's frozen threshold.
        fraction: Share of loci the oracle is taken to route.

    Returns:
        ``[N]`` boolean target.
    """
    margin = np.abs(np.asarray(pb_llr, dtype=np.float64) - float(pb_threshold))
    return margin <= np.quantile(margin, fraction)


def fit_logistic_router(counts: np.ndarray, target: np.ndarray, use_binomial: bool,
                        binomial_threshold: float | None = None,
                        seed: int = 0) -> tuple[LogisticRouter, dict]:
    """Fit the router on a training split.

    Args:
        counts: ``[N, 10]`` training count matrix.
        target: ``[N]`` boolean distillation target.
        use_binomial: Append the cheap binomial LLR and its margin.
        binomial_threshold: Required when ``use_binomial``.
        seed: Passed to the solver, so per-seed replication is meaningful.

    Returns:
        ``(router, fit_diagnostics)``.
    """
    from sklearn.linear_model import LogisticRegression

    matrix = cheap_features(counts).astype(np.float64)
    names = list(CHEAP_FEATURE_NAMES)
    if use_binomial:
        values = binomial_llr(counts)
        margin = -np.abs(values - float(binomial_threshold))
        matrix = np.column_stack([matrix, np.clip(values, -50.0, 50.0) / 10.0,
                                  np.clip(margin, -50.0, 0.0) / 10.0])
        names += ["binomial_llr", "binomial_margin"]

    mean = matrix.mean(axis=0)
    scale = np.where(matrix.std(axis=0) > 0, matrix.std(axis=0), 1.0)
    standardized = (matrix - mean) / scale

    start = time.perf_counter()
    model = LogisticRegression(max_iter=1000, C=1.0, random_state=seed)
    model.fit(standardized, np.asarray(target, dtype=bool))
    seconds = time.perf_counter() - start

    router = LogisticRouter(
        weights=model.coef_.reshape(-1).copy(), intercept=float(model.intercept_[0]),
        feature_names=tuple(names), use_binomial=use_binomial, mean=mean, scale=scale)
    diagnostics = {
        "fit_seconds": seconds, "seed": seed, "loci": int(counts.shape[0]),
        "positives": int(np.sum(target)), "parameter_count": router.parameter_count,
        "coefficients": {name: float(weight)
                         for name, weight in zip(names, router.weights)},
        "intercept": router.intercept,
    }
    return router, diagnostics


def measure_throughput(function, counts: np.ndarray, repeats: int = 3) -> dict:
    """Wall-clock cost of a routing stage.

    Args:
        function: Callable taking the count matrix.
        counts: ``[N, 10]`` count matrix.
        repeats: Timed repetitions; the median is reported.

    Returns:
        Seconds and loci/second.
    """
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        function(counts)
        samples.append(time.perf_counter() - start)
    seconds = float(np.median(samples))
    return {"seconds": seconds, "loci": int(counts.shape[0]),
            "loci_per_second": float(counts.shape[0] / seconds),
            "samples": [float(s) for s in samples]}


__all__ = [
    "CHEAP_FEATURE_NAMES", "ORACLE_TARGET_FRACTION", "cheap_features", "binomial_llr",
    "binomial_margin_uncertainty", "alt_count_uncertainty", "depth_uncertainty",
    "LogisticRouter", "oracle_target", "fit_logistic_router", "measure_throughput",
]
