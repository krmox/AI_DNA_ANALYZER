"""Bit-equality tests of the V2 scalar kernels against the frozen V1 (NumPy/SciPy) oracle.

Run:  python native_v2/tests/test_kernels.py     (needs `make kernels`)
Writes one row per check to $V2_CORRECTNESS_CSV (default: stdout only).
"""
import ctypes
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.special import gammaln

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "oracle"))
import v1_oracle as O  # noqa: E402

LIB = os.environ.get("V2_KERNEL_LIB", "/mnt/archive/AI_DNA_ANALYZER_v2/build/libdnav2_kernels.so")
lib = ctypes.CDLL(LIB)
P = lambda a: a.ctypes.data_as(ctypes.c_void_p)  # noqa: E731
rows = []


def record(name, n, n_bad, detail=""):
    status = "PASS" if n_bad == 0 else "FAIL"
    rows.append((name, n, n_bad, status, detail))
    print(f"[{status}] {name}: n={n} mismatches={n_bad} {detail}", flush=True)


def bits_equal(a, b):
    return np.asarray(a, dtype=np.float64).view(np.uint64) == np.asarray(b, dtype=np.float64).view(np.uint64)


rng = np.random.default_rng(20260924)

# 1. gammaln: every integer argument that can occur (count+1 up to 10001) + random reals
x = np.concatenate([np.arange(1, 10002, dtype=np.float64), rng.uniform(1, 3000, 200000)])
out = np.empty_like(x)
lib.v2_gammaln(P(x), P(out), ctypes.c_int64(x.size))
ok = bits_equal(out, gammaln(x))
record("gammaln_bit_equal", x.size, int((~ok).sum()), f"max_ulp_bad_examples={x[~ok][:3].tolist()}")


# 1b. log1p_cephes == scipy.special.log1p on the arguments that occur (-p) and a wide sweep
from scipy.special import log1p as sp_log1p  # noqa: E402
y = np.concatenate([-np.array([0.01 / 3, 0.5 * 0.99 + 0.5 * (0.01 / 3), 1 - 0.01 + 0.01 / 3]),
                    rng.uniform(-0.999999, 5, 300000), rng.uniform(-1e-3, 1e-3, 100000)])
out = np.empty_like(y)
lib.v2_log1p_cephes(P(y), P(out), ctypes.c_int64(y.size))
record("log1p_cephes_bit_equal_scipy", y.size, int((~bits_equal(out, sp_log1p(y))).sum()))

# 2. Binomial: synthetic count matrices, incl. ref=N, ties, zero depth, extreme depth
def rand_counts(n, max_depth):
    depth = rng.integers(0, max_depth + 1, n)
    c = np.zeros((n, 10))
    for i in range(n):
        d = depth[i]
        if d:
            mode = rng.integers(0, 4)
            vaf = [0.0, rng.uniform(0, 0.05), rng.uniform(0.05, 0.6), rng.uniform(0.6, 1.0)][mode]
            alt = rng.integers(0, 4)
            k = rng.binomial(d, vaf)
            probs = np.full(4, 0.0)
            ref = rng.integers(0, 4)
            probs[ref] = 1.0
            rest = rng.multinomial(d - k, probs)
            c[i, :4] = rest
            c[i, alt] += k
    c[:, 4] = rng.integers(0, 5, n)
    c[:, 5] = rng.integers(0, 5, n)
    c[:, 6] = c[:, :4].sum(1) + rng.integers(0, 10, n)
    c[:, 7] = c[:, :4].sum(1) * 33
    c[:, 8] = c[:, 6] * 60
    c[:, 9] = rng.integers(0, 5, n)
    return c


for max_depth in (5, 60, 400, 1500, 9000):
    c = rand_counts(60000 if max_depth < 2000 else 8000, max_depth)
    # force ties
    c[:200, :4] = np.array([3, 3, 3, 3])
    llr = np.empty(len(c)); alt = np.empty(len(c), np.int32); k = np.empty(len(c)); n = np.empty(len(c))
    lib.v2_binomial(P(c), ctypes.c_int64(len(c)), P(llr), P(alt), P(k), P(n))
    rl, ra, rk, rn = O.binomial(c)
    bad = (~bits_equal(llr, rl)) | (alt != ra) | (~bits_equal(k, rk)) | (~bits_equal(n, rn))
    record(f"binomial_bit_equal_maxdepth{max_depth}", len(c), int(bad.sum()))
    routed_c = np.abs(llr - 7.0) <= 5.411872376933351
    routed_r = np.abs(rl - 7.0) <= 5.411872376933351
    record(f"router_mask_equal_maxdepth{max_depth}", len(c), int((routed_c != routed_r).sum()))


# 3. PB on synthetic V1-tensor rows
def synth_reads(n_loci, n_reads, vaf_kind, bq_kind):
    reads = np.zeros((n_loci, 48, 8), np.uint8)
    reads[:, :, 0] = 255
    counts = np.zeros((n_loci, 10))
    for i in range(n_loci):
        m = n_reads if n_reads else rng.integers(0, 49)
        ref = rng.integers(0, 4)
        alt = (ref + rng.integers(1, 4)) % 4
        vaf = {"zero": 0.0, "low": rng.uniform(0, 0.15), "het": rng.uniform(0.3, 0.7),
               "hom": rng.uniform(0.9, 1.0), "any": rng.uniform(0, 1)}[vaf_kind]
        bqs = {"q30": rng.integers(25, 40, m), "wide": rng.integers(0, 94, m),
               "extreme": rng.choice([0, 1, 2, 12, 13, 60, 93, 255], m),
               "low": rng.integers(13, 22, m)}[bq_kind]
        for j in range(m):
            isalt = rng.random() < vaf
            code = alt if isalt else ref
            if rng.random() < 0.05: code = 5           # non-ACGT
            q = int(bqs[j])
            counted = code <= 3 and q >= 13
            reads[i, j] = [code, q, 60, 0, 10, 10, 1 | (2 if counted else 0), 0]
            if counted: counts[i, code] += 1
        counts[i, 6] = m
        counts[i, 9] = ref + 1
        # the depth>48 subsample effect: full counts can exceed the tensor
        extra = rng.integers(0, 3) * 100 * (n_reads == 48)
        counts[i, alt] += extra
    return reads, counts


for n_reads in (8, 16, 24, 32, 48, 0):
    for vaf_kind in ("zero", "low", "het", "hom", "any"):
        for bq_kind in ("q30", "wide", "extreme", "low"):
            reads, counts = synth_reads(150, n_reads, vaf_kind, bq_kind)
            _, alt, k, n = O.binomial(counts)
            ref = O.pb(reads, counts)
            res = {}
            for mode in (0, 1):
                out = np.empty(len(reads))
                lib.v2_pb(P(reads), P(alt.astype(np.int32)), ctypes.c_int64(len(reads)), mode, P(out))
                res[mode] = out
            bad0 = int((~bits_equal(res[0], ref)).sum())
            bad1 = int((~bits_equal(res[1], ref)).sum())
            dec = int(((res[1] >= 10.5) != (ref >= 10.5)).sum())
            if bad0 or bad1 or dec:
                record(f"pb_synth_R{n_reads}_{vaf_kind}_{bq_kind}", len(reads), bad0 + bad1 + dec,
                       f"ref_bad={bad0} native_bad={bad1} decision_bad={dec}")
tot = 0
record("pb_synthetic_grid_all_cells", 6 * 5 * 4 * 150, 0 if all(r[3] == "PASS" for r in rows if r[0].startswith("pb_synth")) else 1,
       "reference and native vs V1 poisson_binomial_llr(alt_index), bit-level + decision-level "
       "(a failing cell would have printed its own row above)")

# 4. Method C vs verbatim V1
kk = np.concatenate([rng.integers(0, 40, 50000), rng.integers(0, 400, 50000)]).astype(float)
nn = kk + np.concatenate([rng.integers(0, 60, 50000), rng.integers(0, 600, 50000)]).astype(float)
nn[:1000] = kk[:1000]  # hom-like
recs = [dict(k=a, n=b) for a, b in zip(kk, nn)]
gt_o, gq_o, forced_o = O.optimized_method_c(recs)
gt = np.empty(len(kk), np.int32); gq = np.empty(len(kk), np.int32); fo = np.empty(len(kk), np.int32)
lib.v2_methodc(P(kk), P(nn), ctypes.c_int64(len(kk)), P(gt), P(gq), P(fo))
gts = np.where(gt == 2, "1/1", "0/1")
gq_ref = np.array([int(round(float(g))) for g in gq_o])
record("methodc_gt_equal", len(kk), int((gts != gt_o).sum()))
record("methodc_gq_equal", len(kk), int((gq != gq_ref).sum()))
record("methodc_forced_equal", len(kk), int((fo.astype(bool) != forced_o).sum()))

csv = os.environ.get("V2_CORRECTNESS_CSV")
if csv:
    with open(csv, "a") as fh:
        for r in rows:
            fh.write(f"kernel_tests,{time.strftime('%F')},{r[0]},{r[1]},{r[2]},{r[3]},\"{r[4]}\"\n")
sys.exit(1 if any(r[3] == "FAIL" for r in rows) else 0)
