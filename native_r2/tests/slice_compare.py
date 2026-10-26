"""Bit-level comparison of dumped Binomial/PB LLR arrays (window order) against the V1 caches.

slice_compare.py PREFIX FIRST_WINDOW N_WINDOWS [label]
PREFIX.binom.f64 / PREFIX.pb.f64 are what `dnav2 --dump-llr-prefix` writes; the V1 caches are sliced by
window index (cache rows are in window order, 64 loci per window).
"""
import os, sys, time
import numpy as np
prefix, w0, nw = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
label = sys.argv[4] if len(sys.argv) > 4 else "slice"
CACHE = "/mnt/archive/AI_DNA_ANALYZER_benchmark/ai_dna_analyzer_run/cache"
lo, hi = w0 * 64, (w0 + nw) * 64
def sl(f, key):
    with np.load(f"{CACHE}/{f}") as z:
        return z[key][lo:hi] if key != "counts" else z[key][lo:hi]
pb1 = sl("chr20_300x_pb_llr.npz", "pb_llr")
counts = sl("chr20_300x_counts.npz", "counts")
labels = sl("chr20_300x_counts.npz", "labels")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "oracle"))
import v1_oracle as O
b1 = O.v1_binomial_llr(counts)
b2 = np.fromfile(prefix + ".binom.f64", dtype=np.float64)
p2 = np.fromfile(prefix + ".pb.f64", dtype=np.float64)
u = lambda a: a.view(np.uint64)  # noqa: E731
res = [("binomial_llr_bit_exact_all_slice_loci", b1.size, int((u(b1) != u(b2)).sum())),
       ("pb_llr_bit_exact_all_slice_loci", pb1.size, int((u(pb1) != u(p2)).sum()))]
fr = (labels == 0) | (labels == 1)
res.append(("pb_call_exact_frame", int(fr.sum()), int(((pb1 >= 10.5) != (p2 >= 10.5))[fr].sum())))
res.append(("router_mask_exact_frame", int(fr.sum()), int((((np.abs(b1 - 7) <= O.FROZEN_ROUTER_CUTOFF) != (np.abs(b2 - 7) <= O.FROZEN_ROUTER_CUTOFF)) & fr).sum())))
bad = np.nonzero(u(pb1) != u(p2))[0]
for name, n, b in res: print(f"[{'PASS' if b == 0 else 'FAIL'}] {name} ({label}): n={n} mismatches={b}")
if bad.size: print("first mismatching loci (index in slice, V1, V2):", [(int(i), float(pb1[i]), float(p2[i])) for i in bad[:8]])
if os.environ.get("V2_CORRECTNESS_CSV"):
    with open(os.environ["V2_CORRECTNESS_CSV"], "a") as f:
        for name, n, b in res: f.write(f"{label},{time.strftime('%F')},{name},{n},{b},{'PASS' if b==0 else 'FAIL'},\"windows {w0}..{w0+nw}\"\n")
sys.exit(1 if any(b for *_, b in res) else 0)
