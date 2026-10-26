"""Freeze the V1 side of the 3,491-locus correctness gate (read-only over V1 caches).

Output: $WORK/gate/gate_positions.txt (0-based) and $WORK/gate/v1_gate.npz with, per locus,
the V1 counts row, label, V1 Binomial (frozen functions), cached PB LLR, the V1 read-tensor rows
(48x8 uint8) and the V1 cascade / PB-only calls.
"""
import sys, zipfile, time
from pathlib import Path
import numpy as np
import pandas as pd
from numpy.lib import format as npf

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "oracle"))
import v1_oracle as O  # noqa: E402

WORK = Path("/mnt/archive/AI_DNA_ANALYZER_v2/work/gate"); WORK.mkdir(parents=True, exist_ok=True)
CACHE = "/mnt/archive/AI_DNA_ANALYZER_benchmark/ai_dna_analyzer_run/cache"
TSV = str(Path(__file__).resolve().parents[2]) + "/research/benchmark_chr20_300x/routing_analysis/loci_of_interest.tsv"
t0 = time.time()
d = pd.read_csv(TSV, sep="\t")
pos0 = np.sort(d.pos1.values.astype(np.int64) - 1)
np.savetxt(WORK / "gate_positions.txt", pos0, fmt="%d")
cz = np.load(f"{CACHE}/chr20_300x_counts.npz")
positions = cz["positions"]
idx = np.searchsorted(positions, pos0)
assert np.array_equal(positions[idx], pos0)
counts = cz["counts"][idx]; labels = cz["labels"][idx]
pb_cached = np.load(f"{CACHE}/chr20_300x_pb_llr.npz")["pb_llr"][idx]
print("counts/labels/pb loaded", time.time() - t0, flush=True)

zf = zipfile.ZipFile(f"{CACHE}/chr20_300x_reads.npz")
fh = zf.open("reads.npy")
ver = npf.read_magic(fh)
shape, fortran, dtype = (npf.read_array_header_1_0 if ver == (1, 0) else npf.read_array_header_2_0)(fh)
assert shape[0] == positions.size and shape[1:] == (48, 8)
reads = np.empty((idx.size, 48, 8), np.uint8)
CH = 1_000_000; j = 0
for s in range(0, shape[0], CH):
    e = min(s + CH, shape[0])
    blk = np.frombuffer(fh.read((e - s) * 384), dtype=np.uint8).reshape(e - s, 48, 8)
    lo, hi = np.searchsorted(idx, [s, e])
    reads[j:j + hi - lo] = blk[idx[lo:hi] - s]; j += hi - lo
assert j == idx.size
print("reads extracted", time.time() - t0, flush=True)

llr, alt, k, n = O.binomial(counts)
frame = (labels == 0) | (labels == 1)
routed = frame & (np.abs(llr - 7.0) <= O.FROZEN_ROUTER_CUTOFF)
bcall = llr >= 7.0; pcall = pb_cached >= 10.5
ccall = np.where(routed, pcall, bcall)
np.savez(WORK / "v1_gate.npz", pos0=pos0, counts=counts, labels=labels, pb_cached=pb_cached, reads=reads,
         binom_llr=llr, alt=alt, k=k, n=n, frame=frame, routed=routed, bcall=bcall, pcall=pcall, ccall=ccall)
print("saved", idx.size, "loci; routed in gate:", int(routed.sum()), time.time() - t0)
