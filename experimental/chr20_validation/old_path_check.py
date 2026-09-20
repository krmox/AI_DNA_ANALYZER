"""Gate: recompute chr20 chunks with the ORIGINAL path (pure-Python counts + pure-Python reads + original O(n) BED scan,
native disabled) and compare every array with the native-path production cache, byte for byte. usage: CHUNK_START"""
import os, sys, json, argparse
os.environ["AI_DNA_ANALYZER_DISABLE_NATIVE_PILEUP"] = "1"
import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[2])  # project root (path-independent)
sys.path.insert(0, _ROOT + "/experimental/load_reads_opt")
from common import *
use_original_window_scan()
import numpy as np, extract_bench_v12 as X, pileup_counts, read_level_pileup
assert not pileup_counts.NATIVE_AVAILABLE and not read_level_pileup.NATIVE_READS_AVAILABLE
c0 = int(sys.argv[1]); c1 = c0 + 500_000; d = CHR20
args = argparse.Namespace(fasta=d["fasta"], bam=d["bam"], vcf=d["vcf"], bed=d["bed"], contig="chr20", features="none", legacy_pb=False)
r = X.process_chunk(args, c0, c1)["blocks"]
g = np.load("cache/chr20_validation/chr20_15x.npz"); gp = g["positions"]
lo, hi = np.searchsorted(gp, r["positions"][0]), np.searchsorted(gp, r["positions"][-1], side="right")
res = {"chunk": [c0, c1], "loci": int(r["positions"].size)}
for k, gk in [("counts", "counts"), ("labels", "labels"), ("positions", "positions"), ("pb_llr", "pb_llr"), ("binomial_llr", "binomial_llr"), ("depth_", "depth"), ("k_full", "k_full"), ("k_retained", "k_retained"), ("n_counted", "n_counted")]:
    a, b = r[k], g[gk][lo:hi]
    res[k] = bool(a.shape == b.shape and a.dtype == b.dtype and np.array_equal(a, b))
res["all_bytewise_equal"] = all(v for k, v in res.items() if k not in ("chunk", "loci"))
print(json.dumps(res))
