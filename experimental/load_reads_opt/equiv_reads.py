"""Exact equivalence: reference (pure-Python window_reads, original BED scan) vs native load_reads. Byte-level."""
import sys, time, hashlib, json
from common import *
import numpy as np
import read_level_pileup as R
def run(d, region, native):
    R.NATIVE_READS_AVAILABLE = native
    t = time.perf_counter()
    r, l = R.load_reads(d["fasta"], d["bam"], d["vcf"], d["bed"], region, contig=d["contig"])
    return r, l, time.perf_counter() - t
if __name__ == "__main__":
    name, s, e = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); d = DATASETS[name]
    use_original_window_scan()
    ref, rl, tr = run(d, (s, e), False)
    new, nl, tn = run(d, (s, e), True)
    ok = ref.shape == new.shape and ref.dtype == new.dtype
    print(json.dumps(dict(dataset=name, region=[s, e], loci=int(rl.size), valid_read_rows=int((ref[..., 6] & 1).sum()),
        max_valid_rows_per_locus=int((ref[..., 6] & 1).sum(axis=1).max()), shape=list(ref.shape), shape_dtype_equal=ok,
        tensor_cell_mismatches=int((ref != new).sum()) if ok else -1, labels_equal=bool(np.array_equal(rl, nl)),
        sha_ref=hashlib.sha256(ref.tobytes()).hexdigest()[:16], sha_new=hashlib.sha256(new.tobytes()).hexdigest()[:16],
        t_ref=round(tr, 2), t_native=round(tn, 2), speedup=round(tr / tn, 1))))
