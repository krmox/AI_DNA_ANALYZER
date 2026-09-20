"""I/O-bound hypothesis: evict the BAM from the page cache (posix_fadvise DONTNEED, no root) and time native and python load_reads cold vs warm."""
import os, sys, time, json
from common import *
import read_level_pileup as R
d = HG005; region = (1_000_000, 1_100_000)
def evict():
    fd = os.open(d["bam"], os.O_RDONLY); os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED); os.close(fd)
def t(native):
    R.NATIVE_READS_AVAILABLE = native; s = time.perf_counter(); R.load_reads(d["fasta"], d["bam"], d["vcf"], d["bed"], region, contig=d["contig"]); return round(time.perf_counter() - s, 3)
out = {}
for name, nat in (("native", True), ("python", False)):
    evict(); cold = t(nat); warm = t(nat); evict(); cold2 = t(nat)
    out[name] = dict(cold_s=cold, warm_s=warm, cold_again_s=cold2)
print(json.dumps(out))
