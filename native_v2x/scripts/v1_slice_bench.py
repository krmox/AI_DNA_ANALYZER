"""Run the REAL V1 pipeline (frozen modules, native counts/reads backends, NumPy PB) on the benchmark slice
and time each stage, under the same BAM / reference / windows / CPU as the V2 slice benchmark.

V1 code is imported from the untouched main checkout; nothing there is modified.  Output VCFs go to the HDD and
are compared byte-for-byte with V2's.  Stages:  counts | read-tensor | PB LLR (chunked, as the V1 300x script did)
| binomial+router+calls | genotype+VCF.
"""
import json, os, resource, sys, time
import numpy as np

V1 = os.environ.get("AIDNA_V1_ROOT", os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, V1)
# V1 modules FIRST, from the untouched main checkout (it holds the compiled native backends)
from pileup_counts import load_counts, NATIVE_AVAILABLE  # noqa: E402
from read_level_pileup import load_reads, NATIVE_READS_AVAILABLE  # noqa: E402
from quality_error_model import candidate_alt, extract_quality_evidence, poisson_binomial_llr  # noqa: E402
assert NATIVE_AVAILABLE and NATIVE_READS_AVAILABLE, "V1 native backends not active: this would not be the V1 baseline"
sys.path.insert(1, "/mnt/archive/AI_DNA_ANALYZER_v2/repo/native_v2/oracle")
import v1_oracle as O  # noqa: E402  (Method C / VCF text copied verbatim from V1; imports already cached from V1 tree)
assert __import__("pileup_counts").__file__.startswith(V1)

A = "/mnt/archive/AI_DNA_ANALYZER_benchmark"
BAM = f"{A}/HG002_GRCh38_chr20/HG002.GRCh38.300x_chr20.bam"
REF = f"{A}/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
TRUTH = f"{A}/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
BED = f"{V1}/data/giab_hg002_chr20/chr20_highconf.bed"
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
REGION = (31818624, 32398272)            # windows 400000..408191 of the chr20 window list
t = {}; T0 = time.time()

def stage(name, f):
    t0 = time.time(); r = f(); t[name] = time.time() - t0
    print(f"{name}: {t[name]:.2f}s  (peak RSS so far {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.0f} MB)", flush=True)
    return r

counts, labels, positions = stage("counts_extraction", lambda: load_counts(REF, BAM, TRUTH, BED, REGION, contig="chr20", return_positions=True))
reads, _ = stage("read_tensor_extraction", lambda: load_reads(REF, BAM, TRUTH, BED, REGION, contig="chr20", max_reads=48))
n = positions.size
assert reads.shape[0] == n, (reads.shape, n)

def pb_all():
    out = np.zeros(n)
    ref_idx = counts[:, 9].astype(int)
    alt, k, _ = candidate_alt(counts[:, 0:4], ref_idx)
    CH = 65536
    for s in range(0, n, CH):
        e = min(s + CH, n)
        ev = extract_quality_evidence(reads[s:e], ref_idx[s:e])
        out[s:e] = poisson_binomial_llr(ev, k[s:e], alt_index=alt[s:e])["llr"]
    return out
pb = stage("pb_llr", pb_all)

def route():
    b = O.v1_binomial_llr(counts)
    frame = (labels == 0) | (labels == 1)
    pf, cf, bf, pbf = positions[frame], counts[frame], b[frame], pb[frame]
    B = bf >= 7.0; P = pbf >= 10.5
    routed = np.abs(bf - 7.0) <= O.FROZEN_ROUTER_CUTOFF
    return pf, cf, B, P, routed, np.where(routed, P, B)
pf, cf, B, P, routed, casc = stage("binomial_router", route)

def geno():
    for mask, label, name in ((P, "pb_only", "ai_pb_only_chr20_300x.vcf"), (casc, "cascade", "ai_cascade_chr20_300x.vcf")):
        recs = O.extract_records(mask, cf, pf)
        gt, gq, _ = O.optimized_method_c(recs)
        O.write_vcf_c(recs, gt, gq, f"{OUT}/{name}", label)
stage("genotype_vcf", geno)
total = time.time() - T0
summ = {"region": REGION, "n_loci": int(n), "n_routed": int(routed.sum()), "stage_seconds": t, "total_seconds": total,
        "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
        "reads_tensor_bytes": int(reads.nbytes), "counts_bytes": int(counts.nbytes)}
json.dump(summ, open(f"{OUT}/v1_slice_summary.json", "w"), indent=1); print(json.dumps(summ, indent=1))
