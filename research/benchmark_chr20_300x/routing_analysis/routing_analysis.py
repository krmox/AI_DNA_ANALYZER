"""chr20 300x routing / rescue / loss analysis (V2 programme, items 10-14).

Read-only over the frozen V1 artefacts. No threshold, model or output of V1 is changed.
Category definitions: CATEGORY_PREREGISTRATION.md (written before this script ran).

Stage 1 (this file): per-locus table for every locus of interest, plus population-level
routing-dependence tables over all 56.24M frame loci.
"""
import gzip
import json
import os
import sys
import tempfile
import time
import zipfile

import numpy as np
import pysam
from numpy.lib import format as npf

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
# Originally run from the opt-chr20 worktree of this repository; V1 modules are identical.
REPO = os.environ.get("AIDNA_V1_ROOT", _ROOT)
sys.path.insert(0, REPO)
from cheap_router import binomial_llr as compute_binomial_llr  # noqa: E402
from quality_error_model import (candidate_alt, extract_quality_evidence,  # noqa: E402
                                 poisson_binomial_llr)
from read_level_pileup import FLAG_COUNTED, FLAG_VALID  # noqa: E402

ARCH = "/mnt/archive/AI_DNA_ANALYZER_benchmark"
CACHE = f"{ARCH}/ai_dna_analyzer_run/cache"
OUTDIR = _HERE
HAPPY = os.path.join(_HERE, "..", "happy")
TRUTH_VCF = f"{ARCH}/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
BED = f"{ARCH}/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
STRAT = os.path.join(_ROOT, "data", "strat_v31")
CONTIG, CONTIG_LEN = "chr20", 64444167

T_B, T_PB, CUT = 7.0, 10.5, 5.411872376933351
BASES = "ACGT"
REF_TOKENS = ("N", "A", "C", "G", "T")
N_RANDOM, SEED = 1000, 20260924

t0 = time.time()


def log(msg):
    print(f"[{time.time() - t0:7.1f}s] {msg}", flush=True)


# ---------------------------------------------------------------- load -------
cz = np.load(f"{CACHE}/chr20_300x_counts.npz")
counts, labels, positions = cz["counts"], cz["labels"], cz["positions"]
pz = np.load(f"{CACHE}/chr20_300x_pb_llr.npz")
pb_llr = pz["pb_llr"]
assert np.array_equal(pz["positions"], positions) and np.array_equal(pz["labels"], labels)
log(f"loaded {positions.size} loci")

binom = compute_binomial_llr(counts)
frame = (labels == 0) | (labels == 1)
truth = labels == 1
B = binom >= T_B
P = pb_llr >= T_PB
routed = frame & (np.abs(binom - T_B) <= CUT)
casc = np.where(routed, P, B)

bc = counts[:, 0:4]
ref_idx = counts[:, 9].astype(np.int64)
depth = counts[:, 6]
n_cnt = bc.sum(1)
alt_i, k_full, _ = candidate_alt(bc, ref_idx)
vaf = np.where(n_cnt > 0, k_full / np.maximum(n_cnt, 1), np.nan)
mean_bq = np.where(n_cnt > 0, counts[:, 7] / np.maximum(n_cnt, 1), np.nan)
mean_mq = np.where(depth > 0, counts[:, 8] / np.maximum(depth, 1), np.nan)

# sanity vs frozen meta
meta = json.load(open(f"{ARCH}/ai_dna_analyzer_run/pipeline_meta_chr20_300x.json"))
assert int(frame.sum()) == meta["routing"]["n_scored"]
assert int(routed.sum()) == meta["routing"]["n_routed_pb"]
rescuable = frame & (P == truth) & (B != truth)
rescued = rescuable & routed
assert int(rescuable.sum()) == meta["rescue"]["pb_rescuable_loci"]
assert int(rescued.sum()) == meta["rescue"]["rescued_by_router"]
log("frozen meta reproduced (frame, routed, rescuable, rescued)")

# --------------------------------------------------------- stratification ----
def interval_mask(path, gz):
    m = np.zeros(CONTIG_LEN + 1, dtype=bool)
    op = gzip.open if gz else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c, s, e = line.split("\t")[:3]
            if c == CONTIG:
                m[int(s):int(e)] = True
    return m[positions]


in_bed = interval_mask(BED, False)
strat = {name: interval_mask(f"{STRAT}/{name}.bed.gz", True)
         for name in ("alldifficult", "lowmap_segdup", "tandemrepeats")}
log("stratifications loaded")

# ------------------------------------------------------------ populations ---
rng = np.random.default_rng(SEED)
frame_idx = np.nonzero(frame)[0]
random_idx = np.sort(rng.choice(frame_idx, N_RANDOM, replace=False))
disagree = frame & (P != B)
loss_truth = frame & truth & (~casc | ~P)
sel = np.zeros(positions.size, dtype=bool)
sel[random_idx] = True
sel |= routed | rescuable | loss_truth | (frame & (casc != P))
pop = {
    "random1000": np.isin(np.arange(positions.size), random_idx),
    "routed": routed, "rescuable": rescuable, "rescued": rescued,
    "truth_snp_loss": loss_truth, "pb_binom_disagree": disagree,
    "cascade_pbonly_disagree": frame & (casc != P),
}
sel_idx = np.nonzero(sel)[0]
log(f"selected {sel_idx.size} loci for read-level extraction "
    f"(disagree P!=B in frame: {int(disagree.sum())})")
if disagree.sum() > 200000:
    raise SystemExit("disagreement set unexpectedly large; refusing to continue blindly")
sel |= disagree
sel_idx = np.nonzero(sel)[0]

# ------------------------------------------------- stream the read tensor ----
zf = zipfile.ZipFile(f"{CACHE}/chr20_300x_reads.npz")
fh = zf.open("reads.npy")
ver = npf.read_magic(fh)
shape, fortran, dtype = (npf.read_array_header_1_0 if ver == (1, 0) else npf.read_array_header_2_0)(fh)
assert shape[0] == positions.size and shape[1:] == (48, 8) and not fortran
row_bytes = 48 * 8
reads_sel = np.empty((sel_idx.size, 48, 8), dtype=np.uint8)
CH = 1_000_000
j = 0
for start in range(0, shape[0], CH):
    stop = min(start + CH, shape[0])
    buf = fh.read((stop - start) * row_bytes)
    block = np.frombuffer(buf, dtype=np.uint8).reshape(stop - start, 48, 8)
    lo, hi = np.searchsorted(sel_idx, [start, stop])
    reads_sel[j:j + hi - lo] = block[sel_idx[lo:hi] - start]
    j += hi - lo
    if start // CH % 10 == 0:
        log(f"  reads {stop}/{shape[0]}")
assert j == sel_idx.size
np.save(os.path.join(tempfile.gettempdir(), "reads_sel.npy"), reads_sel)
np.save(os.path.join(tempfile.gettempdir(), "sel_idx.npy"), sel_idx)
log("read tensor rows extracted")

# ------------------------------- recompute PB with the frozen functions ------
ev = extract_quality_evidence(reads_sel, ref_idx[sel_idx])
pb_alt = poisson_binomial_llr(ev, k_full[sel_idx], alt_index=alt_i[sel_idx])
pb_min = poisson_binomial_llr(ev, k_full[sel_idx])
cached = pb_llr[sel_idx]
eq_alt = np.array_equal(pb_alt["llr"], cached)
eq_min = np.array_equal(pb_min["llr"], cached)
pb_recompute = {
    "n": int(sel_idx.size),
    "bit_identical_alt_index_variant": bool(eq_alt),
    "bit_identical_min_k_n_variant": bool(eq_min),
    "max_abs_diff_alt_index_variant": float(np.nanmax(np.abs(pb_alt["llr"] - cached))),
    "max_abs_diff_min_k_n_variant": float(np.nanmax(np.abs(pb_min["llr"] - cached))),
    "n_mismatch_alt_index_variant": int(np.sum(pb_alt["llr"] != cached)),
    "n_mismatch_min_k_n_variant": int(np.sum(pb_min["llr"] != cached)),
}
log(f"PB recompute check: {pb_recompute}")
use = pb_alt if eq_alt or not eq_min else pb_min
k_tensor = use["k_effective"]

# ----------------------------------------------------- read-level stats -----
valid = (reads_sel[..., 6] & FLAG_VALID) > 0
counted = (reads_sel[..., 6] & FLAG_COUNTED) > 0
bq = reads_sel[..., 1].astype(float)
mq = reads_sel[..., 2].astype(float)
is_alt = counted & (reads_sel[..., 0] == alt_i[sel_idx][:, None])
is_ref = counted & (reads_sel[..., 0] == (ref_idx[sel_idx] - 1)[:, None])


def mstat(v, m, fn):
    out = np.full(v.shape[0], np.nan)
    for i in range(v.shape[0]):
        x = v[i][m[i]]
        if x.size:
            out[i] = fn(x)
    return out


# ------------------------------------------------------------ truth alleles ---
truth_alt = {}
with pysam.VariantFile(TRUTH_VCF) as vf:
    for r in vf.fetch(CONTIG):
        if len(r.ref) == 1 and all(len(a) == 1 for a in r.alts):
            gt = r.samples[0]["GT"]
            truth_alt[r.pos - 1] = (",".join(r.alts), "/".join(map(str, gt)))


def happy_status(path):
    st = {}
    with pysam.VariantFile(path) as vf:
        for r in vf.fetch(CONTIG):
            s = r.samples
            t, q = s[0]["BD"], s[1]["BD"]
            st[r.pos - 1] = (t or ".", q or ".")
    return st


hp_c = happy_status(f"{HAPPY}/ai_dna_analyzer/ai_cascade.vcf.gz")
hp_p = happy_status(f"{HAPPY}/ai_pb_only/ai_pb_only.vcf.gz")


def vcf_gt(path):
    g = {}
    with open(path) as fh2:
        for line in fh2:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            g[int(f[1]) - 1] = (f[4], f[9].split(":")[0])
    return g


gt_c = vcf_gt(f"{ARCH}/ai_dna_analyzer_run/ai_cascade_chr20_300x.vcf")
gt_p = vcf_gt(f"{ARCH}/ai_dna_analyzer_run/ai_pb_only_chr20_300x.vcf")
log("truth / hap.py / VCF genotypes loaded")


# --------------------------------------------------------- classification ----
def outcome(i):
    if not routed[i]:
        return "NOT_ROUTED"
    b_ok, p_ok = B[i] == truth[i], P[i] == truth[i]
    return {(False, True): "RESCUE", (True, False): "HARMFUL",
            (True, True): "NEUTRAL_BOTH_CORRECT", (False, False): "NEUTRAL_BOTH_WRONG"}[(b_ok, p_ok)]


def subtype(i):
    if B[i] == P[i]:
        return ""
    if truth[i]:
        return "FN_AVOIDED" if P[i] else "FN_CREATED"
    return "FP_AVOIDED" if not P[i] else "FP_CREATED"


def loss_group(i):
    if not (frame[i] and truth[i]):
        return ""
    if rescued[i]:
        g = "C_RESCUED_BY_PB"
    elif not casc[i] and not P[i]:
        g = "A_LOST_BOTH"
    elif not casc[i] and P[i]:
        g = "B_LOST_CASCADE_ONLY"
    elif casc[i] and not P[i]:
        g = "E_LOST_PB_ONLY"
    else:
        g = ""
    if routed[i] and B[i] and not P[i]:
        g = (g + "+" if g else "") + "D_FALSE_RESCUE_LOSS"
    return g


bq_all_med = mstat(bq, counted, np.median)
bq_alt_mean = mstat(bq, is_alt, np.mean)
bq_ref_mean = mstat(bq, is_ref, np.mean)
mq_valid_mean = mstat(mq, valid, np.mean)
mq_alt_mean = mstat(mq, is_alt & valid, np.mean)
mq_valid_min = mstat(mq, valid, np.min)
n_valid = valid.sum(1)
n_counted_t = counted.sum(1)
n_lowbq = (valid & ~counted & (reads_sel[..., 0] <= 3)).sum(1)

cols = ["chrom", "pos1", "ref", "alt_candidate", "label", "in_confident_bed",
        "truth_alt", "truth_gt", "binomial_llr", "pb_llr", "pb_llr_recomputed",
        "routed", "binomial_call", "pb_call", "cascade_call", "outcome", "subtype",
        "loss_group", "near_pb_boundary", "near_router_edge", "depth_capped",
        "depth", "n_counted_full", "k_full", "vaf_full",
        "reads_seen_tensor", "counted_in_tensor", "k_in_tensor", "vaf_tensor",
        "lowbq_rows_in_tensor", "mean_bq_full", "median_bq_tensor", "mean_bq_alt_tensor",
        "mean_bq_ref_tensor", "mean_mapq_full", "mean_mapq_tensor", "min_mapq_tensor",
        "mean_mapq_alt_tensor", "cascade_gt", "pbonly_gt", "genotype_diff",
        "happy_cascade_truth_bd", "happy_cascade_query_bd", "happy_pbonly_query_bd",
        "alldifficult", "lowmap_segdup", "tandemrepeats"] + [f"in_{p}" for p in pop]

rows = []
for j, i in enumerate(sel_idx):
    p0 = int(positions[i])
    ref_b = REF_TOKENS[ref_idx[i]] if 0 <= ref_idx[i] < 5 else "N"
    ta, tg = truth_alt.get(p0, (".", "."))
    cg = gt_c.get(p0, (".", "."))[1]
    pg = gt_p.get(p0, (".", "."))[1]
    hc, hp = hp_c.get(p0, (".", ".")), hp_p.get(p0, (".", "."))
    kt = int(k_tensor[j])
    rows.append([
        CONTIG, p0 + 1, ref_b, BASES[alt_i[i]], int(labels[i]), int(in_bed[i]), ta, tg,
        f"{binom[i]:.6f}", f"{pb_llr[i]:.6f}", f"{use['llr'][j]:.6f}",
        int(routed[i]), int(B[i]), int(P[i]), int(casc[i]), outcome(i), subtype(i),
        loss_group(i), int(abs(pb_llr[i] - T_PB) < 1.0),
        int(routed[i] and abs(binom[i] - T_B) > 0.9 * CUT), int(depth[i] > 48),
        int(depth[i]), int(n_cnt[i]), int(k_full[i]), f"{vaf[i]:.4f}",
        int(n_valid[j]), int(n_counted_t[j]), kt,
        f"{kt / n_counted_t[j]:.4f}" if n_counted_t[j] else "nan", int(n_lowbq[j]),
        f"{mean_bq[i]:.2f}", f"{bq_all_med[j]:.1f}", f"{bq_alt_mean[j]:.2f}",
        f"{bq_ref_mean[j]:.2f}", f"{mean_mq[i]:.2f}", f"{mq_valid_mean[j]:.2f}",
        f"{mq_valid_min[j]:.0f}", f"{mq_alt_mean[j]:.2f}", cg, pg,
        int(cg != "." and pg != "." and cg != pg), hc[0], hc[1], hp[1],
        int(strat["alldifficult"][i]), int(strat["lowmap_segdup"][i]),
        int(strat["tandemrepeats"][i])] + [int(pop[p][i]) for p in pop])

with open(f"{OUTDIR}/loci_of_interest.tsv", "w") as out:
    out.write("\t".join(cols) + "\n")
    for r in rows:
        out.write("\t".join(map(str, r)) + "\n")
log(f"wrote loci_of_interest.tsv ({len(rows)} rows)")

# ------------------------------------------- population routing dependence --
def bin_table(name, x, edges):
    lines = [f"# routing probability vs {name} (population: frame loci, n={int(frame.sum())})",
             "bin_lo\tbin_hi\tn_loci\tn_routed\tp_routed\tn_rescue\tn_harmful\t"
             "n_neutral_correct\tn_neutral_wrong\tpower"]
    out_cls = {
        "rescue": rescued, "harmful": routed & (B == truth) & (P != truth),
        "nc": routed & (B == truth) & (P == truth), "nw": routed & (B != truth) & (P != truth)}
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = frame & (x >= lo) & (x < hi)
        n, r = int(m.sum()), int((m & routed).sum())
        vals = [int((m & out_cls[c]).sum()) for c in ("rescue", "harmful", "nc", "nw")]
        lines.append(f"{lo}\t{hi}\t{n}\t{r}\t{(r / n if n else float('nan')):.3e}\t"
                     + "\t".join(map(str, vals))
                     + f"\t{'UNDERPOWERED' if r < 30 else 'ok'}")
    return "\n".join(lines) + "\n"


with open(f"{OUTDIR}/routing_dependence.tsv", "w") as out:
    out.write(bin_table("depth", depth, [0, 50, 100, 150, 200, 250, 300, 350, 400, 500, 1e9]))
    out.write(bin_table("VAF (k/n, full pileup)", vaf,
                        [0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.65, 0.9, 1.01]))
    out.write(bin_table("mean BQ (counted bases)", mean_bq, [0, 20, 25, 30, 32, 34, 36, 38, 40, 99]))
    out.write(bin_table("mean MAPQ (depth reads)", mean_mq, [0, 20, 30, 40, 50, 55, 58, 60, 61, 256]))
    for s in strat:
        out.write(bin_table(f"stratum {s} (0=outside,1=inside)", strat[s].astype(float), [0, 0.5, 1.5]))
log("wrote routing_dependence.tsv")

summary = {
    "pb_recompute_check": pb_recompute,
    "populations": {p: int((pop[p]).sum()) for p in pop},
    "frame_in_confident_bed": int((frame & in_bed).sum()),
    "frame_total": int(frame.sum()),
    "rows_written": len(rows),
}
json.dump(summary, open(f"{OUTDIR}/stage1_summary.json", "w"), indent=2)
log(json.dumps(summary, indent=2))
