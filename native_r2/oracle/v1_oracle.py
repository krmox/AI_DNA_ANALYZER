"""V1 reference oracle: the frozen V1 semantics, exposed as plain functions.

Nothing here is an improved or re-derived version.  Binomial and PB call the frozen V1 modules
(cheap_router, quality_error_model) unmodified; Method C and the VCF writer are copied verbatim
from the V1 chr20 300x adapter (run_chr20_300x_pipeline.py, whose worktree no longer exists;
its text is preserved in this repository's history of the V2 session and its output is checked
byte-for-byte against the archived V1 VCFs by tests/test_reconstructed_v1_vcf.py).
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import gammaln

REPO = str(Path(__file__).resolve().parents[2])
if REPO not in sys.path:
    sys.path.insert(0, REPO)
from cheap_router import binomial_llr as v1_binomial_llr  # noqa: E402
from quality_error_model import (candidate_alt, extract_quality_evidence,  # noqa: E402
                                 poisson_binomial_llr)

FROZEN_BINOMIAL_THRESHOLD = 7.0
FROZEN_PB_THRESHOLD = 10.5
FROZEN_ROUTER_CUTOFF = 5.411872376933351
EPS = 0.01
GQ_CAP = 99.0
REF_TOKENS = ("N", "A", "C", "G", "T")
CONTIG, CONTIG_LEN, SAMPLE, TAG = "chr20", "64444167", "HG002", "chr20_300x"


def binomial(counts):
    """counts [N,10] float64 -> (llr, alt_index, k, n) exactly as V1."""
    llr = v1_binomial_llr(counts)
    alt, k, n = candidate_alt(counts[:, 0:4], counts[:, 9].astype(int))
    return llr, alt, k, n


def pb(reads, counts, alt=None):
    """V1 PB with alt_index (the frozen contract)."""
    ref_idx = counts[:, 9].astype(int)
    a, k, _ = candidate_alt(counts[:, 0:4], ref_idx)
    alt = a if alt is None else alt
    ev = extract_quality_evidence(reads, ref_idx)
    return poisson_binomial_llr(ev, k, alt_index=alt)["llr"]


# --- verbatim from run_chr20_300x_pipeline.py -------------------------------------------------
def alt_allele(count_row, ref_base):
    base_counts = {"A": count_row[0], "C": count_row[1], "G": count_row[2], "T": count_row[3]}
    base_counts.pop(ref_base, None)
    if not base_counts:
        return None
    best_base, best_count = max(base_counts.items(), key=lambda kv: kv[1])
    if best_count <= 0:
        return None
    return best_base


def extract_records(call_mask, counts_f, positions_f):
    idxs = np.nonzero(call_mask)[0]
    records = []
    for i in idxs:
        pos0 = int(positions_f[i])
        ref_idx = int(counts_f[i][9])
        ref_base = REF_TOKENS[ref_idx] if 0 <= ref_idx < 5 else "N"
        if ref_base == "N":
            continue
        alt = alt_allele(counts_f[i], ref_base)
        if alt is None:
            continue
        base_counts = {"A": counts_f[i][0], "C": counts_f[i][1], "G": counts_f[i][2], "T": counts_f[i][3]}
        k = float(base_counts[alt])
        ref_n = float(base_counts[ref_base])
        n = k + ref_n
        d = int(counts_f[i][6])
        records.append(dict(idx=i, pos0=pos0, pos1=pos0 + 1, ref=ref_base, alt=alt, depth=d, k=k, n=n))
    return records


def optimized_method_c(records):
    k = np.array([r["k"] for r in records], dtype=np.float64)
    n = np.array([r["n"] for r in records], dtype=np.float64)
    p = np.array([
        np.clip(0.0 * (1 - EPS) + 1.0 * EPS, 1e-12, 1 - 1e-12),
        np.clip(0.5 * (1 - EPS) + 0.5 * EPS, 1e-12, 1 - 1e-12),
        np.clip(1.0 * (1 - EPS) + 0.0 * EPS, 1e-12, 1 - 1e-12),
    ])
    log_binom_coeff = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    logp = (log_binom_coeff[:, None] + k[:, None] * np.log(p)[None, :]
            + (n - k)[:, None] * np.log(1 - p)[None, :])
    zero_n = (n == 0)
    logp[zero_n, :] = 0.0
    order_ = np.argsort(-logp, axis=1)
    best_idx = order_[:, 0]
    second_idx = order_[:, 1]
    best_logp = logp[np.arange(len(records)), best_idx]
    second_logp = logp[np.arange(len(records)), second_idx]
    gq_raw = np.minimum(10.0 * (best_logp - second_logp) / np.log(10.0), GQ_CAP)
    gq_raw = np.maximum(gq_raw, 0.0)
    classes = np.array(["0/0", "0/1", "1/1"])
    best_class = classes[best_idx]
    is_theta0 = best_class == "0/0"
    gt_out = np.where(is_theta0, "0/1", best_class)
    gq_capped = np.where(is_theta0, np.minimum(gq_raw, 5.0), gq_raw)
    gq_out = np.round(gq_capped, 1)
    return gt_out, gq_out, is_theta0


HEADER_TMPL = """##fileformat=VCFv4.2
##source=ai_dna_analyzer_{tag}_{label}_C
##contig=<ID={contig},length={contig_len}>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth (passing-filter reads)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (frozen Method C: binomial genotype posterior argmax over {{0/0,0/1,1/1}}, eps=0.01 fixed)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality (10*log10(P_best/P_second), capped at 99)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""


def write_vcf_c(records, gt, gq, out_path, label):
    order_ = np.argsort([r["pos1"] for r in records])
    lines = []
    for i in order_:
        r = records[i]
        lines.append(f"{CONTIG}\t{r['pos1']}\t.\t{r['ref']}\t{r['alt']}\t.\tPASS\t"
                     f"DP={r['depth']}\tGT:GQ\t{gt[i]}:{int(round(float(gq[i])))}")
    with open(out_path, "w") as fh:
        fh.write(HEADER_TMPL.format(label=label, tag=TAG, contig=CONTIG, contig_len=CONTIG_LEN, sample=SAMPLE))
        fh.write("\n".join(lines))
        if lines:
            fh.write("\n")
