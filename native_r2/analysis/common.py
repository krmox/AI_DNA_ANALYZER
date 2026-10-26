"""Shared helpers for the V2.x analysis (Python = analysis only; the production path is C++)."""
import numpy as np, pysam, gzip
from scipy.stats import binom

A = "/mnt/archive/AI_DNA_ANALYZER_benchmark"
TRUTH = f"{A}/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
BED = str(__import__("pathlib").Path(__file__).resolve().parents[2]) + "/data/giab_hg002_chr20/chr20_highconf.bed"
BLOCK = 2_000_000
BASES = "ACGT"

def split_of(pos0):
    return np.where((np.asarray(pos0) // BLOCK) % 2 == 0, "dev", "holdout")

def load_bed():
    iv = []
    for l in open(BED):
        f = l.split()
        if f[0] == "chr20": iv.append((int(f[1]), int(f[2])))
    return np.array(sorted(iv))

def in_bed(pos0, iv=None):
    iv = load_bed() if iv is None else iv
    pos0 = np.asarray(pos0)
    i = np.searchsorted(iv[:, 0], pos0, side="right") - 1
    ok = i >= 0
    out = np.zeros(len(pos0), bool)
    out[ok] = pos0[ok] < iv[i[ok], 1]
    return out

def load_truth_snps():
    """dict pos0 -> (ref, tuple(alts), gtclass)  for SNP records in chr20 (any genotype with an alt)."""
    vf = pysam.VariantFile(TRUTH)
    d = {}
    for r in vf.fetch("chr20"):
        alts = r.alts or ()
        if len(r.ref) == 1 and all(len(a) == 1 for a in alts) and alts:
            gt = r.samples[0]["GT"]
            gtc = "1/1" if (gt[0] == gt[1] and gt[0] != 0) else "0/1"
            d[r.pos - 1] = (r.ref, tuple(alts), gtc, gt)
    return d

EPS = 0.01
_P = np.clip(np.array([EPS, 0.5 * (1 - EPS) + 0.5 * EPS, 1 - EPS]), 1e-12, 1 - 1e-12)
def methodc_gt(k, n):
    """frozen Method C class: '1/1' if hom posterior is argmax else '0/1' (0/0 forced to 0/1)."""
    k = np.asarray(k, float); n = np.asarray(n, float)
    lp = np.stack([binom.logpmf(k, n, p) for p in _P], 1)
    return np.where(lp.argmax(1) == 2, "1/1", "0/1")

def evaluate(calls, truth, in_region):
    """calls: DataFrame-like arrays (pos0, alt(int), gt(str)); truth: dict.  Returns counts inside `in_region(pos0)`.
    TP: truth SNP at pos, call alt in truth alts, genotype class equal.  Else FP.  FN: truth SNPs not TP."""
    tp = set(); fp = 0; fpgt = 0; fpal = 0
    for pos, alt, gt in zip(calls["pos0"], calls["alt"], calls["gt"]):
        if not in_region(pos): continue
        t = truth.get(int(pos))
        if t is not None and BASES[int(alt)] in t[1]:
            if t[2] == gt: tp.add(int(pos)); continue
            fpgt += 1
        elif t is not None: fpal += 1
        fp += 1
    ntruth = sum(1 for p in truth if in_region(p))
    return dict(TP=len(tp), FP=fp, FN=ntruth - len(tp), FPgt=fpgt, FPal=fpal, tp_pos=tp)
