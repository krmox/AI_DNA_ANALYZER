"""chr20 (HG002, 15x) holdout validation: frozen cascade + frozen Method C, applied UNMODIFIED.
Derived from ``experimental/stress_test/run_hg005_pipeline.py`` (the same script the HG005 post-fix
run used, which differed from it only in paths): only the paths, contig name, sample name and VCF
header contig length are parameterised (env PIPE_*); the frame filter, the three frozen constants,
route/cascade logic, Method C (``optimized_method_c``), alt-allele selection and VCF writer are the
same statements. Adds stage timers and writes the per-locus call arrays used by the analysis script.
"""
import glob
import json
import os
import time

import numpy as np
from scipy.special import gammaln

import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[2])  # project root (path-independent)
REPO = _ROOT
OUT = os.environ.get("PIPE_OUT", f"{REPO}/experimental/chr20_validation")
CACHE_GLOB = os.environ.get("PIPE_CACHE", f"{REPO}/cache/chr20_validation/chr20_15x.npz")
CONTIG = os.environ.get("PIPE_CONTIG", "chr20")
CONTIG_LEN = os.environ.get("PIPE_CONTIG_LEN", "64444167")
SAMPLE = os.environ.get("PIPE_SAMPLE", "HG002")
TAG = os.environ.get("PIPE_TAG", "chr20")

FROZEN_BINOMIAL_THRESHOLD = 7.0
FROZEN_PB_THRESHOLD = 10.5
FROZEN_ROUTER_CUTOFF = 5.411872376933351
LABEL_SNP = 1

EPS = 0.01
GQ_CAP = 99.0
REF_TOKENS = ("N", "A", "C", "G", "T")

t_start = time.time()

# --- 1. merge chunks, sorted by position -----------------------------------
files = sorted(glob.glob(CACHE_GLOB))
print(f"merging {len(files)} chunk files")
positions, counts, labels, binomial_llr, pb_llr, depth = [], [], [], [], [], []
for fp in files:
    d = np.load(fp)
    positions.append(d["positions"])
    counts.append(d["counts"])
    labels.append(d["labels"])
    binomial_llr.append(d["binomial_llr"])
    pb_llr.append(d["pb_llr"])
    depth.append(d["depth"])

positions = np.concatenate(positions)
counts = np.concatenate(counts)
labels = np.concatenate(labels)
binomial_llr = np.concatenate(binomial_llr)
pb_llr = np.concatenate(pb_llr)
depth = np.concatenate(depth)

order = np.argsort(positions, kind="stable")
positions, counts, labels = positions[order], counts[order], labels[order]
binomial_llr, pb_llr, depth = binomial_llr[order], pb_llr[order], depth[order]

print(f"merged: {positions.size} loci total, region {positions.min()}-{positions.max()}")

# --- 2. frame filter (verbatim robustness_benchmark.py evaluate_arms) ------
frame = (labels == 0) | (labels == LABEL_SNP)
positions_f = positions[frame]
counts_f = counts[frame]
labels_f = labels[frame]
binomial_llr_f = binomial_llr[frame]
pb_llr_f = pb_llr[frame]
depth_f = depth[frame]
snp_truth = labels_f == LABEL_SNP

print(f"frame: {positions_f.size} loci scored, {int(snp_truth.sum())} truth SNP")

# --- 3. frozen calls (verbatim robustness_benchmark.py route_mask/evaluate_arms)
binomial_calls = binomial_llr_f >= FROZEN_BINOMIAL_THRESHOLD
pb_calls = pb_llr_f >= FROZEN_PB_THRESHOLD
routed = np.abs(binomial_llr_f - FROZEN_BINOMIAL_THRESHOLD) <= FROZEN_ROUTER_CUTOFF
cascade_calls = np.where(routed, pb_calls, binomial_calls)

n_routed = int(routed.sum())
print(f"PB-routed: {n_routed} / {positions_f.size} ({100*n_routed/positions_f.size:.3f}%)")

# --- 4. Method C genotyping (verbatim method_c_regression.py optimized_method_c)
def alt_allele(count_row, ref_base):
    base_counts = {"A": count_row[0], "C": count_row[1], "G": count_row[2], "T": count_row[3]}
    base_counts.pop(ref_base, None)
    if not base_counts:
        return None
    best_base, best_count = max(base_counts.items(), key=lambda kv: kv[1])
    if best_count <= 0:
        return None
    return best_base


def extract_records(call_mask, label):
    idxs = np.nonzero(call_mask)[0]
    records = []
    n_skipped = 0
    for i in idxs:
        pos0 = int(positions_f[i])
        ref_idx = int(counts_f[i][9])
        ref_base = REF_TOKENS[ref_idx] if 0 <= ref_idx < 5 else "N"
        if ref_base == "N":
            n_skipped += 1
            continue
        alt = alt_allele(counts_f[i], ref_base)
        if alt is None:
            n_skipped += 1
            continue
        base_counts = {"A": counts_f[i][0], "C": counts_f[i][1], "G": counts_f[i][2], "T": counts_f[i][3]}
        k = float(base_counts[alt])
        ref_n = float(base_counts[ref_base])
        n = k + ref_n
        d = int(counts_f[i][6])
        records.append(dict(idx=i, pos0=pos0, pos1=pos0 + 1, ref=ref_base, alt=alt, depth=d, k=k, n=n))
    print(f"[{label}] {len(records)} records ({n_skipped} skipped: ref=N or no alt support)")
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
    print(f"wrote {out_path}: {len(lines)} records")


results_meta = {"routing": {"n_scored": int(positions_f.size), "n_routed_pb": n_routed,
                             "pb_routed_fraction": n_routed / positions_f.size}}

for mask, label, out_name in [
    (pb_calls, "pb_only", f"ai_pb_only_{TAG}.vcf"),
    (cascade_calls, "cascade", f"ai_cascade_{TAG}.vcf"),
]:
    _t = time.time()
    records = extract_records(mask, label)
    _t_rec = time.time() - _t
    _t = time.time()
    gt, gq, forced = optimized_method_c(records)
    _t_gt = time.time() - _t
    _t = time.time()
    write_vcf_c(records, gt, gq, f"{OUT}/{out_name}", label)
    results_meta.setdefault("stage_seconds", {})[label] = {
        "record_extraction": _t_rec, "method_c_genotype": _t_gt, "vcf_write": time.time() - _t}
    results_meta[label] = {"n_calls": len(records), "n_forced_00_to_01": int(forced.sum())}

# --- 5. rescue analysis (devlog 14-17 definition, reused verbatim) ---------
# PB-rescuable: PB-only call is correct AND binomial-only call would have been wrong.
pb_correct = pb_calls == snp_truth
binom_correct = binomial_calls == snp_truth
rescuable = pb_correct & ~binom_correct
rescued_by_router = rescuable & routed
rescue_recall = float(rescued_by_router.sum() / rescuable.sum()) if rescuable.sum() > 0 else float("nan")

# false rescue: router sent to PB where binomial was already correct AND PB is wrong
false_rescue = routed & binom_correct & ~pb_correct
false_rescue_rate = float(false_rescue.sum() / routed.sum()) if routed.sum() > 0 else float("nan")
unnecessary_pb = routed & binom_correct  # PB used but binomial alone was already right
unnecessary_pb_rate = float(unnecessary_pb.sum() / routed.sum()) if routed.sum() > 0 else float("nan")

results_meta["rescue"] = {
    "pb_rescuable_loci": int(rescuable.sum()),
    "rescued_by_router": int(rescued_by_router.sum()),
    "rescue_recall": rescue_recall,
    "false_rescue_rate": false_rescue_rate,
    "unnecessary_pb_rate": unnecessary_pb_rate,
}

# cascade vs pb-only vs binomial-only raw agreement (detection layer, pre-hap.py)
cascade_correct = cascade_calls == snp_truth
results_meta["detection_layer_raw"] = {
    "binomial_only_accuracy": float(binom_correct.mean()),
    "pb_only_accuracy": float(pb_correct.mean()),
    "cascade_accuracy": float(cascade_correct.mean()),
}

with open(f"{OUT}/pipeline_meta_{TAG}.json", "w") as fh:
    json.dump(results_meta, fh, indent=2)

print(json.dumps(results_meta, indent=2))
print(f"pipeline wall time: {time.time() - t_start:.1f}s")

# --- 6. per-locus call arrays for the analysis script (additive; no decision logic) ---
np.savez_compressed(f"{OUT}/{TAG}_calls.npz", positions=positions_f, labels=labels_f,
                    binomial_llr=binomial_llr_f, pb_llr=pb_llr_f, depth=depth_f,
                    binomial_calls=binomial_calls, pb_calls=pb_calls, routed=routed,
                    cascade_calls=cascade_calls, counts=counts_f)
