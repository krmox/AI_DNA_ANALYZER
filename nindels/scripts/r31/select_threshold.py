#!/usr/bin/env python3
"""Choose the R3.1 switching thresholds from development runs ONLY (HG004 chr20/chr21), following
docs/R3.1_PROTOCOL.md. Refuses to read any run whose dataset role is not 'dev'.

usage: select_threshold.py OUT_DIR
Reads  <work_root>/runs/dev/<dataset>/<depthkey>/{A,AB}/happy.vcf.gz
Writes OUT_DIR/selection.json and OUT_DIR/selection.md
"""
import gzip
import json
import random
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import records
from common import ENV, SAMTOOLS, registry

BIN_EDGES = [0, 10, 15, 20, 25, 30, 40, 50, 70, 10_000]   # local depth bins [lo, hi)
MIN_TRUTH_PER_BIN = 300
N_BOOT, SEED, ALPHA = 2000, 31, 0.10
DEV = ["HG004_chr20", "HG004_chr21"]


# ---------------------------------------------------------------- pure functions (unit-tested)
def merge_sparse_bins(edges, truth_counts, minimum):
    """edges: list of boundaries (len n+1); truth_counts: per-bin truth totals. Merge bins below `minimum`
    into the neighbour with fewer truth records until all are >= minimum or one bin remains."""
    edges, cnt = list(edges), list(truth_counts)
    while len(cnt) > 1 and min(cnt) < minimum:
        i = cnt.index(min(cnt))
        j = i - 1 if i > 0 and (i == len(cnt) - 1 or cnt[i - 1] <= cnt[i + 1]) else i + 1
        lo, hi = sorted((i, j))
        cnt[lo:hi + 1] = [cnt[lo] + cnt[hi]]
        del edges[hi]
    return edges, cnt


def apply_rule(items):
    """items: list of (label, delta, ci_lo) ordered by increasing depth. Returns the index of the first item such
    that it and every later item has delta > 0 and ci_lo > 0, or None."""
    ok = [d is not None and d > 0 and lo is not None and lo > 0 for _, d, lo in items]
    for i in range(len(items)):
        if all(ok[i:]):
            return i
    return None


def percentile(xs, q):
    xs = sorted(xs)
    k = (len(xs) - 1) * q
    f = int(k)
    return xs[f] + (xs[min(f + 1, len(xs) - 1)] - xs[f]) * (k - f)


def bootstrap_delta(units, by_cfg, label_keys, rng):
    """units: list of unit ids. by_cfg[cfg][label][unit] = counts dict. Returns {label: (delta, lo, hi, p_pos)}
    for dF1 = F1(AB) - F1(A); the same resample is used for both configurations (paired)."""
    def total(cfg, label, sample):
        c = {"TPt": 0, "FN": 0, "TPq": 0, "FP": 0}
        for u in sample:
            x = by_cfg[cfg][label].get(u)
            if x:
                for k in c:
                    c[k] += x[k]
        return c
    out = {}
    pts = {l: records.f1(total("AB", l, units)) - records.f1(total("A", l, units)) for l in label_keys}
    draws = {l: [] for l in label_keys}
    for _ in range(N_BOOT):
        s = [units[rng.randrange(len(units))] for _ in units]
        for l in label_keys:
            draws[l].append(records.f1(total("AB", l, s)) - records.f1(total("A", l, s)))
    for l in label_keys:
        d = draws[l]
        out[l] = (pts[l], percentile(d, ALPHA / 2), percentile(d, 1 - ALPHA / 2), sum(x > 0 for x in d) / len(d))
    return out


# ---------------------------------------------------------------- data
def bin_of(depth, edges):
    for i in range(len(edges) - 1):
        if edges[i] <= depth < edges[i + 1]:
            return i
    return len(edges) - 2


def load_runs(reg, work):
    """{(ds, dk, cfg): rows}; also DP calibration pairs from the A runs."""
    data, dp_pairs = {}, []
    for ds in DEV:
        assert reg["datasets"][ds]["role"] == "dev", f"{ds} is not a development dataset"
        chrom = reg["datasets"][ds]["chrom"]
        lad = json.loads((work / ds / "ladder.json").read_text())["ladder"]
        for dk, info in lad.items():
            for cfg in ("A", "AB"):
                run = work / "runs" / "dev" / ds / dk / cfg
                if not (run / "metrics.json").exists():
                    print(f"missing {ds} {dk} {cfg}", file=sys.stderr)
                    continue
                cache = run / "records.json.gz"
                if cache.exists():
                    rows = json.load(gzip.open(cache, "rt"))
                else:
                    rows = records.load(run / "happy.vcf.gz", info["bam"], chrom)
                    json.dump(rows, gzip.open(cache, "wt"))
                data[(ds, dk, cfg)] = rows
                if cfg == "A":
                    dp_pairs += dp_calibration(run / "calls.pass.vcf.gz", info["bam"], chrom)
    # paired analysis: keep only (dataset, depth) points where both A and A+B finished
    full = {(ds, dk) for (ds, dk, c) in data if c == "A"} & {(ds, dk) for (ds, dk, c) in data if c == "AB"}
    dropped = sorted({(ds, dk) for (ds, dk, _) in data} - full)
    if dropped:
        print("incomplete pairs ignored:", dropped, file=sys.stderr)
    return {k: v for k, v in data.items() if (k[0], k[1]) in full}, dp_pairs


def dp_calibration(vcf, bam, chrom):
    pos_dp = []
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            for kv in f[7].split(";"):
                if kv.startswith("DP="):
                    pos_dp.append((int(f[1]), int(kv[3:])))
    dep = records.local_depths(bam, chrom, [p for p, _ in pos_dp])
    return [(dep.get(p, 0), d) for p, d in pos_dp]


def main():
    out_dir = Path(sys.argv[1]); out_dir.mkdir(parents=True, exist_ok=True)
    reg = registry(); work = Path(reg["work_root"])
    data, dp_pairs = load_runs(reg, work)
    result = {"protocol": "docs/R3.1_PROTOCOL.md", "n_boot": N_BOOT, "alpha": ALPHA, "seed": SEED}

    def analyse(datasets, tag):
        keys = [k for k in data if k[0] in datasets]
        units = sorted({(ds, r["block"]) for (ds, dk, cfg) in keys for r in data[(ds, dk, cfg)]})
        # ---- local depth bins
        truth_per_bin = [0] * (len(BIN_EDGES) - 1)
        for (ds, dk, cfg) in keys:
            if cfg == "A":
                for r in data[(ds, dk, cfg)]:
                    if r["t"]:
                        truth_per_bin[bin_of(r["depth"], BIN_EDGES)] += 1
        edges, tcounts = merge_sparse_bins(BIN_EDGES, truth_per_bin, MIN_TRUTH_PER_BIN)
        by_cfg = {c: defaultdict(lambda: defaultdict(lambda: {"TPt": 0, "FN": 0, "TPq": 0, "FP": 0})) for c in ("A", "AB")}
        for (ds, dk, cfg) in keys:
            for r in data[(ds, dk, cfg)]:
                u = (ds, r["block"])
                lb = f"{edges[bin_of(r['depth'], edges)]}-{edges[bin_of(r['depth'], edges) + 1]}"
                gb = dk
                for lab_space, lab in (("local", lb), ("ladder", gb), (f"ladder:{ds}", gb)):
                    c = by_cfg[cfg][(lab_space, lab)][u]
                    if r["t"] == "TP": c["TPt"] += 1
                    elif r["t"] == "FN": c["FN"] += 1
                    if r["q"] == "TP": c["TPq"] += 1
                    elif r["q"] == "FP": c["FP"] += 1
        labs_local = [("local", f"{edges[i]}-{edges[i + 1]}") for i in range(len(edges) - 1)]
        dks = sorted({dk for (_, dk, _) in keys}, key=lambda s: 1e9 if s == "native" else int(s[1:]))
        labs_ladder = [("ladder", dk) for dk in dks]
        rng = random.Random(SEED)
        local = bootstrap_delta(units, by_cfg, labs_local, rng)
        ladder = bootstrap_delta(units, by_cfg, labs_ladder, rng)
        per_ds = {}
        for ds in datasets:
            labs = [(f"ladder:{ds}", dk) for dk in dks if (f"ladder:{ds}", dk) in by_cfg["A"]]
            u_ds = [u for u in units if u[0] == ds]
            per_ds[ds] = bootstrap_delta(u_ds, by_cfg, labs, rng) if labs else {}
        # ---- rules
        loc_items = [(l[1], local[l][0], local[l][1]) for l in labs_local]
        li = apply_rule(loc_items)
        gi_items = [(l[1], ladder[l][0], ladder[l][1]) for l in labs_ladder]
        both_pos = []
        for (lab, dk) in labs_ladder:
            ok = all(per_ds[ds].get((f"ladder:{ds}", dk), (None,))[0] is not None and per_ds[ds][(f"ladder:{ds}", dk)][0] > 0 for ds in datasets)
            both_pos.append(ok)
        gi_items2 = [(lbl, d if both_pos[i] else None, lo) for i, (lbl, d, lo) in enumerate(gi_items)]
        gi = apply_rule(gi_items2)
        res = {
            "datasets": datasets, "units": len(units), "bin_edges": edges, "truth_per_bin": tcounts,
            "local": {l[1]: {"delta_f1": local[l][0], "ci_lo": local[l][1], "ci_hi": local[l][2], "p_pos": local[l][3]} for l in labs_local},
            "ladder": {l[1]: {"delta_f1": ladder[l][0], "ci_lo": ladder[l][1], "ci_hi": ladder[l][2], "p_pos": ladder[l][3],
                              "per_dataset": {ds: per_ds[ds].get((f"ladder:{ds}", l[1]), (None,))[0] for ds in datasets}} for l in labs_ladder},
            "T_local_lower_edge": None if li is None else edges[li],
            "T_global_ladder_key": None if gi is None else dks[gi],
        }
        return res

    result["all_dev"] = analyse(DEV, "all")
    result["loco"] = {ds: analyse([ds], ds) for ds in DEV}   # each chromosome alone (leave-one-chromosome-out view)
    # DP calibration: engine DP vs samtools local depth, in the 10-80x range where the engine cap (250) is irrelevant
    ratios = [dp / d for d, dp in dp_pairs if 10 <= d <= 80]
    result["ladder_targets"] = {dk: (float(dk[1:]) if dk != "native" else min(json.loads((work / ds / "ladder.json").read_text())["native_median_depth"] for ds in DEV) - 1)
                                for dk in {k[1] for k in data}}
    result["dp_calibration"] = {"n_calls": len(ratios), "median_dp_over_depth": statistics.median(ratios) if ratios else None}
    (out_dir / "selection.json").write_text(json.dumps(result, indent=1))
    print(json.dumps({k: result["all_dev"][k] for k in ("T_local_lower_edge", "T_global_ladder_key")}, indent=1))


if __name__ == "__main__":
    main()
