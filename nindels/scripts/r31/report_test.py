#!/usr/bin/env python3
"""Final-test report: baseline / A / A+B / adaptive R3.1 per sample, chromosome and depth.

usage: report_test.py PHASE OUT_MD OUT_JSON [--datasets A,B] [--bootstrap N]
Reads <work_root>/runs/PHASE/<dataset>/<depth>/<config>/{metrics.json,happy.vcf.gz} and writes
overall metrics, breakdowns (local depth, homopolymer, tandem repeat, insertion/deletion and length, genotype)
and paired dF1 with 2-Mb block-bootstrap intervals. Local depth bins are the ones fixed in the protocol.
"""
import argparse
import gzip
import json
import random
from collections import defaultdict
from pathlib import Path

import records
from common import registry
from select_threshold import BIN_EDGES, bin_of, percentile

CONFIGS = ["baseline", "A", "AB", "R31"]
LABEL = {"baseline": "Baseline", "A": "A", "AB": "A+B", "R31": "Adaptive R3.1"}


def key_depth(k):
    return 1e9 if k == "native" else int(k[1:])


def load_rows(run, bam, chrom):
    cache = run / "records.json.gz"
    if cache.exists():
        return json.load(gzip.open(cache, "rt"))
    rows = records.load(run / "happy.vcf.gz", bam, chrom)
    json.dump(rows, gzip.open(cache, "wt"))
    return rows


def metrics(rows):
    c = records.counts(rows)
    p = c["TPq"] / (c["TPq"] + c["FP"]) if c["TPq"] + c["FP"] else float("nan")
    r = c["TPt"] / (c["TPt"] + c["FN"]) if c["TPt"] + c["FN"] else float("nan")
    gc, gerr = records.gt_concordance(rows)
    return {"P": p, "R": r, "F1": records.f1(c), "FP": c["FP"], "FN": c["FN"], "TP": c["TPt"], "truth": c["TPt"] + c["FN"],
            "GTconc": gc, "GTerr": gerr}


def fmt(x, n=4):
    return "n/a" if x is None or x != x else f"{x:.{n}f}"


def paired_delta(rows_a, rows_b, n_boot, rng):
    """dF1 = F1(b) - F1(a) with a block bootstrap over 2-Mb blocks (same resample for both)."""
    blocks = sorted({r["block"] for r in rows_a} | {r["block"] for r in rows_b})
    ca, cb = defaultdict(lambda: dict(TPt=0, FN=0, TPq=0, FP=0)), defaultdict(lambda: dict(TPt=0, FN=0, TPq=0, FP=0))
    for rows, acc in ((rows_a, ca), (rows_b, cb)):
        for r in rows:
            c = acc[r["block"]]
            if r["t"] == "TP": c["TPt"] += 1
            elif r["t"] == "FN": c["FN"] += 1
            if r["q"] == "TP": c["TPq"] += 1
            elif r["q"] == "FP": c["FP"] += 1
    def tot(acc, sample):
        t = dict(TPt=0, FN=0, TPq=0, FP=0)
        for b in sample:
            for k in t:
                t[k] += acc[b][k] if b in acc else 0
        return t
    pt = records.f1(tot(cb, blocks)) - records.f1(tot(ca, blocks))
    d = []
    for _ in range(n_boot):
        s = [blocks[rng.randrange(len(blocks))] for _ in blocks]
        d.append(records.f1(tot(cb, s)) - records.f1(tot(ca, s)))
    return pt, percentile(d, 0.05), percentile(d, 0.95)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("phase"); ap.add_argument("out_md"); ap.add_argument("out_json")
    ap.add_argument("--datasets"); ap.add_argument("--bootstrap", type=int, default=1000)
    a = ap.parse_args()
    reg = registry(); work = Path(reg["work_root"])
    base = work / "runs" / a.phase
    datasets = a.datasets.split(",") if a.datasets else sorted(p.name for p in base.iterdir() if p.is_dir())
    rng = random.Random(31)
    out, md = {}, []
    for ds in datasets:
        d = reg["datasets"][ds]
        lad = json.loads((work / ds / "ladder.json").read_text())
        nat = lad["native_median_depth"]
        out[ds] = {"role": d["role"], "platform": d["platform"], "native_median_depth": nat, "depths": {}}
        md.append(f"\n## {ds}  (role: {d['role']}; {d['platform']}; native median depth {nat:.0f}x)\n")
        md.append("| Depth | Config | Precision | Recall | F1 | GT conc. | FP | FN | Wall s | CPU s | Peak RSS MB | Chosen |")
        md.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
        dks = sorted((p.name for p in (base / ds).iterdir() if p.is_dir()), key=key_depth)
        for dk in dks:
            rows_by = {}
            for cfg in CONFIGS:
                run = base / ds / dk / cfg
                if not (run / "metrics.json").exists():
                    continue
                m = json.loads((run / "metrics.json").read_text())
                bam = lad["ladder"][dk]["bam"]
                rows = load_rows(run, bam, d["chrom"])
                rows_by[cfg] = rows
                mt = metrics(rows)
                chosen = ""
                if cfg == "R31" and (run / "summary.json").exists():
                    s = json.loads((run / "summary.json").read_text())
                    chosen = ",".join(c["decision"]["config"] for c in s["contigs"])
                    mt["chosen"] = chosen
                mt.update(wall_s=m["wall_s"], cpu_s=m["cpu_s"], rss_mb=m["peak_rss_mb"], reps=m["reps"])
                out[ds]["depths"].setdefault(dk, {})[cfg] = {"overall": mt}
                md.append(f"| {dk} | {LABEL[cfg]} | {fmt(mt['P'])} | {fmt(mt['R'])} | {fmt(mt['F1'])} | {fmt(mt['GTconc'])} ({mt['GTerr']}) | "
                          f"{mt['FP']} | {mt['FN']} | {m['wall_s']:.1f} (n={m['reps']}) | {m['cpu_s']:.0f} | {m['peak_rss_mb']:.0f} | {chosen} |")
                # breakdowns
                def sub(f):
                    return metrics([r for r in rows if f(r)])
                bd = {"local_depth": {}, "context": {}, "kind": {}, "length": {}, "genotype": {}}
                edges = BIN_EDGES
                for i in range(len(edges) - 1):
                    bd["local_depth"][f"{edges[i]}-{edges[i + 1]}"] = sub(lambda r, i=i: bin_of(r["depth"], edges) == i)
                for c in ("homopolymer_ge7", "tandem_repeat_le50", "tandem_repeat_51to200", "other"):
                    bd["context"][c] = sub(lambda r, c=c: r["ctx"] == c)
                for k in ("INS", "DEL", "MULTI"):
                    bd["kind"][k] = sub(lambda r, k=k: r["kind"] == k)
                for k in ("1-5", "6-15", "16+"):
                    bd["length"][k] = sub(lambda r, k=k: r["len"] == k)
                for k in ("het", "homalt", "hetalt"):
                    bd["genotype"][k] = sub(lambda r, k=k: r["gt"] == k)
                out[ds]["depths"][dk][cfg]["breakdown"] = bd
            # paired differences
            pairs = [("A", "baseline"), ("AB", "A"), ("R31", "A"), ("R31", "AB")]
            dl = {}
            for b, a_ in pairs:
                if b in rows_by and a_ in rows_by:
                    dl[f"{LABEL[b]} - {LABEL[a_]}"] = paired_delta(rows_by[a_], rows_by[b], a.bootstrap, rng)
            out[ds]["depths"].setdefault(dk, {})["paired_dF1"] = dl
        md.append("\nPaired dF1 (90% block-bootstrap interval, 2-Mb blocks):\n")
        md.append("| Depth | Comparison | dF1 | 90% interval |\n|---|---|---|---|")
        for dk in dks:
            for name, (pt, lo, hi) in out[ds]["depths"].get(dk, {}).get("paired_dF1", {}).items():
                md.append(f"| {dk} | {name} | {pt:+.4f} | [{lo:+.4f}, {hi:+.4f}] |")
        # breakdown tables
        for dim, title in (("context", "Repeat context"), ("kind", "Insertions / deletions"), ("length", "Indel length (bp)"),
                           ("genotype", "Genotype"), ("local_depth", "Local depth (reads at the position)")):
            md.append(f"\n### {ds}: {title} (F1 [truth n]; recall / precision)\n")
            cols = list(next(iter(out[ds]["depths"].values()))[next(c for c in CONFIGS if c in next(iter(out[ds]["depths"].values()))) ]["breakdown"][dim])
            md.append("| Depth | Config | " + " | ".join(cols) + " |")
            md.append("|---|---|" + "---|" * len(cols))
            for dk in dks:
                for cfg in CONFIGS:
                    e = out[ds]["depths"].get(dk, {}).get(cfg)
                    if not e:
                        continue
                    cells = []
                    for c in cols:
                        m = e["breakdown"][dim][c]
                        cells.append("-" if m["truth"] < 30 else f"{fmt(m['F1'], 3)} [{m['truth']}]; {fmt(m['R'], 3)}/{fmt(m['P'], 3)}")
                    md.append(f"| {dk} | {LABEL[cfg]} | " + " | ".join(cells) + " |")
    Path(a.out_md).write_text("\n".join(md) + "\n")
    Path(a.out_json).write_text(json.dumps(out, indent=1))
    print("\n".join(md[:30]))


if __name__ == "__main__":
    main()
