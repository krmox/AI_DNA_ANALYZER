#!/usr/bin/env python3
"""Benchmark tables and error taxonomy from hap.py outputs (glue only).

usage: analyze.py RUNS_DIR TRUTH_VCF CONF_BED OUT_MD tag [tag ...]

For each tag reads RUNS_DIR/<tag>/{time_*.json, happy.summary.csv,
happy.extended.csv, happy.vcf.gz, calls.vcf.gz} and writes markdown tables:
  1. caller table (median wall/CPU/RSS, SNP/INDEL F1, precision, recall)
  2. INDEL PASS metrics by length class, genotype, repeat stratum
  3. for nindels tags: error taxonomy and concrete FP/FN examples
"""
import bisect
import csv
import glob
import gzip
import json
import statistics
import sys
from collections import Counter, defaultdict


def median_time(d):
    recs = [json.load(open(f)) for f in sorted(glob.glob(f"{d}/time_*.json"))]
    if not recs:
        return None
    med = lambda k: statistics.median(r[k] for r in recs)
    return {"n": len(recs), "wall": med("wall_s"), "cpu": med("cpu_s"), "rss": max(r["peak_rss_mb"] for r in recs)}


def summary(d):
    out = {}
    for r in csv.DictReader(open(f"{d}/happy.summary.csv")):
        if r["Filter"] == "PASS":
            out[r["Type"]] = r
    return out


def extended(d):
    return [r for r in csv.DictReader(open(f"{d}/happy.extended.csv")) if r["Filter"] == "PASS" and r["Type"] == "INDEL"]


def f(x, nd=4):
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return "n/a"


def vcf_records(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            yield line.rstrip("\n").split("\t")


def info_dict(s):
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        elif kv:
            d[kv] = True
    return d


def load_bed(path, chrom):
    iv = []
    for line in open(path):
        p = line.split()
        if p[0] == chrom:
            iv.append((int(p[1]), int(p[2])))
    iv.sort()
    return iv


def in_bed(iv, pos0):
    i = bisect.bisect_right(iv, (pos0, float("inf"))) - 1
    return i >= 0 and iv[i][0] <= pos0 < iv[i][1]


def genotype_concordance(d):
    """TP / (TP + truth indels whose allele was found but genotype differs), confident regions."""
    tp = gt_err = 0
    for r in vcf_records(f"{d}/happy.vcf.gz"):
        if "CONF" not in info_dict(r[7]).get("Regions", ""):
            continue
        fmt = r[8].split(":")
        T = dict(zip(fmt, r[9].split(":")))
        Q = dict(zip(fmt, r[10].split(":")))
        if T.get("BVT") != "INDEL":
            continue
        if T.get("BD") == "TP":
            tp += 1
        elif T.get("BD") == "FN" and Q.get("BK") == "am" and r[6] in (".", "PASS"):
            gt_err += 1
    return tp, gt_err


def taxonomy(d, truth_vcf, conf, chrom, lines):
    # nindels records (all filters) for FN classification
    calls = defaultdict(list)
    for c in vcf_records(f"{d}/calls.vcf.gz"):
        calls[int(c[1])].append(c)
    call_pos = sorted(calls)

    # truth: clustered = another truth variant within 10 bp
    truth_pos = []
    for t in vcf_records(truth_vcf):
        if t[0] == chrom:
            truth_pos.append(int(t[1]))
    truth_pos.sort()

    def clustered(p):
        i = bisect.bisect_left(truth_pos, p)
        return any(0 < abs(truth_pos[j] - p) <= 10 for j in (i - 1, i + 1) if 0 <= j < len(truth_pos))

    def matching_pass(p, ref, alts, w=15):
        """Nearest PASS nindels call with the same length change (hap.py re-normalises POS)."""
        dl = {len(a) - len(ref) for a in alts}
        best = None
        for c in near_calls(p, w):
            if c[6] != "PASS":
                continue
            if any(len(a) - len(c[3]) in dl for a in c[4].split(",")):
                if best is None or abs(int(c[1]) - p) < abs(int(best[1]) - p):
                    best = c
        return [best] if best else []

    def near_calls(p, w=10):
        i = bisect.bisect_left(call_pos, p - w)
        out = []
        while i < len(call_pos) and call_pos[i] <= p + w:
            out += calls[call_pos[i]]
            i += 1
        return out

    fn_cls, fp_cls = Counter(), Counter()
    fn_ex, fp_ex = defaultdict(list), defaultdict(list)
    clus = Counter()
    pure_fp_feats = defaultdict(list)
    tp_feats = defaultdict(list)
    long_truth = Counter()
    fp_gt = Counter()
    for r in vcf_records(f"{d}/happy.vcf.gz"):
        fmt = r[8].split(":")
        T = dict(zip(fmt, r[9].split(":")))
        Q = dict(zip(fmt, r[10].split(":")))
        pos = int(r[1])
        regions = info_dict(r[7]).get("Regions", "")
        ctx = "HP>=7" if "homopolymer_ge7" in regions else ("TR<=50" if "tandem_repeat_le50" in regions else
              ("TR51-200" if "tandem_repeat_51to200" in regions else "other"))
        is_indel_t = T.get("BVT") == "INDEL"
        is_indel_q = Q.get("BVT") == "INDEL"
        alts = r[4].split(",")
        maxlen = max(abs(len(a) - len(r[3])) for a in alts)
        if is_indel_t and T.get("BD") in ("TP", "FN"):
            clus[(clustered(pos), T["BD"])] += 1
            if maxlen > 50:
                long_truth[T["BD"]] += 1
        if is_indel_t and T.get("BD") == "FN":
            near = near_calls(pos)
            if Q.get("BD") == "FP" and Q.get("BK") == "am" and q_pass:
                k = "genotype wrong (allele found, GT differs)"
            elif Q.get("BD") == "FP" and Q.get("BK") == "lm" and q_pass:
                k = "allele/representation differs (nearby PASS call)"
            elif not near:
                k = "never a candidate (discovery miss)"
            elif not any(c[6] == "PASS" for c in near):
                k = "candidate found but filtered (LowProb/UNCERTAIN)"
            elif T.get("BLT") == "hetalt":
                k = "multi-allelic site (truth 1/2), alleles genotyped separately"
            else:
                k = "different allele/representation called nearby"
            if maxlen > 50:
                k += " [>50bp]"
            fn_cls[k] += 1
            if len(fn_ex[k]) < 4:
                ndesc = "; ".join(f"{c[3]}>{c[4]} {c[6]} {c[7][:70]}" for c in near[:2]) or "no candidate"
                fn_ex[k].append(f"{chrom}:{pos} {r[3]}>{r[4]} truth GT {T.get('GT')} [{ctx}] -> {ndesc}")
        q_pass = r[6] in (".", "PASS")  # hap.py copies the query FILTER; only PASS records are calls
        if "CONF" not in regions:
            continue
        if is_indel_q and Q.get("BD") == "FP" and q_pass:
            if Q.get("BK") == "am":
                k = "genotype error (allele right, GT wrong)"
            elif Q.get("BK") == "lm":
                k = "allele/representation mismatch near a truth variant"
            else:
                k = "no truth variant nearby (pure FP)"
            k2 = f"{k} | {ctx}"
            fp_cls[k2] += 1
            c = matching_pass(pos, r[3], alts)
            inf = info_dict(c[0][7]) if c else {}
            fp_gt[(k.split(" (")[0], Q.get("GT"))] += 1
            if k.startswith("no truth"):
                for key in ("NALT", "NREF", "F", "HQ", "QMIN"):
                    if key in inf:
                        pure_fp_feats[key].append(float(inf[key]))
            if len(fp_ex[k2]) < 3:
                fp_ex[k2].append(f"{chrom}:{pos} {r[3]}>{r[4]} called {Q.get('GT')} truth {T.get('GT')} "
                                 f"NALT={inf.get('NALT')} NREF={inf.get('NREF')} F={inf.get('F')} HQ={inf.get('HQ')}")
        if is_indel_q and Q.get("BD") == "TP" and q_pass:
            c = matching_pass(pos, r[3], alts)
            if c:
                inf = info_dict(c[0][7])
                for key in ("NALT", "NREF", "F", "HQ", "QMIN"):
                    if key in inf:
                        tp_feats[key].append(float(inf[key]))

    lines.append("\n**False negatives by cause**\n\n| Cause | Count |\n|---|---|")
    for k, v in fn_cls.most_common():
        lines.append(f"| {k} | {v} |")
    lines.append("\n**False positives by cause and repeat context**\n\n| Cause \\| context | Count |\n|---|---|")
    for k, v in fp_cls.most_common():
        lines.append(f"| {k.replace('|', '\\|')} | {v} |")
    lines.append("\n**FP called genotype**\n")
    for (k, g), v in sorted(fp_gt.items(), key=lambda x: -x[1]):
        lines.append(f"- {k}, called {g}: {v}")
    lines.append("\n**Recall: clustered (another truth variant within 10 bp) vs isolated**\n")
    for c in (False, True):
        tp, fn = clus[(c, "TP")], clus[(c, "FN")]
        if tp + fn:
            lines.append(f"- {'clustered' if c else 'isolated'}: {tp}/{tp + fn} = {tp / (tp + fn):.4f}")
    if long_truth:
        lines.append(f"- truth indels > 50 bp: TP {long_truth['TP']}, FN {long_truth['FN']}")
    lines.append("\n**Median features: TP calls vs pure-FP calls**\n\n| Feature | TP | pure FP |\n|---|---|---|")
    for key in ("NALT", "NREF", "F", "HQ", "QMIN"):
        a, b = tp_feats.get(key), pure_fp_feats.get(key)
        lines.append(f"| {key} | {f(statistics.median(a), 2) if a else 'n/a'} | {f(statistics.median(b), 2) if b else 'n/a'} |")
    if pure_fp_feats.get("NALT"):
        n2 = sum(1 for x in pure_fp_feats["NALT"] if x <= 2)
        lines.append(f"\nPure FPs with NALT <= 2: {n2}/{len(pure_fp_feats['NALT'])}")
    lines.append("\n**Examples**\n")
    for k, ex in list(fn_ex.items()) + list(fp_ex.items()):
        lines.append(f"- *{k}*")
        for e in ex:
            lines.append(f"  - `{e}`")


def gc_str(d):
    tp, e = genotype_concordance(d)
    return f"{tp / (tp + e):.4f} ({e} GT err)" if tp + e else "n/a"


def main():
    runs, truth, conf, out_md = sys.argv[1:5]
    tags = sys.argv[5:]
    chrom = "chr20"
    lines = ["## Benchmark: HG004 chr20, GIAB v4.2.1 confident regions, hap.py (xcmp)\n",
             "| Caller | Wall (s) | CPU (s) | Peak RSS (MB) | SNP F1 | INDEL F1 | INDEL Precision | INDEL Recall | INDEL TP | INDEL FP | INDEL FN | GT concordance |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for t in tags:
        d = f"{runs}/{t}"
        tm = median_time(d)
        s = summary(d)
        snp, ind = s.get("SNP", {}), s.get("INDEL", {})
        snp_f1 = f(snp.get("METRIC.F1_Score")) if snp.get("QUERY.TOTAL", "0") != "0" else "N/A (no SNP calls)"
        wall = f"{tm['wall']:.1f} (n={tm['n']})" if tm else "unknown"
        lines.append(f"| {t} | {wall} | {f(tm['cpu'], 1) if tm else 'unknown'} | {f(tm['rss'], 0) if tm else 'unknown'} | "
                     f"{snp_f1} | {f(ind.get('METRIC.F1_Score'))} | {f(ind.get('METRIC.Precision'))} | "
                     f"{f(ind.get('METRIC.Recall'))} | {ind.get('TRUTH.TP')} | {ind.get('QUERY.FP')} | {ind.get('TRUTH.FN')} | "
                     f"{gc_str(d)} |")

    lines.append("\n## INDEL (PASS) by length class, genotype and repeat context\n")
    lines.append("| Caller | Stratum | Truth | Recall | Precision | F1 |\n|---|---|---|---|---|---|")
    for t in tags:
        for r in extended(f"{runs}/{t}"):
            st, sub, gt = r["Subtype"], r["Subset"], r["Genotype"]
            keep = (sub == "*" and gt == "*" and st != "*") or (sub == "*" and st == "*" and gt != "*") or \
                   (st == "*" and gt == "*" and sub not in ("*", "TS_boundary", "TS_contained"))
            if not keep or r["TRUTH.TOTAL"] in ("0", "0.0"):
                continue
            name = st if st != "*" else (f"GT {gt}" if gt != "*" else sub)
            lines.append(f"| {t} | {name} | {r['TRUTH.TOTAL']} | {f(r['METRIC.Recall'])} | {f(r['METRIC.Precision'])} | {f(r['METRIC.F1_Score'])} |")

    for t in tags:
        if t.startswith("nindels"):
            lines.append(f"\n## Error taxonomy: {t}\n")
            taxonomy(f"{runs}/{t}", truth, conf, chrom, lines)

    open(out_md, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
