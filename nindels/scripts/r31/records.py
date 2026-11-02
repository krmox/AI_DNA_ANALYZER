"""Per-record outcome table from a hap.py run, with local depth taken from the BAM that produced the calls.

One record = one row of happy.vcf.gz inside the confident regions with BVT=INDEL. A row gives
  truth side: TP or FN   (BD of the TRUTH column)
  query side: TP or FP   (BD of the QUERY column)
plus 2-Mb block, local depth (reads overlapping POS), genotype class (BLT), indel kind/length, repeat context.
F1 is computed from counts as 2PR/(P+R) with P = TPq/(TPq+FP), R = TPt/(TPt+FN).
"""
import gzip
import subprocess
import tempfile

from common import ENV, SAMTOOLS

BLOCK = 2_000_000


def parse_info(s):
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def length_class(ref, alts):
    n = max(abs(len(a) - len(ref)) for a in alts if a not in ("*", "."))
    return "1-5" if n <= 5 else ("6-15" if n <= 15 else "16+")


def context(regions):
    for name in ("homopolymer_ge7", "tandem_repeat_le50", "tandem_repeat_51to200"):
        if name in regions:
            return name
    return "other"


def local_depths(bam, chrom, positions):
    """{pos: reads overlapping pos} (primary, non-duplicate, MAPQ >= 20), from samtools depth."""
    out = {}
    if not positions:
        return out
    with tempfile.NamedTemporaryFile("w", suffix=".bed") as bed:
        for p in sorted(set(positions)):
            bed.write(f"{chrom}\t{p - 1}\t{p}\n")
        bed.flush()
        r = subprocess.run([str(SAMTOOLS), "depth", "-a", "-Q", "20", "-b", bed.name, str(bam)],
                           capture_output=True, text=True, check=True, env=ENV)
    for line in r.stdout.splitlines():
        c, p, d = line.split("\t")
        out[int(p)] = int(d)
    return out


def load(happy_vcf, bam, chrom):
    rows = []
    with gzip.open(happy_vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = parse_info(f[7])
            if "CONF" not in info.get("Regions", ""):
                continue
            fmt = f[8].split(":")
            T, Q = dict(zip(fmt, f[9].split(":"))), dict(zip(fmt, f[10].split(":")))
            if T.get("BVT") != "INDEL" and Q.get("BVT") != "INDEL":
                continue
            bt = T.get("BD") if T.get("BD") in ("TP", "FN") else None
            bq = Q.get("BD") if Q.get("BD") in ("TP", "FP") else None
            if bt is None and bq is None:
                continue
            alts = [a for a in f[4].split(",") if a not in ("*", ".")]
            if not alts:
                continue
            kind = ("INS" if len(alts[0]) > len(f[3]) else "DEL") if len(alts) == 1 else "MULTI"
            gt = (T.get("BLT") if bt else Q.get("BLT")) or "NA"
            rows.append({"pos": int(f[1]), "t": bt, "q": bq, "kind": kind, "len": length_class(f[3], alts),
                         "ctx": context(info.get("Regions", "")), "gt": gt, "block": int(f[1]) // BLOCK,
                         "bk": Q.get("BK", ".")})
    dep = local_depths(bam, chrom, [r["pos"] for r in rows])
    for r in rows:
        r["depth"] = dep.get(r["pos"], 0)
    return rows


def counts(rows):
    c = {"TPt": 0, "FN": 0, "TPq": 0, "FP": 0}
    for r in rows:
        if r["t"] == "TP": c["TPt"] += 1
        elif r["t"] == "FN": c["FN"] += 1
        if r["q"] == "TP": c["TPq"] += 1
        elif r["q"] == "FP": c["FP"] += 1
    return c


def gt_concordance(rows):
    """TP / (TP + truth indels whose allele was found but called with the wrong genotype)."""
    tp = sum(1 for r in rows if r["t"] == "TP")
    err = sum(1 for r in rows if r["t"] == "FN" and r.get("bk") == "am")
    return (tp / (tp + err) if tp + err else None), err


def f1(c):
    p = c["TPq"] / (c["TPq"] + c["FP"]) if c["TPq"] + c["FP"] else 0.0
    r = c["TPt"] / (c["TPt"] + c["FN"]) if c["TPt"] + c["FN"] else 0.0
    return 2 * p * r / (p + r) if p + r else 0.0
