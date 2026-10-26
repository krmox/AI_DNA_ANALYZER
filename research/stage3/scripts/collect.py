#!/usr/bin/env python3
"""Aggregate raw Stage 3 runs into RESULTS.tsv and STRATIFIED_RESULTS.tsv (reads only; never alters raw outputs)."""
import csv, gzip, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import stage3 as s3  # noqa: E402

S3 = s3.S3
OUT = Path("/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/stage3")
SPARSE_TRUTH_SNPS = 100
CALLERS = {"varion": "Varion R1.1", "bcftools": "bcftools 1.23.1 (sharded x8)"}
CALLER_CONFIG = {
    "varion": "dnav2 e83d87fd; R1.model 63816707; cascade; --threads 8; edge-windows truth-free (defaults)",
    "bcftools": "mpileup -q20 -Q13 -d10000 -a FORMAT/AD | call -mv | view -v snps -i QUAL>=30; 32 shards, 8 parallel",
}
EXPOSURE = {
    ("HG002", "chr22"): "new contig, exposed sample (region-level only)", ("HG003", "chr22"): "new contig, exposed sample (region-level only)",
    ("HG004", "chr22"): "new contig, exposed sample (region-level only)", ("HG002", "chr18"): "new contig, exposed sample (region-level only)",
    ("HG003", "chr18"): "new contig, exposed sample (region-level only)", ("HG004", "chr18"): "new contig, exposed sample (region-level only)",
    ("HG005", "chr22"): "new sample (unrelated) + new contig", ("HG005", "chr18"): "new sample (unrelated) + new contig",
    ("HG004", "chr20"): "exposed (Stage 2 confirmatory; chr20 = R1 training contig on HG002)", ("HG004", "chr21"): "exposed (Stage 2 confirmatory)",
}


def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def overlap_len(a, b):
    i = j = total = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            total += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def read_bed(path, contig):
    opener = gzip.open if str(path).endswith(".gz") else open
    iv = []
    with opener(path, "rt") as f:
        for line in f:
            p = line.split("\t")
            if p[0] == contig:
                iv.append((int(p[1]), int(p[2])))
    return merge(iv)


_eval_cache = {}


def stratum_bases(sample, contig):
    """Bases of (stratum ∩ confident BED ∩ contig) for every selected stratum."""
    key = (sample, contig)
    if key in _eval_cache:
        return _eval_cache[key]
    conf = read_bed(s3.contig_bed(sample, contig), contig)
    out = {"*": sum(e - s for s, e in conf)}
    for line in open(S3 / "strat" / "strata.tsv"):
        name, rel = line.rstrip("\n").split("\t")
        out[name] = overlap_len(read_bed(S3 / "strat" / rel, contig), conf)
    _eval_cache[key] = out
    return out


def f1_counts(tp, fp, fn):
    return 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else float("nan")


def load_extended(path):
    return list(csv.DictReader(open(path)))


def num(x):
    return float(x) if x not in ("", "nan") else float("nan")


def parse_cond(cond):
    parts = cond.split("_")
    sample, contig, label = parts[0], parts[1], "_".join(parts[2:])
    nominal, seed = "", ""
    if label.endswith("x") or "x_s" in label:
        nominal = label.split("x")[0]
        seed = label.split("_s")[-1] if "_s" in label else ""
    return sample, contig, label, nominal, seed


def main():
    results_rows, strat_rows = [], []
    for run_dir in sorted((S3 / "runs").iterdir()):
        if not run_dir.is_dir() or not (run_dir / "calls.json").exists():
            continue
        cond = run_dir.name
        calls = json.load(open(run_dir / "calls.json"))
        sample, contig, label, nominal, seed = parse_cond(cond)
        cov_path = S3 / "meta" / f"{cond}_cov.json"
        observed = json.load(open(cov_path))["bed_mean_depth"] if cov_path.exists() else ""
        for caller in ("varion", "bcftools"):
            key = "varion" if caller == "varion" else "bcftools8"
            happy_dir = run_dir / f"happy_{caller}"
            ext = happy_dir / "h.extended.csv"
            call_rec = calls.get(key, {})
            row = dict(sample=sample, contig=contig, condition=label, nominal_coverage=nominal, seed=seed, observed_bed_coverage=observed,
                       caller=CALLERS[caller], caller_config=CALLER_CONFIG[caller], exposure=EXPOSURE.get((sample, contig), ""),
                       wall_s=call_rec.get("wall_s", ""), cpu_s=call_rec.get("cpu_s", ""),
                       peak_rss_mb=round(call_rec["peak_rss_kb"] / 1024, 1) if "peak_rss_kb" in call_rec else "",
                       threads=8, input_bytes=calls.get("bam_bytes", ""), loadavg_before=call_rec.get("loadavg_before", ""),
                       bam_sha256=calls.get("bam_sha256", ""), vcf_sha256=call_rec.get("vcf_sha256", ""))
            if not ext.exists():
                row.update(status="NO_HAPPY_OUTPUT", notes="hap.py missing or failed; see happy.log")
                results_rows.append(row)
                continue
            rows = load_extended(ext)
            bases = stratum_bases(sample, contig)
            for r in rows:
                if r["Filter"] != "PASS":
                    continue
                if r["Type"] == "INDEL" and r["Subset"] == "*":
                    row.update(indel_truth_total=r["TRUTH.TOTAL"], indel_tp=r["TRUTH.TP"], indel_fn=r["TRUTH.FN"], indel_fp=r["QUERY.FP"])
                if r["Type"] != "SNP":
                    continue
                tp, fn, fp = int(r["TRUTH.TP"]), int(r["TRUTH.FN"]), int(r["QUERY.FP"])
                qtp = int(r["QUERY.TP"])
                prec = qtp / (qtp + fp) if (qtp + fp) else float("nan")
                rec = tp / (tp + fn) if (tp + fn) else float("nan")
                if r["Subset"] == "*":
                    row.update(TP=tp, FP=fp, FN=fn, QUERY_TP=qtp, truth_total=r["TRUTH.TOTAL"], query_total=r["QUERY.TOTAL"],
                               precision=f"{prec:.6f}", recall=f"{rec:.6f}", F1_counts=f"{f1_counts(tp, fp, fn):.6f}",
                               F1_happy=r["METRIC.F1_Score"], FP_genotype=r["FP.gt"], FP_allele=r["FP.al"], status="OK", notes="")
                if r["Subset"] in bases:
                    strat_rows.append(dict(sample=sample, contig=contig, condition=label, nominal_coverage=nominal, seed=seed,
                                           caller=CALLERS[caller], stratum=r["Subset"], eval_bases=bases[r["Subset"]], truth_snps=r["TRUTH.TOTAL"],
                                           TP=tp, FN=fn, FP=fp, QUERY_TP=qtp, precision=f"{prec:.6f}", recall=f"{rec:.6f}",
                                           F1_counts=f"{f1_counts(tp, fp, fn):.6f}", FP_genotype=r["FP.gt"],
                                           sparse="yes" if int(r["TRUTH.TOTAL"]) < SPARSE_TRUTH_SNPS else "",
                                           stratification="GIAB genome-stratifications v3.6 GRCh38@all", truth="GIAB v4.2.1 GRCh38"))
            results_rows.append(row)
    cols = ["sample", "contig", "condition", "nominal_coverage", "seed", "observed_bed_coverage", "caller", "caller_config", "exposure",
            "TP", "FP", "FN", "QUERY_TP", "truth_total", "query_total", "precision", "recall", "F1_counts", "F1_happy",
            "FP_genotype", "FP_allele", "indel_truth_total", "indel_tp", "indel_fn", "indel_fp",
            "wall_s", "cpu_s", "peak_rss_mb", "threads", "input_bytes", "loadavg_before", "bam_sha256", "vcf_sha256", "status", "notes"]
    with open(OUT / "RESULTS.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(results_rows)
    scols = ["sample", "contig", "condition", "nominal_coverage", "seed", "caller", "stratum", "eval_bases", "truth_snps", "TP", "FN", "FP", "QUERY_TP",
             "precision", "recall", "F1_counts", "FP_genotype", "sparse", "stratification", "truth"]
    with open(OUT / "STRATIFIED_RESULTS.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, scols, delimiter="\t")
        w.writeheader()
        w.writerows(strat_rows)
    print(len(results_rows), "result rows,", len(strat_rows), "stratified rows")


if __name__ == "__main__":
    main()
