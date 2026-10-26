#!/usr/bin/env python3
"""Stage 3 evaluation driver. Does not modify the frozen caller; it only runs it and the comparator.

Subcommands
  derive  SAMPLE CONTIG TARGET SEED     make a downsampled BAM from the native BAM (deterministic samtools -s)
  call    SAMPLE CONTIG COND TAG        run Varion R1.1 and bcftools on bam/<COND>.bam, timed, VCFs to runs/<TAG>/
  happy   SAMPLE CONTIG TAG CALLER      hap.py against the matching GIAB truth, with stratifications
  cov     SAMPLE CONTIG COND            observed coverage inside the confident BED (samtools bedcov)
"""
import argparse, hashlib, json, os, shutil, subprocess, sys, time
from pathlib import Path

S3 = Path("/mnt/archive/AI_DNA_ANALYZER_stage3")
REF_DIR = Path("/mnt/archive/AI_DNA_ANALYZER_benchmark/reference")
REF = REF_DIR / "GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
DNAV2 = S3 / "build" / "dnav2"
MODEL = Path("/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/stage3/frozen/R1.model")
HAPPY_IMAGE = "quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0"
THREADS = 8
BCFTOOLS_SHARDS = 32

TRUTH = {
    "HG002": (Path("/mnt/archive/AI_DNA_ANALYZER_benchmark/truth_HG002_GRCh38_v4.2.1"), "HG002_GRCh38_1_22_v4.2.1_benchmark", "_noinconsistent.bed"),
    "HG003": (Path("/mnt/archive/AI_DNA_ANALYZER_final_eval/truth"), "HG003_GRCh38_1_22_v4.2.1_benchmark", "_noinconsistent.bed"),
    "HG004": (Path("/mnt/archive/AI_DNA_ANALYZER_stage2/truth"), "HG004_GRCh38_1_22_v4.2.1_benchmark", "_noinconsistent.bed"),
    "HG005": (S3 / "truth", "HG005_GRCh38_1_22_v4.2.1_benchmark", ".bed"),
    "HG006": (S3 / "truth", "HG006_GRCh38_1_22_v4.2.1_benchmark", ".bed"),
}


def truth_files(sample):
    d, stem, bed_suffix = TRUTH[sample]
    return d, f"{stem}.vcf.gz", f"{stem}{bed_suffix}" if bed_suffix.startswith("_") else f"{stem}.bed"


def contig_length(contig):
    for line in open(str(REF) + ".fai"):
        f = line.split("\t")
        if f[0] == contig:
            return int(f[1])
    raise SystemExit(f"contig {contig} not in reference index")


def contig_bed(sample, contig):
    """Confident-region BED restricted to one contig (cached)."""
    d, _, bed = truth_files(sample)
    out = S3 / "meta" / f"{sample}_{contig}.bed"
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w") as w:
            for line in open(d / bed):
                if line.split("\t", 1)[0] == contig:
                    w.write(line)
    return out


def timed(cmd, time_file, shell=False):
    """Run under GNU time -v; returns the parsed resource record."""
    full = ["/usr/bin/time", "-v", "-o", str(time_file)] + (["bash", "-c", cmd] if shell else cmd)
    load_before = open("/proc/loadavg").read().split()[0]
    t0 = time.time()
    p = subprocess.run(full, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    rec = parse_time(time_file)
    rec.update(rc=p.returncode, loadavg_before=float(load_before), stderr_tail=p.stderr[-300:], stdout=p.stdout)
    return rec


def parse_time(path):
    out = {}
    for line in open(path):
        line = line.strip()
        if line.startswith("Elapsed (wall clock)"):
            parts = line.rsplit(": ", 1)[1].split(":")
            secs = 0.0
            for x in parts:
                secs = secs * 60 + float(x)
            out["wall_s"] = secs
        elif line.startswith("User time"):
            out["user_s"] = float(line.rsplit(": ", 1)[1])
        elif line.startswith("System time"):
            out["sys_s"] = float(line.rsplit(": ", 1)[1])
        elif line.startswith("Maximum resident set size"):
            out["peak_rss_kb"] = int(line.rsplit(": ", 1)[1])
    out["cpu_s"] = out.get("user_s", 0.0) + out.get("sys_s", 0.0)
    return out


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def cmd_derive(a):
    native = S3 / "bam" / f"{a.sample}_{a.contig}_native.bam"
    meta = json.load(open(S3 / "meta" / f"{a.sample}_{a.contig}_native_cov.json"))
    fraction = a.target / meta["bed_mean_depth"]
    if fraction >= 1.0:
        raise SystemExit(f"target {a.target}x exceeds native {meta['bed_mean_depth']:.1f}x")
    sflag = f"{a.seed}.{round(fraction * 10000):04d}"
    out = S3 / "bam" / f"{a.sample}_{a.contig}_{a.target}x_s{a.seed}.bam"
    cmd = ["samtools", "view", "-b", "-s", sflag, "-@", "4", "-o", str(out), str(native)]
    subprocess.run(cmd, check=True)
    subprocess.run(["samtools", "index", str(out)], check=True)
    rec = dict(sample=a.sample, contig=a.contig, target=a.target, seed=a.seed, samtools_s=sflag,
               fraction=fraction, native_bed_depth=meta["bed_mean_depth"], source=str(native),
               command=" ".join(cmd), bam=str(out), bytes=out.stat().st_size, sha256=sha256(out))
    with open(S3 / "logs" / "derive.jsonl", "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec))


def cmd_cov(a):
    bam = S3 / "bam" / f"{a.cond}.bam"
    bed = contig_bed(a.sample, a.contig)
    total = 0
    bases = 0
    for line in subprocess.run(["samtools", "bedcov", str(bed), str(bam)], capture_output=True, text=True, check=True).stdout.splitlines():
        f = line.split("\t")
        total += int(f[2]) - int(f[1])
        bases += int(f[3])
    rec = dict(cond=a.cond, bed_bases=total, bed_mean_depth=bases / total)
    (S3 / "meta").mkdir(exist_ok=True)
    json.dump(rec, open(S3 / "meta" / f"{a.cond}_cov.json", "w"))
    print(json.dumps(rec))


def run_varion(a, out):
    d = out / "varion"
    d.mkdir(parents=True, exist_ok=True)
    cmd = [str(DNAV2), "run", "--bam", str(S3 / "bam" / f"{a.cond}.bam"), "--ref", str(REF), "--bed", str(contig_bed(a.sample, a.contig)),
           "--out-dir", str(d), "--contig", a.contig, "--contig-len", str(contig_length(a.contig)), "--sample", a.sample,
           "--tag", a.tag, "--threads", str(THREADS), "--mode", "cascade", "--model", str(MODEL)]
    rec = timed(cmd, d / "time.txt")
    rec["command"] = " ".join(cmd)
    (d / "stdout.json").write_text(rec.pop("stdout"))
    vcf = d / f"ai_cascade_{a.tag}.vcf"
    if rec["rc"] == 0:
        subprocess.run(["bcftools", "view", "-Oz", "-o", str(out / "varion.vcf.gz"), str(vcf)], check=True)
        subprocess.run(["bcftools", "index", "-f", "-t", str(out / "varion.vcf.gz")], check=True)
        rec["vcf_sha256"] = sha256(vcf)
        vcf.unlink()
        for extra in d.glob("routed_loci_*.tsv"):
            extra.unlink()
    return rec


def run_bcftools(a, out, parallel):
    d = out / ("bcftools" if parallel > 1 else "bcftools_1t")
    d.mkdir(parents=True, exist_ok=True)
    shards = BCFTOOLS_SHARDS if parallel > 1 else 1
    script = Path(__file__).with_name("bcftools_sharded.sh")
    cmd = [str(script), str(S3 / "bam" / f"{a.cond}.bam"), str(REF), a.contig, str(contig_length(a.contig)), str(d), str(shards), str(parallel)]
    rec = timed(cmd, d / "time.txt")
    rec["command"] = " ".join(cmd)
    rec.pop("stdout")
    if rec["rc"] == 0:
        rec["vcf_sha256"] = sha256(d / "bcftools.vcf.gz")
    return rec


def cmd_call(a):
    out = S3 / "runs" / a.tag
    out.mkdir(parents=True, exist_ok=True)
    result = dict(sample=a.sample, contig=a.contig, cond=a.cond, tag=a.tag,
                  bam_sha256=sha256(S3 / "bam" / f"{a.cond}.bam") if a.hash_bam else None,
                  bam_bytes=(S3 / "bam" / f"{a.cond}.bam").stat().st_size)
    for caller in a.callers.split(","):
        if caller == "varion":
            result["varion"] = run_varion(a, out)
        elif caller == "bcftools8":
            result["bcftools8"] = run_bcftools(a, out, THREADS)
        elif caller == "bcftools1":
            result["bcftools1"] = run_bcftools(a, out, 1)
        result_path = out / "calls.json"
        existing = json.load(open(result_path)) if result_path.exists() else {}
        existing.update(result)
        json.dump(existing, open(result_path, "w"), indent=1)
    print(json.dumps({k: v for k, v in result.items() if isinstance(v, dict) and "wall_s" in v}, indent=1))


def cmd_happy(a):
    out = S3 / "runs" / a.tag
    link = out / "bcftools.vcf.gz"
    if a.caller == "bcftools" and not link.exists():
        for suffix in ("", ".tbi"):
            (out / f"bcftools.vcf.gz{suffix}").symlink_to(f"bcftools/bcftools.vcf.gz{suffix}")
    vcf = out / ("varion.vcf.gz" if a.caller == "varion" else "bcftools.vcf.gz")
    d, tvcf, tbed = truth_files(a.sample)
    hd = out / f"happy_{a.caller}"
    hd.mkdir(exist_ok=True)
    log = hd / "happy.log"
    cmd = ["podman", "--root", "/dev/shm/podman_root", "--storage-driver", "vfs", "run", "--rm", "--security-opt", "label=disable",
           "-v", f"{d}:/truth:ro", "-v", f"{REF_DIR}:/ref:ro", "-v", f"{out}:/q:ro", "-v", f"{hd}:/o", "-v", f"{S3/'strat'}:/strat:ro",
           HAPPY_IMAGE, "/usr/local/bin/hap.py", f"/truth/{tvcf}", f"/q/{vcf.name}", "-f", f"/truth/{tbed}",
           "-r", f"/ref/{REF.name}", "-o", "/o/h", "-l", a.contig, "--threads", str(a.threads),
           "--stratification", "/strat/strata.tsv"]
    t0 = time.time()
    p = subprocess.run(cmd, stdout=open(log, "w"), stderr=subprocess.STDOUT)
    json.dump(dict(rc=p.returncode, wall_s=time.time() - t0, command=" ".join(cmd)), open(hd / "happy_run.json", "w"))
    print(a.tag, a.caller, "rc", p.returncode, f"{time.time()-t0:.0f}s")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("derive"); p.add_argument("sample"); p.add_argument("contig"); p.add_argument("target", type=int); p.add_argument("seed", type=int); p.set_defaults(f=cmd_derive)
    p = sub.add_parser("cov"); p.add_argument("sample"); p.add_argument("contig"); p.add_argument("cond"); p.set_defaults(f=cmd_cov)
    p = sub.add_parser("call"); p.add_argument("sample"); p.add_argument("contig"); p.add_argument("cond"); p.add_argument("tag")
    p.add_argument("--callers", default="varion,bcftools8"); p.add_argument("--hash-bam", action="store_true"); p.set_defaults(f=cmd_call)
    p = sub.add_parser("happy"); p.add_argument("sample"); p.add_argument("contig"); p.add_argument("tag"); p.add_argument("caller")
    p.add_argument("--threads", type=int, default=4); p.set_defaults(f=cmd_happy)
    a = ap.parse_args()
    a.f(a)


if __name__ == "__main__":
    main()
