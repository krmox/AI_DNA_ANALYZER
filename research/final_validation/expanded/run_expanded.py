#!/usr/bin/env python3
"""EXPANDED benchmark driver (HG001 or HG007 via env SAMPLE; both PREVIOUSLY EVALUATED, not blind). Derived from independent_genome/run_campaign.py with only sample/path parametrisation; frozen callers untouched.

Protocol = research/stage3 + stage4 (unchanged):
  * downsample native contig BAM to TARGET x by BED-mean depth: samtools view -s SEED.FRAC(4 digits)
  * callers (same BAM, same BED): R2 frozen (dnav2 + frozen_R2.model --r2), R1.1 frozen (dnav2 R1.1 + R1.model),
    bcftools baseline (research/stage3/scripts/bcftools_sharded.sh, 32 shards)
  * scoring: hap.py 0.3.15 (xcmp), confident-region BED, SNP PASS
Nothing here modifies the frozen source, binaries, models or thresholds. Run: SAMPLE=HG001 run_expanded.py CONTIG TARGET SEED [--callers r2,r11,bcftools]
TARGET may be 'native'.  Appends one JSON line per (condition, caller) to RESULTS.jsonl (incl. per-2Mb-block SNP TP/FP/FN for the preregistered bootstrap); skips finished ones.
"""
import csv, hashlib, json, os, resource, subprocess, sys, time
from pathlib import Path

HOME = Path.home()
VT = HOME / "varion_mamba/envs/vt"
HP = HOME / "varion_mamba/envs/hp"
REPO = Path("/home/kramkost/Projects/Varion_DNA/AI_DNA_ANALYZER")
SAMPLE = os.environ["SAMPLE"]  # HG001 | HG007
D = Path("/mnt/hdd/AI_DNA_ANALYZER_hg001" if SAMPLE == "HG001" else "/mnt/hdd/AI_DNA_ANALYZER_hg007")
BAMDIR = Path(os.environ.get("BAMDIR", str(D / "bam")))  # native + downsampled BAMs
REF = Path("/mnt/hdd/AI_DNA_ANALYZER_benchmark/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna")
TRUTH_VCF = D / f"truth/{SAMPLE}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED = D / f"truth/{SAMPLE}_GRCh38_1_22_v4.2.1_benchmark.bed"
R2_BIN = "/mnt/hdd/AI_DNA_ANALYZER_r2/build/dnav2"
R2_MODEL = REPO / "research/stage4/frozen_R2.model"
R11_BIN = "/mnt/hdd/AI_DNA_ANALYZER_stage3/build/dnav2"
R11_MODEL = REPO / "research/stage3/frozen/R1.model"
BCF_SH = REPO / "research/stage3/scripts/bcftools_sharded.sh"
OUT = REPO / "research/final_validation/expanded"
THREADS = 8
ENV = dict(os.environ, PATH=f"{VT}/bin:{os.environ['PATH']}")
ENV_DNA = dict(ENV, LD_LIBRARY_PATH=f"{VT}/lib")
ENV_HAP = dict(os.environ, PATH=f"{HP}/bin:{VT}/bin:{os.environ['PATH']}")


def run(cmd, env, log=None):
    """Run, return dict(rc, wall_s, user_s, sys_s, cpu_s, peak_rss_kb) from wait4 rusage of the child tree."""
    t0 = time.time(); la = os.getloadavg()[0]
    with open(log or os.devnull, "w") as lf:
        p = subprocess.Popen(cmd, env=env, stdout=lf, stderr=subprocess.STDOUT)
        _, status, ru = os.wait4(p.pid, 0)
    rc = os.waitstatus_to_exitcode(status)
    return dict(rc=rc, loadavg_before=la, wall_s=time.time() - t0, user_s=ru.ru_utime, sys_s=ru.ru_stime, cpu_s=ru.ru_utime + ru.ru_stime, peak_rss_kb=ru.ru_maxrss)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 24), b""):
            h.update(c)
    return h.hexdigest()


def block_counts(vcf, block=2_000_000):
    """Per 2-Mb block SNP TP/FP (query) and FN (truth) from the hap.py annotated VCF (BD/BVT) — preregistered bootstrap input."""
    text = subprocess.run(["bcftools", "query", "-f", "%POS\t[%BD,%BVT\t]\n", str(vcf)], env=ENV, capture_output=True, text=True, check=True).stdout
    blocks = {}
    for line in text.splitlines():
        pos, truth, query = line.rstrip("\t").split("\t")
        acc = blocks.setdefault((int(pos) - 1) // block, [0, 0, 0])
        bd, vt = query.split(",")
        if vt == "SNP" and bd == "TP": acc[0] += 1
        elif vt == "SNP" and bd == "FP": acc[1] += 1
        bd, vt = truth.split(",")
        if vt == "SNP" and bd == "FN": acc[2] += 1
    return [[k] + v for k, v in sorted(blocks.items())]


def contig_bed(c):
    out = D / "meta" / f"{SAMPLE}_{c}.bed"
    out.parent.mkdir(exist_ok=True)
    if not out.exists():
        with open(out, "w") as w:
            for line in open(TRUTH_BED):
                if line.split("\t", 1)[0] == c:
                    w.write(line)
    return out


def contig_len(c):
    for line in open(str(REF) + ".fai"):
        f = line.split("\t")
        if f[0] == c:
            return int(f[1])


def bed_depth(bam, bed):
    tot = bases = 0
    o = subprocess.run(["samtools", "bedcov", str(bed), str(bam)], env=ENV, capture_output=True, text=True, check=True).stdout
    for line in o.splitlines():
        f = line.split("\t")
        tot += int(f[2]) - int(f[1]); bases += int(f[3])
    return bases / tot


def main():
    c, target, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
    callers = (sys.argv[sys.argv.index("--callers") + 1] if "--callers" in sys.argv else "r2,r11,bcftools").split(",")
    native = BAMDIR / f"{SAMPLE}_{c}_native.bam"
    bed = contig_bed(c)
    OUT.mkdir(parents=True, exist_ok=True)
    meta = D / "meta" / f"{SAMPLE}_{c}_native_cov.json"
    if not meta.exists():
        meta.write_text(json.dumps(dict(bed_mean_depth=bed_depth(native, bed))))
    nat_depth = json.load(open(meta))["bed_mean_depth"]
    cond = f"{SAMPLE}_{c}_native" if target == "native" else f"{SAMPLE}_{c}_{target}x_s{seed}"
    rec_base = dict(sample=SAMPLE, contig=c, region=f"{c} (confident BED, {bed})", target=target, seed=None if target == "native" else seed,
                    native_bed_depth=nat_depth, cond=cond)
    if target == "native":
        bam = native
    else:
        frac = int(target) / nat_depth
        sflag = f"{seed}.{round(frac * 10000):04d}"
        bam = BAMDIR / f"{cond}.bam"
        if not bam.exists():
            subprocess.run(["samtools", "view", "-b", "-s", sflag, "-@", "4", "-o", str(bam), str(native)], env=ENV, check=True)
            subprocess.run(["samtools", "index", str(bam)], env=ENV, check=True)
        rec_base.update(samtools_s=sflag, fraction=frac, bam_sha256=sha256(bam), achieved_bed_depth=bed_depth(bam, bed))
    done = set()
    rf = OUT / "RESULTS.jsonl"
    if rf.exists():
        done = {(json.loads(l)["cond"], json.loads(l)["caller"]) for l in open(rf)}
    for caller in callers:
        if (cond, caller) in done:
            continue
        w = OUT / "runs" / cond / caller
        w.mkdir(parents=True, exist_ok=True)
        rec = dict(rec_base, caller=caller)
        if caller == "bcftools":
            cmd = [str(BCF_SH), str(bam), str(REF), c, str(contig_len(c)), str(w), "32", str(THREADS)]
            rec.update(run(cmd, ENV, w / "call.log")); vcf = w / "bcftools.vcf.gz"
        elif caller == "deepvariant":
            # DeepVariant 1.10.0 (google/deepvariant, manifest sha256:962e5a83...) unpacked in user space, executed with bubblewrap (no root/Docker)
            DV_ROOT = Path.home() / "varion_data/dv/rootfs"
            (w / "tmp").mkdir(exist_ok=True)
            cmd = ["bwrap", "--bind", str(DV_ROOT), "/", "--dev", "/dev", "--proc", "/proc", "--bind", str(w), "/out", "--ro-bind", str(REF.parent), "/ref",
                   "--ro-bind", str(bam.parent), "/bam", "--setenv", "TMPDIR", "/out/tmp", "--chdir", "/out",
                   "/opt/deepvariant/bin/run_deepvariant", "--model_type=WGS", f"--ref=/ref/{REF.name}", f"--reads=/bam/{bam.name}", f"--regions={c}",
                   "--output_vcf=/out/dv.vcf.gz", "--intermediate_results_dir=/out/tmp", f"--num_shards={THREADS}"]
            rec.update(run(cmd, ENV, w / "call.log")); vcf = w / "dv.vcf.gz"
            subprocess.run(["rm", "-rf", str(w / "tmp")])  # intermediate examples (large); logs/VCF kept
        else:
            bin_, model, extra = (R2_BIN, R2_MODEL, ["--r2"]) if caller == "r2" else (R11_BIN, R11_MODEL, [])
            cmd = [bin_, "run", "--bam", str(bam), "--ref", str(REF), "--bed", str(bed), "--out-dir", str(w), "--contig", c,
                   "--contig-len", str(contig_len(c)), "--sample", SAMPLE, "--tag", cond, "--threads", str(THREADS), "--mode", "cascade",
                   "--model", str(model)] + extra
            rec.update(run(cmd, ENV_DNA, w / "call.log"))
            if rec["rc"] == 0:
                v = w / f"ai_cascade_{cond}.vcf"
                vcf = w / "varion.vcf.gz"
                subprocess.run(["bcftools", "view", "-Oz", "-o", str(vcf), str(v)], env=ENV, check=True)
                subprocess.run(["bcftools", "index", "-f", "-t", str(vcf)], env=ENV, check=True)
                rec["vcf_sha256"] = sha256(v)
                for f in w.glob("routed_loci_*.tsv"):
                    f.unlink()
        rec["command"] = " ".join(cmd)
        if rec["rc"] != 0:
            rec["status"] = "CALLER_FAILED"; rec["stderr_tail"] = open(w / "call.log").read()[-400:]
        else:
            hd = w / "happy"; hd.mkdir(exist_ok=True)
            hcmd = [str(HP / "bin/python"), str(HP / "bin/hap.py"), str(TRUTH_VCF), str(vcf), "-f", str(bed), "-r", str(REF), "-o", str(hd / "h"),
                    "-l", c, "--threads", "4"]
            h = run(hcmd, ENV_HAP, hd / "happy.log")
            rec["happy_rc"] = h["rc"]; rec["happy_wall_s"] = h["wall_s"]; rec["happy_command"] = " ".join(hcmd)
            try:
                for r in csv.DictReader(open(hd / "h.summary.csv")):
                    if r["Filter"] == "PASS" and r["Type"] in ("SNP", "INDEL"):
                        # hap.py summary has no QUERY.TP column: QTP = QUERY.TOTAL - QUERY.FP - QUERY.UNK
                        rec[r["Type"] + "_QTP"] = int(r["QUERY.TOTAL"]) - int(r["QUERY.FP"]) - int(r["QUERY.UNK"]); rec[r["Type"] + "_TP"] = int(r["TRUTH.TP"]); rec[r["Type"] + "_FN"] = int(r["TRUTH.FN"]); rec[r["Type"] + "_FP"] = int(r["QUERY.FP"])
                        fl = lambda x: float(x) if x not in ("", None) else None  # empty when a class has no calls (R2 is SNP-only)
                        rec[r["Type"] + "_precision"] = fl(r["METRIC.Precision"]); rec[r["Type"] + "_recall"] = fl(r["METRIC.Recall"])
                        rec[r["Type"] + "_F1"] = fl(r["METRIC.F1_Score"])
                rec["snp_blocks_2Mb"] = block_counts(hd / "h.vcf.gz")
                for f in hd.glob("h.*"):
                    if f.name not in ("h.summary.csv", "h.metrics.json", "h.runinfo.json"):
                        f.unlink()  # annotated VCF/extended/roc are large; block counts + summary kept
                rec["status"] = "OK"
            except Exception as e:  # noqa: BLE001
                rec["status"] = f"HAPPY_FAILED: {e}"
        with open(rf, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(cond, caller, rec.get("status"), rec.get("SNP_F1"), f"wall={rec['wall_s']:.1f}s cpu={rec['cpu_s']:.1f}s rss={rec['peak_rss_kb']/1024:.0f}MB", flush=True)
    if target != "native" and "--keep-bam" not in sys.argv:
        for p in (bam, Path(str(bam) + ".bai")):
            p.unlink(missing_ok=True)  # generated, recoverable from native + recorded sha256/flag


if __name__ == "__main__":
    main()
