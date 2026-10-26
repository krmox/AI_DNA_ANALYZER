#!/usr/bin/env python3
"""Run a caller variant on one BAM, convert to VCF.gz, score with hap.py (same image/options as Stage 3, optional stratifications).

usage: run_r2.py OUTROOT SAMPLE CONTIG BAM TAG VARIANT [--strat] [--hap-threads N]
VARIANT: r11 (frozen R1.1 binary + R1.model) | path/to/model.model (R2 binary, --r2)  | modelpath@ood=off ...
Writes OUTROOT/TAG/{calls.json,varion.vcf.gz,happy/*}.  Never touches Stage 3 outputs.
"""
import json, subprocess, sys, time
from pathlib import Path

sys.path.insert(0, "/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/stage3/scripts")
import stage3 as s3  # noqa: E402

R2_BIN = "/mnt/archive/AI_DNA_ANALYZER_r2/build/dnav2"
R11_BIN = str(s3.DNAV2)
R1_MODEL = str(s3.MODEL)


def main():
    root, sample, contig, bam, tag, variant = sys.argv[1:7]
    strat = "--strat" in sys.argv
    hap_threads = int(sys.argv[sys.argv.index("--hap-threads") + 1]) if "--hap-threads" in sys.argv else 4
    extra = [a for a in sys.argv[7:] if a.startswith("--r2-arg=")]
    out = Path(root) / tag
    out.mkdir(parents=True, exist_ok=True)
    d = out / "run"
    d.mkdir(exist_ok=True)
    base = ["run", "--bam", bam, "--ref", str(s3.REF), "--bed", str(s3.contig_bed(sample, contig)), "--out-dir", str(d),
            "--contig", contig, "--contig-len", str(s3.contig_length(contig)), "--sample", sample, "--tag", tag,
            "--threads", str(s3.THREADS), "--mode", "cascade"]
    if variant == "bcftools":
        cmd = []
    elif variant == "r11":
        cmd = [R11_BIN] + base + ["--model", R1_MODEL]
    else:
        cmd = [R2_BIN] + base + ["--model", variant, "--r2"] + [a.split("=", 1)[1] for a in extra]
    if variant == "bcftools":
        a = type("A", (), dict(sample=sample, contig=contig, cond=Path(bam).stem, tag=tag))()
        rec = s3.run_bcftools(a, out, s3.THREADS)
        if rec["rc"] != 0:
            print(tag, "BCFTOOLS FAILED", rec["stderr_tail"]); json.dump(rec, open(out / "calls.json", "w"), indent=1); sys.exit(1)
        for sfx in ("", ".tbi"):
            link = out / f"varion.vcf.gz{sfx}"
            if not link.exists():
                link.symlink_to(f"bcftools/bcftools.vcf.gz{sfx}")
        json.dump(rec, open(out / "calls.json", "w"), indent=1)
    else:
        rec = s3.timed(cmd, d / "time.txt")
        rec["command"] = " ".join(cmd)
        (d / "stdout.json").write_text(rec.pop("stdout"))
        if rec["rc"] != 0:
            print(tag, "CALLER FAILED", rec["stderr_tail"])
            json.dump(rec, open(out / "calls.json", "w"), indent=1)
            sys.exit(1)
        vcf = d / f"ai_cascade_{tag}.vcf"
        subprocess.run(["bcftools", "view", "-Oz", "-o", str(out / "varion.vcf.gz"), str(vcf)], check=True)
        subprocess.run(["bcftools", "index", "-f", "-t", str(out / "varion.vcf.gz")], check=True)
        rec["vcf_sha256"] = s3.sha256(vcf)
        json.dump(rec, open(out / "calls.json", "w"), indent=1)
    truth_dir, tvcf, tbed = s3.truth_files(sample)
    hd = out / "happy"
    hd.mkdir(exist_ok=True)
    cmd = ["podman", "--root", "/dev/shm/podman_root", "--storage-driver", "vfs", "run", "--rm", "--security-opt", "label=disable",
           "-v", f"{truth_dir}:/truth:ro", "-v", f"{s3.REF_DIR}:/ref:ro", "-v", f"{out}:/q:ro", "-v", f"{hd}:/o", "-v", f"{s3.S3/'strat'}:/strat:ro",
           s3.HAPPY_IMAGE, "/usr/local/bin/hap.py", f"/truth/{tvcf}", "/q/varion.vcf.gz", "-f", f"/truth/{tbed}",
           "-r", f"/ref/{s3.REF.name}", "-o", "/o/h", "-l", contig, "--threads", str(hap_threads)]
    if strat:
        cmd += ["--stratification", "/strat/strata.tsv"]
    t0 = time.time()
    p = subprocess.run(cmd, stdout=open(hd / "happy.log", "w"), stderr=subprocess.STDOUT)
    json.dump(dict(rc=p.returncode, wall_s=time.time() - t0, command=" ".join(cmd)), open(hd / "happy_run.json", "w"))
    import csv
    for r in csv.DictReader(open(hd / "h.summary.csv")):
        if r["Type"] == "SNP" and r["Filter"] == "PASS":
            tp, fn, fp = int(r["TRUTH.TP"]), int(r["TRUTH.FN"]), int(r["QUERY.FP"])
            print(f"{tag} {variant.split('/')[-1]} TP {tp} FP {fp} FN {fn} P {1 - fp / max(1, int(r['QUERY.TOTAL'])):.4f} F1 {2 * tp / (2 * tp + fp + fn):.5f} wall {rec['wall_s']:.1f}s cpu {rec['cpu_s']:.0f}s rss {rec['peak_rss_kb']/1024:.0f}MB")


if __name__ == "__main__":
    main()
