#!/usr/bin/env python3
"""BLIND-NEW run: HG006 chr16 with R1.1, R2 (frozen), bcftools baseline, same BAMs, hap.py. Requires R2_LOCK.md to exist (freeze).
usage: run_blind.py COND [COND ...]     COND = native | <depth>x_s<seed>   (e.g. 30x_s101)
Per condition: derive (deterministic samtools -s), run the three callers, score, delete the derived BAM (sha256 logged).
"""
import json, os, subprocess, sys
from pathlib import Path
sys.path.insert(0, "/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/stage3/scripts")
import stage3 as s3
assert Path("/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/stage4/R2_LOCK.md").exists(), "freeze R2 first"
R2 = Path("/mnt/archive/AI_DNA_ANALYZER_r2")
MODEL = "/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/stage4/frozen_R2.model"
HERE = Path(__file__).parent
native = R2 / "blind" / "HG006_chr16_native.bam"
link = s3.S3 / "bam" / "HG006_chr16_native.bam"
for sfx in ("", ".bai"):
    if not Path(str(link) + sfx).exists():
        os.symlink(str(native) + sfx, str(link) + sfx)
covp = s3.S3 / "meta" / "HG006_chr16_native_cov.json"
if not covp.exists():
    s3.cmd_cov(type("A", (), dict(sample="HG006", contig="chr16", cond="HG006_chr16_native"))())
print("native BED depth", json.load(open(covp))["bed_mean_depth"], flush=True)
for cond in sys.argv[1:]:
    name = f"HG006_chr16_{cond}"
    bam = s3.S3 / "bam" / f"{name}.bam"
    derived = cond != "native"
    if derived and not bam.exists():
        t, sd = cond.split("x_s")
        s3.cmd_derive(type("A", (), dict(sample="HG006", contig="chr16", target=int(t), seed=int(sd)))())
    # observed coverage (needed for the table)
    if not (s3.S3 / "meta" / f"{name}_cov.json").exists():
        s3.cmd_cov(type("A", (), dict(sample="HG006", contig="chr16", cond=name))())
    for tag, variant in ((f"b_{cond}_r11", "r11"), (f"b_{cond}_R2", MODEL), (f"b_{cond}_bcftools", "bcftools")):
        if (R2 / "blind_runs" / tag / "happy" / "happy_run.json").exists():
            continue
        subprocess.run([sys.executable, str(HERE / "run_r2.py"), str(R2 / "blind_runs"), "HG006", "chr16", str(bam), tag, variant, "--hap-threads", "3"])
    if derived:
        bam.unlink(); Path(str(bam) + ".bai").unlink()
