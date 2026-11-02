#!/usr/bin/env python3
"""Run one caller configuration on one (dataset, depth) BAM, time it, score it with hap.py.

usage: run_one.py PHASE DATASET DEPTHKEY CONFIG [--reps N] [--extra "nindels args"]
  PHASE     dev | test | seen        (output tree /home/kramkost/r31/runs/PHASE/...; checked against the registry role)
  DEPTHKEY  d10 d15 ... native       (from <work_root>/<DATASET>/ladder.json)
  CONFIG    baseline | A | B | AB | AB_gated:<T> | R31 (via `varion indel --mode auto`)
B-model estimation is timed and cached per BAM; its time is added to the config's runtime.
"""
import argparse
import json
import shutil
import statistics
import sys
from pathlib import Path

from common import HP, NINDELS_BIN, REPO, TABIX, BGZIP, BCFTOOLS, ENV, registry, run, timed
import subprocess


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("phase"); ap.add_argument("dataset"); ap.add_argument("depth"); ap.add_argument("config")
    ap.add_argument("--reps", type=int, default=1)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    if a.phase == "dev" and a.config == "baseline" and a.depth not in ("d30", "native"):
        return  # the dev ladder needs A and A+B only; baseline is kept at 30x and native for reference
    reg = registry()
    d = reg["datasets"][a.dataset]
    allowed = {"dev": {"dev"}, "test": {"test"}, "seen": {"seen"}}[a.phase]
    if d["role"] not in allowed:
        sys.exit(f"refusing: dataset {a.dataset} has role '{d['role']}', not allowed in phase '{a.phase}'")
    root = Path(reg["work_root"]) / a.dataset
    bam = json.loads((root / "ladder.json").read_text())["ladder"][a.depth]["bam"]
    chrom = d["chrom"]
    ref = root / f"ref_{chrom}.fa"
    out = Path(reg["work_root"]) / "runs" / a.phase / a.dataset / a.depth / a.config.replace(":", "_")
    if (out / "metrics.json").exists() and not a.force:
        print("exists", out); return
    out.mkdir(parents=True, exist_ok=True)

    cfg = a.config
    str_tsv = Path(reg["work_root"]) / "runs" / a.phase / a.dataset / a.depth / "str.tsv"
    str_time = {"wall_s": 0.0, "cpu_s": 0.0, "peak_rss_mb": 0.0}
    def need_str():
        nonlocal str_time
        tj = str_tsv.with_suffix(".time.json")
        if not str_tsv.exists():
            t = timed([NINDELS_BIN, "estimate-str", "--bam", bam, "--ref", ref, "--region", chrom,
                       "--threads", a.threads, "--out", str_tsv], log_path=str(str_tsv.with_suffix(".log")))
            tj.write_text(json.dumps({k: t[k] for k in ("wall_s", "cpu_s", "peak_rss_mb")}))
        str_time = json.loads(tj.read_text())

    base = [NINDELS_BIN, "call", "--bam", bam, "--ref", ref, "--region", chrom, "--threads", a.threads,
            "--force-reassess", "--out-vcf", out / "calls.vcf", "--summary-json", out / "summary.json"]
    if cfg == "baseline": extra = []
    elif cfg == "A": extra = ["--genotyper", "multi"]
    elif cfg == "B": need_str(); extra = ["--str-model", str_tsv]
    elif cfg == "AB": need_str(); extra = ["--genotyper", "multi", "--str-model", str_tsv]
    elif cfg.startswith("AB_gated:"):
        need_str(); extra = ["--genotyper", "multi", "--str-model", str_tsv, "--str-min-depth", cfg.split(":")[1]]
    elif cfg == "R31":
        extra = None   # runs through the product CLI below
    else:
        sys.exit(f"unknown config {cfg}")

    times = []
    if cfg == "R31":
        # the product path: depth detection + decision + (estimate-str) + call, exactly as a user runs it
        for r in range(a.reps):
            shutil.rmtree(out / "varion", ignore_errors=True)
            cmd = [sys.executable, "-m", "varion", "indel", "--bam", bam, "--ref", ref, "--region", chrom,
                   "--threads", a.threads, "--out-dir", out / "varion", "--prefix", "r31", "--force", "--quiet",
                   "--engine", NINDELS_BIN]
            times.append(timed(cmd, log_path=str(out / f"run{r + 1}.log"), cwd=REPO))
        shutil.copy(out / "varion/r31.indels.vcf.gz", out / "calls.vcf.gz")
        shutil.copy(out / "varion/r31.indels.summary.json", out / "summary.json")
    else:
        for r in range(a.reps):
            t = timed(base + extra, log_path=str(out / f"run{r + 1}.log"))
            times.append(t)
        run([BGZIP, "-f", out / "calls.vcf"])
    run([TABIX, "-f", "-p", "vcf", out / "calls.vcf.gz"])
    # score exactly what a user receives: PASS records only (Amendment 1)
    run([BCFTOOLS, "view", "-f", "PASS", "-Oz", "-o", out / "calls.pass.vcf.gz", out / "calls.vcf.gz"])
    run([TABIX, "-f", "-p", "vcf", out / "calls.pass.vcf.gz"])

    # hap.py (xcmp), PASS records, stratified. The truth is cut to the contig once (same records, much faster to load).
    truth = root / f"truth_{chrom}.vcf.gz"
    if not truth.exists():
        tmp = root / f"truth_{chrom}.tmp.vcf.gz"
        run([BCFTOOLS, "view", "-r", chrom, "-Oz", "-o", tmp, d["truth"]])
        run([TABIX, "-f", "-p", "vcf", tmp])
        Path(str(tmp) + ".tbi").rename(str(truth) + ".tbi")
        tmp.rename(truth)
    env = dict(ENV, PATH=f"{HP}/bin:" + ENV["PATH"], HGREF=str(ref))
    subprocess.run([str(HP / "bin/python"), str(HP / "bin/hap.py"), str(truth), str(out / "calls.pass.vcf.gz"), "-r", str(ref),
                    "-f", str(root / f"conf_{chrom}.bed"), "-l", chrom, "--stratification", str(root / "strata.tsv"),
                    "--threads", str(a.threads), "-o", str(out / "happy")],
                   check=True, env=env, stdout=open(out / "happy.log", "w"), stderr=subprocess.STDOUT)

    m = {"dataset": a.dataset, "depth": a.depth, "config": cfg, "phase": a.phase, "role": d["role"], "bam": bam,
         "reps": a.reps,
         "wall_s": statistics.median(t["wall_s"] for t in times) + str_time["wall_s"],
         "cpu_s": statistics.median(t["cpu_s"] for t in times) + str_time["cpu_s"],
         "peak_rss_mb": max([t["peak_rss_mb"] for t in times] + [str_time["peak_rss_mb"]]),
         "str_estimation": str_time, "call_wall_s_all_reps": [t["wall_s"] for t in times]}
    (out / "metrics.json").write_text(json.dumps(m, indent=1))
    print(json.dumps(m))


if __name__ == "__main__":
    main()
