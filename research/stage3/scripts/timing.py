#!/usr/bin/env python3
"""Phase 3 timing campaign: interleaved repetitions, sequential, nothing else started by this script.

usage: timing.py COND [COND ...] [--reps 5] [--bcftools1]
COND = existing BAM name in stage3/bam without .bam, e.g. HG002_chr22_30x_s101. Appends to timing/timing.jsonl and
removes each repetition's VCF after hashing it (outputs must be byte-identical across repetitions).
"""
import argparse, json, shutil, subprocess, sys, time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import stage3 as s3  # noqa: E402

S3 = s3.S3


def one(sample, contig, cond, caller, rep):
    tag = f"timing_{cond}_{caller}_rep{rep}"
    out = S3 / "runs_timing" / tag
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    a = type("A", (), dict(sample=sample, contig=contig, cond=cond, tag=tag))()
    if caller == "varion":
        rec = s3.run_varion(a, out)
    elif caller == "bcftools8":
        rec = s3.run_bcftools(a, out, s3.THREADS)
    else:
        rec = s3.run_bcftools(a, out, 1)
    rec.update(cond=cond, caller=caller, rep=rep, loadavg_after=float(open("/proc/loadavg").read().split()[0]),
               bam_bytes=(S3 / "bam" / f"{cond}.bam").stat().st_size)
    shutil.rmtree(out, ignore_errors=True)
    (S3 / "timing").mkdir(exist_ok=True)
    with open(S3 / "timing" / "timing.jsonl", "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(cond, caller, rep, f"wall {rec.get('wall_s')} cpu {rec.get('cpu_s'):.1f} rss {rec.get('peak_rss_kb')} rc {rec['rc']} load {rec['loadavg_before']}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("conds", nargs="+")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--bcftools1", action="store_true")
    args = ap.parse_args()
    for cond in args.conds:
        sample, contig = cond.split("_")[:2]
        subprocess.run(["cat", str(S3 / "bam" / f"{cond}.bam")], stdout=subprocess.DEVNULL, check=True)  # warm the page cache
        callers = ["varion", "bcftools8"] + (["bcftools1"] if args.bcftools1 else [])
        for rep in range(1, args.reps + 1):
            for caller in callers:
                one(sample, contig, cond, caller, rep)


if __name__ == "__main__":
    main()
