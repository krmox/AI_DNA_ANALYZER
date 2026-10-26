#!/usr/bin/env python3
"""Stage 3 pipeline: fetch -> depth ladder -> callers -> hap.py. Idempotent; re-run to resume.

usage: run_all.py UNIT [UNIT ...]       UNIT in {HG002_chr22, HG005_chr22, anchor_HG004_chr21, ...}
Failures are logged and skipped (never silently dropped: see logs/run_all.log and FAILED markers).
"""
import json, os, shutil, subprocess, sys, time, traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import stage3 as s3  # noqa: E402

S3 = s3.S3
HERE = Path(__file__).parent
LOG = S3 / "logs" / "run_all.log"
SEEDS = (101, 102, 103)
LADDER = [(30, SEEDS), (20, SEEDS), (15, SEEDS), (10, SEEDS), (50, (101,))]
HG005_BASE_SEED_FRACTION = "100.3600"  # base subsample of the 300x library; deeper than 100x so 100x can be derived
HAPPY_WORKERS = 3
HAPPY_THREADS = 2
MIN_FREE_GB = 7


def log(msg):
    line = f"{time.strftime('%F %T')} {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def free_gb():
    st = os.statvfs(S3)
    return st.f_bavail * st.f_frsize / 1e9


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def ensure_native(sample, contig):
    native = S3 / "bam" / f"{sample}_{contig}_native.bam"
    if native.exists():
        return native
    while free_gb() < MIN_FREE_GB + 3:
        log(f"waiting for disk space ({free_gb():.1f} GB free)")
        time.sleep(60)
    if sample != "HG005":
        run([str(HERE / "fetch_contig.sh"), sample, contig])
    else:
        url = ("https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/"
               "NHGRI_Illumina300X_Chinesetrio_novoalign_bams/HG005.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam")
        bai = Path("/home/mark/Documents/Projects/AI_DNA_ANALYZER/HG005.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam.bai")
        cmd = ["samtools", "view", "-X", "-b", "-s", HG005_BASE_SEED_FRACTION, "-o", str(native) + ".part", url, str(bai), contig]
        run(cmd)
        os.rename(str(native) + ".part", native)
        run(["samtools", "index", str(native)])
        with open(S3 / "logs" / "fetch.log", "a") as f:
            f.write(f"{sample} {contig} {url} base-subsample -s {HG005_BASE_SEED_FRACTION} {native.stat().st_size} {s3.sha256(native)}\n")
    return native


def measure_cov(sample, contig, cond):
    path = S3 / "meta" / f"{cond}_cov.json"
    if not path.exists():
        s3.cmd_cov(type("A", (), dict(sample=sample, contig=contig, cond=cond))())
    return json.load(open(path))


def call_condition(sample, contig, cond):
    calls = S3 / "runs" / cond / "calls.json"
    if calls.exists() and all(k in json.load(open(calls)) for k in ("varion", "bcftools8")):
        return
    a = type("A", (), dict(sample=sample, contig=contig, cond=cond, tag=cond, callers="varion,bcftools8", hash_bam=True))()
    s3.cmd_call(a)


def happy_job(sample, contig, cond, caller):
    done = S3 / "runs" / cond / f"happy_{caller}" / "happy_run.json"
    if done.exists() and json.load(open(done))["rc"] == 0:
        return
    a = type("A", (), dict(sample=sample, contig=contig, tag=cond, caller=caller, threads=HAPPY_THREADS))()
    s3.cmd_happy(a)


def process_unit(sample, contig, ladder, pool, anchors=None):
    futures = []
    if anchors is None:
        native = ensure_native(sample, contig)
        cov = measure_cov(sample, contig, f"{sample}_{contig}_native")
        log(f"{sample} {contig} native BED depth {cov['bed_mean_depth']:.2f}x")
        conds = [(f"{sample}_{contig}_native", None, None)]
        for target, seeds in ladder:
            conds += [(f"{sample}_{contig}_{target}x_s{seed}", target, seed) for seed in seeds]
    else:
        conds = []
        for label, src in anchors.items():
            link = S3 / "bam" / f"{sample}_{contig}_{label}.bam"
            for suffix in ("", ".bai"):
                if not Path(str(link) + suffix).exists():
                    os.symlink(str(src) + suffix, str(link) + suffix)
            conds.append((f"{sample}_{contig}_{label}", None, None))
    cov_pool = ThreadPoolExecutor(max_workers=1)
    for cond, target, seed in conds:
        try:
            derived = target is not None
            if derived and not (S3 / "bam" / f"{cond}.bam").exists():
                s3.cmd_derive(type("A", (), dict(sample=sample, contig=contig, target=target, seed=seed))())
            cov_future = cov_pool.submit(measure_cov, sample, contig, cond)
            call_condition(sample, contig, cond)
            cov_future.result()
            for caller in ("varion", "bcftools"):
                futures.append(pool.submit(happy_job, sample, contig, cond, caller))
            if derived:
                (S3 / "bam" / f"{cond}.bam").unlink()
                (S3 / "bam" / f"{cond}.bam.bai").unlink()
            log(f"called {cond}")
        except Exception:
            log(f"FAILED {cond}\n{traceback.format_exc()}")
            (S3 / "runs").mkdir(exist_ok=True)
            (S3 / "runs" / f"FAILED_{cond}").write_text(traceback.format_exc())
    return futures


UNITS = {}
for smp in ("HG002", "HG003", "HG004"):
    for ctg in ("chr22", "chr18"):
        UNITS[f"{smp}_{ctg}"] = (smp, ctg, LADDER)
UNITS["HG005_chr22"] = ("HG005", "chr22", [(100, SEEDS), (50, (101,))] + LADDER[:4])
UNITS["HG005_chr18"] = ("HG005", "chr18", [(100, (101,)), (30, SEEDS), (10, (101,))])
STAGE2 = Path("/mnt/archive/AI_DNA_ANALYZER_stage2/bam")
ANCHORS = {c: {"native": STAGE2 / f"HG004_{c}_full.bam", "30x_s21": STAGE2 / f"HG004_{c}_f053_s21.bam", "6x_s21": STAGE2 / f"HG004_{c}_f010_s21.bam"}
           for c in ("chr20", "chr21")}


def main():
    pool = ThreadPoolExecutor(max_workers=HAPPY_WORKERS)
    futures = []
    for name in sys.argv[1:]:
        if name.startswith("anchor_"):
            ctg = name.split("_")[-1]
            futures += process_unit("HG004", ctg, None, pool, anchors=ANCHORS[ctg])
        else:
            smp, ctg, ladder = UNITS[name]
            futures += process_unit(smp, ctg, ladder, pool)
            if os.environ.get("KEEP_NATIVE") != "1":
                for suffix in ("", ".bai"):
                    p = S3 / "bam" / f"{smp}_{ctg}_native.bam{suffix}"
                    if p.exists():
                        p.unlink()
                log(f"removed native BAM {smp} {ctg} (sha256 in logs/fetch.log; re-derivable)")
    for f in futures:
        f.result()
    log(f"done: {' '.join(sys.argv[1:])}")


if __name__ == "__main__":
    main()
