"""`varion indel`: depth-adaptive indel calling (R3.1) on top of the nindels engine."""
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from varion import __version__, RELEASE
from varion.adaptive import Thresholds, choose_config
from varion.tools import (UserError, bam_contigs, check_bam, check_reference, default_contigs, fai_contigs,
                          find_nindels)

SCHEMA = "varion.indel.summary/1"


def _run(cmd, log, what):
    log.write("$ " + " ".join(map(str, cmd)) + "\n")
    log.flush()
    t0 = time.monotonic()
    r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
    log.write(r.stdout + r.stderr)
    log.flush()
    if r.returncode != 0:
        tail = "\n".join((r.stderr or r.stdout).strip().splitlines()[-5:])
        raise UserError(f"{what} failed (exit {r.returncode}); see the log file.\n{tail}")
    return r, round(time.monotonic() - t0, 2)


def _bgzip():
    for c in (os.environ.get("VARION_BGZIP"), shutil.which("bgzip")):
        if c:
            return c
    return None


def run(args):
    engine = find_nindels(args.engine)
    bam = check_bam(args.bam)
    ref = check_reference(args.ref)
    thresholds = Thresholds.load(args.thresholds)
    out_dir = Path(args.out_dir)

    bam_ctgs, bam_sample = bam_contigs(engine, bam, ref)
    ref_ctgs = fai_contigs(ref)
    regions = args.region or default_contigs(bam_ctgs, ref_ctgs)
    if not regions:
        raise UserError("no primary chromosomes with reads found in the BAM that are also in the reference; "
                        "pass --region chrN explicitly")
    for c in regions:
        if c not in bam_ctgs:
            raise UserError(f"contig {c!r} is not in the BAM header (BAM contigs include: "
                            f"{', '.join(list(bam_ctgs)[:5])}...)")
        if c not in ref_ctgs:
            raise UserError(f"contig {c!r} is not in the reference FASTA index")
        if bam_ctgs[c][1] == 0:
            raise UserError(f"contig {c!r} has no mapped reads in the BAM")
    sample = args.sample or (bam_sample if bam_sample and bam_sample != "SAMPLE" else Path(args.bam).name.split(".")[0])
    prefix = args.prefix or sample
    vcf_plain = out_dir / f"{prefix}.indels.vcf"
    vcf_gz = out_dir / f"{prefix}.indels.vcf.gz"
    summary_path = out_dir / f"{prefix}.indels.summary.json"
    log_path = out_dir / f"{prefix}.indels.log"
    if not args.force and (vcf_plain.exists() or vcf_gz.exists() or summary_path.exists()):
        raise UserError(f"outputs for prefix {prefix!r} already exist in {out_dir}; use --force to overwrite or --prefix")

    plan = [{"contig": c, "length": bam_ctgs[c][0], "mapped_reads": bam_ctgs[c][1]} for c in regions]
    if args.dry_run:
        print(f"[dry-run] sample={sample} mode={args.mode} engine={engine}\n[dry-run] contigs: " +
              ", ".join(p["contig"] for p in plan) + f"\n[dry-run] outputs: {vcf_gz} (or .vcf), {summary_path}, {log_path}")
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    work = out_dir / f"{prefix}.work"
    work.mkdir(exist_ok=True)
    t_start = time.monotonic()
    contigs_out, bodies, header = [], [], None
    with open(log_path, "w") as log:
        log.write(f"varion {__version__} ({RELEASE}) indel\nargv: {' '.join(sys.argv)}\n")
        for p in plan:
            c = p["contig"]
            _info(args, f"[{c}] measuring sequencing depth")
            r, t_depth = _run([engine, "depth", "--bam", bam, "--region", c, "--threads", args.threads,
                               "--out", work / f"{c}.depth.json"], log, "depth detection")
            prof = json.loads((work / f"{c}.depth.json").read_text())
            dec = choose_config(args.mode, prof["median"] if prof["windows_sampled"] > prof["windows_empty"] else None, thresholds)
            _info(args, f"[{c}] depth median {prof['median']:.1f}x (p10 {prof['p10']:.1f}, p90 {prof['p90']:.1f}) -> "
                        f"configuration {dec.config}: {dec.reason}")
            rec = {"contig": c, "depth_profile": prof, "decision": dec.as_dict(), "seconds": {"depth": t_depth}}
            call_cmd = [engine, "call", "--bam", bam, "--ref", ref, "--region", c, "--threads", args.threads,
                        "--force-reassess", "--genotyper", "multi", "--min-qual", args.min_qual,
                        "--sample", sample, "--out-vcf", work / f"{c}.vcf", "--summary-json", work / f"{c}.call.json"]
            if dec.use_repeat_model:
                _info(args, f"[{c}] estimating the repeat error model from the BAM (no truth used)")
                _, t_str = _run([engine, "estimate-str", "--bam", bam, "--ref", ref, "--region", c,
                                 "--threads", args.threads, "--out", work / f"{c}.str.tsv"], log, "repeat-model estimation")
                rec["seconds"]["estimate_str"] = t_str
                call_cmd += ["--str-model", work / f"{c}.str.tsv"]
                if dec.str_min_depth > 0:
                    call_cmd += ["--str-min-depth", dec.str_min_depth]
            _info(args, f"[{c}] calling indels")
            _, t_call = _run(call_cmd, log, "indel calling")
            rec["seconds"]["call"] = t_call
            rec["call"] = json.loads((work / f"{c}.call.json").read_text())
            lines = (work / f"{c}.vcf").read_text().splitlines()
            h = [l for l in lines if l.startswith("#")]
            header = header or h
            bodies += [l for l in lines if not l.startswith("#")]
            contigs_out.append(rec)
            if not args.keep_work:
                for f in (work / f"{c}.vcf",):
                    f.unlink(missing_ok=True)

    meta = [f"##varion_version={__version__}", f"##varion_release={RELEASE}", f"##varion_mode={args.mode}"]
    for rec in contigs_out:
        d = rec["decision"]
        meta.append(f"##varion_config_{rec['contig']}={d['config']};global_depth={rec['depth_profile']['median']};"
                    f"str_min_depth={d['str_min_depth']};source={d['source']}")
    idx = next(i for i, l in enumerate(header) if l.startswith("#CHROM"))
    final = header[:idx] + meta + header[idx:]
    vcf_plain.write_text("\n".join(final + bodies) + "\n")
    bg = _bgzip()
    out_vcf = vcf_plain
    if bg:
        subprocess.run([bg, "-f", str(vcf_plain)], check=True)
        out_vcf = vcf_gz
        tabix = shutil.which("tabix") or str(Path(bg).with_name("tabix"))
        if Path(tabix).exists():
            subprocess.run([tabix, "-f", "-p", "vcf", str(vcf_gz)], check=False)
    else:
        print("note: bgzip not found; wrote an uncompressed VCF", file=sys.stderr)
    if not args.keep_work:
        shutil.rmtree(work, ignore_errors=True)

    summary = {
        "schema": SCHEMA, "varion_version": __version__, "release": RELEASE,
        "sample": sample, "mode_requested": args.mode, "inputs": {"bam": str(bam), "reference": str(ref)},
        "thresholds": {"global_min_depth": thresholds.global_min_depth, "local_min_depth": thresholds.local_min_depth,
                       "status": thresholds.status, "provenance": thresholds.provenance},
        "engine": engine, "threads": args.threads, "min_qual": args.min_qual,
        "contigs": contigs_out, "outputs": {"vcf": str(out_vcf), "log": str(log_path)},
        "totals": {"loci": sum(c["call"]["loci"] for c in contigs_out), "pass": sum(c["call"]["pass"] for c in contigs_out),
                   "wall_seconds": round(time.monotonic() - t_start, 2)},
        "host": {"platform": sys.platform, "python": ".".join(map(str, sys.version_info[:3]))},
        "scope_note": "Indels only; research use. Accuracy was measured on GIAB samples (see README); not validated for clinical use.",
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    _info(args, f"done: {summary['totals']['pass']} PASS indel loci -> {out_vcf}\n      summary: {summary_path}")
    return 0


def _info(args, msg):
    if not args.quiet:
        print(msg, file=sys.stderr)
