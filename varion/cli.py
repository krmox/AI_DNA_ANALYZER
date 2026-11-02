"""varion: one entry point for the SNP engine (R2) and the indel engine (R3.1)."""
import argparse
import json
import os
import subprocess
import sys

from varion import __version__, RELEASE
from varion.tools import UserError

EPILOG = """\
examples:
  varion indel --bam sample.bam --ref GRCh38.fa --out-dir out/            # automatic configuration
  varion indel --bam s.bam --ref r.fa --out-dir out/ --region chr20 --mode A+B   # force A+B on one chromosome
  varion depth --bam sample.bam --region chr20                            # what depth does the BAM have?
  varion doctor                                                           # check the installation

legacy commands (nindels-call, dnav2 run) still work; see docs/COMMAND_MIGRATION.md.
Research software: not for clinical use."""


def build_parser():
    p = argparse.ArgumentParser(prog="varion", description=f"Varion {RELEASE}: SNP (R2) and indel (R3.1) calling for Illumina short reads, GRCh38.",
                                epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"varion {__version__} ({RELEASE})")
    sub = p.add_subparsers(dest="command", metavar="COMMAND")

    i = sub.add_parser("indel", help="call indels; picks the configuration from sequencing depth (R3.1)",
                       description="Call germline indels. The sequencing depth of every contig is measured, then either "
                                   "A (joint genotyping) or A+B (A plus a repeat error model) is used. Writes a VCF and a JSON summary.",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    i.add_argument("--bam", required=True, help="coordinate-sorted, indexed BAM/CRAM")
    i.add_argument("--ref", required=True, help="reference FASTA with .fai (GRCh38)")
    i.add_argument("--out-dir", required=True, help="output directory (created if missing)")
    i.add_argument("--sample", help="sample name for the VCF (default: read group SM, else file name)")
    i.add_argument("--prefix", help="output file prefix (default: sample name)")
    i.add_argument("--region", action="append", metavar="CONTIG",
                   help="contig to call, repeatable (default: chr1-22, X, Y that have reads)")
    i.add_argument("--mode", choices=["auto", "A", "A+B"], default="auto",
                   help="auto: choose from depth; A: joint genotyping only; A+B: add the repeat error model everywhere")
    i.add_argument("--threads", type=int, default=os.cpu_count() or 1, help="worker threads")
    i.add_argument("--min-qual", type=float, default=20.0, help="PASS requires phred QUAL >= this (frozen default)")
    i.add_argument("--thresholds", help="alternative thresholds.json (advanced)")
    i.add_argument("--engine", help="path to nindels-call (default: auto-detect)")
    i.add_argument("--force", action="store_true", help="overwrite existing outputs")
    i.add_argument("--keep-work", action="store_true", help="keep per-contig intermediate files")
    i.add_argument("--dry-run", action="store_true", help="validate inputs and print the plan, run nothing")
    i.add_argument("--quiet", action="store_true", help="no progress messages")
    i.set_defaults(func=cmd_indel)

    d = sub.add_parser("depth", help="report sequencing depth of a BAM (sampled windows)",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    d.add_argument("--bam", required=True)
    d.add_argument("--region", required=True, metavar="CONTIG[:START-END]")
    d.add_argument("--window", type=int, default=10000, help="window size (bp)")
    d.add_argument("--samples", type=int, default=400, help="number of random windows")
    d.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    d.add_argument("--engine", help="path to nindels-call")
    d.set_defaults(func=cmd_depth)

    s = sub.add_parser("snp", help="call SNPs with the R2 engine (see README for scope and validation)",
                       description="Thin wrapper over `dnav2 run --r2`. The SNP engine needs a confident-regions BED "
                                   "that defines which positions are scored; it is not validated outside GIAB-style benchmarks.",
                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    s.add_argument("--bam", required=True)
    s.add_argument("--ref", required=True)
    s.add_argument("--bed", required=True, help="regions to call (e.g. GIAB confident regions)")
    s.add_argument("--out-dir", required=True)
    s.add_argument("--contig", required=True, help="contig, e.g. chr20")
    s.add_argument("--sample", default="SAMPLE")
    s.add_argument("--model", help="R2 model file (default: research/stage4/frozen_R2.model)")
    s.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    s.add_argument("--engine", help="path to dnav2")
    s.add_argument("--dry-run", action="store_true")
    s.set_defaults(func=cmd_snp)

    dr = sub.add_parser("doctor", help="check that the engines and tools are installed")
    dr.add_argument("--engine")
    dr.set_defaults(func=cmd_doctor)
    return p


def cmd_indel(a):
    from varion import indel
    if a.threads < 1:
        raise UserError("--threads must be >= 1")
    return indel.run(a)


def cmd_depth(a):
    from varion.tools import check_bam, find_nindels
    eng = find_nindels(a.engine)
    check_bam(a.bam)
    r = subprocess.run([eng, "depth", "--bam", a.bam, "--region", a.region, "--window", str(a.window),
                        "--samples", str(a.samples), "--threads", str(a.threads)], capture_output=True, text=True)
    if r.returncode:
        raise UserError(r.stderr.strip() or "depth failed")
    print(json.dumps(json.loads(r.stdout), indent=2))
    return 0


def cmd_snp(a):
    from pathlib import Path
    from varion.tools import REPO, check_bam, check_reference, find_dnav2, need_file
    eng = find_dnav2(a.engine)
    check_bam(a.bam); check_reference(a.ref); need_file(a.bed, "BED")
    model = a.model or str(REPO / "research/stage4/frozen_R2.model")
    need_file(model, "R2 model")
    cmd = [eng, "run", "--bam", a.bam, "--ref", a.ref, "--bed", a.bed, "--out-dir", a.out_dir, "--contig", a.contig,
           "--sample", a.sample, "--threads", str(a.threads), "--mode", "cascade", "--r2", "--ood", "auto",
           "--truth-free", "--model", model, "--tag", a.contig]
    from pathlib import Path as P
    contig_len = [int(l.split("\t")[1]) for l in open(a.ref + ".fai") if l.split("\t")[0] == a.contig]
    if not contig_len:
        raise UserError(f"contig {a.contig!r} not in the reference index")
    cmd += ["--contig-len", str(contig_len[0])]
    if a.dry_run:
        print(" ".join(cmd)); return 0
    P(a.out_dir).mkdir(parents=True, exist_ok=True)
    return subprocess.call(cmd)


def cmd_doctor(a):
    import shutil
    from varion.tools import find_dnav2, find_nindels
    ok = True
    for name, finder in (("indel engine (nindels-call)", lambda: find_nindels(a.engine)), ("SNP engine (dnav2)", find_dnav2)):
        try:
            print(f"ok       {name}: {finder()}")
        except UserError as e:
            ok = False if "indel" in name else ok
            print(f"MISSING  {name}: {str(e).splitlines()[0]}")
    for t in ("bgzip", "tabix", "samtools"):
        print(f"{'ok      ' if shutil.which(t) else 'optional'} {t}: {shutil.which(t) or 'not on PATH (only needed for compressed/indexed VCF and indexing)'}")
    return 0 if ok else 1


def main(argv=None):
    p = build_parser()
    a = p.parse_args(argv)
    if not getattr(a, "func", None):
        p.print_help()
        return 2
    try:
        return a.func(a)
    except UserError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        return 130
