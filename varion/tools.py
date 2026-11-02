"""Locating native binaries and validating inputs. Errors raised here are UserError: the CLI prints them
without a traceback and exits with status 2."""
import os
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


class UserError(Exception):
    pass


def find_nindels(explicit=None):
    cands = [explicit, os.environ.get("VARION_NINDELS"), REPO / "nindels/build/nindels-call", shutil.which("nindels-call")]
    for c in cands:
        if c and Path(c).is_file() and os.access(c, os.X_OK):
            return str(c)
    raise UserError("indel engine 'nindels-call' not found. Build it with:\n"
                    "    cmake -S nindels -B nindels/build [-DHTSLIB_ROOT=/path/to/htslib] && cmake --build nindels/build -j\n"
                    "or point to it with --engine / $VARION_NINDELS.")


def find_dnav2(explicit=None):
    cands = [explicit, os.environ.get("VARION_DNAV2"), REPO / "native_r2/build/dnav2", REPO / "native_v2x/build/dnav2",
             shutil.which("dnav2")]
    for c in cands:
        if c and Path(c).is_file() and os.access(c, os.X_OK):
            return str(c)
    raise UserError("SNP engine 'dnav2' not found. Build it with: make -C native_r2 engine\n"
                    "or point to it with --engine / $VARION_DNAV2.")


def need_file(path, what):
    p = Path(path)
    if not p.is_file():
        raise UserError(f"{what} not found: {path}")
    return p


def check_bam(path):
    p = need_file(path, "BAM/CRAM")
    idx = [Path(str(p) + s) for s in (".bai", ".csi", ".crai")] + [p.with_suffix(".bai")]
    if not any(i.is_file() for i in idx):
        raise UserError(f"{path} has no index (.bai/.csi/.crai). Create one with: samtools index {path}")
    return p


def check_reference(path):
    p = need_file(path, "reference FASTA")
    if not Path(str(p) + ".fai").is_file():
        raise UserError(f"{path} has no .fai index. Create one with: samtools faidx {path}")
    return p


def fai_contigs(ref):
    out = {}
    with open(str(ref) + ".fai") as fh:
        for line in fh:
            f = line.split("\t")
            out[f[0]] = int(f[1])
    return out


def bam_contigs(engine, bam, ref=None):
    """({contig: (length, mapped_reads)}, sample) from the BAM index via the engine."""
    cmd = [engine, "contigs", "--bam", str(bam)]
    if ref:
        cmd += ["--ref", str(ref)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise UserError(f"cannot read {bam}: {r.stderr.strip() or 'engine failed'}")
    sample, out = None, {}
    for line in r.stdout.splitlines():
        f = line.split("\t")
        if line.startswith("#sample"):
            sample = f[1]
        else:
            out[f[0]] = (int(f[1]), int(f[2]))
    return out, sample


def default_contigs(bam_ctgs, ref_ctgs):
    """Primary chromosomes (chr1-22,X,Y or 1-22,X,Y) that have mapped reads and are in the reference."""
    names = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]] + [str(i) for i in list(range(1, 23)) + ["X", "Y"]]
    return [n for n in names if n in bam_ctgs and n in ref_ctgs and bam_ctgs[n][1] > 0]
