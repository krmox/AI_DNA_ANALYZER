"""Shared helpers for the R3.1 experiment scripts (glue only; no statistics here)."""
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
NINDELS_BIN = Path(os.environ.get("NINDELS_BIN", REPO / "nindels/build/nindels-call"))
VT = Path(os.environ.get("VT_ENV", Path.home() / "varion_mamba/envs/vt"))
HP = Path(os.environ.get("HAPPY_ENV", Path.home() / "varion_mamba/envs/hp"))
SAMTOOLS, BGZIP, TABIX, BCFTOOLS = VT / "bin/samtools", VT / "bin/bgzip", VT / "bin/tabix", VT / "bin/bcftools"

ENV = dict(os.environ, LD_LIBRARY_PATH=f"{VT}/lib:" + os.environ.get("LD_LIBRARY_PATH", ""), VARION_BGZIP=str(VT / "bin/bgzip"))


def registry():
    return json.loads((HERE / "datasets.json").read_text())


def run(cmd, **kw):
    return subprocess.run([str(c) for c in cmd], check=True, env=ENV, **kw)


def timed(cmd, log_path=None, cwd=None):
    """Run cmd, return dict(wall_s, cpu_s, peak_rss_mb) from RUSAGE_CHILDREN deltas of a fresh child."""
    t0 = time.monotonic()
    r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    # peak RSS of the child alone: run through a tiny wrapper process so ru_maxrss is not inherited
    wrapper = ("import resource,subprocess,sys,json;"
               "rc=subprocess.call(sys.argv[1:]);"
               "r=resource.getrusage(resource.RUSAGE_CHILDREN);"
               "print(json.dumps([rc,r.ru_maxrss,r.ru_utime+r.ru_stime]),file=sys.stderr)")
    with open(log_path or os.devnull, "w") as lg:
        p = subprocess.run([sys.executable, "-I", "-c", wrapper] + [str(c) for c in cmd], env=ENV, cwd=cwd,
                           stderr=subprocess.PIPE, stdout=lg, text=True)
    wall = time.monotonic() - t0
    lines = p.stderr.strip().splitlines()
    rc, maxrss_kb, cpu = json.loads(lines[-1])
    if log_path:
        Path(log_path).with_suffix(".stderr").write_text("\n".join(lines[:-1]) + "\n")
    if rc != 0:
        raise RuntimeError(f"command failed rc={rc}: {' '.join(map(str, cmd))}\n" + "\n".join(lines[-10:]))
    return {"wall_s": round(wall, 2), "cpu_s": round(cpu, 2), "peak_rss_mb": round(maxrss_kb / 1024, 1),
            "stderr": "\n".join(lines[:-1])}
