#!/usr/bin/env python3
"""Prepare per-dataset inputs: contig reference, confident BED, seeded depth ladder.

usage: prep.py DATASET DEPTH [DEPTH ...]     (DEPTH = target mean depth, or 'native')
Downsampling: samtools view -s 21.FRAC (pairs kept together), FRAC = target / native median depth.
"""
import json
import sys
from pathlib import Path

from common import NINDELS_BIN, SAMTOOLS, registry, run, VT
import subprocess


def native_depth(bam, chrom):
    out = subprocess.run([str(NINDELS_BIN), "depth", "--bam", bam, "--region", chrom, "--threads", "8"],
                         check=True, capture_output=True, text=True, env={"LD_LIBRARY_PATH": f"{VT}/lib"}).stdout
    return json.loads(out)["median"]


def main():
    reg = registry()
    name, depths = sys.argv[1], sys.argv[2:]
    d = reg["datasets"][name]
    root = Path(reg["work_root"]) / name
    root.mkdir(parents=True, exist_ok=True)
    ref = root / f"ref_{d['chrom']}.fa"
    if not ref.with_suffix(".fa.fai").exists():
        with open(ref, "w") as fh:
            run([SAMTOOLS, "faidx", reg["ref_fasta"], d["chrom"]], stdout=fh)
        run([SAMTOOLS, "faidx", ref])
    conf = root / f"conf_{d['chrom']}.bed"
    conf.write_text("".join(l for l in open(d["bed"]) if l.split("\t")[0] == d["chrom"]))
    (root / "strata.tsv").write_text("".join(f"{k}\t{v}\n" for k, v in reg["strata"].items()))
    nat = native_depth(d["bam"], d["chrom"])
    meta = {"native_median_depth": nat, "ladder": {}}
    for dep in depths:
        if dep == "native":
            link = root / "native.bam"
            for suf in ("", ".bai"):
                if not Path(str(link) + suf).exists():
                    Path(str(link) + suf).symlink_to(d["bam"] + suf)
            meta["ladder"]["native"] = {"bam": str(link), "target": nat, "fraction": 1.0}
            continue
        t = float(dep)
        if t >= nat * 0.97:
            print(f"skip {dep}x: native depth is {nat:.1f}x")
            continue
        out = root / f"d{int(t)}.bam"
        frac = f"{t / nat:.4f}".split(".")[1]
        if not Path(str(out) + ".bai").exists():
            run([SAMTOOLS, "view", "-@4", "-b", "-s", f"21.{frac}", "-o", out, d["bam"], d["chrom"]])
            run([SAMTOOLS, "index", out])
        meta["ladder"][f"d{int(t)}"] = {"bam": str(out), "target": t, "fraction": float("0." + frac)}
    old = {}
    if (root / "ladder.json").exists():
        old = json.loads((root / "ladder.json").read_text())
    old.update({"native_median_depth": nat}); old.setdefault("ladder", {}).update(meta["ladder"])
    (root / "ladder.json").write_text(json.dumps(old, indent=1))
    print(json.dumps(old, indent=1))


if __name__ == "__main__":
    main()
