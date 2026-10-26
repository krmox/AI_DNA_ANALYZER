#!/usr/bin/env python3
"""Per-record hap.py decisions (SNP only) -> TSV: pos(1-based) ref alt truth_gt query_gt truth_bd query_bd truth_bk query_bk regions."""
import subprocess, sys
vcf, out = sys.argv[1], sys.argv[2]
fmt = "%POS\t%REF\t%ALT\t[%GT\t%BD\t%BK\t%BVT\t]%INFO/Regions\n"
txt = subprocess.run(["bcftools", "query", "-f", fmt, vcf], capture_output=True, text=True, check=True).stdout
with open(out, "w") as w:
    w.write("pos\tref\talt\ttruth_gt\ttruth_bd\ttruth_bk\ttruth_vt\tquery_gt\tquery_bd\tquery_bk\tquery_vt\tregions\n")
    for line in txt.splitlines():
        f = line.split("\t")
        if "SNP" not in (f[6], f[10]):
            continue
        w.write("\t".join(f[:12]) + "\n")
