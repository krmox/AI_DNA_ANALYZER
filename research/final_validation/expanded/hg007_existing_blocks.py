#!/usr/bin/env python3
"""Analysis-only: per-2Mb SNP TP/FP/FN blocks for the EXISTING HG007 runs, from the retained hap.py annotated VCFs
(same function as run_campaign.block_counts). Also re-derives TP/FP/FN and checks them against the hap.py summary in RESULTS.jsonl."""
import json, subprocess, sys, os
from pathlib import Path
VT = Path.home() / "varion_mamba/envs/vt"
ENV = dict(os.environ, PATH=f"{VT}/bin:{os.environ['PATH']}")
R = Path("/home/kramkost/Projects/Varion_DNA/AI_DNA_ANALYZER/research/final_validation/hg007")

def block_counts(vcf, block=2_000_000):
    text = subprocess.run(["bcftools", "query", "-f", "%POS\t[%BD,%BVT\t]\n", str(vcf)], env=ENV, capture_output=True, text=True, check=True).stdout
    blocks = {}
    for line in text.splitlines():
        pos, truth, query = line.rstrip("\t").split("\t")
        acc = blocks.setdefault((int(pos) - 1) // block, [0, 0, 0])
        bd, vt = query.split(",")
        if vt == "SNP" and bd == "TP": acc[0] += 1
        elif vt == "SNP" and bd == "FP": acc[1] += 1
        bd, vt = truth.split(",")
        if vt == "SNP" and bd == "FN": acc[2] += 1
    return [[k] + v for k, v in sorted(blocks.items())]

bad = 0
with open(Path(__file__).parent / "hg007_existing_blocks.jsonl", "w") as out:
    for l in open(R / "RESULTS.jsonl"):
        r = json.loads(l)
        v = R / "runs" / r["cond"] / r["caller"] / "happy" / "h.vcf.gz"
        b = block_counts(v)
        tp, fp, fn = (sum(x[i] for x in b) for i in (1, 2, 3))
        ok = (fp == r["SNP_FP"] and fn == r["SNP_FN"])
        bad += (not ok)
        out.write(json.dumps(dict(cond=r["cond"], caller=r["caller"], snp_blocks_2Mb=b, blocks_match_summary=ok, sum_tp=tp, sum_fp=fp, sum_fn=fn, summary_tp=r["SNP_TP"], summary_fp=r["SNP_FP"], summary_fn=r["SNP_FN"])) + "\n")
        print(r["cond"], r["caller"], ok, tp, r["SNP_TP"], fp, r["SNP_FP"], fn, r["SNP_FN"], flush=True)
print("MISMATCHES", bad)
