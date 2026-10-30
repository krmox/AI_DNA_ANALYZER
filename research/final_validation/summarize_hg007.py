#!/usr/bin/env python3
"""Aggregate hg007/RESULTS.jsonl -> hg007/SUMMARY.md (SNP PASS, hap.py). Point estimates; seed range only (3 seeds on one BAM are not independent samples)."""
import json, statistics as st, collections
R=[json.loads(l) for l in open('hg007/RESULTS.jsonl')]
G=collections.defaultdict(dict)
for r in R: G[(r['contig'],r['target'],r['seed'])][r['caller']]=r
def fmt(x): return f"{x:.5f}"
out=["# HG007 fresh-data results (frozen R2) — hap.py 0.3.15, SNP PASS, confident-region BED, one contig at a time","",
"Per-seed rows are in RESULTS.jsonl. 'min–max' is the range over seeds 101-103 (same native BAM, different downsampling) — NOT a confidence interval.",""]
for c in sorted({k[0] for k in G}):
    out+= [f"## {c}","","| target | n | R2 F1 mean (min–max) | R1.1 F1 | bcftools F1 | R2−R1.1 | R2−bcftools | R2 P / R | R2 TP / FP / FN (sum or single) |","|---|---|---|---|---|---|---|---|---|"]
    for t in ['10','20','30','native']:
        ks=[k for k in G if k[0]==c and k[1]==t]
        rows=[G[k] for k in ks if all(x in G[k] and G[k][x].get('status')=='OK' for x in ('r2','r11','bcftools'))]
        if not rows: continue
        f=lambda cal:[r[cal]['SNP_F1'] for r in rows]
        r2=f('r2'); r11=f('r11'); bc=f('bcftools')
        d1=[a-b for a,b in zip(r2,r11)]; d2=[a-b for a,b in zip(r2,bc)]
        mm=lambda v:f"{fmt(st.mean(v))} ({fmt(min(v))}–{fmt(max(v))})" if len(v)>1 else fmt(v[0])
        p=st.mean(r['r2']['SNP_precision'] for r in rows); rc=st.mean(r['r2']['SNP_recall'] for r in rows)
        tp=" / ".join(str(round(st.mean(r['r2'][k] for r in rows))) for k in ('SNP_TP','SNP_FP','SNP_FN'))
        out.append(f"| {t if t=='native' else t+'x'} | {len(rows)} | {mm(r2)} | {mm(r11)} | {mm(bc)} | {st.mean(d1):+.5f} | {st.mean(d2):+.5f} (min {min(d2):+.5f}) | {p:.5f} / {rc:.5f} | {tp} (mean) |")
    out.append("")
out+=["## Resources (R2, 8 threads, wait4 rusage; machine otherwise idle except chr17 download during early chr19 runs)","","| contig | target | caller | wall s | CPU s | peak RSS MB |","|---|---|---|---|---|---|"]
for c in sorted({k[0] for k in G}):
    for t in ['10','20','30','native']:
        rr=[r for r in R if r['contig']==c and r['target']==t and r.get('status')=='OK']
        for cal in ('r2','r11','bcftools'):
            v=[r for r in rr if r['caller']==cal]
            if v: out.append(f"| {c} | {t} | {cal} | {st.mean(x['wall_s'] for x in v):.1f} | {st.mean(x['cpu_s'] for x in v):.1f} | {max(x['peak_rss_kb'] for x in v)/1024:.0f} |")
bad=[r for r in R if r.get('status')!='OK']; out+=["",f"Failures: {len(bad)}"]+[f"- {r['cond']} {r['caller']}: {r.get('status')}" for r in bad]
open('hg007/SUMMARY.md','w').write("\n".join(out)+"\n"); print("\n".join(out))
