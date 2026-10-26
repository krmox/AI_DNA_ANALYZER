"""Paired 2-Mb block bootstrap of dF1 (candidate minus control) from hap.py per-record output. 4000 resamples, seed 20261025."""
import subprocess, sys, numpy as np
H="/mnt/archive/AI_DNA_ANALYZER_depth_robust/happy/out/"
def blocks(tag, odd_only):
    out=subprocess.run(["bcftools","query","-f","%CHROM\t%POS\t[%BD,%BVT\t]\n",H+tag+".vcf.gz"],capture_output=True,text=True,check=True).stdout.splitlines()
    d={}
    for l in out:
        c,pos,tr,qu=l.rstrip("\t").split("\t"); b=(int(pos)-1)//2_000_000
        if odd_only and b%2==0: continue
        a=d.setdefault(b,[0,0,0])
        (bd,vt)=qu.split(",")
        if vt=="SNP" and bd=="TP": a[0]+=1
        if vt=="SNP" and bd=="FP": a[1]+=1
        (bd,vt)=tr.split(",")
        if vt=="SNP" and bd=="FN": a[2]+=1
    return d
def f1(a): tp,fp,fn=a.sum(0); return 2*tp/(2*tp+fp+fn)
rng=np.random.default_rng(20261025)
for ds,odd in (("t1_chr8_30x",False),("t2_chr20_30x",True)):
    ctrl=blocks(f"{ds}_control",odd)
    for cand in ("R1","R2","oldAC"):
        c=blocks(f"{ds}_{cand}",odd); ks=sorted(set(ctrl)|set(c))
        A=np.array([ctrl.get(k,[0,0,0]) for k in ks]); Bm=np.array([c.get(k,[0,0,0]) for k in ks])
        d0=f1(Bm)-f1(A); n=len(ks); ds_=[]
        for _ in range(4000):
            i=rng.integers(0,n,n); ds_.append(f1(Bm[i])-f1(A[i]))
        lo,hi=np.percentile(ds_,[2.5,97.5]); print(f"{ds:14} {cand:6} blocks={n:3} F1 ctrl={f1(A):.5f} cand={f1(Bm):.5f} dF1={d0:+.5f} 95%CI [{lo:+.5f},{hi:+.5f}]")
