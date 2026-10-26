"""Paired 2-Mb block bootstrap (candidate minus baseline) of dF1, dRecall, dPrecision from hap.py per-record output.
4000 resamples, seed 20261026. Baseline = R1 --truth-free (H2); candidate = H1+H2 (primary), also H1 vs B0 (reference)."""
import subprocess, numpy as np
H="/mnt/archive/AI_DNA_ANALYZER_stage2/happy/out/"
def blocks(tag):
    out=subprocess.run(["bcftools","query","-f","%POS\t[%BD,%BVT\t]\n",H+tag+".vcf.gz"],capture_output=True,text=True,check=True).stdout.splitlines()
    d={}
    for l in out:
        pos,tr,qu=l.rstrip("\t").split("\t"); b=(int(pos)-1)//2_000_000
        a=d.setdefault(b,[0,0,0])
        bd,vt=qu.split(",")
        if vt=="SNP" and bd=="TP": a[0]+=1
        if vt=="SNP" and bd=="FP": a[1]+=1
        bd,vt=tr.split(",")
        if vt=="SNP" and bd=="FN": a[2]+=1
    return d
def metrics(a):
    tp,fp,fn=a.sum(0); return np.array([2*tp/(2*tp+fp+fn), tp/(tp+fn), tp/(tp+fp)])
rng=np.random.default_rng(20261026)
print("contig depth    cand vs base | n_blk | dF1 [95% CI] | dRecall [95% CI] | dPrecision [95% CI]")
for c in ("chr21","chr20"):
    for d in ("f053_s21","full","f010_s21"):
        for cand,base in (("H1H2","H2"),("H1","B0")):
            A=blocks(f"c_{c}_{d}_{base}"); B=blocks(f"c_{c}_{d}_{cand}"); ks=sorted(set(A)|set(B))
            X=np.array([A.get(k,[0,0,0]) for k in ks]); Y=np.array([B.get(k,[0,0,0]) for k in ks])
            d0=metrics(Y)-metrics(X); n=len(ks)
            bs=np.array([metrics(Y[i])-metrics(X[i]) for i in (rng.integers(0,n,n) for _ in range(4000))])
            ci=np.percentile(bs,[2.5,97.5],axis=0)
            print(f"{c} {d:9s} {cand}-{base:3s} | {n:3d} | "+" | ".join(f"{d0[j]:+.5f} [{ci[0,j]:+.5f},{ci[1,j]:+.5f}]" for j in range(3)))
