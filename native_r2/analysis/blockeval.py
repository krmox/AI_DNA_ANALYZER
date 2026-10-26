import numpy as np, pandas as pd
from lib import *
def vcf_calls(path):
    vc=[l.split("\t") for l in open(path) if not l.startswith("#")]
    return pd.DataFrame(dict(pos0=[int(x[1])-1 for x in vc],alt=[BASES.index(x[4]) for x in vc],gt=[x[9].split(":")[0] for x in vc]))
def per_block(path):
    """TP/FP/FN per 2-Mb block from a VCF using the hap.py-equivalent locus evaluator (validated == hap.py on control)."""
    t=truth(); calls=vcf_calls(path)
    nb=64444167//BLOCK+1
    TP=np.zeros(nb,int); FP=np.zeros(nb,int); FN=np.zeros(nb,int)
    tp=set()
    for p,a,g in zip(calls["pos0"],calls["alt"],calls["gt"]):
        if not BM[p]: continue
        b=p//BLOCK; tr=t.get(int(p))
        if tr is not None and BASES[a] in tr[1] and tr[2]==g: TP[b]+=1; tp.add(int(p))
        else: FP[b]+=1
    for p in t:
        if BM[p]: FN[p//BLOCK]+=1
    FN-=TP
    return pd.DataFrame(dict(block=np.arange(nb),TP=TP,FP=FP,FN=FN))
def f1(tp,fp,fn):
    p=tp/(tp+fp) if tp+fp else 0.0; r=tp/(tp+fn) if tp+fn else 0.0
    return p,r,(2*p*r/(p+r) if p+r else 0.0)
def boot_delta(bm,bc,split,nboot=4000,seed=7):
    """block bootstrap of (model - control): resample 2-Mb blocks of `split` with replacement."""
    rng=np.random.default_rng(seed)
    sel=[b for b in bm.block if (split=="full" or (b%2==0)==(split=="dev"))]
    m=bm.set_index("block").loc[sel]; c=bc.set_index("block").loc[sel]
    n=len(sel); out=[]
    for _ in range(nboot):
        ix=rng.integers(0,n,n)
        a=m.iloc[ix][["TP","FP","FN"]].sum(); b=c.iloc[ix][["TP","FP","FN"]].sum()
        pa,ra,fa=f1(*a); pb,rb,fb=f1(*b)
        out.append((fa-fb,pa-pb,ra-rb,a.FP-b.FP,a.FN-b.FN))
    o=np.array(out); lo,hi=np.percentile(o,[2.5,97.5],axis=0)
    a=m[["TP","FP","FN"]].sum(); b=c[["TP","FP","FN"]].sum()
    pa,ra,fa=f1(*a); pb,rb,fb=f1(*b)
    pt=(fa-fb,pa-pb,ra-rb,int(a.FP-b.FP),int(a.FN-b.FN))
    return pt,lo,hi
