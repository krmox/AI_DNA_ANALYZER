import numpy as np, pandas as pd
from common import *
W="/mnt/archive/AI_DNA_ANALYZER_v2x/work"
_truth=None
def truth():
    global _truth
    if _truth is None: _truth=load_truth_snps()
    return _truth
def bedmask():
    m=np.zeros(64444168,bool)
    for a,b in load_bed(): m[a:b]=True
    return m
BM=bedmask()

BASE=["depth","n","k","nref","k2","vaf","bllr"]
A_FEATS=["alt_bq_mean","ref_bq_mean","bq_diff","alt_bq_hi_frac","alt_lowbq","alt_mq_mean","ref_mq_mean","mq_diff","alt_mq_lt40_frac",
 "alt_fwd_frac","ref_fwd_frac","strand_diff","alt_min_strand","alt_end_mean","ref_end_mean","alt_near_end_frac",
 "alt_nm_mean","ref_nm_mean","alt_clip_frac","ref_clip_frac","alt_uniq_frac","gap_frac","ins_frac","del_next_frac"]
D_CTX=["nc2","nc5","nc10","nc25","nearest","lmr25","lgap10","depth_ratio25","neigh_max_vaf25"]
D_REF=["hp_len","left_is_alt","right_is_alt","dinuc10","is_ts","ref_cg","cpg"]

def load():
    df=pd.read_pickle(f"{W}/feat/features_labeled.pkl")
    df=df[df.inbed].reset_index(drop=True)          # evaluation population = high-confidence BED
    df["y"]=df.truth_gt_ok.astype(int)
    df["block"]=(df.pos0//BLOCK).astype(int)
    df["ctrl"]=np.where(df.routed==1, df.pb48>=10.5, df.bllr>=7.0).astype(bool)
    return df

def truth_counts():
    t=truth(); out={}
    pos=np.array([p for p in t if BM[p]])
    for s in ("dev","holdout"): out[s]=int((split_of(pos)==s).sum())
    return out
TC=truth_counts()

def metrics(df, call, split):
    """offline hap.py-equivalent (validated exactly on the control) for candidate-table decisions in `split`."""
    m=(df.split==split).values
    c=np.asarray(call)&m
    tp=int((c&(df.y==1).values).sum()); fp=int(c.sum())-tp
    fn=TC[split]-tp
    p=tp/(tp+fp) if tp+fp else 0; r=tp/(tp+fn); f1=2*p*r/(p+r) if p+r else 0
    return dict(TP=tp,FP=fp,FN=fn,precision=p,recall=r,F1=f1)

def extra_calls():
    """control calls that are not in the candidate table (k<3): carried unchanged by every variant."""
    vc=[l.split("\t") for l in open(f"{W}/runs/V2_CONTROL_t8_warm/ai_cascade_chr20_300x.vcf") if not l.startswith("#")]
    df=pd.read_pickle(f"{W}/feat/features_labeled.pkl")
    S=set(df.pos0.values); t=truth(); ex=[]
    for x in vc:
        p=int(x[1])-1
        if p not in S and BM[p]:
            ok=(p in t and x[4] in t[p][1] and t[p][2]==x[9].split(":")[0])
            ex.append((p,bool(ok)))
    return ex
EX=extra_calls()
def metrics(df, call, split):
    m=(df.split==split).values
    c=np.asarray(call)&m
    tp=int((c&(df.y==1).values).sum()); fp=int(c.sum())-tp
    for p,ok in EX:
        if split_of([p])[0]==split: tp+=ok; fp+=(not ok)
    fn=TC[split]-tp
    pr=tp/(tp+fp) if tp+fp else 0; r=tp/(tp+fn); f1=2*pr*r/(pr+r) if pr+r else 0
    return dict(TP=tp,FP=fp,FN=fn,precision=pr,recall=r,F1=f1)
