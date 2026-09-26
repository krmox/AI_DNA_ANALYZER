"""Mechanism D ablation: does local context add information beyond the per-locus model?  (A+C base; add window families)
All fits on DEVELOPMENT only (block-CV); holdout applied frozen."""
import sys; sys.path.insert(0,'.')
from lib import *
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import GroupKFold
LOGF=("depth","n","k","nref","k2","alt_lowbq")
def prep(X):
    X=X.copy()
    for c in LOGF:
        if c in X: X[c]=np.log1p(X[c])
    return X
df=load(); df["routed_f"]=df.routed.astype(float); isdev=(df.split=="dev").values
def run(feats,C=0.1):
    X=prep(df[feats]).values; y=df.y.values
    Xd,yd,gd=X[isdev],y[isdev],df.block.values[isdev]; dd=df[isdev]
    mk=lambda: make_pipeline(StandardScaler(),LogisticRegression(C=C,max_iter=3000))
    p=np.zeros(len(yd))
    for tr,te in GroupKFold(8).split(Xd,yd,gd): p[te]=mk().fit(Xd[tr],yd[tr]).predict_proba(Xd[te])[:,1]
    best=max(((metrics(dd,p>=t,"dev")["F1"],t) for t in np.linspace(0.05,0.95,91)))
    m=mk().fit(Xd,yd); ph=m.predict_proba(X)[:,1]
    h=metrics(df,ph>=best[1],"holdout"); d=metrics(dd,p>=best[1],"dev")
    return d,h,best[1]
base=BASE+["routed_f"]+A_FEATS
fams={"A+C (no context)":base,
 "+D1 (±2bp: nc2)":base+["nc2"],
 "+D2 (±5bp: nc2,nc5)":base+["nc2","nc5"],
 "+D3 (±10bp: nc2,nc5,nc10)":base+["nc2","nc5","nc10"],
 "+D4 (±25bp: nc2..nc25,nearest)":base+["nc2","nc5","nc10","nc25","nearest"],
 "+D-other (lmr25,lgap10,depth_ratio25,neigh_max_vaf25)":base+["lmr25","lgap10","depth_ratio25","neigh_max_vaf25"],
 "+D-ref (hp_len,slippage,dinuc,ts,CpG)":base+D_REF,
 "+D all":base+D_CTX+D_REF}
rows=[]
for k,f in fams.items():
    d,h,t=run(f)
    rows.append(dict(model=k,nfeat=len(f),thr=round(t,2),dev_TP=d["TP"],dev_FP=d["FP"],dev_FN=d["FN"],dev_F1=round(d["F1"],5),
                     hold_TP=h["TP"],hold_FP=h["FP"],hold_FN=h["FN"],hold_F1=round(h["F1"],5)))
    print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(f"{W}/context_window_ablation.csv",index=False)
