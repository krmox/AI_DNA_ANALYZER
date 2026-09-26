"""Offline ablation grid over mechanisms A,B,C,D.  Every fitted quantity is fitted on DEVELOPMENT loci only
(block-CV inside dev for C/threshold selection); HOLDOUT is only *applied* (frozen) and reported."""
import sys, itertools, json, hashlib; sys.path.insert(0,'.')
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
df=load(); df["routed_f"]=df.routed.astype(float)
dev=df[df.split=="dev"]
PBV=["pb24","pb48","pb64","pb96","pb128","pb0","pbblk_med","pbblk_min","pbblk_max"]

def routed_call(pbcol): return np.where(df.routed==1, df[pbcol]>=10.5, df.bllr>=7.0)
def dev_f1(call): return metrics(df,call,"dev")["F1"]

def choose_B():
    res={c:dev_f1(routed_call(c)) for c in PBV}
    base=res["pb48"]; best=max(res,key=res.get)
    return best,res

def veto_search(feats, base_call, min_gain=0.0):
    """best single hard veto  call &= ~(feat >/< thr)  on dev (one stump: 2 fitted numbers)."""
    d=dev; bc=base_call[(df.split=="dev").values]; best=(dev_f1(base_call),None)
    for f in feats:
        v=d[f].values
        for q in np.unique(np.quantile(v[np.isfinite(v)],np.linspace(0.02,0.98,49))):
            for dirn in (">=","<="):
                veto=(df[f].values>=q) if dirn==">=" else (df[f].values<=q)
                f1=dev_f1(base_call&~veto)
                if f1>best[0]+min_gain: best=(f1,(f,dirn,float(q)))
    return best[1]
def apply_veto(call,rule):
    if rule is None: return call
    f,dirn,q=rule; v=(df[f].values>=q) if dirn==">=" else (df[f].values<=q)
    return call&~v

def fit_lr(feats):
    X=prep(df[feats]).values; y=df.y.values; isdev=(df.split=="dev").values
    Xd,yd,gd=X[isdev],y[isdev],df.block.values[isdev]
    best=None
    for C in (0.03,0.1,0.3,1.0,3.0):
        mk=lambda: make_pipeline(StandardScaler(),LogisticRegression(C=C,max_iter=3000))
        p=np.zeros(len(yd))
        for tr,te in GroupKFold(8).split(Xd,yd,gd):
            p[te]=mk().fit(Xd[tr],yd[tr]).predict_proba(Xd[te])[:,1]
        dd=df[isdev]
        for t in np.linspace(0.05,0.95,91):
            f1=metrics(dd,p>=t,"dev")["F1"]
            if best is None or f1>best[0]: best=(f1,C,t,p.copy())
    f1,C,t,poof=best
    m=make_pipeline(StandardScaler(),LogisticRegression(C=C,max_iter=3000)).fit(Xd,yd)
    pall=m.predict_proba(X)[:,1]
    pdev_oof=pall.copy(); pdev_oof[isdev]=poof     # dev uses out-of-fold probabilities
    return dict(model=m,C=C,thr=t,p=pdev_oof,p_full=pall)

def fit_band(p,thr,pbcol):
    """uncertainty band around thr: uncertain -> PB decides (frozen 10.5) when PB is available."""
    best=(dev_f1(p>=thr),0.0)
    for w in (0.02,0.05,0.1,0.15,0.2,0.3,0.4):
        unc=(np.abs(p-thr)<=w)&df[pbcol].notna().values
        call=np.where(unc, df[pbcol].fillna(-1e9).values>=10.5, p>=thr)
        f1=dev_f1(call)
        if f1>best[0]+1e-6: best=(f1,w)
    return best[1]

rows=[]; ENG={}
PBIDX={"pb24":0,"pb48":1,"pb64":2,"pb96":3,"pb128":4,"pb0":5,"pbblk_med":6,"pbblk_min":7,"pbblk_max":8}
def write_model(name,kind,feats=None,L=None,w=0.0,B=False,pbcol="pb48",vetos=()):
    fn=name.replace("+","_")
    lines=[f"name {fn}",f"kind {kind}","kmin 3",f"routed_variant {PBIDX[pbcol]}"]
    if kind=="logistic":
        sc=L["model"].named_steps["standardscaler"]; lr=L["model"].named_steps["logisticregression"]
        t=L["thr"]; lg=lambda q: float(np.log(q/(1-q)))
        lines+= [f"intercept {float(lr.intercept_[0])!r}",f"zthr {lg(t)!r}"]
        if B: lines+=[f"zlo {lg(t-w)!r}",f"zhi {lg(t+w)!r}",f"band_variant {PBIDX[pbcol]}"] if w>0 else []
        for j,f in enumerate(feats): lines.append(f"feat {f} {float(sc.mean_[j])!r} {float(sc.scale_[j])!r} {float(lr.coef_[0][j])!r}")
    for r in vetos:
        if r: lines.append(f"veto {r[0]} {'ge' if r[1]=='>=' else 'le'} {r[2]!r}")
    open(f"../models/{fn}.model","w").write("\n".join(lines)+"\n")
def record(name,mech,call,notes="",pb_eval=0):
    for s in ("dev","holdout"):
        m=metrics(df,call,s); m.update(model=name,mechanism=mech,split=s,notes=notes,pb_evaluations_offline=pb_eval)
        rows.append(m)
Bbest,Bres=choose_B(); print("B dev F1 per PB variant:",{k:round(v,5) for k,v in Bres.items()},"-> best",Bbest)
record("V2_CONTROL","control",df.ctrl.values)
for c in PBV:
    record(f"B[{c}]","B",routed_call(c),notes="PB variant on routed loci, frozen thr 10.5")
cfgs={}
for A,B,C,D in itertools.product((0,1),repeat=4):
    if not (A or B or C or D): continue
    name="+".join(n for n,f in zip("ABCD",(A,B,C,D)) if f)
    feats=BASE+["routed_f"]+(A_FEATS if A and C else [])+((D_CTX+D_REF) if D and C else [])
    pbcol=Bbest if B else "pb48"
    if C:
        L=fit_lr(feats); p=L["p"]; thr=L["thr"]
        call=p>=thr; note=f"logreg C={L['C']} thr={thr:.2f} nfeat={len(feats)}"; pbev=0
        if B:
            w=fit_band(p,thr,pbcol); unc=(np.abs(p-thr)<=w)&df[pbcol].notna().values
            call=np.where(unc, df[pbcol].fillna(-1e9).values>=10.5, p>=thr); note+=f"; PB exception band ±{w} via {pbcol}"; pbev=int(unc.sum())
        # holdout probs come from the dev-only-fitted model
        # (p already = dev OOF on dev, full-dev model on holdout)
        cfgs[name]=dict(L=L,feats=feats)
        # engine-equivalent decision: full-dev model applied everywhere
        pf=L["p_full"]; ce=pf>=thr
        if B and w>0:
            unc=(np.abs(pf-thr)<=w)&df[pbcol].notna().values; ce=np.where(unc, df[pbcol].fillna(-1e9).values>=10.5, pf>=thr)
        ENG[name]=df.pos0.values[ce]
        write_model(name,"logistic",feats,L,(w if B else 0.0),B,pbcol)
    else:
        call=routed_call(pbcol) if B else df.ctrl.values.copy(); note=f"PB variant {pbcol}" if B else ""; pbev=int((df.routed==1).sum()) if B else 0
        rules=[]
        for flag,group in ((A,A_FEATS),(D,D_CTX+D_REF)):
            if flag:
                r=veto_search(group,call); rules.append(r); call=apply_veto(call,r)
        if rules: note+=" veto:"+json.dumps(rules)
        ENG[name]=df.pos0.values[call]
        write_model(name,"stump",vetos=[r for r in rules if r],pbcol=pbcol)
    record(name,name,call,note,pbev)
import pickle; pickle.dump(ENG,open(f"{W}/eng_expected_calls.pkl","wb"))
out=pd.DataFrame(rows)
out.to_csv(f"{W}/ablation_offline.csv",index=False)
pd.set_option("display.width",250)
print(out[out.split=="dev"][["model","TP","FP","FN","precision","recall","F1"]].round(5).to_string())
print(out[out.split=="holdout"][["model","TP","FP","FN","precision","recall","F1"]].round(5).to_string())
print(out[["model","notes"]].drop_duplicates("model").to_string())
