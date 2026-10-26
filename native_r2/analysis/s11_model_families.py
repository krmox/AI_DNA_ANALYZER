"""Mechanism C: candidate model families (logistic regression = calibrated linear score; shallow decision trees; GBT only as an
information CEILING, not a production candidate).  Fit on DEVELOPMENT only; block-CV for threshold; holdout applied once."""
import sys; sys.path.insert(0,'.')
from lib import *
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import GroupKFold
LOGF=("depth","n","k","nref","k2","alt_lowbq")
def prep(X):
    X=X.copy()
    for c in LOGF:
        if c in X: X[c]=np.log1p(X[c])
    return X
df=load(); df["routed_f"]=df.routed.astype(float); isdev=(df.split=="dev").values; dd=df[isdev]
fs={"C (Binomial/router outputs)":BASE+["routed_f"],"A+C":BASE+["routed_f"]+A_FEATS,"C+D":BASE+["routed_f"]+D_CTX+D_REF,"A+C+D":BASE+["routed_f"]+A_FEATS+D_CTX+D_REF}
mods={"logreg(C=0.1)":lambda: make_pipeline(StandardScaler(),LogisticRegression(C=0.1,max_iter=3000)),
      "tree depth3":lambda: DecisionTreeClassifier(max_depth=3,min_samples_leaf=50,random_state=0),
      "tree depth5":lambda: DecisionTreeClassifier(max_depth=5,min_samples_leaf=50,random_state=0),
      "GBT depth4 (ceiling only)":lambda: HistGradientBoostingClassifier(max_depth=4,max_iter=200,learning_rate=0.1,random_state=0)}
rows=[]
for fn,f in fs.items():
    X=prep(df[f]).values; y=df.y.values
    for mn,mk in mods.items():
        p=np.zeros(isdev.sum())
        for tr,te in GroupKFold(8).split(X[isdev],y[isdev],dd.block):
            p[te]=mk().fit(X[isdev][tr],y[isdev][tr]).predict_proba(X[isdev][te])[:,1]
        best=max((metrics(dd,p>=t,"dev")["F1"],t) for t in np.linspace(0.05,0.95,91))
        m=mk().fit(X[isdev],y[isdev]); ph=m.predict_proba(X)[:,1]
        d=metrics(dd,p>=best[1],"dev"); h=metrics(df,ph>=best[1],"holdout")
        rows.append(dict(features=fn,model=mn,thr=round(best[1],2),dev_oof_F1=round(d["F1"],5),dev_FP=d["FP"],dev_FN=d["FN"],hold_F1=round(h["F1"],5),hold_TP=h["TP"],hold_FP=h["FP"],hold_FN=h["FN"]))
        print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(f"{W}/model_families.csv",index=False)
