import sys; sys.path.insert(0,'.')
from lib import *
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import GroupKFold
df=load(); dev=df[df.split=="dev"].reset_index(drop=True)
print("EX",EX, "ctrl dev",metrics(df,df.ctrl.values,"dev"))
def prep(X):
    X=X.copy()
    for c in ("depth","n","k","nref","k2","alt_lowbq"):
        if c in X: X[c]=np.log1p(X[c])
    return X
def oof(feats, mk, d=dev, nf=8):
    X=prep(d[feats]).values; y=d.y.values; p=np.zeros(len(d))
    for tr,te in GroupKFold(nf).split(X,y,d.block):
        m=mk().fit(X[tr],y[tr]); p[te]=m.predict_proba(X[te])[:,1]
    return p
def best_thr(d,p):
    best=(-1,0)
    for t in np.linspace(0.05,0.95,91):
        f=metrics(d,p>=t,"dev")["F1"]
        if f>best[0]: best=(f,t)
    return best
lr=lambda: make_pipeline(StandardScaler(),LogisticRegression(C=1.0,max_iter=2000))
gb=lambda: HistGradientBoostingClassifier(max_depth=4,max_iter=200,learning_rate=0.1)
dt=lambda: DecisionTreeClassifier(max_depth=5,min_samples_leaf=50)
sets={"C: base":BASE,"A+C":BASE+A_FEATS,"D+C":BASE+D_CTX+D_REF,"A+D+C":BASE+A_FEATS+D_CTX+D_REF}
print("control dev",metrics(dev,dev.ctrl.values,"dev"))
for name,f in sets.items():
    for mn,mk in (("logreg",lr),("tree5",dt),("GBT(ceiling)",gb)):
        p=oof(f,mk); f1,t=best_thr(dev,p)
        m=metrics(dev,p>=t,"dev")
        print(f"{name:8s} {mn:13s} thr={t:.2f} F1={m['F1']:.5f} P={m['precision']:.4f} R={m['recall']:.4f} TP={m['TP']} FP={m['FP']} FN={m['FN']}")
