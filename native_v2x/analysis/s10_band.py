"""C2/C3 question: is PB useful as a last-resort path for the loci a cheap classifier is unsure about?
Model = the frozen dev-fitted A+C model file (parsed independently of sklearn).  Band = |p - thr| <= w."""
import sys; sys.path.insert(0,'.')
from lib import *
def load_model(fn):
    m=dict(feat=[]); 
    for l in open(fn):
        t=l.split()
        if t[0]=="feat": m["feat"].append((t[1],float(t[2]),float(t[3]),float(t[4])))
        elif t[0] in ("intercept","zthr"): m[t[0]]=float(t[1])
    return m
M=load_model("../models/A_C.model")
df=load(); df["routed_f"]=df.routed.astype(float)
z=np.full(len(df),M["intercept"])
for n,mu,sc,co in M["feat"]:
    x=df[n].values.astype(float)
    if n in ("depth","n","k","nref","k2","alt_lowbq"): x=np.log1p(x)
    z+=co*((x-mu)/sc)
p=1/(1+np.exp(-z)); thr=1/(1+np.exp(-M["zthr"]))
rows=[]
for var in ("pb24","pb48","pbblk_min","pb0"):
    for w in (0.0,0.02,0.05,0.1,0.2,0.3,0.5):
        unc=(np.abs(p-thr)<=w)&df[var].notna().values
        call=np.where(unc,df[var].fillna(-1e9).values>=10.5,p>=thr)
        for s in ("dev","holdout"):
            m=metrics(df,call,s)
            sub=(df.split==s).values&unc
            rows.append(dict(pb_variant=var,band_w=w,split=s,n_uncertain_with_PB=int(sub.sum()),F1=round(m["F1"],5),TP=m["TP"],FP=m["FP"],FN=m["FN"],
                clf_correct_in_band=int((((p>=thr)==(df.y==1).values)&sub).sum()),pb_correct_in_band=int((((df[var].fillna(-1e9).values>=10.5)==(df.y==1).values)&sub).sum())))
out=pd.DataFrame(rows); out.to_csv(f"{W}/pb_exception_band.csv",index=False)
pd.set_option("display.width",250); print(out[out.pb_variant.isin(["pb48","pbblk_min"])].to_string())
