import sys,os; sys.path.insert(0,'.')
from lib import *
from blockeval import per_block, boot_delta
def vp(k): return f"{W}/runs/M_{k}_t8_warm/ai_cascade_chr20_300x.vcf"
pairs=[("A_C","A"),("A_C","C"),("C","A"),("A_B_C_D","A_C"),("A_C_D","A_C"),("C_D","C"),("A_B","A"),("A_C","A_B_C")]
cache={}
def pb(k):
    if k not in cache: cache[k]=per_block(vp(k)); return cache[k]
    return cache[k]
rows=[]
for a,b in pairs:
    for split in ("dev","holdout","full"):
        pt,lo,hi=boot_delta(pb(a),pb(b),split)
        rows.append(dict(model=a.replace("_","+"),vs=b.replace("_","+"),split=split,dF1=pt[0],lo=lo[0],hi=hi[0],dFP=pt[3],dFN=pt[4],significant=bool(lo[0]>0 or hi[0]<0)))
o=pd.DataFrame(rows); o.to_csv(f"{W}/bootstrap_pairs.csv",index=False); print(o.round(5).to_string())
