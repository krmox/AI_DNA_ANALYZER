"""Paired block bootstrap (resampling whole 2-Mb blocks; loci within a block are correlated) of every finalist vs V2 control,
plus a cross-check of the locus evaluator against real hap.py TP/FP/FN."""
import sys, os; sys.path.insert(0,'.')
from lib import *
from blockeval import per_block, boot_delta, f1
import glob
runs=sorted(glob.glob(f"{W}/runs/M_*_t8_warm"))
NAMES={"A":"A","B":"B","C":"C","D":"D","A_B":"A+B","A_C":"A+C","A_D":"A+D","B_C":"B+C","B_D":"B+D","C_D":"C+D","A_B_C":"A+B+C","A_B_D":"A+B+D","A_C_D":"A+C+D","B_C_D":"B+C+D","A_B_C_D":"A+B+C+D"}
bc=per_block(f"{W}/runs/V2_CONTROL_t8_warm/ai_cascade_chr20_300x.vcf")
rows=[]; xc=[]
for d in runs:
    key=os.path.basename(d)[2:-8]
    if key not in NAMES: continue
    bm=per_block(f"{d}/ai_cascade_chr20_300x.vcf")
    for split in ("dev","holdout","full"):
        pt,lo,hi=boot_delta(bm,bc,split)
        rows.append(dict(model=NAMES[key],split=split,dF1=pt[0],dF1_lo=lo[0],dF1_hi=hi[0],dPrecision=pt[1],dP_lo=lo[1],dP_hi=hi[1],dRecall=pt[2],dR_lo=lo[2],dR_hi=hi[2],dFP=pt[3],dFP_lo=lo[3],dFP_hi=hi[3],dFN=pt[4],dFN_lo=lo[4],dFN_hi=hi[4],
                         significant_F1=bool(lo[0]>0 or hi[0]<0)))
    tot=bm[["TP","FP","FN"]].sum()
    xc.append(dict(model=NAMES[key],eval_TP=int(tot.TP),eval_FP=int(tot.FP),eval_FN=int(tot.FN)))
out=pd.DataFrame(rows); out.to_csv(f"{W}/bootstrap_deltas.csv",index=False)
pd.DataFrame(xc).to_csv(f"{W}/evaluator_totals.csv",index=False)
pd.set_option("display.width",250)
print(out[out.split=="holdout"].round(5).to_string())
print(pd.DataFrame(xc).to_string())
