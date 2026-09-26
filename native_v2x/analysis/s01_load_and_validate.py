import numpy as np, pandas as pd, pickle, struct
from common import *
W="/mnt/archive/AI_DNA_ANALYZER_v2x/work"
df=pd.read_csv(f"{W}/feat/features_chr20.tsv",sep="\t",float_precision="round_trip")
print(df.shape)
truth=load_truth_snps(); iv=load_bed()
df["split"]=split_of(df.pos0.values)
df["inbed"]=in_bed(df.pos0.values,iv)
df["gtc"]=methodc_gt(df.k,df.k+df.nref)
tr_alt=np.array([ (BASES[a] in truth[p][1]) if p in truth else False for p,a in zip(df.pos0,df.alt)])
df["tp_pos"]=tr_alt
df["truth_gt_ok"]=[ (p in truth and BASES[a] in truth[p][1] and truth[p][2]==g) for p,a,g in zip(df.pos0,df.alt,df.gtc)]
print("candidates",len(df),"in truth (alt match)",int(tr_alt.sum()))
# --- PB48 bit-exactness vs engine routed table (control) ---
rt=pd.read_csv(f"{W}/runs/V2_CONTROL_t8_warm/routed_loci_chr20_300x.tsv",sep="\t",float_precision="round_trip")
rt["pb"]=rt.pb_llr_hex.map(float.fromhex); rt["bl"]=rt.binomial_llr_hex.map(float.fromhex)
m=rt.merge(df[["pos0","pb48","bllr"]],on="pos0",how="left")
print("routed",len(rt),"found in cand table",m.pb48.notna().sum())
mm=m.dropna()
print("PB48 bit-identical:",int((mm.pb.values==mm.pb48.values).sum()),"/",len(mm)," max|d|",np.abs(mm.pb.values-mm.pb48.values).max())
print("bllr max|d|",np.abs(mm.bl.values-mm.bllr.values).max())
# --- evaluator vs hap.py on the V2 control VCF ---
vc=[l.split("\t") for l in open(f"{W}/runs/V2_CONTROL_t8_warm/ai_cascade_chr20_300x.vcf") if not l.startswith("#")]
calls=pd.DataFrame(dict(pos0=[int(x[1])-1 for x in vc],alt=[BASES.index(x[4]) for x in vc],gt=[x[9].split(":")[0] for x in vc]))
r=evaluate(calls,truth,lambda p: True) # region filter below
ib=lambda p: bool(in_bed(np.array([p]),iv)[0])
# fast in_bed via set
bed_mask=np.zeros(64444167+1,bool)
for a,b in iv: bed_mask[a:b]=True
r=evaluate(calls,truth,lambda p: bed_mask[int(p)])
print("my evaluator on V2 control:",{k:v for k,v in r.items() if k!='tp_pos'},"  hap.py: TP 70012 FP 2649 FN 1321 FPgt 38 FPal 159")
df.to_pickle(f"{W}/feat/features_labeled.pkl")
