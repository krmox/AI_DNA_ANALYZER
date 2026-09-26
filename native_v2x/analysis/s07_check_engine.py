"""Engine-vs-offline gate: the C++ engine's calls must equal the calls the dev-fitted model implies offline
(independent numpy/sklearn implementation), for candidate loci; non-candidate loci must equal the V2 control."""
import sys, pickle; sys.path.insert(0,'.')
from lib import *
from blockeval import vcf_calls
label,name=sys.argv[1],sys.argv[2]
ENG=pickle.load(open(f"{W}/eng_expected_calls.pkl","rb"))
df_all=pd.read_pickle(f"{W}/feat/features_labeled.pkl")
cand=set(df_all.pos0.values)
ctrl=vcf_calls(f"{W}/runs/V2_CONTROL_t8_warm/ai_cascade_chr20_300x.vcf")
non_cand_ctrl=set(p for p in ctrl.pos0 if p not in cand)
# ENG holds only in-BED candidates; candidates outside the BED are not evaluated offline -> compare inside BED
got=vcf_calls(f"{W}/runs/{label}_t8_warm/ai_cascade_chr20_300x.vcf")
gs=set(p for p in got.pos0 if BM[p]); 
exp=set(int(p) for p in ENG[name])|set(p for p in non_cand_ctrl if BM[p])
print(f"{label}: engine calls in BED {len(gs)}  expected {len(exp)}  only-engine {len(gs-exp)}  only-expected {len(exp-gs)}")
# outside-BED non-candidates must equal control
gn=set(p for p in got.pos0 if p not in cand); print("non-candidate calls: engine",len(gn),"control",len(non_cand_ctrl),"equal",gn==non_cand_ctrl)
