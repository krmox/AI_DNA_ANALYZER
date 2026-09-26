import sys; sys.path.insert(0,'.')
from lib import *
from sklearn.metrics import roc_auc_score
df=load(); dev=df[df.split=="dev"]
grp={}
for f in BASE: grp[f]="base(existing)"
for f in A_FEATS: grp[f]="A"
for f in D_CTX: grp[f]="D-context"
for f in D_REF: grp[f]="D-refctx"
grp["lowmq_frac"]="A"
cost={  # qualitative, from the code (features.cpp / engine.cpp compute_ctx)
 "alt_nm_mean":("bam_aux_get('NM') per ALT read (tag scan ~1e2 ns)","high(rel.)"),
 "ref_nm_mean":("bam_aux_get('NM') on <=32 REF reads","high(rel.)"),
 "alt_clip_frac":("cigar first/last op, ALT reads","low"),"ref_clip_frac":("cigar first/last op, <=32 REF reads","low"),
 "alt_uniq_frac":("sort of ALT (start,strand) keys","low"),
}
rows=[]
cand_ctrl=df[df.ctrl]
for f in BASE+A_FEATS+["lowmq_frac"]+D_CTX+D_REF:
    x=dev[f].values
    def auc(d,target,name):
        v=d[f].values; y=d[target].values
        if len(np.unique(y))<2 or np.nanstd(v)==0: return np.nan
        return roc_auc_score(y,v)
    c=dev[dev.ctrl]
    a_tpfp=auc(c,"y","")             # among V2 calls: does the feature separate TP from FP?
    a_all=auc(dev,"y","")            # among all candidates
    miss=dev[(dev.y==1)]
    miss=miss.assign(missed=(~miss.ctrl).astype(int))
    a_fn=auc(miss,"missed","")       # among true SNPs: does it separate V2-missed (FN) from V2-found?
    src="pileup(existing count)" if f in BASE else ("neighbour arrays" if f in D_CTX else ("reference bases" if f in D_REF else "same pileup column, candidate loci only"))
    c1,c2=cost.get(f,("integer accumulate over the column's reads" if grp[f]=="A" else ("O(1..50) array reads" if f not in BASE else "0 (already computed)"),"low" if f not in BASE else "none"))
    rows.append(dict(feature=f,group=grp[f],source=src,extra_bam_traversal="no",compute_cost=c1,cost_class=c2,
        extra_memory="0 (scalar accumulators)" if f not in D_CTX else "0 (task-local arrays that already exist)",
        auc_TPvsFP_among_V2_calls_dev=a_tpfp,auc_true_vs_false_all_candidates_dev=a_all,auc_V2FN_vs_found_among_true_dev=a_fn))
out=pd.DataFrame(rows).round(4)
out.to_csv(f"{W}/feature_audit.csv",index=False)
pd.set_option("display.width",220,"display.max_rows",200)
print(out[["feature","group","cost_class","auc_TPvsFP_among_V2_calls_dev","auc_true_vs_false_all_candidates_dev","auc_V2FN_vs_found_among_true_dev"]].to_string())
