import sys,pickle; sys.path.insert(0,'.')
from lib import *
ENG=pickle.load(open(f"{W}/eng_expected_calls.pkl","rb"))
df=load(); t=truth()
tp_all=np.array(sorted(p for p in t if BM[p]))
rows=[]
for name in ("V2_CONTROL","A+C","A+B+C+D"):
    calls=df.pos0[df.ctrl].values if name=="V2_CONTROL" else np.array(sorted(ENG[name]))
    for s in ("dev","holdout"):
        d=df[df.split==s]
        truth_pos=tp_all[split_of(tp_all)==s]
        good=np.intersect1d(calls,d.pos0[d.y==1].values)
        fn=np.setdiff1d(truth_pos,good)
        iscand=np.isin(fn,d.pos0.values); wascalled=np.isin(fn,calls)
        cin=np.intersect1d(calls,d.pos0.values); fp=np.setdiff1d(cin,d.pos0[d.y==1].values)
        rows.append(dict(model=name,split=s,TP=len(good),FN_total=len(fn),FN_not_a_candidate_k_lt_3_or_out_of_frame=int((~iscand).sum()),
            FN_candidate_not_called=int((iscand&~wascalled).sum()),FN_called_but_wrong_GT=int((iscand&wascalled).sum()),FP_in_candidates=len(fp)))
print(pd.DataFrame(rows).to_string())
pd.DataFrame(rows).to_csv(f"{W}/error_decomposition.csv",index=False)
