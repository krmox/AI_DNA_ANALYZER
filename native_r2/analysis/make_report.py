import pandas as pd, numpy as np
OUT=str(__import__("pathlib").Path(__file__).resolve().parents[2]/"research"); V2W="/mnt/archive/AI_DNA_ANALYZER_v2x/work"
B=pd.read_csv(f"{OUT}/V2X_ACCURACY_PERFORMANCE_BENCHMARK.csv"); B["reg"]=B.development_or_holdout.map(lambda s:"dev" if s.startswith("dev") else ("holdout" if s.startswith("holdout") else "full"))
BS=pd.read_csv(f"{OUT}/V2X_BOOTSTRAP_DELTAS_VS_V2.csv"); T=pd.read_csv(f"{OUT}/V2X_TIMING_RESULTS.csv")
W=T[T.plan=="warm8"]; g=W.groupby("model").agg(wall=("wall_s","median"),lo=("wall_s","min"),hi=("wall_s","max"),cpu=("user_s","median"),sysc=("sys_s","median"),rss=("peak_rss_mb","max"))
g["cpus"]=g.cpu+g.sysc; c0=g.loc["V2_CONTROL"]
def row(m,reg):
    key="V2_CONTROL" if m=="V2_CONTROL" else "V2.x_"+m
    return B[(B.model==key)&(B.reg==reg)].iloc[0]
def dF(m,reg):
    if m=="V2_CONTROL": return "—"
    r=BS[(BS.model==m)&(BS.split==reg)].iloc[0]; return f"{r.dF1:+.5f} [{r.dF1_lo:+.5f}; {r.dF1_hi:+.5f}]"
def rt(m):
    k="V2_CONTROL" if m=="V2_CONTROL" else m
    if k in g.index:
        r=g.loc[k]; return f"{r.wall:.1f}", f"{(r.wall/c0.wall-1)*100:+.1f}% wall / {(r.cpus/c0.cpus-1)*100:+.1f}% CPU", f"{r.rss:.0f}"
    return "n/t", "≈ same path", f"{row(m,'holdout').peak_rss_mb:.0f}"
order=["V2_CONTROL","A","B","C","D","A+B","A+C","A+D","B+C","B+D","C+D","A+B+C","A+B+D","A+C+D","B+C+D","A+B+C+D"]
def main_table(reg):
    L=["| Model | F1 | Precision | Recall | FP | FN | Runtime, s (8 thr, warm) | Δtime vs V2 | RAM, MB | PB calls | ΔF1 vs V2 [95% block-CI] |","|---|---|---|---|---|---|---|---|---|---|---|"]
    for m in order:
        r=row(m,reg); a,b,c=rt(m)
        nm="**V2 CONTROL**" if m=="V2_CONTROL" else ("V2.x "+m+(" (ALL)" if m=="A+B+C+D" else ""))
        L.append(f"| {nm} | {r.F1:.5f} | {r.precision:.5f} | {r.recall:.5f} | {r.FP} | {r.FN} | {a} | {b} | {c} | {r.pb_evaluations} | {dF(m,reg)} |")
    return "\n".join(L)
tab_hold=main_table("holdout"); tab_full=main_table("full"); tab_dev=main_table("dev")
# accuracy detail for key models
def acc_table(models):
    L=["| Model | Region | TP | FP | FN | Precision | Recall | F1 | FP.gt | FP.al | Ti/Tv (query) | het/hom (query) | QUERY.TOTAL | QUERY.UNK |","|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for m in models:
        for reg,nm in (("holdout","holdout"),("dev","dev (in-sample)"),("full","full chr20")):
            r=row(m,reg); L.append(f"| {'V2 CONTROL' if m=='V2_CONTROL' else 'V2.x '+m} | {nm} | {r.TP} | {r.FP} | {r.FN} | {r.precision:.6f} | {r.recall:.6f} | {r.F1:.6f} | {r.FP_gt} | {r.FP_al} | {r.TiTv:.4f} | {r.het_hom:.4f} | {r.query_total} | {r.query_unk} |")
    return "\n".join(L)
acc=acc_table(["V2_CONTROL","A","C","A+C","A+B+C+D"])
# speed table
S=T[T.plan.isin(["warm8","cold8","scaling"])]
def speed():
    L=["| Model | Threads | Cache | Runs | Wall s (median; min–max) | CPU s (user+sys) | CPU util | loci/s | Peak RSS MB |","|---|---|---|---|---|---|---|---|---|"]
    for (m,th,ca),d in S.groupby(["model","threads","cache"],sort=False):
        L.append(f"| {m} | {th} | {ca} | {len(d)} | {d.wall_s.median():.1f} ({d.wall_s.min():.1f}–{d.wall_s.max():.1f}) | {(d.user_s+d.sys_s).median():.0f} | {d.cpu_pct.astype(float).mean():.0f}% | {int(d.loci_per_s.median()):,} | {d.peak_rss_mb.max():.0f} |")
    return "\n".join(L)
sp=speed()
Bf=pd.read_csv(f"{OUT}/V2X_B_MAX_READS_RESULTS.csv")
def btab():
    L=["| PB input for routed loci (threshold 10.5 frozen) | holdout F1 | holdout FP | holdout FN | full-chr20 F1 | full FP.al |","|---|---|---|---|---|---|"]
    for v,d in Bf.groupby("variant",sort=False):
        h=d[d.split=="holdout"].iloc[0]; f=d[d.split=="full"].iloc[0]
        L.append(f"| {v} | {h.F1:.5f} | {h.FP} | {h.FN} | {f.F1:.5f} | {f.FP_al} |")
    return "\n".join(L)
bt=btab()
open(f"{V2W}/report_tables.md","w").write("## HOLD\n"+tab_hold+"\n## FULL\n"+tab_full+"\n## DEV\n"+tab_dev+"\n## ACC\n"+acc+"\n## SPEED\n"+sp+"\n## B\n"+bt+"\n")
print(open(f"{V2W}/report_tables.md").read())
