"""Assemble the deliverable CSVs from raw hap.py summaries, engine run metadata and the timing harness."""
import os, json, glob, re, sys
import pandas as pd, numpy as np
V2W="/mnt/archive/AI_DNA_ANALYZER_v2x/work"; OUT=str(__import__("pathlib").Path(__file__).resolve().parents[2]/"research")
NFRAME=56242697; NLOCI=56266816
NAMES={"V2_CONTROL":"V2_CONTROL","A":"A","B":"B","C":"C","D":"D","A_B":"A+B","A_C":"A+C","A_D":"A+D","B_C":"B+C","B_D":"B+D","C_D":"C+D",
       "A_B_C":"A+B+C","A_B_D":"A+B+D","A_C_D":"A+C+D","B_C_D":"B+C+D","A_B_C_D":"A+B+C+D"}
DESC={"V2_CONTROL":"frozen V2: Binomial>=7.0; router cutoff 5.4119 -> PB(48 reads)>=10.5","A":"V2 + 1 dev-fitted hard veto on a cheap read-evidence feature",
 "B":"V2 with PB on 48-read hash-blocks (block-min) instead of first-48 for routed loci; thr 10.5 frozen","C":"logistic classifier on existing Binomial/router outputs only",
 "D":"V2 + 1 dev-fitted hard veto on a local-context feature","A+C":"logistic on Binomial outputs + cheap read evidence","C+D":"logistic on Binomial outputs + local context",
 "A+B+C+D":"V2.x_ALL: logistic(A+D features) -> confident call; uncertain band -> PB block-min (frozen 10.5)"}
def timing():
    p=f"{V2W}/timing.csv"
    return pd.read_csv(p) if os.path.exists(p) else pd.DataFrame()
T=timing()
def run_meta(label):
    d=f"{V2W}/runs/{label}_t8_warm"
    if os.path.exists(f"{d}/SAME_VCF_AS"): 
        src=open(f"{d}/SAME_VCF_AS").read().strip(); hp=src
    else: hp=label
    j=json.load(open(f"{d}/stdout.json")); tt=open(f"{d}/time.txt").read()
    el=re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)",tt).group(1); p=[float(x) for x in el.split(":")]
    wall=p[-1]+60*p[-2]+(3600*p[-3] if len(p)>2 else 0)
    rss=int(re.search(r"Maximum resident set size \(kbytes\): (\d+)",tt).group(1))/1024
    return j,wall,rss,hp
def happy(label,region):
    f=f"{V2W}/happy/{label}/{region}/ai_{label}.summary.csv"
    if not os.path.exists(f): return None
    d=pd.read_csv(f); r=d[(d.Type=="SNP")&(d.Filter=="ALL")].iloc[0]
    return dict(TP=int(r["TRUTH.TP"]),FP=int(r["QUERY.FP"]),FN=int(r["TRUTH.FN"]),precision=r["METRIC.Precision"],recall=r["METRIC.Recall"],F1=r["METRIC.F1_Score"],
                FP_gt=int(r["FP.gt"]),FP_al=int(r["FP.al"]),TiTv=r["QUERY.TOTAL.TiTv_ratio"],het_hom=r["QUERY.TOTAL.het_hom_ratio"],query_total=int(r["QUERY.TOTAL"]),query_unk=int(r["QUERY.UNK"]))
REG={"dev":"development (fitted on this region; hap.py on dev BED)","holdout":"holdout (never used for fitting; hap.py on holdout BED)","full":"full_chr20 (dev+holdout mixed: engineering benchmark, NOT independent)"}
rows=[]
for fn,name in NAMES.items():
    label="V2_CONTROL" if fn=="V2_CONTROL" else f"M_{fn}"
    if not os.path.exists(f"{V2W}/runs/{label}_t8_warm/stdout.json"): continue
    j,wall,rss,hp=run_meta(label)
    tm=T[(T.model==name)&(T.plan=="warm8")] if len(T) else pd.DataFrame()
    if len(tm): wall_r=float(tm.wall_s.median()); rss_r=float(tm.peak_rss_mb.max()); tnote=f"median of {len(tm)} warm 8-thread reps (interleaved)"
    else: wall_r=float("nan"); rss_r=rss; tnote="NOT separately timed (single pipeline run 118-169 s for the same binary depending on page-cache state: not comparable); compute path identical to a timed model"
    for reg in ("dev","holdout","full"):
        h=happy(hp,reg)
        if h is None: continue
        mech="control" if name=="V2_CONTROL" else name
        rows.append(dict(model="V2_CONTROL" if name=="V2_CONTROL" else f"V2.x_{name}",mechanism=mech,development_or_holdout=REG[reg],**h,
            runtime_seconds=(round(wall_r,1) if wall_r==wall_r else float("nan")),peak_rss_mb=round(rss_r,1),threads=8,pb_evaluations=j["n_pb_evaluated"],pb_fraction=j["n_pb_evaluated"]/NFRAME,
            loci_per_second=(round(NLOCI/wall_r) if wall_r==wall_r else float("nan")),notes=f"{DESC.get(name,'')}; runtime={tnote}; pb_evaluations counted on the full chr20 run"+("" if hp==label else f"; hap.py reused from {hp} (byte-identical VCF)")))
df=pd.DataFrame(rows)
df.to_csv(f"{OUT}/V2X_ACCURACY_PERFORMANCE_BENCHMARK.csv",index=False)
ab=df.drop(columns=["runtime_seconds","peak_rss_mb","threads","loci_per_second"]).copy()
ab.insert(2,"A",ab.mechanism.map(lambda m: int("A" in m and m!="control"))); ab.insert(3,"B",ab.mechanism.map(lambda m: int("B" in m and m!="control")))
ab.insert(4,"C",ab.mechanism.map(lambda m: int("C" in m and m!="control"))); ab.insert(5,"D",ab.mechanism.map(lambda m: int("D" in m and m!="control")))
ab.to_csv(f"{OUT}/V2X_ABLATION_RESULTS.csv",index=False)
pd.set_option("display.width",250)
print(df[["model","development_or_holdout","TP","FP","FN","precision","recall","F1","FP_gt","FP_al","pb_evaluations","runtime_seconds"]].round(5).to_string())

# ---- B family: MAX_READS / PB-variant table (real hap.py) ----
BF={"Bcap24":"B1: cap 24","Bcap48":"B2: cap 48 (= V2 control PB input)","Bcap64":"B3: cap 64","Bcap96":"B4: cap 96","Bcap128":"B5: cap 128","Bcapall":"all admitted reads (no cap)","Bcapblkmed":"B6a: median PB over disjoint 48-read blocks","Bcapblkmax":"B6c: max PB over blocks","B":"B6b: min PB over blocks (dev-selected)"}
rows=[]
for fn,desc in BF.items():
    lab=f"M_{fn}"
    for reg in ("dev","holdout","full"):
        h=happy(lab,reg) if fn!="B" else happy("M_B",reg)
        if h: rows.append(dict(variant=desc,split=reg,**h))
pd.DataFrame(rows).to_csv(f"{OUT}/V2X_B_MAX_READS_RESULTS.csv",index=False)
print(pd.DataFrame(rows)[["variant","split","TP","FP","FN","F1","FP_gt","FP_al"]].round(5).to_string())
