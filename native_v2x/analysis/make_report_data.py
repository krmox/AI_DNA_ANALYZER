"""Verify evaluator==hap.py, write V2X_CORRECTNESS_RESULTS.csv, copy supplementary CSVs into research/."""
import shutil, hashlib, glob, os
import pandas as pd, numpy as np
V2W="/mnt/archive/AI_DNA_ANALYZER_v2x/work"; OUT=str(__import__("pathlib").Path(__file__).resolve().parents[2]/"research")
ev=pd.read_csv(f"{V2W}/evaluator_totals.csv"); ab=pd.read_csv(f"{OUT}/V2X_ABLATION_RESULTS.csv")
full=ab[ab.development_or_holdout.str.startswith("full")].copy(); full["m"]=full.model.str.replace("V2.x_","")
mm=ev.merge(full[["m","TP","FP","FN"]],left_on="model",right_on="m")
ok=((mm.eval_TP==mm.TP)&(mm.eval_FP==mm.FP)&(mm.eval_FN==mm.FN)).sum()
print("evaluator == hap.py (TP,FP,FN) on full chr20:",ok,"/",len(mm))
rows=[]
def add(check,scope,result,detail): rows.append(dict(check=check,scope=scope,result=result,detail=detail))
add("V2 control reproduced by V2.x binary (no model): cascade VCF byte-identical to frozen V1 VCF","full chr20, 56,266,816 loci","PASS","cmp of VCF with feature layer ON (REGRESS run) and OFF (control run)")
add("3,491-locus gate, V2.x binary without model","3,491 loci","PASS","17/17 field checks, 0 mismatches incl. 1,340,544 V1 tensor bytes; PB LLR bit-exact (max|d|=0)")
add("Feature dump unchanged by PbStore refactor","163,705 candidates x 60 columns","PASS","cmp of TSV before/after refactor: byte-identical")
add("PB(48) from the new PB path == V2 routed PB","674 routed candidates","PASS","bit-identical doubles, max|d|=0")
add("Model path B(cap48) (stump, no veto) == V2 control","full chr20","PASS","VCF body byte-identical to V2 control")
add("C++ unit tests (pb_llr_n vs pb_native_v2 vs pb_reference; PB variants; model parse/logit/veto)","20,000 random loci + hand-built cases","PASS","tests/test_v2x.cpp; logit matches Python to 1e-12")
add("Cheap read-evidence features vs independent pysam re-implementation","300 candidate loci x 16 features (4,800 comparisons)","PASS","0 mismatches (tests/test_features_pysam.py)")
add("Local-context features vs independent re-implementation","80 loci x 16 features","PASS (with documented edge effect)","76 loci exact; 4 loci differ only in lmr25/depth_ratio25/dinuc10 and all lie 3-14 bp from a task start (context is truncated at task edges; deterministic, task partition is thread-independent)")
add("Engine calls == offline (independent numpy/sklearn) decisions","15 models, full chr20, in-BED calls","PASS","0 only-engine, 0 only-expected calls for all 15 (logs: work/gate_engine_vs_offline.log)")
add("Non-candidate loci (k<3) decisions == V2 control","15 models","PASS","the one non-candidate control call is identical in every model")
add("Determinism: VCF body md5 identical across 1/2/4/8 threads, cold/warm, 3 repetitions","5 models timed (control, A, C, A+C, A+B+C+D)","PASS","1 unique md5 per model over all 28 timing runs")
add("Locus evaluator == real hap.py 0.3.15 (TP, FP, FN)","15 models, full chr20",("PASS" if ok==len(mm) else "FAIL"),f"{ok}/{len(mm)} identical; used only for bootstrap and offline ablation, all reported metrics are real hap.py")
add("hap.py reproduces V2 control metrics (podman instead of docker)","full chr20","PASS","TP 70012 FP 2649 FN 1321 F1 0.972429 FP.gt 38 FP.al 159 = archived V2 hap.py run")
pd.DataFrame(rows).to_csv(f"{OUT}/V2X_CORRECTNESS_RESULTS.csv",index=False)
cp={f"{V2W}/timing.csv":"V2X_TIMING_RESULTS.csv",f"{V2W}/bootstrap_deltas.csv":"V2X_BOOTSTRAP_DELTAS_VS_V2.csv",f"{V2W}/bootstrap_pairs.csv":"V2X_BOOTSTRAP_PAIRS.csv",
    f"{V2W}/feature_audit.csv":"V2X_FEATURE_COST_AUDIT.csv",f"{V2W}/error_decomposition.csv":"V2X_ERROR_DECOMPOSITION.csv",f"{V2W}/context_window_ablation.csv":"V2X_D_CONTEXT_WINDOW_ABLATION.csv",
    f"{V2W}/model_families.csv":"V2X_C_MODEL_FAMILIES.csv",f"{V2W}/pb_exception_band.csv":"V2X_PB_EXCEPTION_BAND.csv",f"{V2W}/ablation_offline.csv":"V2X_ABLATION_OFFLINE_DEV_OOF.csv"}
for a,b in cp.items(): shutil.copy(a,f"{OUT}/{b}")
sh=[dict(model_file=os.path.basename(f),sha256=hashlib.sha256(open(f,'rb').read()).hexdigest()) for f in sorted(glob.glob("/mnt/archive/AI_DNA_ANALYZER_v2x/repo/native_v2x/models/*.model"))]
pd.DataFrame(sh).to_csv(f"{OUT}/V2X_MODEL_FILES_SHA256.csv",index=False)
print("written")
