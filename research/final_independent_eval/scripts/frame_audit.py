#!/usr/bin/env python3
"""M-1 window-frame audit for chr8 (descriptive, after the run): rebuild the V1/V2 window list from the confident BED
(global 64-bp tiling, whole windows fully inside the BED) and check it against the engine's n_windows before using it."""
import json, sys, numpy as np, pysam
L=145138636; CH=500000; W=64
bed=[(int(a[1]),int(a[2])) for a in (l.split() for l in open("truth/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")) if a[0]=="chr8"]
cover=np.zeros(L,dtype=np.uint8)
for s,e in bed: cover[s:e]=1
cs=np.concatenate([[0],np.cumsum(cover,dtype=np.int64)])
wins=[]
for s in range(0,L-W+1,W):   # global 64-bp tiling as in native_v2x/src/frame.cpp (window fully inside the BED)
    if cs[s+W]-cs[s]==W: wins.append(s)
wins=np.array(wins); n_win=len(wins)
eng=json.load(open("runs/A_C_chr8/stdout.json"))
print("reconstructed windows",n_win,"engine n_windows",eng["n_windows"])
ok = n_win==eng["n_windows"]
inwin=np.zeros(L,dtype=bool)
for s in wins: inwin[s:s+W]=True
vf=pysam.VariantFile("happy/out/ai_A_C_chr8.vcf.gz")
tot=fn=fn_out=tot_out=0; snp_in_bed=0
for rec in vf:
    s=rec.samples["TRUTH"]
    if s.get("BVT")!="SNP": continue
    bd=s.get("BD"); p=rec.pos-1
    if bd in ("TP","FN"):
        tot+=1
        out = not inwin[p]
        tot_out+=out
        if bd=="FN":
            fn+=1; fn_out+=out
res=dict(reconstruction_matches_engine_window_count=ok, n_windows_reconstructed=int(n_win), n_windows_engine=int(eng["n_windows"]),
         n_loci_engine=int(eng["n_loci"]), truth_snp_total_happy=tot, truth_snp_outside_retained_windows=int(tot_out),
         fraction_outside=tot_out/tot, fn_total=fn, fn_outside_retained_windows=int(fn_out), fraction_of_fn_outside=fn_out/fn if fn else None,
         recall_upper_bound_from_window_filter=(tot-tot_out)/tot)
json.dump(res,open("results_chr8_frame_audit.json","w"),indent=1); print(json.dumps(res,indent=1))
