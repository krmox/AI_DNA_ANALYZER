# Read-only descriptive FN taxonomy for a hap.py run on chr8 (HG003). Frame rebuilt as in native_v2x/src/frame.cpp.
import sys, json, numpy as np, pysam
from collections import Counter
tag = sys.argv[1]
E="/mnt/archive/AI_DNA_ANALYZER_final_eval"; R="/mnt/archive/AI_DNA_ANALYZER_depth_robust"
REF="/mnt/archive/AI_DNA_ANALYZER_benchmark/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
L=145138636; W=64
bed=sorted((int(a[1]),int(a[2])) for a in (l.split() for l in open(f"{E}/truth/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed")) if a[0]=="chr8")
cover=np.zeros(L,np.uint8)
for s,e in bed: cover[s:e]=1
fa=pysam.FastaFile(REF); seq=np.frombuffer(fa.fetch("chr8",0,L).upper().encode(),dtype=np.uint8)
isN=(seq==ord('N'))
cs=np.concatenate([[0],np.cumsum(cover,dtype=np.int64)]); cn=np.concatenate([[0],np.cumsum(isN,dtype=np.int64)])
win=np.zeros(L,np.uint8)   # 1 = retained window
partial=np.zeros(L,np.uint8) # window touching BED but not fully inside
for s in range(0,L-W+1,W):
    c=cs[s+W]-cs[s]
    if c==0: continue
    if c<W: partial[s:s+W]=1; continue
    if cn[s+W]-cn[s]>W//2: continue
    win[s:s+W]=1
# distance of each BED base to nearest BED edge (for boundary-distance stat)
tv=pysam.VariantFile(f"{E}/truth/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
indel_lab=np.zeros(L,np.uint8)
for r in tv.fetch("chr8"):
    g=r.samples[0].get("GT")
    if g is None or all(a in (0,None) for a in g): continue
    for alt in (r.alts or ()):
        if alt.startswith("<"): continue
        if len(alt)!=len(r.ref):
            off=0 if len(alt)>len(r.ref) else 1; span=1 if len(alt)>len(r.ref) else max(1,len(r.ref)-len(alt))
            indel_lab[r.start+off:r.start+off+span]=1
vf=pysam.VariantFile(f"{R}/happy/out/{tag}.vcf.gz")
c=Counter(); tot=0
import bisect
starts=[b[0] for b in bed]
for rec in vf:
    s=rec.samples["TRUTH"]
    if s.get("BVT")!="SNP" or s.get("BD") not in("TP","FN"): continue
    p=rec.pos-1; fn=s.get("BD")=="FN"; tot+=1
    if win[p]: cat="in_window_indel_label" if indel_lab[p] else "in_window_frame"
    elif partial[p]: cat="outside_partial_bed_window"
    elif isN[(p//W)*W:(p//W)*W+W].sum()>W//2: cat="outside_N_filter"
    else: cat="outside_other(tail/uncovered)"
    c[("FN" if fn else "TP",cat)]+=1
out={}
for cat in sorted({k[1] for k in c}):
    out[cat]={"TP":c[("TP",cat)],"FN":c[("FN",cat)]}
print(tag,"truth SNP",tot,json.dumps(out,indent=1))
