"""Independent pysam/pyfaidx-free re-implementation of the local-context (D) features for sampled candidates."""
import sys, numpy as np, pandas as pd, pysam
A="/mnt/archive/AI_DNA_ANALYZER_benchmark"; BAM=f"{A}/HG002_GRCh38_chr20/HG002.GRCh38.300x_chr20.bam"
REF=f"{A}/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
W="/mnt/archive/AI_DNA_ANALYZER_v2x/work"
df=pd.read_pickle(f"{W}/feat/features_labeled.pkl")
s=df.sample(80,random_state=11)
bam=pysam.AlignmentFile(BAM); fa=pysam.FastaFile(REF)
B="ACGT"
def counts(lo,hi):
    """per-position (A,C,G,T,gap,ins,depth) with the engine's admission rules, positions lo..hi-1"""
    c=np.zeros((hi-lo,7),int)
    for col in bam.pileup("chr20",lo,hi,truncate=True,stepper="samtools",max_depth=8000,min_base_quality=0,min_mapping_quality=0,ignore_overlaps=True,ignore_orphans=True):
        i=col.reference_pos-lo
        for p in col.pileups:
            a=p.alignment
            if a.mapping_quality<20 or a.is_duplicate or a.is_qcfail or a.is_secondary or a.is_supplementary: continue
            c[i,6]+=1
            if p.is_del or p.is_refskip: c[i,4]+=1; continue
            q=a.query_qualities[p.query_position]
            b=a.query_sequence[p.query_position]
            if b not in B or q<13: continue
            c[i,B.index(b)]+=1
            if p.indel>0: c[i,5]+=1
    return c
bad={};chk=0
for _,r in s.iterrows():
    pos=int(r.pos0); lo=pos-30; hi=pos+31
    c=counts(lo,hi); ref=[ (B.index(x)+1 if x in B else 0) for x in fa.fetch("chr20",lo,hi).upper()]
    n=c[:,:4].sum(1); k=np.zeros(len(n)); 
    for j in range(len(n)):
        rj=ref[j]; m=c[j,:4].astype(float).copy()
        if rj>=1: m[rj-1]=-1
        k[j]=c[j,:4][int(np.argmax(m))] if m.max()>=0 else 0
    nonref=np.array([n[j]-(c[j,ref[j]-1] if ref[j]>=1 else 0) for j in range(len(n))],float)
    cand=lambda j: k[j]>=3 and n[j]>0 and k[j]/n[j]>=0.05
    i=30; ex={}
    ws=lambda w:[j for j in range(max(0,i-w),min(len(n)-1,i+w)+1) if j!=i]
    for w,nm in ((2,"nc2"),(5,"nc5"),(10,"nc10"),(25,"nc25")): ex[nm]=sum(cand(j) for j in ws(w))
    d=[abs(j-i) for j in ws(25) if cand(j)]; ex["nearest"]=min(d) if d else 26
    w25=ws(25); ex["lmr25"]=nonref[w25].sum()/n[w25].sum() if n[w25].sum()>0 else 0
    w10=ws(10); ex["lgap10"]=(c[w10,4]+c[w10,5]).sum()/c[w10,6].sum() if c[w10,6].sum()>0 else 0
    ex["depth_ratio25"]=c[i,6]/c[w25,6].mean() if c[w25,6].mean()>0 else 1.0
    ex["neigh_max_vaf25"]=max([k[j]/n[j] for j in w25 if cand(j)] or [0])
    rr=ref[i]; a=i; b=i
    while a>0 and ref[a-1]==rr and i-a<30: a-=1
    while b+1<len(ref) and ref[b+1]==rr and b-i<30: b+=1
    ex["hp_len"]=b-a+1; alt=int(r.alt)
    ex["left_is_alt"]=int(ref[i-1]==alt+1); ex["right_is_alt"]=int(ref[i+1]==alt+1)
    ex["dinuc10"]=sum(1 for j in range(i-10,i+9) if ref[j]!=0 and ref[j]==ref[j+2])
    ex["is_ts"]=int(rr>=1 and ((rr-1)^alt)==2); ex["ref_cg"]=int(rr in (2,3)); ex["cpg"]=int((rr==2 and ref[i+1]==3) or (rr==3 and ref[i-1]==2))
    chk+=1
    for kf,v in ex.items():
        if not np.isclose(r[kf],v,rtol=1e-9,atol=1e-9): bad.setdefault(kf,[]).append((pos,float(r[kf]),float(v)))
print("loci checked",chk,"features per locus",len(ex))
for kf in ex: print(f"  {kf:16s} mismatches {len(bad.get(kf,[]))}")
print("TOTAL",sum(len(v) for v in bad.values()))
if bad: print({k:v[:4] for k,v in bad.items()})
