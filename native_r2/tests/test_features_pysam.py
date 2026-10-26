"""Independent re-implementation (pysam, pure Python) of the cheap read-evidence features for sampled candidate loci,
compared with the C++ engine's dump.  Same read admission as the engine: MAPQ>=20, not dup/qcfail/secondary/supplementary,
proper-pair, BQ>=13, base in ACGT, htslib overlap handling (pysam default)."""
import sys, numpy as np, pandas as pd, pysam
sys.path.insert(0,"analysis")
A="/mnt/archive/AI_DNA_ANALYZER_benchmark"; BAM=f"{A}/HG002_GRCh38_chr20/HG002.GRCh38.300x_chr20.bam"
W="/mnt/archive/AI_DNA_ANALYZER_v2x/work"
df=pd.read_pickle(f"{W}/feat/features_labeled.pkl")
rng=np.random.default_rng(20260925)
s=pd.concat([df[df.tp_pos].sample(120,random_state=1),df[~df.tp_pos].sample(120,random_state=2),df[df.routed==1].sample(60,random_state=3)]).drop_duplicates("pos0")
bam=pysam.AlignmentFile(BAM)
BASES="ACGT"
bad={}; n=0
for _,r in s.iterrows():
    pos=int(r.pos0); alt=BASES[int(r.alt)]; ref=BASES[int(r.ref_idx)-1]
    alt_bq=[];ref_bq=[];alt_mq=[];ref_mq=[];afw=0;arv=0;aend=[];anm=[];aclip=0;starts=set();rfw=0;rrv=0;lowbq=0;hi=0;lt40=0;near=0
    nplp=0; nadm=0
    for col in bam.pileup("chr20",pos,pos+1,truncate=True,stepper="samtools",max_depth=8000,min_base_quality=0,min_mapping_quality=0,ignore_overlaps=True,ignore_orphans=True):
        nplp=col.nsegments if False else len(col.pileups)
        for p in col.pileups:
            a=p.alignment
            if a.mapping_quality<20 or a.is_duplicate or a.is_qcfail or a.is_secondary or a.is_supplementary: 
                continue
            nadm+=1
            if p.is_del or p.is_refskip: continue
            q=a.query_qualities[p.query_position] if a.query_qualities is not None else 0
            b=a.query_sequence[p.query_position]
            if b not in BASES: continue
            if q<13:
                if b==alt: lowbq+=1
                continue
            rev=a.is_reverse
            L=a.query_length; e=min(p.query_position,max(L-1-p.query_position,0))
            cig=a.cigartuples
            if b==alt:
                alt_bq.append(q);alt_mq.append(a.mapping_quality);(arv if rev else afw).__class__
                if rev: arv+=1
                else: afw+=1
                aend.append(e);near+=e<10;hi+=q>=20;lt40+=a.mapping_quality<40
                anm.append(a.get_tag("NM") if a.has_tag("NM") else 0)
                aclip+=(cig[0][0]==4 or cig[-1][0]==4); starts.add((a.reference_start<<1)|int(rev))
            elif b==ref:
                ref_bq.append(q);ref_mq.append(a.mapping_quality)
                if rev: rrv+=1
                else: rfw+=1
    k=len(alt_bq); exp=dict(k=k,nref=len(ref_bq),alt_bq_mean=np.mean(alt_bq) if k else 0,ref_bq_mean=np.mean(ref_bq) if ref_bq else 0,
        alt_mq_mean=np.mean(alt_mq) if k else 0,alt_fwd_frac=afw/k if k else 0,alt_min_strand=min(afw,arv),alt_end_mean=np.mean(aend) if k else 0,
        alt_near_end_frac=near/k if k else 0,alt_nm_mean=np.mean(anm) if k else 0,alt_clip_frac=aclip/k if k else 0,alt_uniq_frac=len(starts)/k if k else 0,
        alt_bq_hi_frac=hi/k if k else 0,alt_lowbq=lowbq,alt_mq_lt40_frac=lt40/k if k else 0,lowmq_frac=(nplp-nadm)/nplp if nplp else 0)
    n+=1
    for kf,v in exp.items():
        if not np.isclose(r[kf],v,rtol=1e-9,atol=1e-9): bad.setdefault(kf,[]).append((pos,float(r[kf]),float(v)))
print(f"loci checked: {n}; features compared per locus: {len(exp)}")
for kf in exp: print(f"  {kf:20s} mismatches: {len(bad.get(kf,[]))}")
tot=sum(len(v) for v in bad.values()); print("TOTAL MISMATCHES",tot)
if bad: print({k:v[:3] for k,v in bad.items()})
sys.exit(1 if tot else 0)
