"""Timing / memory / determinism harness (orchestration only).  Runs the C++ engine, parses /usr/bin/time -v,
appends rows to work/timing.csv.  Usage: bench_timing.py PLAN   (PLAN = warm8 | cold8 | scaling)"""
import subprocess, sys, os, re, json, hashlib, csv, time, itertools
V2W="/mnt/archive/AI_DNA_ANALYZER_v2x/work"; R="/mnt/archive/AI_DNA_ANALYZER_v2x/repo/native_v2x"
A="/mnt/archive/AI_DNA_ANALYZER_benchmark"
BAM=f"{A}/HG002_GRCh38_chr20/HG002.GRCh38.300x_chr20.bam"
REF=f"{A}/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
TRUTH=f"{A}/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
BED=f"{A}/truth_HG002_GRCh38_v4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
CTRL_BIN="/mnt/archive/AI_DNA_ANALYZER_v2x/build/dnav2"; X_BIN="/mnt/archive/AI_DNA_ANALYZER_v2x/build_x/dnav2"
MODELS={"V2_CONTROL":(CTRL_BIN,None),"A":(X_BIN,"A"),"C":(X_BIN,"C"),"A+C":(X_BIN,"A_C"),"A+B+C+D":(X_BIN,"A_B_C_D"),"B":(X_BIN,"B"),"D":(X_BIN,"D")}
CSV=f"{V2W}/timing.csv"
def run(name,threads,cold,rep,tag):
    binp,mf=MODELS[name]
    out=f"{V2W}/bench/{tag}_{name.replace('+','_')}_t{threads}_{'cold' if cold else 'warm'}_r{rep}"; os.makedirs(out,exist_ok=True)
    if cold:
        fd=os.open(BAM,os.O_RDONLY); os.posix_fadvise(fd,0,0,os.POSIX_FADV_DONTNEED); os.close(fd)
    cmd=["/usr/bin/time","-v","-o",f"{out}/time.txt",binp,"run","--bam",BAM,"--ref",REF,"--truth",TRUTH,"--bed",BED,"--out-dir",out,"--tag","chr20_300x","--threads",str(threads),"--mode","cascade"]
    if mf: cmd+=["--model",f"{R}/models/{mf}.model"]
    t0=time.time(); subprocess.run(cmd,stdout=open(f"{out}/stdout.json","w"),stderr=open(f"{out}/stderr.txt","w"),check=True); wall_py=time.time()-t0
    tt=open(f"{out}/time.txt").read()
    g=lambda p: re.search(p,tt).group(1)
    el=g(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): (\S+)"); parts=[float(x) for x in el.split(":")]
    wall=parts[-1]+60*parts[-2]+(3600*parts[-3] if len(parts)>2 else 0)
    j=json.load(open(f"{out}/stdout.json"))
    body=b"".join(l.encode() for l in open(f"{out}/ai_cascade_chr20_300x.vcf") if not l.startswith("##"))
    row=dict(plan=tag,model=name,threads=threads,cache="cold" if cold else "warm",rep=rep,wall_s=round(wall,2),user_s=float(g(r"User time \(seconds\): (\S+)")),
             sys_s=float(g(r"System time \(seconds\): (\S+)")),cpu_pct=g(r"Percent of CPU this job got: (\S+)%"),peak_rss_mb=round(int(g(r"Maximum resident set size \(kbytes\): (\d+)"))/1024,1),
             n_loci=j["n_loci"],loci_per_s=round(j["n_loci"]/wall),n_pb_evaluated=j["n_pb_evaluated"],n_candidates=j.get("n_candidates",0),n_calls=j["cascade"]["n_calls"],
             vcf_md5=hashlib.md5(body).hexdigest(),date=time.strftime("%F %T"))
    new=not os.path.exists(CSV)
    with open(CSV,"a",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(row)); 
        if new: w.writeheader()
        w.writerow(row)
    print(row,flush=True)
    # keep only the VCF of the run (drop bulky nothing else) 
plan=sys.argv[1]
if plan=="warm8":
    names=["V2_CONTROL","A","C","A+C","A+B+C+D"]
    run("V2_CONTROL",8,False,0,"warmup")      # page-cache warm-up (recorded, excluded from stats)
    for rep in (1,2,3):
        for i in range(len(names)): run(names[(i+rep)%len(names)],8,False,rep,"warm8")   # rotated order
elif plan=="cold8":
    for n in ("V2_CONTROL","A+C","A+B+C+D"): run(n,8,True,1,"cold8")
elif plan=="scaling":
    for th in (4,2,1):
        for n in ("V2_CONTROL","A+C","A+B+C+D"): run(n,th,False,1,"scaling")
