import csv, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
R="/home/mark/Documents/Projects/AI_DNA_ANALYZER/research/"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7,"axes.linewidth":0.6,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":300})
BLUE="#0072B2";ORG="#D55E00";GRN="#009E73";GREY="#6b7280";PUR="#CC79A7";YEL="#E69F00"

# ---------- Fig 1: overview ----------
fig,ax=plt.subplots(figsize=(6.9,2.3)); ax.set_xlim(0,100); ax.set_ylim(0,34); ax.axis("off")
def box(x,y,w,h,title,body,ec,fc="white",ls="-"):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.1,rounding_size=0.8",ec=ec,fc=fc,lw=1,ls=ls))
    ax.text(x+w/2,y+h-1.6,title,ha="center",va="top",fontsize=5.3,fontweight="bold")
    ax.text(x+w/2,y+h-5.2,body,ha="center",va="top",fontsize=4.9,color="#333",linespacing=1.3)
def arrow(x1,y1,x2,y2,c="#222"):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=6,lw=0.8,color=c))
y=17;h=15.5;w=12.6;gap=1.9
xs=[0.4+i*(w+gap) for i in range(7)]
specs=[("Input","BAM, MAPQ ≥ 20\nreference\nGIAB BED","#9ca3af","#f9fafb"),
("Extraction","64-bp windows\n9 count channels\none pileup pass","#6b7280","#f3f4f6"),
("Binomial screen","ε = 0.01\ncall if LLR ≥ 7.0\nevery locus","#0072B2","#e6f1fa"),
("Router","|LLR − 7.0| ≤ 5.41\n0.04–0.16 % of\nloci (15×)","#E69F00","#fff6e6"),
("PB model","per-read error q\n≤ 48 reads\ncall if LLR ≥ 10.5","#0072B2","#e6f1fa"),
("Allele calls","V1 / V2 output\nMethod C\n(GT, GQ)","#009E73","#e8f7f1"),
("V2.x classifier","candidates, k ≥ 3\n(0.29 % of loci)\nlogistic, 32 inputs","#CC79A7","#faf0f6")]
for (t,bd,ec,fc),x in zip(specs,xs):
    box(x,y,w,h,t,bd,ec,fc,ls="--" if t.startswith("V2.x") else "-")
for i in range(6):
    arrow(xs[i]+w+0.2,y+h/2,xs[i+1]-0.2,y+h/2, PUR if i==5 else "#222")
def bracket(x1,x2,yy,label,c):
    ax.plot([x1,x1,x2,x2],[yy+1.1,yy,yy,yy+1.1],color=c,lw=1.2)
    ax.text((x1+x2)/2 - (3 if x1>80 else 0),yy-0.8,label,ha="center",va="top",fontsize=5.3,color=c,linespacing=1.3)
bracket(xs[1],xs[4]+w,10.5,"V1 (frozen protocol): Python, optional C/htslib extraction.\nConstants 7.0 / 10.5 / 5.41 selected once, never changed.",GREY)
bracket(xs[6],xs[6]+w,10.5,"V2.x: changes calls;\nnot in the V1 protocol",PUR)
bracket(xs[1],xs[5]+w,3.0,"V2: same decision logic as a multi-threaded C++20 engine; VCFs byte-identical to V1 (HG002 chr20, 300×)",BLUE)
ax.text(0.4,33.8,"Unrouted loci keep the binomial call. Not drawn: PB-only control (PB at every locus), used for the V1 paired comparison.",fontsize=4.9,color=GREY,va="bottom")
fig.savefig("fig1_overview.png",bbox_inches="tight",pad_inches=0.03); plt.close(fig)

# ---------- Fig 2: forest (V1) + V2.x deltas ----------
v1=[("v13  HG002 chr21 (13 regions)",0.103,0.0,0.22,"PRESERVED/IMPROVED",GREY),
("v14  HG002, 4×1 Mb\n[experiment NEGATIVE]",0.581,0.075,1.089,"pooled IMPROVED; exp. NEGATIVE",ORG),
("v15  HG002, 4 contigs",2.686,1.968,3.414,"IMPROVED",BLUE),
("v16  HG002 chr13",-0.035,-0.171,0.070,"PRESERVED",GREY),
("v17  HG002 chr17 segdup",1.892,1.250,2.615,"IMPROVED",BLUE),
("v18  HG003 chr13",0.262,0.076,0.483,"IMPROVED",BLUE),
("v19  HG004 chr2/3/5",1.014,0.787,1.254,"IMPROVED",BLUE),
("v20  HG002 chr20 (whole)\n[DEGRADED]",-0.0564,-0.105,-0.0141,"DEGRADED",ORG)]
rows={}
for r in csv.DictReader(open(R+"V2X_BOOTSTRAP_DELTAS_VS_V2.csv")):
    if r["split"]=="holdout": rows[r["model"]]=(float(r["dF1"]),float(r["dF1_lo"]),float(r["dF1_hi"]))
order=[("B","B: PB on 48-read blocks"),("D","D: neighbour context"),("A","A: read-evidence veto"),("C","C: classifier, binomial inputs"),("A+C","A+C (V2.x)"),("A+B+C+D","A+B+C+D")]
fig,(a,b)=plt.subplots(1,2,figsize=(6.9,2.9),gridspec_kw={"width_ratios":[1.15,1],"wspace":0.9})
a.axvspan(-1,1,color="#e5e7eb",zorder=0)
for i,(lab,d,lo,hi,v,c) in enumerate(v1):
    yy=len(v1)-1-i; a.plot([lo,hi],[yy,yy],color=c,lw=1.4); a.plot(d,yy,"o",color=c,ms=3.8)
a.set_yticks(range(len(v1))); a.set_yticklabels([x[0] for x in v1][::-1],fontsize=6)
a.axvline(0,color="k",lw=0.6); a.set_xlim(-1.3,3.8)
a.set_xlabel("ΔF1 (cascade − PB-only), ×10⁻³\ninternal evaluator, depths pooled; 95 % locus-bootstrap CI;\ngrey band = pre-specified margin |ΔF1| < 0.001",fontsize=6.2)
a.set_title("(a) V1: routing vs PB-everywhere",fontsize=7,loc="left",fontweight="bold")
for i,(k,lab) in enumerate(order):
    d,lo,hi=rows[k]; yy=len(order)-1-i; c=PUR if k=="A+C" else BLUE
    b.plot([lo*1e3,hi*1e3],[yy,yy],color=c,lw=1.4); b.plot(d*1e3,yy,"o",color=c,ms=3.8)
b.set_yticks(range(len(order))); b.set_yticklabels([x[1] for x in order][::-1],fontsize=6)
b.axvline(0,color="k",lw=0.6); b.set_xlim(-1,28)
b.set_xlabel("ΔF1 vs V2 control, ×10⁻³\nHG002 chr20 holdout blocks; 95 % block-bootstrap CI",fontsize=6.2)
b.set_title("(b) V2.x: components, holdout",fontsize=7,loc="left",fontweight="bold")
fig.savefig("fig2_results.png",bbox_inches="tight",pad_inches=0.03); plt.close(fig)

# ---------- Fig 3: precision-recall + error counts ----------
import numpy as np
fig,(a,b)=plt.subplots(1,2,figsize=(6.9,2.9),gridspec_kw={"width_ratios":[1.2,1],"wspace":0.32})
rec=np.linspace(0.97,1.0005,300)
for f in (0.975,0.98,0.985,0.99,0.995):
    p=f*rec/(2*rec-f); m=(p<=1.0005)&(p>0.94)
    a.plot(rec[m],p[m],color="#d1d5db",lw=0.6,zorder=0)
    r0=0.9965; p0=f*r0/(2*r0-f)
    if 0.946<p0<1.0: a.text(r0,p0+0.0012,f"F1={f}",fontsize=4.8,color="#9ca3af",ha="right")
pts=[("V1/V2 control",0.981481,0.963543,GREY,"s",True),("DeepVariant 1.10.0",0.995836,0.999395,ORG,"D",True),("GATK 4.6.2.0",0.996117,0.987347,GRN,"^",True),
("V2.x A+C, whole chr20*",0.981117,0.997563,BLUE,"o",True),("V2.x A+C, chr20 holdout",0.981781,0.996362,BLUE,"o",False)]
for n,r_,p_,c,m,fill in pts:
    a.plot(r_,p_,m,color=c,mfc=c if fill else "white",ms=5.5,label=n,mew=1.1)
a.plot(0.980603,0.997567,"*",color=PUR,ms=8,label="V2.x A+C, HG003 chr8 (other data)")
a.set_xlim(0.975,1.0); a.set_ylim(0.945,1.001)
a.set_xlabel("Recall (hap.py, SNP, PASS)"); a.set_ylabel("Precision")
a.legend(fontsize=5.3,loc="upper center",bbox_to_anchor=(0.45,-0.17),ncol=2,frameon=False,handlelength=1.2)
a.set_title("(a) HG002 chr20, 300×: precision vs recall",fontsize=7,loc="left",fontweight="bold")
# panel b: FP and FN
lab=["V2 control\nchr20 holdout","V2.x A+C\nchr20 holdout","V2.x A+C\nHG003 chr8"]
fp=[1323,122,446]; fn=[606,620,3618]
x=np.arange(3); wd=0.36
b.bar(x-wd/2,fp,wd,color=ORG,label="false positives"); b.bar(x+wd/2,fn,wd,color=BLUE,label="false negatives")
for xi,v in zip(x-wd/2,fp): b.text(xi,v+40,f"{v:,}",ha="center",fontsize=5.5)
for xi,v in zip(x+wd/2,fn): b.text(xi,v+40,f"{v:,}",ha="center",fontsize=5.5)
b.set_xticks(x); b.set_xticklabels(lab,fontsize=5.8); b.set_ylabel("SNP count (hap.py)"); b.set_ylim(0,4100)
b.legend(fontsize=5.5,frameon=False,loc="upper left")
b.set_title("(b) error counts (different sizes)",fontsize=7,loc="left",fontweight="bold")
fig.savefig("fig3_pr_errors.png",bbox_inches="tight",pad_inches=0.03); plt.close(fig)
