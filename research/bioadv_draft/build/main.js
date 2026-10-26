const fs = require("fs");
const B = require("./build.js");
const { D, C, runs, base, P, H1, H2, figure, table, caption, afterTable, FONT, PAGE_W, PAGE_H, MARGIN, GAP, FULL, COL, BODY, SMALL, CAP, REF } = B;
const { Document, Packer, Paragraph, TextRun, AlignmentType, Footer, Header, PageNumber, SectionType, BorderStyle } = D;

const OUT = process.argv[2];
const FIG = (n) => `${__dirname}/${n}`;

// ---------- section assembly ----------
const sections = [];
let cur = null;
const open = (cols) => { cur = { cols, children: [] }; sections.push(cur); };
const add = (...xs) => { for (const x of xs.flat()) cur.children.push(x); };
const full = (...xs) => { open(1); add(...xs); open(2); };

// ---------- table definitions ----------
const t1w = [820, 1300, 1000, 1653];
const T1 = table(t1w, [
  ["Version", "What it is", "Changes calls?", "Evidence in this paper"],
  ["V1", "Binomial screen, router, PB on routed loci; Python with optional C/htslib extraction; constants fixed", "Reference behaviour", "Eight experiments against a PB-only control (Sect. 4.1–4.2)"],
  ["V2", "C++20 re-implementation of V1; no read tensor; multi-threaded", "No (byte-identical VCFs, HG002 chr20)", "Equivalence gates and runtime (Sect. 4.3)"],
  ["V2.x", "V2 plus cheap read-level evidence and a logistic classifier (A+C)", "Yes; not part of the V1 protocol", "Development, holdout, final run (Sect. 4.4–4.5)"],
], { leftAll: true });

const t2w = [560, 2150, 1500, 1500, 1050, 3146];
const T2 = table(t2w, [
  ["Exp.", "Sample / region (depth)", "PB-only F1 → cascade F1", "ΔF1 ×10^−3^ [95 % CI]", "Discordant loci (true SNPs lost)", "Verdict (pre-specified rule)"],
  ["v13", "HG002 chr21, 13 regions (mixed)", "0.98162 → 0.98172", "+0.10 [0.00, +0.22]", "15 (0)", "PRESERVED in 12/13 regions, 1 IMPROVED"],
  ["v14", "HG002, 4 × 1 Mb (full/30×/15×)", "0.955726 → 0.956306", "+0.58 [+0.08, +1.09]", "41 (5 events, 4 loci)", "Pooled IMPROVED; **experiment NEGATIVE**"],
  ["v15", "HG002, 4 held-out contigs", "0.925079 → 0.927766", "+2.69 [+1.97, +3.41]", "80 (1)", "IMPROVED"],
  ["v16", "HG002 chr13:70–73 Mb", "0.994600 → 0.994565", "−0.04 [−0.17, +0.07]", "3 (0)", "PRESERVED"],
  ["v17", "HG002 chr17:18–21 Mb (segdup)", "0.948070 → 0.949962", "+1.89 [+1.25, +2.62]", "38 (2)", "IMPROVED"],
  ["v18", "HG003 chr13:70–73 Mb", "0.993751 → 0.994013", "+0.26 [+0.08, +0.48]", "7 (0)", "IMPROVED"],
  ["v19", "HG004 chr2, chr3, chr5 (3 Mb each)", "0.986405 → 0.987419", "+1.01 [+0.79, +1.25]", "87 (1)", "IMPROVED"],
  ["v20", "HG002 chr20, whole (15×)", "0.97809 → 0.97803", "**−0.056** [−0.105, −0.014]", "10 (1)", "**DEGRADED**"],
], { hl: [2, 8] });

const t3w = [1900, 700, 750, 1423];
const T3 = table(t3w, [
  ["Configuration (HG002 chr20, 300×)", "Wall time", "Peak RSS", "Note"],
  ["V1, frozen record (1 thread)", "18,578 s", "8 GB; tensor 21.6 GB", "Another day and cache state"],
  ["V2 pb_all (both VCFs), 1 thread", "1,059.5 s", "114 MB", "≈ 17.5× vs V1, same work"],
  ["V2 pb_all, 8 threads", "182.5 s", "225 MB", "≈ 102× vs V1"],
  ["V2 cascade, 8 threads", "119.6 s", "225 MB", "PB on 707 loci only: different work"],
  ["V2 control, 8 threads (median of 3)", "127.3 s", "219 MB", "Range 118.4–129.9 s"],
  ["V2.x A+C, 8 threads (median of 3)", "128.0 s", "221 MB", "+1.0 % CPU time vs V2 control"],
], { leftAll: true });

const t4w = [1790, 600, 600, 700, 540, 553];
const T4 = table(t4w, [
  ["Model, data (tier)", "Prec.", "Recall", "F1", "FP", "FN"],
  ["V2 control, chr20 holdout", "0.9619", "0.9822", "0.97195", "1,323", "606"],
  ["V2.x A only (ALT-read base-quality veto), holdout", "0.9933", "0.9821", "0.98769", "225", "608"],
  ["V2.x A+C, chr20 even blocks (development, in-sample)", "0.9987", "0.9805", "0.98950", "49", "727"],
  ["**V2.x A+C, chr20 odd blocks (holdout)**", "0.9964", "0.9818", "**0.98902**", "122", "620"],
  ["V2 control, whole chr20 (mixed)", "0.9635", "0.9815", "0.97243", "2,649", "1,321"],
  ["V2.x A+C, whole chr20 (mixed)", "0.9976", "0.9811", "0.98927", "171", "1,347"],
  ["**V2.x A+C, HG003 chr8 (final run; no control)**", "0.9976", "0.9806", "**0.98901**", "446", "3,618"],
], { hl: [4, 7] });

const t5w = [1650, 560, 560, 650, 1353];
const T5 = table(t5w, [
  ["Caller (HG002 chr20, 300×)", "Prec.", "Recall", "F1", "Wall-clock time (unmatched conditions)"],
  ["V1 cascade", "0.9635", "0.9815", "0.97243", "≈ 5 h 10 min (record; 1 thread)"],
  ["V2 control (= V1 output)", "0.9635", "0.9815", "0.97243", "127.3 s, 8 threads"],
  ["V2.x A+C, whole chr20 (mixed)", "0.9976", "0.9811", "0.98927", "128.0 s, 8 threads"],
  ["V2.x A+C, holdout blocks", "0.9964", "0.9818", "0.98902", "as above"],
  ["DeepVariant 1.10.0 (CPU)", "0.9994", "0.9958", "0.99761", "1 h 28 min, 8 shards"],
  ["GATK HC 4.6.2.0", "0.9873", "0.9961", "0.99171", "≈ 3 h 43 min + 20 min RG step"],
], { hl: [3, 4] });

// ---------- front matter (single column) ----------
open(1);
add(
  new Paragraph({ spacing: { after: 60 }, children: [new TextRun({ text: "WORKING DRAFT FOR AUTHOR REVIEW — not submission-ready. Target journal: Bioinformatics Advances (OUP). Original Article.", font: FONT, size: 15, bold: true, color: "B00020" })] }),
  new Paragraph({ spacing: { after: 70, line: 250 }, children: [new TextRun({ text: C.title, font: FONT, size: 34, bold: true })] }),
  new Paragraph({ spacing: { after: 20 }, children: runs("Mark Kostogryzov^1,*^", base(BODY)) }),
  new Paragraph({ spacing: { after: 90 }, border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: "000000", space: 3 } }, children: runs("^1^Independent researcher. [[c0]]^*^Correspondence: m.2researcher141941241@proton.me[[/c0]]. Preprint status of the source manuscript: not peer reviewed. [[c15]]Nothing was re-run for this draft; numbers are taken from the 36-page manuscript and its result files.[[/c15]]", base(SMALL)) }),
  new Paragraph({ spacing: { after: 40 }, children: runs("**Abstract** [[c1]](≤ 200 words)[[/c1]]", base(BODY)) }),
  ...C.abstract.map((t) => P(t, { size: 19, after: 40 })),
  new Paragraph({ spacing: { before: 40, after: 80 }, border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: "000000", space: 3 } },
    children: runs("**Keywords:** " + C.keywords + ".  [[c2]]Layout and style assumptions: see note.[[/c2]]", base(18)) }),
);

// ---------- body (two columns) ----------
open(2);
add(H1("1 Introduction"), ...C.intro.map((t) => P(t)));

add(H1("2 Methods"));
add(H2("2.1 Versions and evidence tiers"));
add(P("Table 1 lists the versions. The chronology is V1, V2, V2.x. V2.x was developed after the V1 protocol had been frozen and archived (release v1.0, 2026-09-20), and no V1 constant was changed for it. Results belong to four evidence tiers, stated in each table: (i) *development* data used to choose constants, features or coefficients; (ii) *holdout* data from the same sample and chromosome, withheld from fitting; (iii) one *final independent evaluation* on another sample and chromosome; (iv) *external caller comparisons*, which are descriptive."));
add(caption("**Table 1.** Versions of AI_DNA_ANALYZER and the evidence used for each."), T1, afterTable(""));

add(H2("2.2 V1: the frozen cascade"));
add(P("At a locus with *n* usable A/C/G/T observations and *k* reads supporting the most frequent non-reference base, a binomial log-likelihood ratio (LLR~bin~) compares heterozygous and homozygous-alternative models against the null, with a fixed error ε = 0.01 (null *p* = ε/3). A locus is called if LLR~bin~ ≥ 7.0. The router sends a locus to the PB model if |LLR~bin~ − 7.0| ≤ 5.4119, that is, LLR~bin~ ∈ [1.588, 12.412]. For routed loci the call is LLR~PB~ ≥ 10.5. The PB model uses the error probability *q*~i~ = 10^−Q~i~/10^ (clipped to [10^−4^, 0.25]) of each of at most 48 reads, ordered by a hash of the read name, and an exact log-space recursion for the count distribution. The router measures uncertainty, not error: a locus at which the binomial test is confidently wrong is never routed."));
add(P("The three constants (7.0, 10.5, 5.4119) were selected once on HG002 chr21:32–40 Mb at about 15× (training, validation and test blocks by index) and never changed. The thresholds maximised validation F1. [[c20]]The router cutoff is the smallest routed fraction (0.1 % of validation loci) whose validation F1 lay within 0.0005 of PB-everywhere; every larger fraction tied at 0.9784.[[/c20]] The constants are depth-specific. Reads are admitted at MAPQ ≥ 20 and base quality ≥ 13. The region is tiled in 64-bp windows, and a window is dropped if any base lies outside the high-confidence BED (the “M-1 window filter”), which removes 1.4–1.8 % of truth SNPs from reach. A genotype layer (Method C, a binomial-likelihood GT/GQ with no fitted parameter) is applied after the allele calls and does not change them. The V1 estimand is ΔF1 = F1(cascade) − F1(PB-only), paired on the same loci, where PB-only evaluates PB at every locus with threshold 10.5."));

add(H2("2.3 V2: execution engine"));
add(P("V2 (`dnav2`) is a C++20 executable with the same decision logic. One streaming htslib pileup (Bonfield *et al.*, 2021) per task of up to 2,048 windows feeds counting, the binomial LLR and the router. Routed loci pass at most 48 (base, quality) pairs to a scalar PB kernel, so no read tensor exists. Worker threads claim tasks from an atomic counter and results are written in genomic order, which makes output independent of thread count. Exactness relies on ports of the SciPy `gammaln` and `log1p` routines, the glibc math library and a build with `-ffp-contract=off`. The acceptance criterion was byte-identical VCFs."));

add(H2("2.4 V2.x: read-evidence classifier"));
add(P("In the V2 control, 91 % of the false positives on HG002 chr20 lay outside the router’s reach, in low-VAF loci whose ALT reads have low base quality. V2.x therefore adds a second stage at *candidate loci* (≥ 3 ALT reads; 163,705 loci, 0.29 % of scored loci). Component **A** collects 24 read-level statistics in the same pileup pass (ALT/REF base and mapping quality, strand balance, distance to read end, edit distance, soft-clipping, distinct start positions, gap/indel fraction). Component **C** is a logistic regression on 32 standardised inputs: the 24 A statistics plus 8 binomial/router outputs (depth, *n*, *k*, reference count, third-base count, VAF, LLR~bin~, routed flag). Coefficients and the decision threshold are stored in a text file. Two further components were tested and not adopted: **B** (PB on more reads or on independent 48-read blocks) and **D** (features from loci within 2–25 bp). Loci with *k* < 3 keep the V2 path.")); cur.children[cur.children.length-1].__marker = "m1";

add(H1("3 Experimental design and datasets"));
add(H2("3.1 Data, truth and scoring"));
add(P("Truth sets are GIAB v4.2.1 GRCh38 benchmark VCFs and high-confidence BEDs (Zook *et al.*, 2019; Wagner *et al.*, 2022). V1 experiments used the NIST Illumina 2×250 NovoAlign BAMs of the HG002–HG004 trio (Zook *et al.*, 2016), at native depth (47–74×) and at 30× and 15× downsamples (`samtools view -s`; Danecek *et al.*, 2021), so depth cells are correlated. V2 and V2.x used HG002 chr20 at 300× (Illumina HiSeq 2500, 2×148 bp, NovoAlign). The final evaluation used HG003 chr8 at about 303× from the same GIAB sequencing project. Two scorers are used and always named. The *internal locus evaluator* compares boolean call and truth masks, ignores genotype and applies the BED per window. *hap.py* 0.3.15 (xcmp engine; Krusche *et al.*, 2019) is allele- and genotype-aware and applies the BED per record. All results are SNP-only."));

add(H2("3.2 V1 experiments and decision rule"));
add(P("Eight experiments (v13–v20) cover other regions, chromosomes and samples of the trio and one whole chromosome (Table 2). None is independent of the selection data in every respect: v13 includes three regions inside the selection span, v14–v17 change only region and chromosome, v18–v19 change sample within one trio, and v20 (whole HG002 chr20 at 15×) changes only the chromosome, with 1 Mb previously scored in v14. ΔF1 carries a 95 % paired locus-bootstrap interval (10,000 resamples). For v19 and v20 a 50-kb block bootstrap (Künsch, 1989) and, for v20, an exact McNemar test (McNemar, 1947) were added. [[c4]]Per cell (≥ 100 SNPs) and pooled, the pre-specified rule is:[[/c4]] PRESERVED if the interval contains 0 and |ΔF1| < 0.001; IMPROVED if the interval excludes 0 with ΔF1 > 0; DEGRADED if the interval excludes 0 with ΔF1 < 0, or |ΔF1| ≥ 0.001. The last clause, as written, also fires for positive differences (a disclosed defect). There is no equivalence test and no multiplicity correction. The condition for an overall CONFIRMED verdict cannot be met, so no global verdict is asserted."));

add(H2("3.3 V2 and V2.x protocols"));
add(P("V2 was checked against V1 on HG002 chr20 (300×) with a 3,491-locus gate of real loci, a full run on all 56,266,816 loci, and comparison of VCFs and of results for 1, 2, 4 and 8 threads. Timings come from one machine (Intel i7-6700, 4 cores/8 threads, 31 GiB, BAM on a 5400-rpm disk) with a desktop session running. For V2.x, HG002 chr20 was divided into 2-Mb blocks. The split was written before any feature table was inspected: even blocks are development, odd blocks are holdout. Features, coefficients, the threshold (maximising out-of-fold development F1) and the configuration were chosen on development blocks only. Each finalist was frozen and hashed before its holdout score was computed. Differences against V2 use a paired block bootstrap (whole 2-Mb blocks, 4,000 resamples). [[c7]]In an early exploratory script, holdout numbers of the B variants were printed next to development numbers, although B was selected on development F1; no model was refitted afterwards.[[/c7]] V2.x runtime is the median of three interleaved warm-cache runs with 8 threads."));
add(P("The final evaluation was pre-registered before any V2.x output existed for the data. It fixed sample (HG003), region (whole chr8, 145 Mb, not used by any V2.x analysis), model, command, scorer and stopping rule: one run, no tuning, and a failure would be reported, not repaired. It deliberately excluded a V1/V2 control, other variants and external callers. DeepVariant 1.10.0 (CPU, WGS model, 8 shards; Poplin *et al.*, 2018) and GATK 4.6.2.0 HaplotypeCaller (6 pair-HMM threads; McKenna *et al.*, 2010; Poplin *et al.*, 2017) were run on the same host and chr20 data. The BAM has no read group, so GATK used a copy with an added read group (about 20 min, not counted)."));

// ---------- Results ----------
add(H1("4 Results"));
add(H2("4.1 V1: routing versus PB-everywhere"));
add(P("The router sent 0.039–0.108 % of loci to PB in the regional experiments and 0.137 % on whole chr20. Pooled ΔF1 ranged from −5.6×10^−5^ to +2.7×10^−3^ (Table 2, Fig. 2a): two experiments were PRESERVED, five IMPROVED and one, whole chr20, DEGRADED. [[c5]]These are the verdicts of the pre-specified rule, with no global verdict.[[/c5]]"));
cur.__afterT2 = true;

// Table 2 and Fig 2 full width, inserted here
full(
  caption("**Table 2.** V1 cascade versus PB-only (internal evaluator, depths pooled; tier: V1 development/validation, see Sect. 3.2 for the relation of each experiment to the selection data). v13 pools 13 chr21 regions, three of them inside the selection span. Intervals are paired locus-bootstrap; v19 block-bootstrap CI [+0.52, +1.60]×10^−3^ and v20 block-bootstrap CI [−0.108, −0.0075]×10^−3^ agree. True SNPs lost are events counted per validation cell."),
  T2, afterTable("Values transcribed from the source manuscript (Table 4 there)."),
  ...figure(FIG("fig2_results.png"), FULL, "**Figure 2.** Main quantitative results. (a) V1: ΔF1 = cascade − PB-only per experiment (internal evaluator, 95 % locus-bootstrap interval). Grey band: pre-specified margin |ΔF1| < 0.001. Orange marks the experiment-level NEGATIVE (v14) and DEGRADED (v20) outcomes. [[c17]]Values from Table 2.[[/c17]] (b) V2.x: ΔF1 against the V2 control for each component on the odd (holdout) 2-Mb blocks of HG002 chr20, 95 % paired block-bootstrap interval (4,000 resamples). Panels use different data, evaluators and scales (×10^−3^) and must not be compared with each other. The figure shows effect sizes with uncertainty. It does not show that routing saves compute, that the V2.x gain transfers to other samples, or any ranking of callers."),
);

add(P("Two results carry a negative verdict, and both are reported at full size. In v14, one 1-Mb cell (HG002 chr20:33–34 Mb, full depth, 765 SNPs) had 26 false positives under PB-only and 28 under the cascade: ΔF1 = −0.001262, interval [−0.003849, +0.001228]. The interval contains zero, but the cell was DEGRADED through the |ΔF1| ≥ 0.001 clause, which made the experiment-level verdict NEGATIVE although the pooled result was IMPROVED. In v20 (70,324 truth SNPs scored), PB-only reached F1 0.97809 and the cascade 0.97803 (ΔF1 = −5.64×10^−5^; locus-bootstrap CI [−1.05, −0.14]×10^−4^). Ten loci disagreed: the cascade was wrong at nine and right at one (McNemar exact *p* = 0.021, resting on ten events), with 1 fewer true positive, 7 more false positives and 1 more false negative. hap.py gave the same direction and size (ΔF1 −5.5×10^−5^). The difference is 18-fold smaller than the 0.001 margin; it is detectable only because 70,324 SNPs give the cell high power. A random router at matched PB coverage lost 0.0186 F1, so the frozen router does real work."));
add(P("The IMPROVED verdicts are driven by loci in segmental duplications (Bailey *et al.*, 2002), where PB produced false positives from low-VAF reads that the binomial screen rejects (26 of 41 v14 departures, 36 of 38 in v17, 83 of 87 in v19). That reading is an interpretation: no read-level phasing was done. On chr20 this class supplies one departure of ten, and the sign of ΔF1 reverses. The pre-specified claim, that accuracy is not lost, is therefore supported only to within small margins, and it failed by a very small amount on one chromosome. HG003 and HG004 reproduce the regional pattern, but the trio was sequenced in one batch, which is weak evidence of independence."));

add(H2("4.2 V1 failure modes and compute"));
add(P("Ten true SNPs were lost across the validation cells (nine before chr20, one on chr20; per-cell counts, so the number of distinct loci is not established). They share VAF 0.11–0.13, high base quality, and an LLR~bin~ that lies 0.02–6.9 units below the router’s lower edge, so they were never routed. The cascade-introduced false positives show the opposite profile (VAF 0.15–0.60, mean base quality 9.6–26.5, LLR~bin~ just above the upper edge). That a fixed-ε binomial test is blind to both classes is an interpretation supported by an ε-shift analysis, not tested at read level. An ε-sensitivity safety layer that escalated flagged loci to PB fixed 5 calls and broke 69 on 9.5 M held-out loci and was rejected. At 300× the router engages on 707 of 56.2 M loci (0.001 %), so the V1 cascade there is close to binomial-only; its accuracy depends on the depth at which the constants were chosen."));
add(P("[[c6]]Evidence for computational savings in V1 is limited.[[/c6]] The benchmark extractor computes PB at every locus because the PB-only control needs it, and router-gated extraction was not implemented. The only measured BAM-to-calls gain from restricting PB to routed loci is 1.65× (one 0.5-Mb HG004 region, *n* = 3). Projected caller-stage speedups (171–396×) are projections, not results. Native C extraction sped single functions up 21–35×, but the pipeline gained 2.7× (Amdahl’s bound 3.06×) because PB then took 93–96 % of the time. Two statements from an earlier native-backend report (65× end-to-end; “17–80× faster” than three established callers) are withdrawn."));

add(H2("4.3 V2: exact re-implementation and runtime"));
add(P("All numeric kernels were bit-identical to the V1 NumPy/SciPy code on the tested grids. The 3,491-locus gate passed 17 of 17 field checks with no mismatch. A full run reproduced binomial and PB LLRs on all 56,266,816 loci, the same 707 routed loci and calls, and cascade and PB-only VCFs byte-identical to the archived V1 files, for 1, 2, 4 and 8 threads. Re-scoring the V2 control with hap.py reproduced V1 exactly (TP 70,012, FP 2,649, FN 1,321). Equivalence was verified on this one BAM only. V2 inherits every V1 accuracy limit, including the window filter."));
add(caption("**Table 3.** Runtime and memory (HG002 chr20, 300×; single runs unless stated). Tier: engineering benchmark. V1 is a frozen record measured under other conditions."), T3);
add(afterTable("Same-work comparison: both VCFs, PB at every locus (pb_all)."));
add(P("[[c18]]On HG002 chr20, the C++ engine processed the PB-everywhere workload in 1,059.5 s on one thread and 182.5 s on eight threads. The archived V1 run took 18,578 s, which corresponds to roughly 17.5-fold and 102-fold differences in elapsed time. Because V1 was measured under different execution and cache conditions, these ratios are indicative rather than controlled speed-up estimates (about ±15 %).[[/c18]] They reflect implementation and parallelism, not the cascade. The cascade row (119.6 s) also skips PB at most loci, an algorithmic saving, but it does different work and cannot produce the PB-only control. A single-thread profile attributed 52 % of time to htslib decoding, 22 % to counting and about 3 % to PB arithmetic. Extra BGZF threads slowed the engine and an AVX2+FMA build was no faster. Peak memory fell from gigabytes to 114–225 MB. Hardware counters were unavailable, so we make no claim about whether the engine is memory- or compute-bound."));

add(H2("4.4 V2.x: development and holdout"));
add(P("On the holdout blocks, A+C raised hap.py F1 from 0.97195 to 0.98902 (ΔF1 = +0.01706, 95 % CI +0.01183 to +0.02326). Precision rose from 0.9619 to 0.9964 and recall fell from 0.9822 to 0.9818: false positives fell from 1,323 to 122 and false negatives rose from 606 to 620 (Table 4, Fig. 2b). The transition/transversion ratio of the calls moved from 2.172 to 2.307 (truth: 2.315), a check that does not depend on matching calls to truth. The in-sample development F1 (0.98950) is optimistic by about 0.0005; the out-of-fold development estimate was 0.98899."));
add(caption("**Table 4.** V2.x accuracy by evidence tier (hap.py 0.3.15, SNP, PASS). “Mixed” means development and holdout blocks together, which is an engineering benchmark and not independent."), T4);
add(afterTable("HG003 chr8 has 186,522 truth SNPs; chr20 whole has 71,333. FP/FN counts across different data are not comparable."));
add(P("**A against A+C.** [[c16]]A single rule, a veto when the mean base quality of the ALT reads is ≤ 25.33, delivers about 92 % of the A+C gain[[/c16]] (holdout ΔF1 +0.01574 against +0.01706; precision 0.9933 against 0.9964; 225 against 122 false positives; 608 against 620 false negatives). The classifier on binomial outputs alone (C) gives +0.01631. The extra 0.0013 F1 of A+C over A is within the interval of either estimate, so this holdout does not separate them. The threshold reflects this data set’s quality distribution and would need refitting for other platforms or quality binning. B did not help: PB on more reads was monotonically worse under the frozen threshold (holdout F1 0.97195 with 48 reads, 0.97056 with 128, 0.96957 with all). D changed holdout F1 by −0.00016 to +0.00003. After the classifier, only 33 PB evaluations remain on chr20 (707 for V2). The extra work costs about 1.0 % CPU time (795 → 803 s) and +2.4 MB of memory; the wall-time difference (+0.5 %) lies inside the run-to-run spread of about ±5 %. Of the 620 holdout false negatives, 558 (90 %) are truth SNPs that the window frame cannot reach or that have fewer than three ALT reads, so the classifier is not what limits recall."));

add(H2("4.5 Final evaluation on HG003 chr8"));
add(P("In the single pre-registered run, A+C reached F1 0.989012 (block-bootstrap 95 % CI 0.98713–0.99048), precision 0.997567 and recall 0.980603 (TP 182,904; FP 446; FN 3,618) in 412.5 s at 326 MB (8 threads; BAM on a 5400-rpm disk, partly cached, so not comparable with the warm chr20 medians). This equals the chr20 holdout value to within 0.00001. [[c8]]No V2 control was run on HG003 chr8, so the size of the gain on this sample is unknown and no number is inferred for it.[[/c8]] Of the 3,618 false negatives, 2,873 (79 %) are truth SNPs outside all retained windows; they cap recall at 0.9846 by construction. The result shows that one frozen model, fitted on one sample and chromosome, gave a similar F1 on another sample and chromosome of the same platform family, aligner and read length. HG003 is the father of the development sample, so relatedness and shared technical properties could favour transfer. This is not an independent biological replication."));

add(H2("4.6 Comparison with established callers"));
add(P("[[c9]]On HG002 chr20 at 300×, DeepVariant and GATK scored higher than V2.x[[/c9]] (Table 5, Fig. 3a): F1 0.99761 and 0.99171, against 0.98927 for V2.x on the whole chromosome and 0.98902 on the holdout blocks. V2.x precision (0.9976) lies between GATK’s (0.9873) and DeepVariant’s (0.9994). Its recall (0.981) is below both (0.996), partly because the V1/V2 window frame cannot reach 1,041 of the 1,321 false negatives of the V2 control. In wall-clock time, V2.x completed in about 128 s, against about 1.5 h for DeepVariant and 3.7 h for GATK. These timings describe the tested configurations and do not establish a general ranking of caller efficiency (see Sect. 6)."));
cur.children[cur.children.length-1].__marker = "m2";
add(caption("**Table 5.** hap.py 0.3.15 (xcmp) SNP results, HG002 chr20, 300× (`TRUTH.TOTAL` = 71,333). Tier: external comparison on development data for V1/V2/V2.x. Descriptive; no ranking."), T5);
add(afterTable(""));

// Fig 1 and Fig 3 full width
// ---------- Discussion etc ----------
add(H1("5 Discussion"));
add(P("**V1.** Routing kept accuracy to within small margins, not exactly. The defensible statement is “preserved to within small margins in the evaluated GIAB/GRCh38 SNP-only setting”. The sign of ΔF1 is set by composition: where PB makes paralogous false positives the cascade gains, and where that class is nearly absent (chr20) the cascade has little to gain and its small losses are not offset. The IMPROVED verdicts are evidence about PB’s weakness in duplications, and the chr20 result is the cleaner test of “accuracy is not lost”. The router works where the binomial test is uncertain and fails where it is confidently wrong. The author’s view is that V1 traded a small amount of accuracy for lower computation, as the chr20 result shows, and that this accuracy may degrade at depths other than the one at which the constants were set. [[c25]]The measured compute gain from routing is limited to 1.65× in one region (Sect. 4.2), so the trade cannot be quantified from these data.[[/c25]]"));
add(P("**V2.** The speed-up comes from implementation and parallelism, not from changed science. The engine removed the read tensor, and PB arithmetic stopped being the bottleneck (about 3 % of the all-locus mode). The comparison with V1 rests on a record from other conditions. For orientation, the V1 pipeline (about 5 h 10 min on chr20) was slower and less accurate than both DeepVariant (1 h 28 min, F1 0.99761) and GATK (about 3 h 43 min, F1 0.99171) on the same data; the V2 engine changed the first of these differences, V2.x the second only in part."));
add(P("**V2.x.** Most of the V1/V2 error was outside the router’s reach, and cheap read-level evidence removed most of it. The gain is a precision gain, and most of it comes from a single base-quality rule. The mechanism is the reason for caution: low ALT-read base quality separates false positives (about 19) from true variants (about 36) in this data, and the classifier’s allele-fraction terms assume a diploid germline model. The model is simple and inspectable; it is not a general improvement of variant calling and should not be expected to hold for somatic or low-frequency variants, other quality binning or other sequencers without refitting. The equality of the chr8 and chr20 holdout F1 to five decimals is consistent with transfer, but a single run cannot say how much of that equality is chance. A V2 control on HG003 would have been needed to measure the gain itself."));
add(P("**External callers.** In the tested HG002 chr20 dataset, V2.x completed the analysis in 128 s with 221 MB of peak memory. This makes the implementation lightweight in this benchmark, although its accuracy remained below that of DeepVariant and GATK. The comparison is not like for like: DeepVariant and GATK also call indels, ran with different thread counts, containers and preprocessing, once each, on data that are the development chromosome of V2.x. Memory was not profiled for them."));

add(P("**Author’s assessment and outlook.** The author regards V2 as the central result of this work: it shows that the V1 decision logic can be executed exactly, and far faster, on modest hardware. Whether the small V1 loss on chr20 is an exception or a trend cannot be decided from one chromosome. It needs a larger set of samples, and the composition explanation in Sect. 4.1 remains a hypothesis. The main engineering lesson was that bottlenecks must be found by measurement. Native code was introduced stepwise, first for extraction in V1 and then for the whole pipeline in V2, each step guided by profiling, while micro-optimisation beyond the measured hot spots was not pursued. Equally important was verifying that every optimisation left the output unchanged (Sect. 4.3). The working hypothesis is that most of the calling work can be done cheaply with little loss of accuracy. The base-quality rule that recovers about 92 % of the V2.x gain is consistent with this, but the present data do not test it more broadly. [[c22]]A further idea, not implemented or tested here, is to use the router as a point at which other detectors, for example for specific mutation patterns, can be added or removed.[[/c22]] Planned work is to refine the formulas and the router, using low-level code (C, and assembly where profiling shows a benefit) to improve accuracy at no cost in speed. Further steps are to evaluate other genomes and chromosomes, to add the missing V2 control on HG003, to compare the classifier with the simple rule on more samples, and to compare the results with established callers to see where accuracy improves and where it does not."));

add(H1("6 Limitations"));
add(P("(1) All results are SNP-only, Illumina short reads, GRCh38 and GIAB high-confidence regions; difficult regions outside the BED were not evaluated, and truth is least reliable in segmental duplications, where cascade and PB differ most. (2) V2.x was fitted on one sample, chromosome, depth and platform. The final evaluation adds one related sample and one chromosome with the same platform family, aligner and read length, run once, without a control. Ancestry, other platforms, depths, somatic or low-VAF variants and polyploid samples were not tested. (3) The window filter bounds recall of every version (0.986 on chr20, 0.9846 on chr8). (4) V1 constants are depth-specific (selected at about 15×), and V1 intervals treat loci as independent outside v19 and v20. The v20 verdict rests on ten discordant loci. (5) V2 equivalence was verified on one BAM; the native counts backend differs from pysam on adversarial synthetic BAMs with repeated read names, and the native/reference gate on chr20 covers two chunks (1.7 % of loci). [[c14]]The original V1 tests were lost; the replacements are reconstructions.[[/c14]] (6) Timings come from one desktop machine, mostly single runs, with ±5 % page-cache variation. Disk space and hardware limited how many additional datasets could be processed. (7) [[c13]]HG005 results are record-only and omitted here.[[/c13]] (8) There is no environment lockfile, and the V2/V2.x code was committed after the analyses (intermediate states are identified by SHA-256 digests). (9) No clinical validation of any kind."));

add(H1("7 Conclusions"));
add(P("The results separate three contributions of AI_DNA_ANALYZER. The V1 experiments show that routing uncertain loci to PB does not guarantee preserved accuracy in every tested region: one chromosome was degraded by a very small amount (ΔF1 −5.6×10^−5^, ten discordant loci in 56.2 M), one experiment-level verdict was negative, and ten true SNPs were lost. The C++ implementation reproduces the archived V1 calls byte for byte on HG002 chr20 and runs far faster than the archived V1 record, although that record was measured under other conditions and routing itself was never shown to shorten a V1 run. The V2.x read-evidence layer raises holdout F1 from 0.97195 to 0.98902, mainly through a base-quality rule, and reaches a similar F1 (0.989012) on HG003 chr8. That evaluation has no V2 control, and both DeepVariant and GATK achieve higher F1 on the chr20 comparison. The next step is to test the contribution of the classifier against the simple base-quality rule on additional samples and sequencing datasets."));

add(H1("Data and code availability"));
add(P("[[c11]]Code, models, the final-evaluation artefacts and the analysis materials are in release v2.0.0[[/c11]] (2026-09-26): https://github.com/krmox/AI_DNA_ANALYZER/releases/tag/v2.0.0, archived at https://doi.org/10.5281/zenodo.22978775. The engine binary and model file are identified by SHA-256 digests (`adb2f4a0…bb7bdb`, `e9942a22…2d51ce`). Input data are public GIAB files; paths and checksums are in the repository. Regional BAMs are not redistributed. Software: hap.py 0.3.15 (container digest recorded), pysam 0.24.0/htslib 1.23.1, gcc 15.3.1 (C++20), Fedora 43.", { align: AlignmentType.LEFT }));
add(P("**Author contributions and AI assistance.** The author formulated the statistical models and hypotheses, designed the structure of the functions and the architecture, set the conditions of the experiments, read and analysed the results, and interpreted them using the literature, mathematical reasoning and discussion. [[c21]]Claude (Anthropic) acted as programmer: it implemented the code, ran the experiments and wrote the result records, and it drafted and edited the manuscript text under the author’s direction.[[/c21]] The author reviewed the results and takes responsibility for the conclusions. **Competing interests:** none declared. **Funding:** this research received no external funding.", { align: AlignmentType.LEFT }));

// ---------- References ----------
add(H1("References"));
const refs = [
  "Bailey,J.A. *et al.* (2002) Recent segmental duplications in the human genome. *Science*, **297**, 1003–1007.",
  "Bonfield,J.K. *et al.* (2021) HTSlib: C library for reading/writing high-throughput sequencing data. *GigaScience*, **10**, giab007.",
  "Cui,M. *et al.* (2024) miniSNV: accurate and fast single nucleotide variant calling from nanopore sequencing data. *Brief. Bioinform.*, **25**, bbae473.",
  "Danecek,P. *et al.* (2021) Twelve years of SAMtools and BCFtools. *GigaScience*, **10**, giab008.",
  "Gu,A. and Dao,T. (2023) Mamba: linear-time sequence modeling with selective state spaces. *arXiv*, 2312.00752. [[c23]](not machine-verified)[[/c23]]",
  "Hong,Y. (2013) On computing the distribution function for the Poisson binomial distribution. *Comput. Stat. Data Anal.*, **59**, 41–51.",
  "Krusche,P. *et al.* (2019) Best practices for benchmarking germline small-variant calls in human genomes. *Nat. Biotechnol.*, **37**, 555–560.",
  "Künsch,H.R. (1989) The jackknife and the bootstrap for general stationary observations. *Ann. Stat.*, **17**, 1217–1241.",
  "Liu,Y. *et al.* (2021) Psi-Caller: a lightweight short read-based variant caller with high speed and accuracy. *Front. Cell Dev. Biol.*, **9**, 731424.",
  "McKenna,A. *et al.* (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res.*, **20**, 1297–1303.",
  "McNemar,Q. (1947) Note on the sampling error of the difference between correlated proportions or percentages. *Psychometrika*, **12**, 153–157.",
  "Poplin,R. *et al.* (2017) Scaling accurate genetic variant discovery to tens of thousands of samples. *bioRxiv*, doi:10.1101/201178.",
  "Poplin,R. *et al.* (2018) A universal SNP and small-indel variant caller using deep neural networks. *Nat. Biotechnol.*, **36**, 983–987.",
  "Viola,P. and Jones,M. (2001) Rapid object detection using a boosted cascade of simple features. In *Proc. IEEE CVPR*, **1**, I-511–I-518.",
  "Wagner,J. *et al.* (2022) Benchmarking challenging small variants with linked and long reads. *Cell Genomics*, **2**, 100128.",
  "Wilm,A. *et al.* (2012) LoFreq: a sequence-quality aware, ultra-sensitive variant caller. *Nucleic Acids Res.*, **40**, 11189–11201.",
  "Zheng,Z. *et al.* (2022) Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. *Nat. Comput. Sci.*, **2**, 797–803. [[c10]](reference check note)[[/c10]]",
  "Zook,J.M. *et al.* (2016) Extensive sequencing of seven human genomes to characterize benchmark reference materials. *Sci. Data*, **3**, 160025.",
  "Zook,J.M. *et al.* (2019) An open resource for accurately benchmarking small variant and reference calls. *Nat. Biotechnol.*, **37**, 561–566.",
];
refs.forEach((r) => add(new Paragraph({ children: runs(r, base(REF)), spacing: { after: 20, line: 196 }, indent: { left: 170, hanging: 170 }, alignment: AlignmentType.LEFT })));

// Figures 1 and 3 are appended at positions via post-processing below
module.exports = { sections };

// Insert full-width figures: Fig 1 after Methods intro (after 2.4), Fig 3 after Section 4.6
function insertFull(afterMarker, blocks) {
  for (let i = 0; i < sections.length; i++) {
    const s = sections[i];
    if (s.cols !== 2) continue;
    const idx = s.children.findIndex((c) => c.__marker === afterMarker);
    if (idx >= 0) {
      const tail = s.children.splice(idx + 1);
      sections.splice(i + 1, 0, { cols: 1, children: blocks }, { cols: 2, children: tail });
      return true;
    }
  }
  return false;
}
const ok1 = insertFull("m1", figure(FIG("fig1_overview.png"), FULL, "**Figure 1.** Method overview and version relationships. V1 is the frozen cascade (binomial screen, router, PB model); the PB-only control, not drawn, evaluates PB at every locus. V2 re-implements the same decision logic as a C++ engine with byte-identical output on HG002 chr20. V2.x adds a logistic classifier at candidate loci (≥ 3 ALT reads) and is not part of the V1 protocol. The schematic shows data flow and frozen constants only; it carries no performance data."));
const ok2 = insertFull("m2", figure(FIG("fig3_pr_errors.png"), FULL, "**Figure 3.** Precision–recall and error analysis (hap.py 0.3.15, SNP, PASS). (a) HG002 chr20 at 300×; grey curves are F1 isolines. The whole-chromosome V2.x point (filled) includes blocks on which V2.x was fitted; the holdout point (open) does not. The HG003 chr8 point (star) is a different sample and chromosome and is shown for reference only. (b) False positives and false negatives for the V2 control and V2.x A+C on the chr20 holdout blocks and for A+C on HG003 chr8 (no control exists for chr8; counts across different data sizes are not comparable). On chr8, 2,873 of the 3,618 false negatives (79 %) lie outside all retained windows. The figure demonstrates that V2.x gains come from fewer false positives and that DeepVariant and GATK score higher on this data. It does not show an advantage of any caller outside the tested setting, nor runtime, memory or indel calling."));
if (!ok1 || !ok2) console.error("marker insertion failed", ok1, ok2);

// ---------- Document ----------
const footer = () => new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [
  new TextRun({ text: "AI_DNA_ANALYZER — working draft for author review — page ", font: FONT, size: 15, color: "555555" }),
  new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 15, color: "555555" }),
  new TextRun({ text: " of ", font: FONT, size: 15, color: "555555" }),
  new TextRun({ children: [PageNumber.TOTAL_PAGES], font: FONT, size: 15, color: "555555" })] })] });
const header = () => new Header({ children: [new Paragraph({ alignment: AlignmentType.RIGHT, children: [
  new TextRun({ text: "Draft — Bioinformatics Advances (Original Article)", font: FONT, size: 15, italics: true, color: "555555" })] })] });

const doc = new Document({
  creator: "Working draft (Claude, for author review)", title: "AI_DNA_ANALYZER working draft",
  comments: { children: Object.entries(C.comments).map(([id, text]) => ({ id: +id, author: "Editorial draft", initials: "ED", date: new Date("2026-10-25T12:00:00Z"),
    children: [new Paragraph({ children: [new TextRun({ text, font: FONT, size: 18 })] })] })) },
  styles: { default: { document: { run: { font: FONT, size: BODY } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true, run: { font: FONT, size: 21, bold: true }, paragraph: { outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true, run: { font: FONT, size: BODY, bold: true, italics: true }, paragraph: { outlineLevel: 1 } } ] },
  sections: sections.map((s, i) => ({
    properties: { type: i === 0 ? SectionType.NEXT_PAGE : SectionType.CONTINUOUS,
      page: { size: { width: PAGE_W, height: PAGE_H }, margin: { top: 1000, bottom: 900, left: MARGIN, right: MARGIN, header: 450, footer: 400 } },
      column: s.cols === 2 ? { count: 2, space: GAP, equalWidth: true } : { count: 1 } },
    headers: { default: header() }, footers: { default: footer() }, children: s.children })),
});
Packer.toBuffer(doc).then((b) => { fs.writeFileSync(OUT, b); console.log("written", OUT, sections.map(s => s.cols + ":" + s.children.length).join(" ")); });
