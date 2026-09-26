#!/usr/bin/env python3
"""Machine check of the new (V2 / V2.x / final evaluation / chr20-300x) numbers of PAPER_MANUSCRIPT_V2.md against the
existing artefacts.  Every entry: (what, value as printed in the manuscript, value computed from an artefact, tolerance).
Run from research/:  python3 audit_manuscript_numbers.py     Exit status 1 if any check fails."""
import csv, json, statistics, sys, os
R = os.path.dirname(os.path.abspath(__file__))
fails = 0; n = 0
def chk(what, printed, actual, tol=0.5e-6):
    global fails, n
    n += 1
    ok = abs(float(printed) - float(actual)) <= tol
    if not ok:
        fails += 1
    print(("PASS " if ok else "FAIL ") + f"{what}: manuscript {printed} vs artefact {actual}")
def rows(p): return list(csv.DictReader(open(os.path.join(R, p))))
acc = rows("V2X_ACCURACY_PERFORMANCE_BENCHMARK.csv")
def A(model, tag):
    for r in acc:
        if r["model"] == model and r["development_or_holdout"].startswith(tag): return r
    raise KeyError((model, tag))
# ---- Table 11 (holdout)
T11 = {"V2_CONTROL": (0.97195, 0.96193, 0.98219, 1323, 606), "V2.x_B": (0.97359, 0.96517, 0.98216, 1206, 607),
       "V2.x_D": (0.97572, 0.97416, 0.97728, 882, 773), "V2.x_A": (0.98769, 0.99331, 0.98213, 225, 608),
       "V2.x_C": (0.98826, 0.99544, 0.98119, 153, 640), "V2.x_A+C": (0.98902, 0.99636, 0.98178, 122, 620),
       "V2.x_A+B+C+D": (0.98890, 0.99586, 0.98204, 139, 611)}
for m, (f1, p, r, fp, fn) in T11.items():
    a = A(m, "holdout")
    chk(f"T11 {m} F1", f1, a["F1"], 0.5e-5); chk(f"T11 {m} P", p, a["precision"], 0.5e-5); chk(f"T11 {m} R", r, a["recall"], 0.5e-5)
    chk(f"T11 {m} FP", fp, a["FP"], 0); chk(f"T11 {m} FN", fn, a["FN"], 0)
a = A("V2.x_A+C", "holdout")
chk("A+C holdout F1 (6 dp)", 0.989017, a["F1"]); chk("A+C holdout P (6 dp)", 0.996362, a["precision"]); chk("A+C holdout R (6 dp)", 0.981781, a["recall"])
chk("A+C holdout TP", 33410, a["TP"], 0)
c = A("V2_CONTROL", "holdout"); chk("control holdout F1 (6 dp)", 0.971953, c["F1"]); chk("control holdout TP", 33424, c["TP"], 0)
full = A("V2.x_A+C", "full"); chk("A+C full F1", 0.989271, full["F1"]); chk("A+C full P", 0.997563, full["precision"]); chk("A+C full R", 0.981117, full["recall"])
chk("A+C full TP", 69986, full["TP"], 0); chk("A+C full FP", 171, full["FP"], 0); chk("A+C full FN", 1347, full["FN"], 0)
cf = A("V2_CONTROL", "full"); chk("control full F1", 0.972429, cf["F1"]); chk("control full TP", 70012, cf["TP"], 0); chk("control full FP", 2649, cf["FP"], 0); chk("control full FN", 1321, cf["FN"], 0)
dv = A("V2.x_A+C", "dev"); chk("A+C dev F1", 0.989503, dv["F1"]); chk("A+C dev TP", 36576, dv["TP"], 0); chk("A+C dev FP", 49, dv["FP"], 0); chk("A+C dev FN", 727, dv["FN"], 0)
chk("A+C dev P", 0.998662, dv["precision"]); chk("A+C dev R", 0.980511, dv["recall"])
chk("Ti/Tv control", 2.172, float(c["TiTv"]), 0.0005); chk("Ti/Tv A+C", 2.307, float(a["TiTv"]), 0.0005)
# ---- bootstrap deltas
bd = {(r["model"], r["split"]): r for r in rows("V2X_BOOTSTRAP_DELTAS_VS_V2.csv")}
for m, (d, lo, hi) in {"A+C": (0.01706, 0.01183, 0.02326), "A": (0.01574, 0.01043, 0.02185), "C": (0.01631, 0.01111, 0.02248),
                       "B": (0.00164, 0.00105, 0.00233), "D": (0.00377, 0.00162, 0.00626), "A+B+C+D": (0.01695, 0.01176, 0.02307)}.items():
    r = bd[(m, "holdout")]; chk(f"dF1 holdout {m}", d, r["dF1"], 0.5e-5); chk(f"dF1 lo {m}", lo, r["dF1_lo"], 0.5e-5); chk(f"dF1 hi {m}", hi, r["dF1_hi"], 0.5e-5)
r = bd[("A+C", "full")]; chk("dF1 full A+C", 0.01684, r["dF1"], 0.5e-5); chk("dF1 full lo", 0.01345, r["dF1_lo"], 0.5e-5); chk("dF1 full hi", 0.02060, r["dF1_hi"], 0.5e-5)
chk("dFP holdout A+C", -1201, float(A("V2.x_A+C","holdout")["FP"]) - float(c["FP"]), 0)
# ---- Table 12 (timing): medians of the three warm8 reps
tim = rows("V2X_TIMING_RESULTS.csv")
def med(model, key):
    xs = [float(r[key]) for r in tim if r["model"] == model and r["plan"] == "warm8"]
    return statistics.median(xs), min(xs), max(xs), len(xs)
for m, (w, lo, hi, cpu, rss, pb) in {"V2_CONTROL": (127.3, 118.4, 129.9, 795, 219, 707), "A": (124.7, 119.3, 130.5, 804, 221, 707), "C": (130.5, 119.5, 131.7, 801, 221, 33),
                                    "A+C": (128.0, 119.3, 129.7, 803, 221, 33), "A+B+C+D": (128.0, 120.4, 131.8, 803, 222, 36)}.items():
    md, mn, mx, k = med(m, "wall_s"); chk(f"T12 {m} wall median (n={k})", w, md, 0.06); chk(f"T12 {m} wall min", lo, mn, 0.06); chk(f"T12 {m} wall max", hi, mx, 0.06)
    cpus = [float(r["user_s"]) + float(r["sys_s"]) for r in tim if r["model"] == m and r["plan"] == "warm8"]; chk(f"T12 {m} CPU", cpu, statistics.median(cpus), 1.0)
    chk(f"T12 {m} RSS", rss, statistics.median([float(r["peak_rss_mb"]) for r in tim if r["model"] == m and r["plan"] == "warm8"]), 1.0)
    chk(f"T12 {m} PB evals", pb, [int(r["n_pb_evaluated"]) for r in tim if r["model"] == m and r["plan"] == "warm8"][0], 0)
# ---- chr8 final evaluation
J = json.load(open(os.path.join(R, "final_independent_eval/results_chr8_analysis.json")))
s = J["happy_summary_SNP"]["PASS"]
chk("chr8 TRUTH.TOTAL", 186522, s["TRUTH.TOTAL"], 0); chk("chr8 TP", 182904, s["TRUTH.TP"], 0); chk("chr8 FN", 3618, s["TRUTH.FN"], 0); chk("chr8 FP", 446, s["QUERY.FP"], 0)
chk("chr8 QUERY.TOTAL", 183350, s["QUERY.TOTAL"], 0); chk("chr8 precision", 0.997567, s["METRIC.Precision"]); chk("chr8 recall", 0.980603, s["METRIC.Recall"]); chk("chr8 F1", 0.989012, s["METRIC.F1_Score"])
chk("chr8 FP.gt", 91, s["FP.gt"], 0); chk("chr8 FP.al", 65, s["FP.al"], 0)
chk("chr8 TiTv query", 1.9322, s["QUERY.TOTAL.TiTv_ratio"], 5e-5); chk("chr8 TiTv truth", 1.9302, s["TRUTH.TOTAL.TiTv_ratio"], 5e-5)
chk("chr8 het/hom query", 1.7008, s["QUERY.TOTAL.het_hom_ratio"], 5e-5); chk("chr8 het/hom truth", 1.6968, s["TRUTH.TOTAL.het_hom_ratio"], 5e-5)
ci = J["block_bootstrap_95ci"]
for k, v in {"f1": (0.98713, 0.99048), "precision": (0.99655, 0.99828), "recall": (0.97770, 0.98299)}.items():
    chk(f"chr8 CI {k} lo", v[0], ci[k][0], 5e-6); chk(f"chr8 CI {k} hi", v[1], ci[k][1], 5e-6)
chk("chr8 blocks", 73, J["n_blocks"], 0)
E = json.load(open(os.path.join(R, "final_independent_eval/run/stdout.json")))
chk("chr8 n_loci", 133987200, E["n_loci"], 0); chk("chr8 n_frame_scored", 133935319, E["n_frame_scored"], 0); chk("chr8 windows", 2093550, E["n_windows"], 0)
chk("chr8 candidates", 338184, E["n_candidates"], 0); chk("chr8 routed", 1018, E["n_routed_pb"], 0); chk("chr8 pb evals", 102, E["n_pb_evaluated"], 0)
chk("chr8 calls", 183350, E["cascade"]["n_calls"], 0); chk("chr8 forced 0/0->0/1", 219, E["cascade"]["n_forced_00_to_01"], 0)
chk("chr8 wall", 412.5, E["wall_seconds"], 0.05); chk("chr8 CPU", 1930.2, E["cpu_user_seconds"] + E["cpu_sys_seconds"], 0.05)
chk("chr8 RSS kB", 333900, E["peak_rss_kb"], 0); chk("chr8 RSS MB", 326, E["peak_rss_kb"] / 1024, 0.5); chk("chr8 loci/s", 324800, E["n_loci"] / E["wall_seconds"], 50)
chk("chr8 candidates % scored", 0.25, 100 * E["n_candidates"] / E["n_frame_scored"], 0.005)
F = json.load(open(os.path.join(R, "final_independent_eval/results_chr8_frame_audit.json")))
chk("frame audit windows match", 1, int(F["reconstruction_matches_engine_window_count"]), 0); chk("truth outside windows", 2873, F["truth_snp_outside_retained_windows"], 0)
chk("outside % of truth", 1.54, 100 * F["fraction_outside"], 0.005); chk("FN outside share", 79.4, 100 * F["fraction_of_fn_outside"], 0.05); chk("recall bound", 0.98460, F["recall_upper_bound_from_window_filter"], 5e-6)
chk("chr8 FN in frame", 745, F["fn_total"] - F["fn_outside_retained_windows"], 0)
cov = open(os.path.join(R, "final_independent_eval/logs/samtools_coverage_chr8.txt")).read().split("\n")[1].split("\t")
chk("chr8 mean depth", 302.7, cov[6], 0.05); chk("chr8 covered %", 99.2, cov[5], 0.05); chk("chr8 mean BQ", 35.5, cov[7], 0.05); chk("chr8 mean MAPQ", 69.5, cov[8], 0.05); chk("chr8 reads", 297503843, cov[3], 0)
# ---- chr20 300x external (Table 16) from hap.py summaries
def happy(path, typ="SNP"):
    for r in csv.DictReader(open(os.path.join(R, path))):
        if r["Type"] == typ and r["Filter"] == "PASS": return r
B = "benchmark_chr20_300x/happy/"
for name, (p, vals) in {"V1 cascade": (B + "ai_dna_analyzer/ai_cascade.summary.csv", (0.963543, 0.981481, 0.972429, 70012, 2649, 1321)),
        "PB-only": (B + "ai_pb_only/ai_pb_only.summary.csv", (0.960849, 0.981579, 0.971104, 70019, 2853, 1314)),
        "DeepVariant": (B + "deepvariant/deepvariant.summary.csv", (0.999395, 0.995836, 0.997613, 71036, 43, 297)),
        "GATK": (B + "gatk/gatk.summary.csv", (0.987347, 0.996117, 0.991712, 71056, 911, 277))}.items():
    r = happy(p); pr, re_, f1, tp, fp, fn = vals
    chk(f"T16 {name} P", pr, r["METRIC.Precision"]); chk(f"T16 {name} R", re_, r["METRIC.Recall"]); chk(f"T16 {name} F1", f1, r["METRIC.F1_Score"])
    chk(f"T16 {name} TP", tp, r["TRUTH.TP"], 0); chk(f"T16 {name} FP", fp, r["QUERY.FP"], 0); chk(f"T16 {name} FN", fn, r["TRUTH.FN"], 0); chk(f"T16 {name} TRUTH.TOTAL", 71333, r["TRUTH.TOTAL"], 0)
chk("DeepVariant SNP QUERY.UNK", 19234, happy(B + "deepvariant/deepvariant.summary.csv")["QUERY.UNK"], 0); chk("GATK SNP QUERY.UNK", 42700, happy(B + "gatk/gatk.summary.csv")["QUERY.UNK"], 0)
p20 = "/mnt/archive/AI_DNA_ANALYZER_v2x/work/happy/M_A_C/full/ai_M_A_C.summary.csv"   # raw hap.py output of the V2.x run (outside the repository)
if os.path.exists(p20): chk("chr20 truth Ti/Tv (hap.py)", 2.315, happy(p20)["TRUTH.TOTAL.TiTv_ratio"], 0.0006)
print(f"\n{n - fails}/{n} checks passed"); sys.exit(1 if fails else 0)
