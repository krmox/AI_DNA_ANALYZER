"""Overlap / disagreement analysis: AI Cascade vs PB-only vs DeepVariant/Clair3/GATK
on the common HG002 chr21:32,000,000-44,000,000 region.

SNP-only, restricted to the GIAB v4.2.1 confident BED. Truth comes from the same
truth VCF used for the hap.py runs. This is a locus-level (not hap.py-graph-based)
comparison, explicitly secondary/simplified relative to the hap.py numbers already
produced for each external caller individually; it exists to let us compare AI
Cascade/PB directly against each external caller's calls at the same loci, which
hap.py alone (two VCFs at a time, no truth-free direct query-vs-query mode) does
not give us.
"""
import gzip
import subprocess
import numpy as np

REGION_CHROM = "chr21"
REGION_START = 32_000_000
REGION_END = 44_000_000

BASE = "/home/mark/Documents/Projects/AI_DNA_ANALYZER"
WORK = "/tmp/head_to_head_work/experimental/head_to_head"

TRUTH_VCF = f"{BASE}/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M.vcf.gz"
CONF_BED = f"{BASE}/data/giab_hg002_chr21_12Mb/hg002_chr21_32_44M_highconf.bed"

EXTERNAL_VCFS = {
    "DeepVariant": f"{WORK}/deepvariant/deepvariant.vcf.gz",
    "Clair3": f"{WORK}/clair3/merge_output.vcf.gz",
    "GATK_HC": f"{WORK}/gatk/gatk_hc.vcf.gz",
}

STRAT_BEDS = {
    "lowmap_segdup": f"{BASE}/data/strat_v31/lowmap_segdup.bed.gz",
    "tandem_repeat": f"{BASE}/data/strat_v31/tandemrepeats.bed.gz",
    "alldifficult": f"{BASE}/data/strat_v31/alldifficult.bed.gz",
}


def load_bed_intervals(path, chrom):
    intervals = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if parts[0] != chrom:
                continue
            s, e = int(parts[1]), int(parts[2])
            if e < REGION_START or s > REGION_END:
                continue
            intervals.append((s, e))
    return intervals


def make_mask(positions, intervals):
    mask = np.zeros(len(positions), dtype=bool)
    for s, e in intervals:
        mask |= (positions >= s) & (positions < e)
    return mask


def load_vcf_snp_positions(path, chrom, start, end):
    """Return a dict pos->(ref,alt,gt) for PASS/unfiltered SNP-only, biallelic records in region."""
    out = {}
    cmd = ["bcftools", "view", "-H", "-r", f"{chrom}:{start}-{end}", path]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        chrom_f, pos, _id, ref, alt = f[0], int(f[1]), f[2], f[3], f[4]
        filt = f[6]
        fmt_fields = f[8].split(":") if len(f) > 8 else []
        sample = f[9].split(":") if len(f) > 9 else []
        gt = None
        if "GT" in fmt_fields and sample:
            gt = sample[fmt_fields.index("GT")]
        if alt == ".":
            continue  # ref-only / no-call record
        alts = alt.split(",")
        if any(len(a) != 1 for a in alts) or len(ref) != 1:
            continue  # SNP-only
        if filt not in (".", "PASS", "RefCall"):
            # keep RefCall too (some callers use it for 0/0 with alt candidate) but exclude other hard filters
            if filt not in ("PASS", "."):
                continue
        if gt in ("0/0", "0|0", None):
            continue
        out[pos] = (ref, alts[0], gt)
    return out


print("Loading truth SNPs (already used for hap.py, reused here for direct locus comparison)...")
truth = load_vcf_snp_positions(TRUTH_VCF, REGION_CHROM, REGION_START, REGION_END)
print(f"truth SNP loci in region: {len(truth)}")

conf_intervals = load_bed_intervals(CONF_BED, REGION_CHROM)


def in_conf(pos):
    # simple linear scan would be too slow; build array mask once instead
    return None


conf_positions_mask_builder = conf_intervals  # used below via make_mask on arrays

external_calls = {}
for name, path in EXTERNAL_VCFS.items():
    external_calls[name] = load_vcf_snp_positions(path, REGION_CHROM, REGION_START, REGION_END)
    print(f"{name}: {len(external_calls[name])} SNP calls with non-ref GT in region")

# Restrict truth and calls to confident BED
def restrict_to_bed(d, intervals):
    positions = np.array(sorted(d.keys())) if d else np.array([], dtype=np.int64)
    if len(positions) == 0:
        return {}
    mask = make_mask(positions, intervals)
    keep = set(positions[mask].tolist())
    return {p: v for p, v in d.items() if p in keep}


truth_conf = restrict_to_bed(truth, conf_intervals)
external_conf = {name: restrict_to_bed(d, conf_intervals) for name, d in external_calls.items()}
print(f"truth SNP loci in confident BED: {len(truth_conf)}")
for name, d in external_conf.items():
    print(f"{name} SNP calls in confident BED: {len(d)}")

# Load AI cascade / PB-only per-locus arrays (already confident-BED restricted upstream
# by the frozen extractor, which only scores loci inside the highconf BED)
ai = dict(np.load(f"{WORK}/ai_cascade_per_locus.npz", allow_pickle=True))
ai_positions = ai["positions"]
ai_snp_truth = ai["snp"].astype(bool)
ai_pb_call = ai["pb_calls"].astype(bool)
ai_cascade_call = ai["cascade_calls"].astype(bool)
ai_binomial_call = ai["binomial_calls"].astype(bool)

pb_call_set = set(ai_positions[ai_pb_call].tolist())
cascade_call_set = set(ai_positions[ai_cascade_call].tolist())
truth_set = set(ai_positions[ai_snp_truth].tolist())

print(f"AI-frame loci: {len(ai_positions)}; AI truth SNP count in-frame: {len(truth_set)}")
print(f"AI truth SNP count in-frame == hap.py-style truth_conf count check: {len(truth_set)} vs {len(truth_conf)}")

results = {}
callers = {"PB-only": pb_call_set, "AI-Cascade": cascade_call_set, **{k: set(v.keys()) for k, v in external_conf.items()}}

for name, call_set in callers.items():
    tp = len(call_set & truth_set)
    fp = len(call_set - truth_set)
    fn = len(truth_set - call_set)
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else float("nan")
    results[name] = dict(tp=tp, fp=fp, fn=fn, precision=precision, recall=recall, f1=f1)
    print(name, results[name])

import json
with open(f"{WORK}/locus_level_secondary_evaluator_results.json", "w") as f:
    json.dump(results, f, indent=2, default=str)

# Pairwise overlap: AI-Cascade vs each external caller, and PB-only vs each external caller
pairwise = {}
for base_name, base_set in [("AI-Cascade", cascade_call_set), ("PB-only", pb_call_set)]:
    for ext_name, ext_set in external_conf.items():
        ext_call_set = set(ext_set.keys())
        shared_tp = len((base_set & ext_call_set) & truth_set)
        base_only_tp = len((base_set - ext_call_set) & truth_set)
        ext_only_tp = len((ext_call_set - base_set) & truth_set)
        base_only_fp = len((base_set - ext_call_set) - truth_set)
        ext_only_fp = len((ext_call_set - base_set) - truth_set)
        key = f"{base_name}_vs_{ext_name}"
        pairwise[key] = dict(
            shared_tp=shared_tp, base_only_tp=base_only_tp, ext_only_tp=ext_only_tp,
            base_only_fp=base_only_fp, ext_only_fp=ext_only_fp,
        )
        print(key, pairwise[key])

with open(f"{WORK}/pairwise_overlap_results.json", "w") as f:
    json.dump(pairwise, f, indent=2, default=str)

# Disagreement table: AI-Cascade vs each external caller, at loci where they differ
strat_masks = {}
for sname, spath in STRAT_BEDS.items():
    ivals = load_bed_intervals(spath, REGION_CHROM)
    strat_masks[sname] = make_mask(ai_positions, ivals)

pos_to_idx = {p: i for i, p in enumerate(ai_positions)}
depth = ai["depth"]
counts = ai["counts"]
binomial_llr = ai["binomial_llr"]
pb_llr = ai["pb_llr"]
routed = ai["routed"].astype(bool)


def vaf_at(idx):
    c = counts[idx]
    # counts layout: [A,C,G,T, ...] first 4 are base counts typically; use depth-based proxy
    d = depth[idx]
    if d <= 0:
        return 0.0
    alt_count = max(c[0], c[1], c[2], c[3]) if len(c) >= 4 else 0
    total = c[0] + c[1] + c[2] + c[3] if len(c) >= 4 else d
    ref_like = max(c[:4]) if len(c) >= 4 else 0
    non_ref = total - ref_like if total else 0
    return float(non_ref / total) if total else 0.0


rows = []
for ext_name, ext_set in external_conf.items():
    ext_call_set = set(ext_set.keys())
    disagree_positions = (cascade_call_set ^ ext_call_set) & set(ai_positions.tolist())
    for p in sorted(disagree_positions):
        idx = pos_to_idx.get(p)
        if idx is None:
            continue
        truth_here = p in truth_set
        ai_call = p in cascade_call_set
        ext_call = p in ext_call_set
        pb_call = p in pb_call_set
        rows.append(dict(
            competitor=ext_name, chrom=REGION_CHROM, position=p,
            truth=truth_here, AI_call=ai_call, competitor_call=ext_call, PB_call=pb_call,
            depth=float(depth[idx]), vaf=vaf_at(idx),
            binomial_LLR=float(binomial_llr[idx]), PB_LLR=float(pb_llr[idx]),
            router_margin=float(abs(binomial_llr[idx] - 5.411872376933351)) if False else float(binomial_llr[idx] - 7.0),
            routed=bool(routed[idx]),
            lowmap_segdup=bool(strat_masks["lowmap_segdup"][idx]),
            tandem_repeat=bool(strat_masks["tandem_repeat"][idx]),
            alldifficult=bool(strat_masks["alldifficult"][idx]),
        ))

import csv
with open(f"{WORK}/ERROR_ANALYSIS.csv", "w", newline="") as f:
    if rows:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    else:
        f.write("no disagreements found\n")

print(f"Wrote {len(rows)} disagreement rows to ERROR_ANALYSIS.csv")
