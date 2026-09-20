"""bench_v20 -- chr20 holdout validation of the FROZEN cascade (HG002, 15x, whole chromosome).

Nothing here fits, selects or modifies a threshold: the frozen constants are imported from
``cascade`` and only read. Scoring reuses, unchanged, the project's existing statistical
decision framework (``bench_v19_unseen.score_cell`` -> ``evaluate_arms`` paired bootstrap,
block bootstrap, McNemar, ``bench_v14_crosschrom.classify``: PRESERVED / DEGRADED / IMPROVED /
UNDERPOWERED). New in this file: (1) the M-1 BED/window audit, (2) the rescue composition split
into true-SNP rescue vs false-positive avoidance, (3) a disagreement taxonomy. Thresholds used
by (3) are fixed here, before the results were opened.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pysam

from bench_v13_robustness import disagreement_records, scoring_frame
from bench_v14_crosschrom import classify, load_cell
from bench_v19_unseen import enrich_disagreements, score_cell
from cascade import FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD, FROZEN_ROUTER_CUTOFF
from config import LABEL_SNP
from providers import GiabAlignmentProvider, load_bed_regions
from robustness_benchmark import evaluate_arms, route_mask

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger("bench_v20")

import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[0])  # project root (path-independent)
M = _ROOT
D = f"{M}/data/giab_hg002_chr20"
CACHE = "cache/chr20_validation/chr20_15x.npz"
OUT = Path("results/bench_v20")
CONTIG, CLEN, CHUNK_BP, SEQ_LEN = "chr20", 64_444_167, 500_000, 64
PREVIOUSLY_SCORED = (33_000_000, 34_000_000)   # chr20:33-34 Mb was a bench_v14 cell (frozen, not fit)

# Taxonomy thresholds -- fixed before opening any chr20 result.
TAX = {"low_vaf_lt": 0.15, "low_depth_lt": 10, "high_depth_ge": 30, "low_bq_lt": 25.0, "low_mapq_lt": 40.0}

SPEC = {"contig": CONTIG, "start": 0, "stop": CLEN, "base": "chr20", "depth_tag": "15x", "role": "holdout",
        "regime": "whole chromosome, GIAB high-confidence windows", "path": CACHE}


def subset(region: dict, mask: np.ndarray, name: str) -> dict:
    out = dict(region)
    for key in ("counts", "binomial_llr", "pb_llr", "labels", "depth", "mean_base_quality", "mean_mapq", "positions"):
        out[key] = region[key][mask]
    out["name"] = name
    return out


def truth_snp_positions() -> tuple[np.ndarray, dict]:
    """0-based positions of truth SNP bases, labelled exactly as the provider labels them."""
    vcf = pysam.VariantFile(f"{D}/chr20.vcf.gz")
    pos, n_records, n_hom_ref, n_other = [], 0, 0, 0
    for rec in vcf.fetch(CONTIG):
        n_records += 1
        if rec.alts is None:
            continue
        if GiabAlignmentProvider._is_homozygous_reference(rec):
            n_hom_ref += 1
            continue
        for alt in rec.alts:
            if alt is None or alt.startswith("<"):
                continue
            label, off, span = GiabAlignmentProvider._classify_allele(rec.ref or "", alt)
            if label == LABEL_SNP:
                pos.extend(rec.start + off + k for k in range(span))
            else:
                n_other += 1
    p = np.unique(np.asarray(pos, dtype=np.int64))
    return p, {"vcf_records": n_records, "hom_ref_records_skipped": n_hom_ref, "non_snp_alleles": n_other}


def m1_audit(positions: np.ndarray, labels: np.ndarray) -> dict:
    """Truth SNPs inside the GIAB BED vs those inside retained windows, with the drop reason."""
    truth, meta = truth_snp_positions()
    beds = load_bed_regions(f"{D}/chr20_highconf.bed", CONTIG)
    starts = np.asarray([s for s, _ in beds]); ends = np.asarray([e for _, e in beds])

    def in_bed(p):
        i = np.searchsorted(starts, p, side="right") - 1
        return (i >= 0) & (p < ends[np.maximum(i, 0)])

    truth_in_bed = truth[in_bed(truth)]
    retained_mask = np.isin(truth_in_bed, positions)
    retained = truth_in_bed[retained_mask]
    lab = labels[np.searchsorted(positions, retained)]
    retained_snp_label = int((lab == LABEL_SNP).sum())
    dropped = truth_in_bed[~retained_mask]

    fa = pysam.FastaFile(f"{M}/data/reference/chr20_full.fa")
    seq_name = fa.references[0]
    reasons = {"tile_remainder_at_chunk_end": 0, "window_partly_outside_bed": 0, "window_mostly_N": 0, "other": 0}
    for p in dropped:
        cs = (int(p) // CHUNK_BP) * CHUNK_BP
        ce = min(cs + CHUNK_BP, CLEN)
        w0 = cs + ((int(p) - cs) // SEQ_LEN) * SEQ_LEN
        if w0 + SEQ_LEN > ce:
            reasons["tile_remainder_at_chunk_end"] += 1
        else:
            j = int(np.searchsorted(starts, w0, side="right")) - 1
            inside = j >= 0 and ends[j] >= w0 + SEQ_LEN and starts[j] <= w0
            if not inside:
                reasons["window_partly_outside_bed"] += 1
            elif fa.fetch(seq_name, w0, w0 + SEQ_LEN).upper().count("N") > SEQ_LEN // 2:
                reasons["window_mostly_N"] += 1
            else:
                reasons["other"] += 1
    n_bed = int(truth_in_bed.size)
    return {
        "definition": "truth SNP bases (provider labelling: biallelic SNPs + MNPs decomposed; 0/0 records skipped)",
        "truth_snp_bases_chr20_total": int(truth.size),
        "N_BED": n_bed,
        "N_retained": int(retained.size),
        "N_retained_and_scored_as_SNP_in_frame": retained_snp_label,
        "N_retained_but_label_overwritten_by_indel": int(retained.size) - retained_snp_label,
        "N_dropped": int(dropped.size),
        "fraction_dropped": float(dropped.size / max(n_bed, 1)),
        "drop_reason_counts": reasons,
        "vcf_meta": meta,
        "bed_intervals": len(beds),
        "note": "hap.py TRUTH.TOTAL uses its own normalisation and will differ slightly from these base counts",
    }


def integrity_gates(region: dict) -> dict:
    """Reference-integrity and coordinate gates (the C-1 lesson): full population, not a sample."""
    fa = pysam.FastaFile(f"{M}/data/reference/chr20_full.fa")
    ref = np.frombuffer(fa.fetch(fa.references[0]).upper().encode(), dtype=np.uint8)
    lut = np.zeros(256, dtype=np.int64)
    for i, b in enumerate("NACGT"):
        lut[ord(b)] = i
    pos = region["positions"]
    mism = int((lut[ref[pos]] != region["counts"][:, 9].astype(np.int64)).sum())
    beds = load_bed_regions(f"{D}/chr20_highconf.bed", CONTIG)
    starts = np.asarray([s for s, _ in beds]); ends = np.asarray([e for _, e in beds])
    j = np.searchsorted(starts, pos, side="right") - 1
    in_bed = (j >= 0) & (pos < ends[np.maximum(j, 0)])
    return {"loci": int(pos.size), "reference_index_mismatches": mism,
            "positions_strictly_increasing": bool((np.diff(pos) > 0).all()),
            "positions_in_bed": bool(in_bed.all()), "positions_min": int(pos.min()), "positions_max": int(pos.max()),
            "all_windows_whole": bool(pos.size % SEQ_LEN == 0)}


def rescue_composition(snp, binomial_calls, pb_calls, routed) -> dict:
    """Router rescue, split by population. 'Rescuable' = PB is right where the binomial is wrong."""
    pb_ok = pb_calls == snp
    bn_ok = binomial_calls == snp
    resc = pb_ok & ~bn_ok
    tp_resc = resc & snp            # PB calls a true SNP the binomial missed
    fp_avoid = resc & ~snp          # binomial false positive that PB rejects
    harm = ~pb_ok & bn_ok           # PB wrong where the binomial is right

    def frac(a, b):
        return float(a / b) if b else None
    return {
        "total_pb_rescuable_loci": int(resc.sum()),
        "true_snp_rescuable_loci": int(tp_resc.sum()),
        "false_positive_avoidance_loci": int(fp_avoid.sum()),
        "true_snp_rescuable_routed": int((tp_resc & routed).sum()),
        "fp_avoidance_routed": int((fp_avoid & routed).sum()),
        "true_snp_rescue_recall (routed true-SNP-rescuable / true-SNP-rescuable)": frac((tp_resc & routed).sum(), tp_resc.sum()),
        "fp_avoidance_capture_rate (routed FP-avoidance / FP-avoidance)": frac((fp_avoid & routed).sum(), fp_avoid.sum()),
        "all_rescuable_capture_rate": frac((resc & routed).sum(), resc.sum()),
        "pb_harm_loci (PB wrong, binomial right)": int(harm.sum()),
        "pb_harm_routed (router sends these to PB: cascade inherits the PB error)": int((harm & routed).sum()),
        "false_rescue_rate (routed & binomial right & PB wrong / routed)": frac((routed & bn_ok & ~pb_ok).sum(), routed.sum()),
        "unnecessary_pb_rate (routed & binomial already right / routed)": frac((routed & bn_ok).sum(), routed.sum()),
        "n_routed": int(routed.sum()),
        "populations": "true-SNP rescue and FP avoidance are different populations; each rate names its denominator",
    }


def taxonomy(records: list[dict]) -> dict:
    counts: dict[str, dict] = {}
    for r in records:
        flags = ["low_vaf" if r["vaf"] < TAX["low_vaf_lt"] else "high_vaf",
                 "low_depth" if r["depth"] < TAX["low_depth_lt"] else "high_depth" if r["depth"] >= TAX["high_depth_ge"] else "mid_depth"]
        if r["mean_base_quality"] < TAX["low_bq_lt"]:
            flags.append("low_bq")
        if r.get("mean_mapq", 99) < TAX["low_mapq_lt"]:
            flags.append("low_mapq")
        flags += list(r.get("annotations", []))          # lowmap_segdup / tandemrepeats / alldifficult
        r["flags"] = flags
        r["outside_router_band"] = not r["routed_to_pb"]
        r["primary_class"] = ("segdup/lowmap" if "lowmap_segdup" in flags else "tandem_repeat" if "tandemrepeats" in flags
                              else "low_bq" if "low_bq" in flags else "low_mapq" if "low_mapq" in flags
                              else "low_vaf_other" if "low_vaf" in flags else "high_depth_other" if "high_depth" in flags else "other")
        b = counts.setdefault(r["effect"], {"n": 0, "flags": {}, "primary_class": {}})
        b["n"] += 1
        for f in flags:
            b["flags"][f] = b["flags"].get(f, 0) + 1
        b["primary_class"][r["primary_class"]] = b["primary_class"].get(r["primary_class"], 0) + 1
    fn = [r for r in records if r["effect"] == "router_fn"]
    mech = [r for r in fn if r["vaf"] < TAX["low_vaf_lt"] and "lowmap_segdup" in r["flags"] and r["outside_router_band"]]
    keys = ("position", "vaf", "depth", "mean_base_quality", "binomial_llr", "pb_llr", "router_margin", "flags")
    return {"thresholds": TAX, "by_effect": counts,
            "known_mechanism_router_fn": {
                "definition": "router_fn AND VAF<0.15 AND in lowmap_segdup AND unrouted (binomial margin outside cutoff)",
                "router_fn_events": len(fn), "matching_mechanism": len(mech),
                "router_fn_not_matching": [{k: r.get(k) for k in keys} for r in fn if r not in mech]}}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    region = load_cell("chr20_15x", SPEC)
    blob = np.load(CACHE, allow_pickle=True)
    region["positions"] = blob["positions"].astype(np.int64)
    log.info("loaded %d loci", region["labels"].size)

    gates = integrity_gates(region)
    log.info("integrity gates: %s", gates)
    assert gates["reference_index_mismatches"] == 0 and gates["positions_in_bed"] and gates["positions_strictly_increasing"]
    assert (FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD, FROZEN_ROUTER_CUTOFF) == (7.0, 10.5, 5.411872376933351)
    frozen = {"binomial_threshold": FROZEN_BINOMIAL_THRESHOLD, "pb_threshold": FROZEN_PB_THRESHOLD,
              "router_cutoff": FROZEN_ROUTER_CUTOFF}

    log.info("M-1 audit")
    m1 = m1_audit(region["positions"], region["labels"])

    log.info("score_cell (existing framework)")
    entry = score_cell(region, None)

    frame = scoring_frame(region)
    snp = region["labels"][frame] == LABEL_SNP
    binomial_llr = region["binomial_llr"][frame]
    pb_llr = region["pb_llr"][frame]
    routed = route_mask(binomial_llr, FROZEN_ROUTER_CUTOFF)
    binomial_calls = binomial_llr >= FROZEN_BINOMIAL_THRESHOLD
    pb_calls = pb_llr >= FROZEN_PB_THRESHOLD
    entry["rescue_composition"] = rescue_composition(snp, binomial_calls, pb_calls, routed)

    bb = entry["block_bootstrap"]
    lo, hi = bb["ci95"]
    contains0 = lo <= 0 <= hi
    entry["verdict_block_bootstrap_rule"] = ("UNDERPOWERED" if entry["snp"] < 100 or bb["ci_half_width"] > 0.01 else
                                              "PRESERVED" if contains0 and abs(bb["delta_f1"]) < 0.001 else
                                              "IMPROVED" if (not contains0 and bb["delta_f1"] > 0) else "DEGRADED")

    keep = ~((region["positions"] >= PREVIOUSLY_SCORED[0]) & (region["positions"] < PREVIOUSLY_SCORED[1]))
    rx = subset(region, keep, "chr20_15x_excl_v14_1Mb")
    ex = evaluate_arms(rx)
    ex["verdict"] = classify(ex)
    ex["excluded_range"] = list(PREVIOUSLY_SCORED)

    dis = disagreement_records(region)
    region["spec"] = SPEC
    dis = enrich_disagreements(region, dis)
    tax = taxonomy(dis["records"])

    report = {"frozen": frozen, "integrity_gates": gates, "m1_audit": m1, "cell": entry,
              "sensitivity_excluding_v14_chr20_33_34Mb": ex, "taxonomy": tax,
              "threshold_fit_span": "chr21:32-40 Mb (router cutoff); chr20 is a different chromosome"}
    (OUT / "chr20_results.json").write_text(json.dumps(report, indent=1, default=float))
    (OUT / "chr20_disagreements.json").write_text(json.dumps(dis, indent=1, default=float))
    log.info("wrote %s", OUT)
    log.info("verdict=%s delta_f1=%.5f pb_f1=%.5f cascade_f1=%.5f",
             entry["verdict"], entry["delta_f1_vs_pb"], entry["pb_only"]["f1"], entry["cheap_router_pb"]["f1"])


if __name__ == "__main__":
    main()
