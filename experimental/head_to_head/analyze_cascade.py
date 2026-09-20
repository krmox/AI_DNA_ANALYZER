"""Head-to-head cascade analysis on the common HG002 chr21:32-44M region.

Reuses frozen cascade constants and the existing robustness_benchmark helpers
(route_mask, accuracy_block) exactly as bench_v12_stage1.py does. Nothing here
redefines a threshold or refits anything.
"""
import json
import sys

import numpy as np

import pathlib as _pl; _ROOT = str(_pl.Path(__file__).resolve().parents[2])  # project root (path-independent)
sys.path.insert(0, _ROOT)

from cascade import FROZEN_BINOMIAL_THRESHOLD, FROZEN_PB_THRESHOLD, FROZEN_ROUTER_CUTOFF
from config import LABEL_SNP
from robustness_benchmark import accuracy_block, route_mask
from residual_metrics import vaf_from_counts
from evaluate_quality_error import paired_bootstrap

blob = np.load(_ROOT + "/cache_h2h/hg002_32_44M_15x.npz", allow_pickle=True)
print("keys:", list(blob.keys()))

counts = blob["counts"]
labels = blob["labels"]
binomial_llr = blob["binomial_llr"].astype(np.float64)
pb_llr = blob["pb_llr"].astype(np.float64)
depth = blob["depth"].astype(np.float64) if "depth" in blob else counts[:, 9].astype(np.float64)
region = [int(x) for x in blob["region"]]

frame = (labels == 0) | (labels == LABEL_SNP)
snp = labels[frame] == LABEL_SNP
b_llr = binomial_llr[frame]
p_llr = pb_llr[frame]
d = depth[frame]

binomial_calls = b_llr >= FROZEN_BINOMIAL_THRESHOLD
pb_calls = p_llr >= FROZEN_PB_THRESHOLD
routed = route_mask(b_llr, FROZEN_ROUTER_CUTOFF)
cascade_calls = np.where(routed, pb_calls, binomial_calls)

result = {
    "region": region,
    "loci_total": int(frame.sum()),
    "snp_total": int(snp.sum()),
    "binomial_only": accuracy_block(binomial_calls, snp),
    "pb_only": accuracy_block(pb_calls, snp),
    "cascade": accuracy_block(cascade_calls, snp),
    "fraction_routed_to_pb": float(routed.mean()),
    "pb_compute_loci": int(routed.sum()),
    "disagreement_cascade_vs_pb": int(np.count_nonzero(cascade_calls != pb_calls)),
    "true_snp_losses_cascade_vs_pb": int(np.count_nonzero((cascade_calls == 0) & (pb_calls == 1) & snp)),
    "cascade_only_fp": int(np.count_nonzero((cascade_calls == 1) & (pb_calls == 0) & (~snp))),
    "pb_only_fp_not_cascade": int(np.count_nonzero((pb_calls == 1) & (cascade_calls == 0) & (~snp))),
}
result["delta_f1_cascade_vs_pb"] = result["cascade"]["f1"] - result["pb_only"]["f1"]
result["vs_pb_paired_bootstrap"] = paired_bootstrap(cascade_calls, pb_calls, snp)

# Rescue analysis
pb_correct = (pb_calls.astype(bool) == snp)
binomial_correct = (binomial_calls.astype(bool) == snp)
pb_rescuable = pb_correct & (~binomial_correct)
rescued = pb_rescuable & routed
missed_rescue = pb_rescuable & (~routed)
binomial_already_correct_and_routed = binomial_correct & routed
# false rescue: binomial correct, routed to PB, and PB does not change the outcome for the better
# (i.e. router sent it to PB unnecessarily since binomial was already right)
false_rescue = binomial_correct & routed

result["rescue"] = {
    "pb_rescuable_total": int(pb_rescuable.sum()),
    "rescued_by_router": int(rescued.sum()),
    "missed_by_router": int(missed_rescue.sum()),
    "rescue_recall": float(rescued.sum() / pb_rescuable.sum()) if pb_rescuable.sum() else None,
    "false_rescue_count": int(false_rescue.sum()),
    "false_rescue_rate_of_routed": float(false_rescue.sum() / routed.sum()) if routed.sum() else None,
    "total_routed": int(routed.sum()),
}

with open(_ROOT + "/experimental/head_to_head/ai_cascade_region_result.json", "w") as f:
    json.dump(result, f, indent=2, default=str)

print(json.dumps(result, indent=2, default=str))

# Save per-locus arrays for the disagreement/overlap analysis against external callers
np.savez(
    _ROOT + "/experimental/head_to_head/ai_cascade_per_locus.npz",
    counts=counts[frame], labels=labels[frame], binomial_llr=b_llr, pb_llr=p_llr,
    depth=d, binomial_calls=binomial_calls, pb_calls=pb_calls, cascade_calls=cascade_calls,
    routed=routed, snp=snp,
)
print("Saved per-locus arrays.")
