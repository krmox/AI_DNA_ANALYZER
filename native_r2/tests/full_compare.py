"""Full-chromosome comparison of a V2 `--mode pb_all --dump-llr-prefix` run against the V1 caches / VCFs.

full_compare.py RUN_DIR [label]     (RUN_DIR holds llr.binom.f64, llr.pb.f64, ai_*_chr20_300x.vcf)
Everything is compared at bit level over all 56,266,816 loci (windows in V1 cache order).
"""
import filecmp, os, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "oracle"))
import v1_oracle as O  # noqa: E402
run, label = Path(sys.argv[1]), (sys.argv[2] if len(sys.argv) > 2 else "full_pb_all")
CACHE = "/mnt/archive/AI_DNA_ANALYZER_benchmark/ai_dna_analyzer_run/cache"; RUN = "/mnt/archive/AI_DNA_ANALYZER_benchmark/ai_dna_analyzer_run"
t0 = time.time()
cd = np.load(f"{CACHE}/chr20_300x_counts.npz"); counts, labels = cd["counts"], cd["labels"]
pb1 = np.load(f"{CACHE}/chr20_300x_pb_llr.npz")["pb_llr"]
b1 = O.v1_binomial_llr(counts)
b2 = np.fromfile(run / "llr.binom.f64", dtype=np.float64); p2 = np.fromfile(run / "llr.pb.f64", dtype=np.float64)
u = lambda a: a.view(np.uint64)  # noqa: E731
fr = (labels == 0) | (labels == 1)
res = [("n_loci", b1.size, int(b1.size != b2.size or pb1.size != p2.size)),
       ("binomial_llr_bit_exact_all_loci", b1.size, int((u(b1) != u(b2)).sum())),
       ("pb_llr_bit_exact_all_loci", pb1.size, int((u(pb1) != u(p2)).sum())),
       ("router_mask_exact_frame", int(fr.sum()), int((((np.abs(b1 - 7) <= O.FROZEN_ROUTER_CUTOFF) != (np.abs(b2 - 7) <= O.FROZEN_ROUTER_CUTOFF)) & fr).sum())),
       ("binomial_call_exact_frame", int(fr.sum()), int(((b1 >= 7) != (b2 >= 7))[fr].sum())),
       ("pb_call_exact_frame", int(fr.sum()), int(((pb1 >= 10.5) != (p2 >= 10.5))[fr].sum()))]
routed = fr & (np.abs(b1 - 7) <= O.FROZEN_ROUTER_CUTOFF)
res.append(("routed_loci_count_equals_707", int(routed.sum()), int(routed.sum() != 707)))
casc1 = np.where(routed, pb1 >= 10.5, b1 >= 7); casc2 = np.where(routed, p2 >= 10.5, b2 >= 7)
res.append(("cascade_call_exact_frame", int(fr.sum()), int((casc1 != casc2)[fr].sum())))
for name in ("ai_cascade_chr20_300x.vcf", "ai_pb_only_chr20_300x.vcf"):
    same = filecmp.cmp(run / name, f"{RUN}/{name}", shallow=False)
    res.append((f"vcf_byte_identical_{name}", os.path.getsize(run / name), 0 if same else 1))
for name, n, bad in res: print(f"[{'PASS' if bad == 0 else 'FAIL'}] {name} ({label}): n={n} mismatches={bad}")
if (bad_idx := np.nonzero(u(pb1) != u(p2))[0]).size:
    print("first PB mismatches:", [(int(i), float(pb1[i]), float(p2[i])) for i in bad_idx[:10]])
print("elapsed", round(time.time() - t0, 1))
if os.environ.get("V2_CORRECTNESS_CSV"):
    with open(os.environ["V2_CORRECTNESS_CSV"], "a") as f:
        for name, n, bad in res: f.write(f"{label},{time.strftime('%F')},{name},{n},{bad},{'PASS' if bad == 0 else 'FAIL'},\n")
sys.exit(1 if any(b for *_, b in res) else 0)
