"""The reconstructed V1 adapter (oracle/v1_oracle.py) must reproduce the archived V1 VCFs byte-for-byte
from the archived caches -- otherwise it is not a valid oracle."""
import filecmp, os, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "oracle"))
import v1_oracle as O  # noqa: E402

CACHE = "/mnt/archive/AI_DNA_ANALYZER_benchmark/ai_dna_analyzer_run/cache"
RUN = "/mnt/archive/AI_DNA_ANALYZER_benchmark/ai_dna_analyzer_run"
OUT = Path("/mnt/archive/AI_DNA_ANALYZER_v2/work/oracle"); OUT.mkdir(parents=True, exist_ok=True)
t0 = time.time()
cd = np.load(f"{CACHE}/chr20_300x_counts.npz"); counts, labels, positions = cd["counts"], cd["labels"], cd["positions"]
pb = np.load(f"{CACHE}/chr20_300x_pb_llr.npz")["pb_llr"]
binom = O.v1_binomial_llr(counts)
frame = (labels == 0) | (labels == 1)
pf, cf, bf, pbf = positions[frame], counts[frame], binom[frame], pb[frame]
B = bf >= 7.0; Pc = pbf >= 10.5
routed = np.abs(bf - 7.0) <= O.FROZEN_ROUTER_CUTOFF
casc = np.where(routed, Pc, B)
res = []
for mask, label, name in ((Pc, "pb_only", "ai_pb_only_chr20_300x.vcf"), (casc, "cascade", "ai_cascade_chr20_300x.vcf")):
    recs = O.extract_records(mask, cf, pf)
    gt, gq, forced = O.optimized_method_c(recs)
    O.write_vcf_c(recs, gt, gq, str(OUT / name), label)
    same = filecmp.cmp(OUT / name, f"{RUN}/{name}", shallow=False)
    res.append((f"oracle_reproduces_archived_{label}_vcf", len(recs), 0 if same else 1))
    print(f"[{'PASS' if same else 'FAIL'}] oracle_reproduces_archived_{label}_vcf: records={len(recs)} forced={int(forced.sum())}")
print("routed", int(routed.sum()), "time", round(time.time() - t0, 1))
if os.environ.get("V2_CORRECTNESS_CSV"):
    with open(os.environ["V2_CORRECTNESS_CSV"], "a") as f:
        for n, c, b in res: f.write(f"oracle,{time.strftime('%F')},{n},{c},{b},{'PASS' if b==0 else 'FAIL'},\n")
sys.exit(1 if any(b for *_, b in res) else 0)
