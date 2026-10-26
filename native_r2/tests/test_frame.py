"""Frame equivalence: V2 window list and truth-derived labels == V1 cache, over ALL 56.27M loci."""
import os, sys, time
import numpy as np
W = os.environ.get("V2W", "/mnt/archive/AI_DNA_ANALYZER_v2/work") + "/setup"
CACHE = "/mnt/archive/AI_DNA_ANALYZER_benchmark/ai_dna_analyzer_run/cache/chr20_300x_counts.npz"
cz = np.load(CACHE)
pos = cz["positions"]; lab = cz["labels"]
ws = np.fromfile(f"{W}/windows_chr20_300x.i64", dtype=np.int64)
lb = np.fromfile(f"{W}/labels_chr20_300x.u8", dtype=np.uint8)
pos_v2 = (ws[:, None] + np.arange(64)[None, :]).ravel()
res = []
res.append(("window_positions_equal_all_loci", pos.size, int((pos_v2.size != pos.size) or (pos_v2 != pos).sum())))
res.append(("labels_equal_all_loci", lab.size, int((lb.astype(np.int64) != lab).sum())))
frame_v1 = ((lab == 0) | (lab == 1)).sum()
res.append(("frame_count_equal", int(frame_v1), int(frame_v1 != ((lb == 0) | (lb == 1)).sum())))
for name, n, bad in res:
    print(f"[{'PASS' if bad == 0 else 'FAIL'}] {name}: n={n} mismatches={bad}")
if os.environ.get("V2_CORRECTNESS_CSV"):
    with open(os.environ["V2_CORRECTNESS_CSV"], "a") as f:
        for name, n, bad in res:
            f.write(f"frame,{time.strftime('%F')},{name},{n},{bad},{'PASS' if bad == 0 else 'FAIL'},\n")
sys.exit(1 if any(b for _, _, b in res) else 0)
