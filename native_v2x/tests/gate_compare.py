"""3,491-locus correctness gate: V2 engine dump vs frozen V1 (field-by-field, bit-level).

Usage: gate_compare.py DUMP.tsv [label]     (DUMP produced by `dnav2 run --loci-file ... --dump-loci`)
"""
import os, sys, struct, time
from pathlib import Path
import numpy as np
import pandas as pd

WORK = Path(os.environ.get("V2W", "/mnt/archive/AI_DNA_ANALYZER_v2/work")) / "gate"
dump, label = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else "gate")
g = np.load(WORK / "v1_gate.npz")
d = pd.read_csv(dump, sep="\t")
assert np.array_equal(d.pos0.values, g["pos0"]), "locus set/order differs"
n = len(d)
hexf = lambda col: np.array([float.fromhex(x) for x in d[col]], dtype=np.float64)  # noqa: E731
u64 = lambda a: np.asarray(a, np.float64).view(np.uint64)  # noqa: E731
res = []
def chk(name, bad, extra=""):
    res.append((name, n, int(bad), extra))

counts = g["counts"]
chk("counts_cols0-8_exact", (d[[f"c{i}" for i in range(9)]].values != counts[:, :9]).any(axis=1).sum())
chk("reference_index_exact", (d.ref_idx.values != counts[:, 9]).sum())
chk("label_exact", (d.label.values != g["labels"]).sum())
chk("frame_flag_exact", (d.frame.values.astype(bool) != g["frame"]).sum())
chk("binomial_llr_bit_exact", (u64(hexf("binomial_llr_hex")) != u64(g["binom_llr"])).sum())
chk("candidate_alt_exact", (d.alt.values != g["alt"]).sum())
chk("k_exact", (d.k.values != g["k"]).sum())
chk("n_exact", (d.n.values != g["n"]).sum())
chk("router_mask_exact", (d.routed.values.astype(bool) != g["routed"]).sum(), f"routed={int(g['routed'].sum())}")
chk("binomial_call_exact", (d.binomial_call.values.astype(bool) != g["bcall"]).sum())
# a locus with no admitted read is never emitted by the window pass; V1 then holds an all-padding row
# (base_code 255, all else 0) and PB = 0.0.  Require exactly that, and that counts also show depth 0.
pad = np.zeros(384, np.uint8); pad[0::8] = 255
unev = d.pb_evaluated.values == 0
chk("unevaluated_loci_are_exactly_V1_padding_rows_with_depth0",
    (~np.all(g["reads"].reshape(n, -1)[unev] == pad, axis=1)).sum() + (d.c6.values[unev] != 0).sum(),
    f"n_unevaluated={int(unev.sum())}")
pbv2 = hexf("pb_llr_hex")
chk("pb_llr_bit_exact_vs_v1_cache", (u64(pbv2) != u64(g["pb_cached"])).sum(),
    f"max_abs_diff={np.max(np.abs(pbv2 - g['pb_cached'])):.3g}")
chk("pb_call_exact", (d.pb_call.values.astype(bool) != g["pcall"]).sum())
chk("cascade_call_exact", (d.cascade_call.values.astype(bool) != g["ccall"])[g["frame"]].sum(), "frame loci only")
# V1 tensor rows byte-for-byte
rows_v1 = g["reads"].reshape(n, -1)
rows_v2 = np.tile(pad, (n, 1))          # V2 has no rows for unemitted loci -> padding by definition
have = d.rows_hex.values != "-"
for i in np.nonzero(have)[0]:
    rows_v2[i] = np.frombuffer(bytes.fromhex(d.rows_hex.values[i]), dtype=np.uint8)
chk("v1_tensor_rows_byte_exact(48x8 per locus)", (rows_v1 != rows_v2).any(axis=1).sum(),
    f"bytes compared={rows_v1.size}; loci with real reads={int(have.sum())}")
# tensor-derived consistency
ncnt_v1 = ((g["reads"][:, :, 6] & 2) > 0).sum(1)
chk("n_counted_tensor_exact", (d.n_counted_tensor.values != ncnt_v1).sum())
kt_v1 = (((g["reads"][:, :, 6] & 2) > 0) & (g["reads"][:, :, 0] == g["alt"][:, None])).sum(1)
chk("k_tensor_exact", (d.k_tensor.values != kt_v1).sum())

bad_total = sum(r[2] for r in res)
for name, nn, bad, extra in res:
    print(f"[{'PASS' if bad == 0 else 'FAIL'}] {name}: n={nn} mismatches={bad} {extra}")
csv = os.environ.get("V2_CORRECTNESS_CSV")
if csv:
    with open(csv, "a") as f:
        for name, nn, bad, extra in res:
            f.write(f"{label},{time.strftime('%F')},{name},{nn},{bad},{'PASS' if bad == 0 else 'FAIL'},\"{extra}\"\n")
sys.exit(1 if bad_total else 0)
