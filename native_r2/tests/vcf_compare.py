"""Byte-level VCF comparison V2 vs frozen V1.

vcf_compare.py V2.vcf V1.vcf [gate_positions.txt]
With a positions file, V1 is restricted to the 64-bp windows containing those loci (the V2 gate run
only processes those windows).  Reports header equality, record counts, and every differing line.
"""
import os, sys, time
import numpy as np
v2, v1 = sys.argv[1], sys.argv[2]
def read(p):
    hdr, recs = [], []
    with open(p) as f:
        for line in f:
            (hdr if line.startswith("#") else recs).append(line)
    return hdr, recs
h2, r2 = read(v2); h1, r1 = read(v1)
scope = "full"
if len(sys.argv) > 3:
    pos = np.loadtxt(sys.argv[3], dtype=np.int64)
    wins = set((pos // 64 * 64).tolist())
    r1 = [l for l in r1 if ((int(l.split("\t")[1]) - 1) // 64 * 64) in wins]
    scope = f"windows of {len(pos)} gate loci"
same_header = h1 == h2
diff = [(a, b) for a, b in zip(r1, r2) if a != b]
n_bad = len(diff) + abs(len(r1) - len(r2)) + (0 if same_header else 1)
res = [("vcf_header_identical", 1, 0 if same_header else 1),
       ("vcf_record_lines_byte_identical", len(r1), len(diff)),
       ("vcf_record_count_equal", len(r1), abs(len(r1) - len(r2)))]
for name, n, bad in res:
    print(f"[{'PASS' if bad == 0 else 'FAIL'}] {name} ({scope}): n={n} mismatches={bad}")
for a, b in diff[:5]: print("V1:", a.strip(), "\nV2:", b.strip())
csv = os.environ.get("V2_CORRECTNESS_CSV")
if csv:
    with open(csv, "a") as f:
        for name, n, bad in res:
            f.write(f"{os.environ.get('V2_LABEL','vcf')},{time.strftime('%F')},{name},{n},{bad},{'PASS' if bad==0 else 'FAIL'},\"{scope}\"\n")
sys.exit(1 if n_bad else 0)
