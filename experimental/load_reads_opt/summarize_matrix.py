import json, statistics as st, collections
rows = [json.loads(l) for l in open("results/matrix_load_reads.jsonl")]
g = collections.OrderedDict()
for r in rows: g.setdefault((r["dataset"], r["label"]), []).append(r)
base, shas = {}, {}
for (ds, lab), rs in g.items():
    if lab.startswith("baseline/serial"): base[ds] = st.median(x["wall_s"] for x in rs); shas[ds] = rs[0]["sha256"]
print("| dataset | config | n | median s | min | max | speedup vs baseline | CPU cores | peak MB | output sha == baseline |\n|---|---|--:|--:|--:|--:|--:|--:|--:|---|")
for (ds, lab), rs in g.items():
    w = [x["wall_s"] for x in rs]
    print(f"| {ds} | {lab} | {len(rs)} | {st.median(w):.2f} | {min(w):.2f} | {max(w):.2f} | {base[ds]/st.median(w):.2f}x | {st.median(x['cpu_util_cores'] for x in rs):.2f} | {int(max(x['peak_rss_mb'] for x in rs))} | {all(x['sha256']==shas[ds] for x in rs)} |")
