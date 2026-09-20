"""ONE candidate/config once in a fresh process -> JSON. usage: bench_candidates.py CAND DATASET START END CHUNK_BP WORKERS  (env BENCH_MODE=serial|threads|procs)"""
import sys, os, time, json, hashlib, threading
from common import *
import numpy as np, psutil
import candidates as C
cand, ds, s, e, chunk_bp, workers = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
mode = os.environ.get("BENCH_MODE", "serial"); d = DATASETS[ds]; fn = getattr(C, cand)
peak = [0]; stop = threading.Event()
def sampler():
    me = psutil.Process()
    while not stop.is_set():
        try: peak[0] = max(peak[0], me.memory_info().rss + sum(c.memory_info().rss for c in me.children(recursive=True)))
        except Exception: pass
        time.sleep(0.2)
threading.Thread(target=sampler, daemon=True).start()
t0 = os.times(); w0 = time.perf_counter()
out = (C.run_serial(fn, d, (s, e), chunk_bp) if mode == "serial" else C.run_threads(fn, d, (s, e), chunk_bp, workers) if mode == "threads"
       else C.run_procs(cand, d, (s, e), chunk_bp, workers))
wall = time.perf_counter() - w0; t1 = os.times(); stop.set()
cpu = (t1.user - t0.user) + (t1.system - t0.system) + (t1.children_user - t0.children_user) + (t1.children_system - t0.children_system)
h = hashlib.sha256(); n = 0
for reads, labels in out: h.update(reads.tobytes()); h.update(labels.tobytes()); n += labels.size
print(json.dumps(dict(cand=cand, mode=mode, workers=workers, dataset=ds, region=[s, e], chunk_bp=chunk_bp, loci=n, wall_s=round(wall, 3),
      cpu_s=round(cpu, 2), cpu_util_cores=round(cpu / wall, 2), peak_rss_mb=round(peak[0] / 2**20), sha256=h.hexdigest())))
