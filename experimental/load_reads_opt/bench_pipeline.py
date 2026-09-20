"""E2E chunk-pipeline benchmark, one mode per process. usage: bench_pipeline.py MODE READS(python|native) START END WORKERS
MODE: serial | async | threads | procs (procs uses extract_bench_v12.process_chunk). env CHUNK_BP (default 100000).
reads=python means: pure-Python load_reads + original BED scan (the OLD code; load_counts stays native, as before)."""
import sys, os, time, json, hashlib, threading, asyncio, argparse
mode, reads_impl, s, e, workers = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
if reads_impl == "python": os.environ["AI_DNA_ANALYZER_DISABLE_NATIVE_READS"] = "1"
from common import *
import numpy as np, psutil
import read_level_pileup as R
import extract_bench_v12 as X
if reads_impl == "python": use_original_window_scan()
d = HG005
args = argparse.Namespace(fasta=d["fasta"], bam=d["bam"], vcf=d["vcf"], bed=d["bed"], contig=d["contig"], features="none",
                          legacy_pb=False, workers=workers, chunk_bp=int(os.environ.get("CHUNK_BP", 100_000)), region=[s, e])
chunks = [(c, min(c + args.chunk_bp, e)) for c in range(s, e, args.chunk_bp)]
peak = [0]; stop = threading.Event()
def sampler():
    me = psutil.Process()
    while not stop.is_set():
        try: peak[0] = max(peak[0], me.memory_info().rss + sum(c.memory_info().rss for c in me.children(recursive=True)))
        except Exception: pass
        time.sleep(0.2)
threading.Thread(target=sampler, daemon=True).start()
def stage_load(c0, c1):
    t = time.perf_counter(); counts, labels, positions = X.load_counts(args.fasta, args.bam, args.vcf, args.bed, (c0, c1), contig=args.contig, return_positions=True); tc = time.perf_counter() - t
    t = time.perf_counter(); reads, rl = X.load_reads(args.fasta, args.bam, args.vcf, args.bed, (c0, c1), contig=args.contig); tr = time.perf_counter() - t
    return counts, labels, positions, reads, dict(counts=tc, reads=tr)
def stage_pb(counts, labels, positions, reads):
    t = time.perf_counter(); pbd = X.poisson_binomial_diagnostics(reads, counts); tp = time.perf_counter() - t
    t = time.perf_counter(); binom = X.BinomialVariantCaller(error_rate=0.01).score_counts(counts[:, 0:4], counts[:, 9].astype(int)).llr; tb = time.perf_counter() - t
    return dict(pb_llr=pbd["pb_llr"], binomial_llr=binom, counts=counts.astype(np.float32), labels=labels, positions=positions), dict(pb=tp, binomial=tb)
stages = dict(counts=0.0, reads=0.0, pb=0.0, binomial=0.0)
def acc(t):
    for k, v in t.items(): stages[k] += v
t0 = os.times(); w0 = time.perf_counter(); outs = []
if mode == "serial":
    for c0, c1 in chunks:
        counts, labels, positions, reads, t = stage_load(c0, c1); acc(t)
        o, t = stage_pb(counts, labels, positions, reads); acc(t); outs.append(o)
elif mode == "async":
    async def run():
        loop = asyncio.get_running_loop()
        from concurrent.futures import ThreadPoolExecutor
        load_ex, pb_ex = ThreadPoolExecutor(1), ThreadPoolExecutor(1)
        nxt = loop.run_in_executor(load_ex, stage_load, *chunks[0])
        for i in range(len(chunks)):
            counts, labels, positions, reads, t = await nxt; acc(t)
            nxt = loop.run_in_executor(load_ex, stage_load, *chunks[i + 1]) if i + 1 < len(chunks) else None   # prefetch next chunk
            o, t = await loop.run_in_executor(pb_ex, stage_pb, counts, labels, positions, reads); acc(t); outs.append(o)
    asyncio.run(run())
elif mode == "threads":
    from concurrent.futures import ThreadPoolExecutor
    def full(c):
        counts, labels, positions, reads, t1 = stage_load(*c); o, t2 = stage_pb(counts, labels, positions, reads); return o, {**t1, **t2}
    with ThreadPoolExecutor(workers) as ex:
        for o, t in ex.map(full, chunks): acc(t); outs.append(o)
elif mode == "procs":
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    with ProcessPoolExecutor(workers, mp_context=mp.get_context("fork")) as ex:
        for r in ex.map(X.process_chunk, [args] * len(chunks), [c[0] for c in chunks], [c[1] for c in chunks]):
            for k in ("counts", "reads", "pb"): stages[k] += r["timing"][k + "_seconds"]
            b = r["blocks"]; outs.append(dict(pb_llr=b["pb_llr"], binomial_llr=b["binomial_llr"], counts=b["counts"], labels=b["labels"], positions=b["positions"]))
wall = time.perf_counter() - w0; t1 = os.times(); stop.set()
cpu = (t1.user - t0.user) + (t1.system - t0.system) + (t1.children_user - t0.children_user) + (t1.children_system - t0.children_system)
h = hashlib.sha256()
for o in outs:
    for k in ("pb_llr", "binomial_llr", "counts", "labels", "positions"): h.update(np.ascontiguousarray(o[k]).tobytes())
print(json.dumps(dict(mode=mode, reads=reads_impl, workers=workers, region=[s, e], chunks=len(chunks), loci=int(sum(o["labels"].size for o in outs)),
    wall_s=round(wall, 2), cpu_s=round(cpu, 1), cpu_util_cores=round(cpu / wall, 2), peak_rss_mb=int(peak[0] / 2**20),
    stage_seconds_summed={k: round(v, 1) for k, v in stages.items()}, sha256=h.hexdigest())))
