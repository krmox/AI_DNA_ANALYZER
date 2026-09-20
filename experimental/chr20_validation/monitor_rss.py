"""Sample total RSS of a process tree every 2 s until the root exits; write peak to JSON."""
import sys, time, json, psutil
root = psutil.Process(int(sys.argv[1])); peak = 0; n = 0; t0 = time.time()
while root.is_running() and root.status() != psutil.STATUS_ZOMBIE:
    try:
        peak = max(peak, sum(p.memory_info().rss for p in [root] + root.children(recursive=True))); n += 1
    except psutil.Error: pass
    time.sleep(2)
json.dump({"peak_rss_mb_total_process_tree": peak / 2**20, "samples": n, "wall_s_observed": time.time() - t0}, open(sys.argv[2], "w"), indent=1)
