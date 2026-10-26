#!/usr/bin/env python3
"""Re-run hap.py for every condition/caller whose output is missing or failed (idempotent)."""
import json, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import stage3 as s3
import run_all

def todo():
    for d in sorted((s3.S3 / "runs").iterdir()):
        if not (d / "calls.json").exists() or not all((d / "calls.json").read_text().count(k) for k in ('"varion"', '"bcftools8"')):
            continue
        sample, contig = d.name.split("_")[:2]
        for caller in ("varion", "bcftools"):
            f = d / f"happy_{caller}" / "happy_run.json"
            if not f.exists() or json.load(open(f))["rc"] != 0:
                yield sample, contig, d.name, caller

jobs = list(todo())
print(len(jobs), "hap.py jobs to (re)run", flush=True)
with ThreadPoolExecutor(max_workers=int(sys.argv[1]) if len(sys.argv) > 1 else 2) as pool:
    list(pool.map(lambda j: run_all.happy_job(*j), jobs))
