#!/usr/bin/env python3
"""Run a shell command and record wall time, CPU time and peak RSS.

Glue only. /usr/bin/time is not installed here, so this uses
resource.getrusage(RUSAGE_CHILDREN). For a pipeline (a | b), peak RSS is that
of the largest single process, not the sum.

usage: timed.py OUT_JSON -- command ...
"""
import json
import resource
import subprocess
import sys
import time


def main():
    out = sys.argv[1]
    cmd = sys.argv[sys.argv.index("--") + 1:]
    t0 = time.monotonic()
    rc = subprocess.call(["bash", "-o", "pipefail", "-c", " ".join(cmd)])
    wall = time.monotonic() - t0
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    rec = {
        "cmd": " ".join(cmd),
        "rc": rc,
        "wall_s": round(wall, 2),
        "user_s": round(ru.ru_utime, 2),
        "sys_s": round(ru.ru_stime, 2),
        "cpu_s": round(ru.ru_utime + ru.ru_stime, 2),
        "peak_rss_mb": round(ru.ru_maxrss / 1024, 1),
    }
    with open(out, "w") as f:
        json.dump(rec, f, indent=1)
    print(json.dumps(rec), file=sys.stderr)
    sys.exit(rc)


if __name__ == "__main__":
    main()
