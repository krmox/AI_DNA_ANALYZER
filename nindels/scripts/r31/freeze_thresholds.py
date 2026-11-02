#!/usr/bin/env python3
"""Turn selection.json (development data only) into varion/thresholds.json. Run once; commit before any test run.

usage: freeze_thresholds.py SELECTION_JSON OUT_THRESHOLDS_JSON
"""
import hashlib
import json
import subprocess
import sys
from datetime import date
from pathlib import Path


def main():
    sel_path, out_path = Path(sys.argv[1]), Path(sys.argv[2])
    sel = json.loads(sel_path.read_text())
    dev = sel["all_dev"]
    tg, tl = dev["T_global_ladder_key"], dev["T_local_lower_edge"]
    ratio = sel["dp_calibration"]["median_dp_over_depth"] or 1.0
    notes = []
    if tg is None:
        g = None
        notes.append("no ladder depth met the rule (every higher depth: dF1>0 with 90% CI lower bound >0 and positive on both chromosomes): auto mode uses A everywhere")
    else:
        g = float(sel["ladder_targets"][tg])
    if tl is None:
        l = None
        notes.append("no local-depth bin met the rule")
    else:
        l = max(1, int(round(tl * ratio)))
        notes.append(f"local gate: lower bin edge {tl}x (samtools depth) x median engine-DP/depth {ratio:.3f} -> {l} reads")
    loco = {ds: (r["T_global_ladder_key"], r["T_local_lower_edge"]) for ds, r in sel["loco"].items()}
    notes.append("leave-one-chromosome-out (global key, local edge): " + json.dumps(loco))
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        commit = "unknown"
    out = {
        "status": "selected" if g is not None else "no depth found at which B reliably helps",
        "global_min_depth": g,
        "local_min_depth": l if g is not None else None,
        "provenance": f"{date.today()}; selection on HG004 chr20+chr21 depth ladder only (docs/R3.1_THRESHOLD_SELECTION.md); "
                      f"code commit {commit}; selection.json sha256 {hashlib.sha256(sel_path.read_bytes()).hexdigest()[:16]}; " + "; ".join(notes),
    }
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
