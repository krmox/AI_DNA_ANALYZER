"""Assemble the machine-readable verdict for the read-level experiment.

The success criterion is encoded here rather than argued in prose, so it cannot
drift from the numbers. The declared bar was: *a reproducible improvement over
the binomial baseline with few additional false positives, preferably
concentrated in the low-depth regime* -- not "F1 went up by 0.001".

Classification
--------------
* **positive** -- beats binomial v1 with a paired CI strictly above zero, at no
  more than a small FP cost, and replicates on a second seed;
* **weak_positive** -- point estimate above binomial but the CI touches zero or
  the effect does not replicate;
* **null** -- statistically indistinguishable from binomial v1;
* **negative** -- significantly worse than binomial v1.

Two further diagnostics decide the *scientific* reading, which is separate from
the leaderboard position:

* ``read_level_information_helps`` -- does the full model beat CONTROL 6, in
  which the base-to-attribute association is destroyed while every aggregate
  count is preserved byte for byte? This is the direct test of the hypothesis.
  If CONTROL 6 is *better*, the read-level association is not merely
  uninformative, it is actively harmful.
* ``architecture_matches_aggregate`` -- does CONTROL 1, which strips the input
  down to what the 14-channel representation already contains, match the
  existing aggregate model? This separates "the representation carries nothing"
  from "this architecture cannot exploit it".
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def classify(results: dict) -> dict:
    """Apply the declared criteria to the assembled test results."""
    primary = results["primary"]
    versus_binomial = results["vs_binomial_paired_bootstrap"]
    versus_aggregate = results["vs_aggregate_mamba_paired_bootstrap"]

    binomial = primary["binomial_v1"]
    model = primary["read_level"]
    delta = versus_binomial["read_level"]
    replicate = versus_binomial.get("readlevel_seed424242")
    decouple = primary.get("readlevel_ctrl6_decouple")
    aggregate_only = primary.get("readlevel_ctrl1_aggregate_only")
    permute = primary.get("readlevel_ctrl2_permute")
    errors = results["error_analysis"]

    evidence = {
        "read_level_f1": model["f1"],
        "binomial_f1": binomial["f1"],
        "delta_f1_vs_binomial": delta["delta_f1"],
        "delta_ci95": delta["ci95"],
        "ci_excludes_zero_positive": delta["ci95"][0] > 0,
        "ci_excludes_zero_negative": delta["ci95"][1] < 0,
        "replicate_delta_f1": replicate["delta_f1"] if replicate else None,
        "replicate_also_worse": bool(replicate and replicate["ci95"][1] < 0),
        "false_positives_binomial": binomial["fp"],
        "false_positives_read_level": model["fp"],
        "binomial_fn_recovered": errors["binomial_only_fn_recovered_by_model"],
        "new_false_positives": errors["model_only_fp"],
        "binomial_tp_overturned": errors["binomial_tp_overturned_by_model"],
        "delta_f1_vs_aggregate_mamba": versus_aggregate["read_level"]["delta_f1"],
    }

    if decouple is not None:
        evidence["control6_decouple_f1"] = decouple["f1"]
        evidence["full_minus_decoupled_f1"] = model["f1"] - decouple["f1"]
        evidence["read_level_information_helps"] = model["f1"] > decouple["f1"]
        evidence["read_level_association_is_harmful"] = decouple["f1"] > model["f1"]
    if aggregate_only is not None:
        evidence["control1_aggregate_only_f1"] = aggregate_only["f1"]
        evidence["full_minus_aggregate_only_f1"] = model["f1"] - aggregate_only["f1"]
    if permute is not None:
        evidence["control2_permutation_identical"] = (
            permute["tp"] == model["tp"] and permute["fp"] == model["fp"]
            and permute["fn"] == model["fn"])

    if delta["ci95"][1] < 0:
        verdict = "negative"
    elif delta["ci95"][0] > 0 and evidence["replicate_also_worse"] is False and replicate \
            and replicate["ci95"][0] > 0 and model["fp"] <= binomial["fp"] + 5:
        verdict = "positive"
    elif delta["delta_f1"] > 0:
        verdict = "weak_positive"
    else:
        verdict = "null"

    return {"verdict": verdict, "evidence": evidence}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", default="results/readlevel/readlevel_experiment_results.json")
    parser.add_argument("--sanity", default="results/readlevel/sanity_test.json")
    parser.add_argument("--training-json", nargs="*", default=[])
    parser.add_argument("--out", default="results/readlevel/readlevel_experiment_summary.json")
    args = parser.parse_args()

    results = json.loads(Path(args.results).read_text())
    sanity = json.loads(Path(args.sanity).read_text())
    verdict = classify(results)

    summary = {
        "experiment": "read_level_pileup_15x",
        "question": ("Can a neural model exploit read-level evidence that is unavailable "
                     "to the aggregate binomial caller?"),
        "representation": {
            "shape": "[loci, 48 reads, 14 per-read features] + locus reference one-hot",
            "max_reads": 48,
            "truncation_rate": sanity["summary"]["truncation_rate"],
            "read_order": "deterministic BLAKE2b hash of read name",
            "counts_reconstructed_exactly": sanity["checks"]["counts_reconstructed_exactly"],
        },
        "sanity_gate_passed": sanity["passed"],
        "primary_metric": "SNP vs Normal on the 15x held-out test region",
        "test_results": results["primary"],
        "paired_delta_vs_binomial": results["vs_binomial_paired_bootstrap"],
        "paired_delta_vs_aggregate_mamba": results["vs_aggregate_mamba_paired_bootstrap"],
        "by_depth": results["by_depth"],
        "error_analysis": results["error_analysis"],
        "oracle_diagnostic": results["oracle_diagnostic"],
        "runtime": results["runtime"],
        "training_runs": {},
        **verdict,
    }
    for path in args.training_json:
        run = json.loads(Path(path).read_text())
        summary["training_runs"][run["tag"]] = {
            "control": run["control"], "drop_groups": run["drop_groups"],
            "parameters": run["parameters"], "selected_epoch": run["selected_epoch"],
            "selected_val_macro_f1": run["selected_val_macro_f1"],
            "frozen_threshold": run["frozen_threshold"],
            "seconds_per_epoch": sum(h["train_seconds"] for h in run["history"])
            / max(len(run["history"]), 1),
        }

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(summary, indent=2))
    print(json.dumps({"verdict": summary["verdict"], "evidence": summary["evidence"]}, indent=2))


if __name__ == "__main__":
    main()
