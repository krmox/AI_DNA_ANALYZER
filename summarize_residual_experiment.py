"""Assemble the machine-readable verdict for the residual-binomial experiment.

Applies the pre-declared decision criteria to the artefacts produced by the
run and writes a single summary JSON. The classification rule is encoded here
rather than asserted in prose so that it cannot drift from the numbers:

* **strong_positive** -- paired dF1 CI strictly above 0, the gain replicates
  across seeds, no unacceptable FP increase, and the residual carries genuine
  locus-specific information (standalone AUC clearly above 0.5);
* **weak_positive** -- a positive point estimate that either fails to
  replicate, has a CI touching 0, or is fully explained by a constant shift;
* **null** -- indistinguishable from the binomial baseline;
* **negative** -- worse than the binomial baseline.

The residual-information test is what separates "the model learned something"
from "the model slid the threshold": a residual whose standalone ROC-AUC is
~0.5 and whose effect is reproduced by adding a scalar to the LLR is not new
evidence, whatever it does to F1.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

#: A residual must beat this standalone AUC to count as carrying real
#: locus-specific information rather than a global bias.
INFORMATIVE_AUC = 0.60


def classify(main: dict, replicate: dict, selection: dict) -> dict:
    """Apply the decision criteria to the assembled results."""
    binomial = main["primary"]["binomial_v1"]
    selected = main["primary"]["residual_selected"]
    trained = main["primary"]["residual_best_trained"]
    delta_selected = main["vs_binomial_paired_bootstrap"]["residual_selected"]
    delta_trained = main["vs_binomial_paired_bootstrap"]["residual_best_trained"]
    delta_replicate = replicate["vs_binomial_paired_bootstrap"]["residual_best_trained"]
    shift = main["constant_shift_analysis"]["residual_best_trained"]

    evidence = {
        "validation_prefers_untrained_model": selection["selected"]["epoch"] == 0,
        "selected_arm_delta_f1": delta_selected["delta_f1"],
        "selected_arm_ci_excludes_zero": delta_selected["ci95"][0] > 0,
        "best_trained_delta_f1": delta_trained["delta_f1"],
        "best_trained_ci_excludes_zero": delta_trained["ci95"][0] > 0,
        "replicate_delta_f1": delta_replicate["delta_f1"],
        "replicate_ci_excludes_zero": delta_replicate["ci95"][0] > 0,
        "gain_replicates_across_seeds": (delta_trained["ci95"][0] > 0
                                         and delta_replicate["ci95"][0] > 0),
        "residual_standalone_auc": shift["residual_alone_as_snp_score"]["roc_auc"],
        "residual_is_informative": (shift["residual_alone_as_snp_score"]["roc_auc"]
                                    > INFORMATIVE_AUC),
        "constant_shift_explains_most_of_gain": (
            shift["constant_shift_null"]["f1"] - binomial["f1"]
            >= 0.5 * (trained["f1"] - binomial["f1"])
            if trained["f1"] > binomial["f1"] else True),
        "false_positives_binomial": binomial["fp"],
        "false_positives_residual": trained["fp"],
    }

    def classify_arm(delta: dict, stats: dict) -> str:
        """Criteria applied to a single arm's paired delta against binomial."""
        if delta["delta_f1"] < 0 and delta["ci95"][1] < 0:
            return "negative"
        if (delta["ci95"][0] > 0 and evidence["gain_replicates_across_seeds"]
                and evidence["residual_is_informative"] and stats["fp"] <= binomial["fp"]):
            return "strong_positive"
        if delta["delta_f1"] > 0 and delta["ci95"][0] > 0:
            return "weak_positive"
        return "null"

    # The headline verdict belongs to the arm validation actually selects. The
    # best-trained arm is reported alongside it, but it cannot carry the
    # headline: it is the answer to "what if training is forced to move",
    # which is a diagnostic question, not the deployment question.
    arms = {
        "residual_selected": classify_arm(delta_selected, selected),
        "residual_best_trained": classify_arm(delta_trained, trained),
    }
    return {"verdict": arms["residual_selected"], "verdict_by_arm": arms,
            "headline_arm": "residual_selected",
            "headline_rationale": (
                "Validation selects the untrained (epoch 0) model at every lambda, so the "
                "deployed residual caller IS the binomial caller. The best-trained arm's "
                "+0.0038 F1 does not replicate across seeds and is reproduced by adding a "
                "scalar to the LLR, so it does not upgrade the headline."),
            "evidence": evidence}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--init-sanity", default="results/residual_init_sanity.json")
    parser.add_argument("--selection", default="results/residual_lambda_selection.json")
    parser.add_argument("--main", default="results/residual_experiment_results.json")
    parser.add_argument("--replicate", default="results/residual_seed424242_test.json")
    parser.add_argument("--out", default="results/residual_experiment_summary.json")
    args = parser.parse_args()

    sanity = json.loads(Path(args.init_sanity).read_text())
    selection = json.loads(Path(args.selection).read_text())
    main_results = json.loads(Path(args.main).read_text())
    replicate = json.loads(Path(args.replicate).read_text())

    verdict = classify(main_results, replicate, selection)
    summary = {
        "experiment": "residual_binomial_mamba_15x",
        "question": ("Can the Mamba architecture extract information from the raw pileup "
                     "that improves on the binomial LLR caller, rather than relearning or "
                     "distorting the same statistic?"),
        "architecture": "final_logit = binomial_LLR_prior + residual(14ch raw pileup), "
                        "residual head zero-initialized",
        "initialization_gate": {
            "passed": sanity["passed"],
            "pearson": sanity["seeds"][sorted(sanity["seeds"])[0]]["pearson_logit_vs_llr"],
            "spearman": sanity["seeds"][sorted(sanity["seeds"])[0]]["spearman_logit_vs_llr"],
            "max_abs_logit_difference":
                sanity["seeds"][sorted(sanity["seeds"])[0]]["max_abs_logit_difference"],
        },
        "lambda_sweep": {"values": [0.0, 1e-4, 1e-3, 1e-2],
                         "selected": selection["selected"],
                         "best_trained": selection["best_trained"],
                         "selected_on": "validation split only"},
        "primary_metric": "SNP vs Normal on the 15x held-out test region",
        "test_results": main_results["primary"],
        "paired_delta_vs_binomial": main_results["vs_binomial_paired_bootstrap"],
        "residual_diagnostics": main_results["residual_diagnostics"],
        "constant_shift_analysis": main_results["constant_shift_analysis"],
        "seed_replicate": {
            "seed": 424242,
            "primary": replicate["primary"]["residual_best_trained"],
            "delta_vs_binomial":
                replicate["vs_binomial_paired_bootstrap"]["residual_best_trained"],
            "residual_diagnostics": replicate["residual_diagnostics"]["residual_best_trained"],
        },
        "by_depth": main_results["by_depth"],
        **verdict,
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(summary, indent=2))
    print(json.dumps({"verdict": summary["verdict"], "evidence": summary["evidence"]}, indent=2))


if __name__ == "__main__":
    main()
