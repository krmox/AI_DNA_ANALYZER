"""Pick the residual penalty lambda and the checkpoint, on validation only.

Reads the per-lambda training JSONs and reports two validation-only choices:

* ``selected`` -- the best over *all* candidates including epoch 0, the
  untrained model. Epoch 0 is a legitimate candidate: it is what early
  stopping would return if no amount of training improves on the prior, and
  excluding it would force a worse model into the comparison to manufacture a
  result;
* ``best_trained`` -- the best over epochs 1..N only, which is what the
  residual diagnostics need in order to have a non-zero residual to analyse.

Ties on validation F1 are broken by validation PR-AUC. This matters here: the
validation split holds ~251 SNPs, so F1 moves in steps of about 0.002 and
several epochs land on exactly the same value. PR-AUC is continuous and, like
F1, is computed on validation only.

The test set is not read by this script.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def candidates(runs: dict, include_epoch_zero: bool):
    """Yield ``(lambda, epoch, val_metrics)`` triples across every run."""
    for name, run in runs.items():
        if include_epoch_zero:
            yield run["residual_lambda"], 0, run["epoch0_untrained_binomial"], name
        for record in run["history"]:
            yield run["residual_lambda"], record["epoch"], record["val"], name


def best(runs: dict, include_epoch_zero: bool) -> dict:
    """Argmax over (val F1, val PR-AUC)."""
    lam, epoch, metrics, name = max(
        candidates(runs, include_epoch_zero),
        key=lambda item: (item[2]["f1"], item[2]["pr_auc"]))
    return {"residual_lambda": lam, "epoch": epoch, "run": name,
            "val_f1": metrics["f1"], "val_pr_auc": metrics["pr_auc"],
            "val_roc_auc": metrics["roc_auc"], "threshold": metrics["threshold"],
            "val_mean_squared_residual": metrics["mean_squared_residual"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-json", nargs="+", required=True)
    parser.add_argument("--out", default="results/residual_lambda_selection.json")
    args = parser.parse_args()

    runs = {}
    for path in args.training_json:
        run = json.loads(Path(path).read_text())
        runs[run["tag"]] = run

    report = {
        "runs": {tag: {"residual_lambda": run["residual_lambda"],
                       "epoch0_val_f1": run["epoch0_untrained_binomial"]["f1"],
                       "best_trained_val_f1": max(r["val"]["f1"] for r in run["history"]),
                       "final_mean_squared_residual":
                           run["history"][-1]["val"]["mean_squared_residual"],
                       "per_epoch": [{"epoch": r["epoch"], "val_f1": r["val"]["f1"],
                                      "val_pr_auc": r["val"]["pr_auc"],
                                      "threshold": r["val"]["threshold"],
                                      "mean_squared_residual":
                                          r["val"]["mean_squared_residual"]}
                                     for r in run["history"]]}
                 for tag, run in runs.items()},
        "selected": best(runs, include_epoch_zero=True),
        "best_trained": best(runs, include_epoch_zero=False),
        "note": "validation split only; the test region was not read by this script",
    }
    report["any_trained_epoch_beats_untrained"] = bool(
        report["best_trained"]["val_f1"] > report["selected"]["val_f1"]
        or report["selected"]["epoch"] > 0)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2))
    print(json.dumps({k: report[k] for k in
                      ("selected", "best_trained", "any_trained_epoch_beats_untrained")},
                     indent=2))


if __name__ == "__main__":
    main()
