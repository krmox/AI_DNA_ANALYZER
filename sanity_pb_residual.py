"""Pre-flight checks for the PB-residual experiment. Run before any GPU time.

Three questions, cheapest first:

1. **Is the initialization guarantee real on real data?** The untrained model
   must reproduce the Poisson-binomial caller's decision score *exactly*, not
   approximately, or "the model improved on PB" is not a meaningful claim.

2. **Is the representation non-degenerate?** Constant, non-finite or absurdly
   scaled channels would waste a training run.

3. **Is there any exploitable signal above PB at all?** A logistic-regression
   residual probe -- ``score = PB_llr + w . features`` -- is fitted on the fit
   split and evaluated on the validation split. This is a linear model with no
   sequence context and ~45 parameters; it takes seconds. If even it cannot
   beat PB on validation, a 6-layer Mamba is unlikely to, and the cheap
   negative result is worth far more than a GPU-hour.

   The probe is a *lower bound*, not a prediction: it cannot use neighbouring
   loci and cannot express interactions, both of which are the backbone's job.

Nothing here touches the test region.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import torch

from config import LABEL_NORMAL, LABEL_SNP
from locus_evidence import FEATURE_GROUPS, FEATURE_NAMES
from model_pb_residual import PoissonBinomialResidualMamba
from residual_metrics import prf
from residual_prior import prior_logits_from_llr, snp_decision_score
from train_raw_pileup import BATCH_SIZE, SEQ_LEN, TRAIN_FRACTION
from train_residual_binomial import select_threshold

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


def split_point(n_loci: int) -> int:
    """The fit/validation boundary, byte-identical to the training scripts."""
    n_windows = n_loci // SEQ_LEN
    n_batches = (n_windows + BATCH_SIZE - 1) // BATCH_SIZE
    return int(n_batches * TRAIN_FRACTION) * BATCH_SIZE * SEQ_LEN


def check_initialization(features: np.ndarray, llr: np.ndarray, device) -> dict:
    """The untrained model must equal the PB caller exactly."""
    torch.manual_seed(0)
    model = PoissonBinomialResidualMamba().to(device).eval()

    sample = features[: 64 * 64]
    sample_llr = llr[: 64 * 64]
    feature_windows = torch.tensor(sample, dtype=torch.float).reshape(-1, SEQ_LEN,
                                                                     features.shape[1])
    prior_windows = torch.tensor(prior_logits_from_llr(sample_llr),
                                 dtype=torch.float).reshape(-1, SEQ_LEN, 4)
    with torch.no_grad():
        logits, residual, _ = model(feature_windows.to(device), prior_windows.to(device),
                                    return_parts=True)
    score = snp_decision_score(logits).reshape(-1).cpu().numpy().astype(np.float64)
    max_deviation = float(np.max(np.abs(score - sample_llr)))
    logger.info("Init check: max |model_score - PB_llr| = %.3e over %d loci",
                max_deviation, sample_llr.size)
    return {"max_abs_deviation": max_deviation,
            "residual_all_zero": bool(torch.all(residual == 0).item()),
            "loci": int(sample_llr.size)}


def check_representation(features: np.ndarray) -> dict:
    """Per-channel ranges and degeneracy check."""
    report = {}
    degenerate = []
    for index, name in enumerate(FEATURE_NAMES):
        column = features[:, index]
        entry = {"min": float(column.min()), "max": float(column.max()),
                 "mean": float(column.mean()), "std": float(column.std())}
        report[name] = entry
        if entry["std"] == 0.0:
            degenerate.append(name)
        if not np.isfinite(column).all():
            raise ValueError(f"channel {name} is not finite")
    logger.info("Representation: %d channels, %d constant, max |value| = %.2f",
                features.shape[1], len(degenerate), float(np.abs(features).max()))
    if degenerate:
        logger.warning("Constant channels (carry no information here): %s", degenerate)
    return {"channels": report, "constant_channels": degenerate}


def fit_logistic_probe(x: np.ndarray, y: np.ndarray, offset: np.ndarray,
                       epochs: int = 300, lr: float = 0.05, device=None) -> torch.Tensor:
    """Fit ``sigmoid(offset + w . x + b)`` by gradient descent on balanced classes.

    The PB LLR enters as a fixed offset, so the weights can only learn what PB
    does *not* already express -- the linear analogue of the residual model.

    Args:
        x: ``[N, F]`` standardized features.
        y: ``[N]`` float labels.
        offset: ``[N]`` PB log-odds, held fixed.
        epochs: Full-batch gradient steps.
        lr: Learning rate.
        device: Torch device.

    Returns:
        The fitted weight vector including its bias as the last element.
    """
    x_t = torch.tensor(x, dtype=torch.float32, device=device)
    y_t = torch.tensor(y, dtype=torch.float32, device=device)
    offset_t = torch.tensor(offset, dtype=torch.float32, device=device)

    # Class-balanced weighting, mirroring the training loop's intent.
    positive = y_t.sum()
    weight = torch.where(y_t > 0, y_t.numel() / (2 * positive.clamp(min=1)),
                         y_t.numel() / (2 * (y_t.numel() - positive).clamp(min=1)))

    parameters = torch.zeros(x.shape[1] + 1, device=device, requires_grad=True)
    optimizer = torch.optim.Adam([parameters], lr=lr)
    for _ in range(epochs):
        optimizer.zero_grad()
        logit = offset_t + x_t @ parameters[:-1] + parameters[-1]
        loss = torch.nn.functional.binary_cross_entropy_with_logits(
            logit, y_t, weight=weight)
        loss.backward()
        optimizer.step()
    return parameters.detach()


def probe_experiment(features, llr, labels, device) -> dict:
    """Fit the residual probe on the fit split, evaluate it on validation."""
    boundary = split_point(labels.size)
    frame = (labels == LABEL_NORMAL) | (labels == LABEL_SNP)

    fit_mask = np.zeros(labels.size, dtype=bool)
    fit_mask[:boundary] = True
    fit_mask &= frame
    val_mask = np.zeros(labels.size, dtype=bool)
    val_mask[boundary:] = True
    val_mask &= frame

    snp = labels == LABEL_SNP
    mean = features[fit_mask].mean(axis=0)
    std = features[fit_mask].std(axis=0)
    std[std == 0] = 1.0

    def standardize(mask):
        return (features[mask] - mean) / std

    weights = fit_logistic_probe(standardize(fit_mask), snp[fit_mask].astype(np.float32),
                                 llr[fit_mask], device=device)

    val_x = standardize(val_mask)
    learned = val_x @ weights[:-1].cpu().numpy() + float(weights[-1])
    pb_score = llr[val_mask]
    probe_score = pb_score + learned
    truth = snp[val_mask]

    results = {}
    for name, score in (("poisson_binomial", pb_score), ("probe", probe_score)):
        # ``labels[val_mask]`` is already the SNP-vs-Normal frame, so
        # select_threshold's internal masking is a no-op here.
        threshold, _ = select_threshold(score, labels[val_mask])
        results[name] = {"threshold": threshold, **prf(score >= threshold, truth)}

    ranking = sorted(zip(FEATURE_NAMES, weights[:-1].cpu().numpy().tolist()),
                     key=lambda pair: -abs(pair[1]))
    results["top_weights"] = [{"feature": name, "weight": value} for name, value in ranking[:12]]
    results["learned_term"] = {
        "mean": float(learned.mean()), "std": float(learned.std()),
        "min": float(learned.min()), "max": float(learned.max()),
    }
    results["group_weight_norm"] = {
        group: float(np.linalg.norm([weights[FEATURE_NAMES.index(name)].item()
                                     for name in names]))
        for group, names in FEATURE_GROUPS.items()
    }
    results["validation_loci"] = int(val_mask.sum())
    results["validation_snp"] = int(truth.sum())
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-evidence", default="cache/train_15x_evidence.npz")
    parser.add_argument("--out", default="results/pb_residual_sanity.json")
    args = parser.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    blob = np.load(args.train_evidence, allow_pickle=True)
    features, labels, llr = blob["features"], blob["labels"], blob["pb_llr"]
    logger.info("Loaded %d loci x %d channels", features.shape[0], features.shape[1])

    report = {
        "initialization": check_initialization(features, llr, device),
        "representation": check_representation(features),
        "probe": probe_experiment(features, llr, labels, device),
    }

    probe = report["probe"]
    logger.info("VALIDATION | PB    thr=%.1f F1=%.4f P=%.4f R=%.4f TP=%d FP=%d FN=%d",
                probe["poisson_binomial"]["threshold"], probe["poisson_binomial"]["f1"],
                probe["poisson_binomial"]["precision"], probe["poisson_binomial"]["recall"],
                probe["poisson_binomial"]["tp"], probe["poisson_binomial"]["fp"],
                probe["poisson_binomial"]["fn"])
    logger.info("VALIDATION | PROBE thr=%.1f F1=%.4f P=%.4f R=%.4f TP=%d FP=%d FN=%d",
                probe["probe"]["threshold"], probe["probe"]["f1"],
                probe["probe"]["precision"], probe["probe"]["recall"],
                probe["probe"]["tp"], probe["probe"]["fp"], probe["probe"]["fn"])
    logger.info("Learned term: mean=%+.3f std=%.3f", probe["learned_term"]["mean"],
                probe["learned_term"]["std"])
    logger.info("Group weight norms: %s", probe["group_weight_norm"])
    logger.info("Top weights: %s", [(entry["feature"], round(entry["weight"], 3))
                                    for entry in probe["top_weights"][:8]])

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2, default=float))
    logger.info("Wrote %s", args.out)


if __name__ == "__main__":
    main()
