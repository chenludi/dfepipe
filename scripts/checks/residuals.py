#!/usr/bin/env python3
"""Fit-quality metrics for a model-vs-data SFS comparison.

Residual convention follows dadi exactly: Anscombe_Poisson_residual computes
    resid = 1.5*(datatrans - modeltrans)/model**(1/6)
and returns -resid, so a POSITIVE residual means the model OVERESTIMATES that
frequency class. Bins where data == 0 are masked, as dadi masks them, and the
first and last entries are the fixed categories and never enter the metrics.

Metrics are reported over all unmasked bins and split into "rare" (allele counts
1..rare_max_count) and "rest", because the first few bins carry most of the SNPs
and swing hardest.
"""
import numpy as np

__all__ = ["anscombe_residual", "residual_metrics", "classify"]


def anscombe_residual(model, data):
    """dadi's Anscombe Poisson residual, masked where data == 0."""
    model = np.asarray(model, dtype=float)
    data = np.asarray(data, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        datatrans = data ** (2.0 / 3) - np.power(data, -1.0 / 3) / 9
        modeltrans = model ** (2.0 / 3) - np.power(model, -1.0 / 3) / 9
        resid = 1.5 * (datatrans - modeltrans) / model ** (1.0 / 6)
    masked = (data == 0) | ~np.isfinite(resid)
    return -resid, masked


def _block(resid):
    if resid.size == 0:
        return {"max_abs_resid": 0.0, "rms_resid": 0.0, "mean_resid": 0.0,
                "longest_same_sign_run": 0}
    signs = np.sign(resid)
    longest = best = 0
    prev = 0
    for s in signs:
        if s != 0 and s == prev:
            best += 1
        else:
            best = 1 if s != 0 else 0
        prev = s
        longest = max(longest, best)
    return {
        "max_abs_resid": float(np.max(np.abs(resid))),
        "rms_resid": float(np.sqrt(np.mean(resid ** 2))),
        "mean_resid": float(np.mean(resid)),
        "longest_same_sign_run": int(longest),
    }


def residual_metrics(model, data, rare_max_count=4):
    """Metrics over all bins, plus the same four split rare vs rest."""
    resid, masked = anscombe_residual(model, data)
    counts = np.arange(len(resid))
    interior = np.zeros(len(resid), dtype=bool)
    interior[1:-1] = True                      # drop the fixed categories
    use = interior & ~masked

    out = _block(resid[use])
    out["n_masked_bins"] = int(np.sum(interior & masked))
    out["n_snps"] = float(np.sum(np.asarray(data, dtype=float)[interior]))

    rare = use & (counts <= rare_max_count)
    rest = use & (counts > rare_max_count)
    for suffix, sel in (("rare", rare), ("rest", rest)):
        for k, v in _block(resid[sel]).items():
            out[f"{k}_{suffix}"] = v
    return out


def classify(metrics, good, warn, bias):
    """PASS / WARN / FAIL plus human-readable reasons.

    The criterion is fit quality, not convergence: convergence is reported
    elsewhere and never decides this verdict.
    """
    reasons = []
    status = "PASS"
    m = metrics.get("max_abs_resid", 0.0)
    if m > warn:
        status = "FAIL"
        reasons.append(f"max|resid| {m:.2f} exceeds {warn}")
    elif m > good:
        status = "WARN"
        reasons.append(f"max|resid| {m:.2f} exceeds {good}")
    mean = metrics.get("mean_resid", 0.0)
    if abs(mean) > bias:
        direction = "overestimates" if mean > 0 else "underestimates"
        if status == "PASS":
            status = "WARN"
        reasons.append(f"mean resid {mean:.2f}: model systematically {direction} the data")
    return status, reasons
