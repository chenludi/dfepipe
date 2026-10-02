#!/usr/bin/env python3
"""Fit-quality metrics for a model-vs-data SFS comparison.

For reference only. This check flags fits worth a second look; it is not an
accept/reject test. A warning means "the residuals here are large", which is
sometimes the expected answer - the neutral null is meant to fit badly on sites
under selection, and it does.

Residual and deviance thresholds only ever warn, however extreme, which is why
the config keys are named `..._LARGE` and `..._BADFIT`. They were once named
`..._WARN` and `..._FAIL`, which implied the second tier stopped a run; it never
did. FAIL is reserved for structural problems - no usable bins, or a model giving
no usable probability where variants are observed - and those do block a chain.

The thresholds were set from 137 gamma and neugamma fits (calibrated 2026-09-29).
The calibration selected runs by `*_model_sfs.fs`, and the `_neu` directories
hold none, so no neutral null contributed to them and there is no past-run range
for that model at all. Messages therefore state the value and the threshold it
crossed, and claim nothing about where previous runs landed. The calibration data
is not shipped with this package.

The metrics kept are the ones that separated the one labelled bad fit from the
137:

- ``max_abs_rare``    worst |residual| over allele counts 1..rare_max_count.
                      0 of 137 accepted fits exceed the bad fit (47.2), and it
                      does not scale with SNP count (r = -0.12).
- ``g2_per_snp``      likelihood-ratio deviance per SNP. Also 0 of 137 worse
                      than the bad fit (0.0395); the good control fit on the
                      same data scores 0.000114.
- ``mean_rest_norm``  bulk bias, normalised by sqrt(n_snps/1e5) because the raw
                      |mean_rest| correlates +0.51 with log10(n_snps).
- ``max_abs_rest``    only as "the worst bin sits outside the rare classes".

``rms``, ``max_abs_all``, raw signed bias, same-sign run length and raw G2 were
all tested and failed to separate (10-23 accepted fits scored worse than the
labelled bad fit), so they are deliberately not computed here.

Conventions:

- Residuals are dadi's Anscombe Poisson residual, with dadi's sign convention:
  dadi computes ``resid = 1.5*(datatrans - modeltrans)/model**(1./6)`` and
  returns ``-resid``, so a POSITIVE residual means the model OVERESTIMATES that
  frequency class.
- Index 0 and the last index are the fixed categories and never enter anything.
- Bins with ``data == 0`` are legitimate (no variants observed): they are
  excluded from the metrics and counted in ``n_masked_data_zero``.
- Bins with ``data > 0`` whose residual is not finite mean the model gives ~0
  probability where variants are observed. They are excluded from the metrics,
  counted in ``n_masked_model_zero``, and always FAIL: they must never be
  silently dropped the way ``data == 0`` bins legitimately are.
"""
import numpy as np

try:                                            # pragma: no cover - env dependent
    from dadi import Inference as _dadi_Inference
    USING_DADI = True
except Exception:                               # pragma: no cover - env dependent
    _dadi_Inference = None
    USING_DADI = False

__all__ = ["anscombe_residual", "residual_metrics", "classify", "USING_DADI"]


def anscombe_residual(model, data):
    """dadi's Anscombe Poisson residual, as a plain float array.

    Prefers ``dadi.Inference.Anscombe_Poisson_residual`` so the pipeline and
    dadi can never drift apart. The local formula is only a fallback for
    development machines where dadi is not installed (it is installed on the
    cluster, where the real fits run); ``USING_DADI`` says which path ran.

    Non-finite entries are left in place: the caller decides what they mean.
    """
    model = np.asarray(model, dtype=float)
    data = np.asarray(data, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        if USING_DADI:
            resid = _dadi_Inference.Anscombe_Poisson_residual(model, data)
        else:
            datatrans = data ** (2.0 / 3) - np.power(data, -1.0 / 3) / 9
            modeltrans = model ** (2.0 / 3) - np.power(model, -1.0 / 3) / 9
            resid = -(1.5 * (datatrans - modeltrans) / model ** (1.0 / 6))
    # dadi hands back a masked array when given Spectrum objects.
    resid = np.ma.filled(np.ma.asarray(resid, dtype=float), np.nan)
    return np.asarray(resid, dtype=float)


def _max_abs(values):
    return float(np.max(np.abs(values))) if values.size else 0.0


def residual_metrics(model, data, rare_max_count=4):
    """The calibrated metric set for one model-vs-data SFS comparison.

    Returns exactly the keys documented in the module docstring, nothing more.
    """
    model = np.asarray(model, dtype=float)
    data = np.asarray(data, dtype=float)
    resid = anscombe_residual(model, data)

    counts = np.arange(data.size)
    interior = np.zeros(data.size, dtype=bool)
    interior[1:-1] = True                       # drop the fixed categories

    data_zero = interior & (data == 0)
    model_zero = interior & (data > 0) & ~np.isfinite(resid)
    use = interior & (data > 0) & np.isfinite(resid)

    rare = resid[use & (counts <= rare_max_count)]
    rest = resid[use & (counts > rare_max_count)]

    n_snps = float(np.sum(data[interior]))
    mean_rest = float(np.mean(rest)) if rest.size else 0.0
    if n_snps > 0:
        mean_rest_norm = abs(mean_rest) / np.sqrt(n_snps / 1e5)
    else:
        mean_rest_norm = 0.0

    d = data[use]
    m = model[use]
    with np.errstate(divide="ignore", invalid="ignore"):
        g2 = float(2.0 * np.sum(d * np.log(d / m) - (d - m))) if d.size else 0.0

    return {
        "max_abs_rare": _max_abs(rare),
        "max_abs_rest": _max_abs(rest),
        "mean_rest": mean_rest,
        "mean_rest_norm": float(mean_rest_norm),
        "g2": g2,
        "g2_per_snp": float(g2 / n_snps) if n_snps > 0 else 0.0,
        "n_snps": n_snps,
        "n_bins_used": int(np.sum(use)),
        "n_masked_data_zero": int(np.sum(data_zero)),
        "n_masked_model_zero": int(np.sum(model_zero)),
    }


def classify(metrics, cfg):
    """PASS / WARN / FAIL plus short human-readable reasons.

    ``cfg`` carries ``max_rare_large``, ``max_rare_badfit``, ``bulk_bias_warn``,
    ``rest_over_rare_min``, ``g2_large`` and ``g2_badfit``. The verdict is about
    fit quality only: convergence is reported elsewhere and never decides it.

    Remember the verdict is advisory (see the module docstring): a warning means
    "this fit is worth looking at", not "wrong", and nothing here blocks.
    """
    max_rare = float(metrics.get("max_abs_rare", 0.0))
    max_rest = float(metrics.get("max_abs_rest", 0.0))
    mean_rest = float(metrics.get("mean_rest", 0.0))
    bias = float(metrics.get("mean_rest_norm", 0.0))
    g2_per_snp = float(metrics.get("g2_per_snp", 0.0))
    n_model_zero = int(metrics.get("n_masked_model_zero", 0))

    fails = []
    if int(metrics.get("n_bins_used", 1)) == 0:
        fails.append("no usable bins: every interior bin was masked, so nothing was measured")
    if n_model_zero > 0:
        # The symptom is a non-finite residual where variants exist. The usual cause is a
        # model probability of ~0 there, but an invalid model value gives the same symptom,
        # so the wording covers both rather than asserting the cause.
        fails.append(
            f"{n_model_zero} bin(s) where the model gives no usable probability "
            "(~zero or invalid) but variants are observed")
    if fails:
        return "FAIL", fails

    # Residual and deviance thresholds only ever warn, which is why the config
    # keys are named ..._LARGE and ..._BADFIT rather than ..._WARN and ..._FAIL:
    # the old names implied the second tier stopped a run, and it never did.
    # Residuals are tricky - they scale with SNP count and the rare bins swing
    # hard for real reasons - so FAIL is reserved above for structural problems.
    warns = []
    # Plain wording, and no claim about past runs. The thresholds come from
    # gamma and neugamma fits only, so "beyond anything seen before" was not true
    # of the neutral null, which has no calibration set at all.
    if max_rare > cfg["max_rare_badfit"]:
        warns.append(
            f"the rare-bin residual is large ({max_rare:.1f}, over "
            f"{cfg['max_rare_badfit']}): the model is not a good fit to the rare "
            "frequency classes; look at the figure")
    elif max_rare > cfg["max_rare_large"]:
        warns.append(
            f"the rare-bin residual is on the large side ({max_rare:.1f}, over "
            f"{cfg['max_rare_large']}): worth a look at the figure")
    if g2_per_snp > cfg["g2_badfit"]:
        warns.append(
            f"the deviance per SNP is large ({g2_per_snp:.4g}, over "
            f"{cfg['g2_badfit']}): the model is not a good fit; look at the "
            "figure")
    elif g2_per_snp > cfg["g2_large"]:
        warns.append(
            f"the deviance per SNP is on the large side ({g2_per_snp:.4g}, over "
            f"{cfg['g2_large']})")
    if bias > cfg["bulk_bias_warn"]:
        direction = "overestimates" if mean_rest > 0 else "underestimates"
        warns.append(
            f"the model systematically {direction} the common frequency classes "
            f"(bulk bias {bias:.2f}, over {cfg['bulk_bias_warn']})")
    if max_rest > max_rare and max_rest > cfg["rest_over_rare_min"]:
        warns.append(
            f"the worst-fitting bin is not a rare one: the common-class residual "
            f"({max_rest:.1f}) is worse than the rare-class residual "
            f"({max_rare:.1f})")
    if warns:
        return "WARN", warns

    return "PASS", []
