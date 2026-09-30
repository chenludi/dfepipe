#!/usr/bin/env python3
"""Checkpoint for the demo stage: does the demographic model fit the control SFS?

FOR REFERENCE ONLY, as everywhere in this pipeline: a WARN means "unusual
compared with 137 previous fits", not "wrong". The figure remains the arbiter.

FAIL is reserved for structural problems, which do block a chain:
    no fits at all                     nothing was inferred
    best_inference.txt missing/unparseable   the two keys later stages grep for
                                       (`max_ll_key:`, `Na:`) are not there
    Na <= 0                            the cache stage derives its gamma bounds
                                       from 2*s*Na, so a non-positive Na makes
                                       every downstream number meaningless
    no usable bins / model gives no probability where variants are observed
                                       (from residuals.classify)

WARN is advisory: the calibrated residual and deviance thresholds, which never
FAIL however extreme, plus poor convergence. Convergence across starting points
is desirable - it makes a single local optimum less likely - but a converged bad
fit is still a bad fit and a scattered good fit is still a good fit, so it only
ever warns.

The residual metrics are not recomputed here: run_bestfit.sh writes the
best-fit model spectrum beside the fits as `bestfit_model.fs`, and this script
compares it with the control spectrum through scripts/checks/residuals.py.
"""
import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checkfile import write_check  # noqa: E402
from conf import conf  # noqa: E402
from residuals import classify, residual_metrics  # noqa: E402

# Top fits within this many log-likelihood units count as converged. One ll unit
# is the conventional "same optimum" tolerance for a nested-model comparison.
CONVERGENCE_LL_SPREAD = 1.0

# Written by run_bestfit.sh; the model spectrum for the maximum-likelihood fit.
MODEL_FS_NAME = "bestfit_model.fs"
BEST_INFERENCE_NAME = "best_inference.txt"
# The control spectrum the fits were made against. The tail of the name carries
# the bed's short name and CTRL_TAG, so it is globbed rather than rebuilt: the
# file the sfs stage wrote is the authority.
CTRL_FS_GLOB = "quiecent_ctrl.*.chrs_1t22.fs"

# best_inference.txt keys that must be present and numeric for the file to count
# as parseable. `max_ll_key` and `Na` are the two the later stages grep for.
_REQUIRED_KEYS = ("max_ll_key", "Na")
_FLOAT_KEYS = ("Na", "theta", "ll", "mutation_ctrl_diploid", "likelyhood_cutoff")


# ------------------------------------------------------------------- verdict --

def assess_demo(metrics, top_lls, cfg, na=None, n_fits=None,
                ll_spread_warn=CONVERGENCE_LL_SPREAD):
    """(status, reasons) for one bed. See the module docstring.

    `metrics` is a residual_metrics() mapping, or None when the spectra could
    not be read. `top_lls` are the best log-likelihoods across fits, `na` the Na
    from best_inference.txt (None when that file is missing or unparseable) and
    `cfg` the threshold mapping residuals.classify takes.
    """
    fails = []
    if n_fits is not None and n_fits < 1:
        fails.append(
            "no fits produced: no *.InferDM.bestfits in the run directory. Either "
            "run_demo.sh has not run, or every array task failed - check the SGE "
            "logs, and that DEMO_P0_DIR holds the p0 files")
    if na is None:
        fails.append(
            "%s is missing or unparseable: the later stages grep it for "
            "'max_ll_key:' and 'Na:', so nothing downstream can run. Re-run "
            "scripts/demo/run_bestfit.sh" % BEST_INFERENCE_NAME)
    elif na <= 0:
        fails.append(
            "Na is %g, not positive: theta or the control mutation count is "
            "wrong, so every 2*s*Na gamma bound downstream would be too" % na)
    if metrics is None:
        fails.append(
            "the model and data spectra could not be read, so the fit was not "
            "measured: check that %s sits beside the fits" % MODEL_FS_NAME)
    if fails:
        return "FAIL", fails

    # Residual and deviance thresholds only ever WARN; classify FAILs only on
    # structural problems (no usable bins, no model probability where variants
    # are observed). A structural FAIL is the whole verdict: convergence is a
    # remark about an inference that at least produced something measurable.
    status, reasons = classify(metrics, cfg)
    if status == "FAIL":
        return status, reasons

    lls = [float(x) for x in top_lls]
    if len(lls) >= 2:
        spread = max(lls) - min(lls)
        if spread > ll_spread_warn:
            if status == "PASS":
                status = "WARN"
            reasons.append(
                "the top %d fits span %.1f log-likelihood units (> %.1f): they "
                "did not converge on one optimum, so the best fit may be local"
                % (len(lls), spread, ll_spread_warn))
    return status, reasons


# ------------------------------------------------------------------- readers --

def parse_best_inference(path):
    """best_inference.txt as a dict, or None if it is absent or unusable.

    The file is `key:value` per line, written by run_bestfit.sh in the
    historical order. Values that should be numbers are converted; a value that
    will not convert leaves the key as the raw string, except for the keys the
    verdict depends on, whose absence makes the whole file unusable.
    """
    try:
        with open(path) as fh:
            lines = fh.read().splitlines()
    except OSError:
        return None

    out = {}
    for line in lines:
        if ":" not in line:
            continue
        key, _, value = line.partition(":")
        out[key.strip()] = value.strip()

    if any(k not in out for k in _REQUIRED_KEYS):
        return None
    for k in _FLOAT_KEYS:
        if k in out:
            try:
                out[k] = float(out[k])
            except ValueError:
                if k in _REQUIRED_KEYS:
                    return None
    if not isinstance(out.get("Na"), float):
        return None
    return out


def read_bestfits_ll(path):
    """The top converged log-likelihood in a dadi-cli .bestfits file, or None.

    Same rule the best-fit picker uses (plot_bestll_demo_args_0421.py:139-180):
    the first data row of the "# Converged results" block. A file with no such
    block did not converge and contributes nothing.
    """
    try:
        with open(path) as fh:
            lines = fh.read().splitlines()
    except OSError:
        return None
    found = False
    for line in lines:
        stripped = line.strip()
        if stripped == "# Converged results":
            found = True
            continue
        if not found:
            continue
        if stripped.startswith("#"):
            if stripped == "#":          # end of the converged block
                return None
            continue
        if not stripped:
            continue
        try:
            return float(stripped.split("\t")[0])
        except ValueError:
            return None
    return None


def read_fs(path):
    """A dadi 1-D .fs as a float list with masked bins zeroed, or None.

    Zeroing is how a masked bin is excluded: residual_metrics drops every bin
    where data == 0, and a folded 160-chromosome spectrum masks its whole upper
    half, which would otherwise be read as real zeros in the model.
    """
    try:
        with open(path) as fh:
            rows = [r for r in fh.read().splitlines() if r.strip()]
    except OSError:
        return None
    rows = [r for r in rows if not r.startswith("#")]
    if len(rows) < 2:
        return None
    try:
        values = [float(x) for x in rows[1].split()]
    except ValueError:
        return None
    if len(rows) >= 3:
        mask = rows[2].split()
        if len(mask) == len(values):
            values = [0.0 if m not in ("0", "0.0", "False") else v
                      for v, m in zip(values, mask)]
    return values


def _bed_base(bed):
    base = os.path.basename(bed)
    return base[:-4] if base.endswith(".bed") else base


# ---------------------------------------------------------------------- rows --

def row_for(demo_dir, ctrl_fs_dir, base, cfg, rare_max_count, num_convergence):
    """The single check row for one bed."""
    fits = sorted(glob.glob(os.path.join(demo_dir, "*.InferDM.bestfits")))
    lls = [ll for ll in (read_bestfits_ll(f) for f in fits) if ll is not None]
    lls.sort(reverse=True)
    top = lls[:num_convergence]

    best = parse_best_inference(os.path.join(demo_dir, BEST_INFERENCE_NAME))
    na = best["Na"] if best else None

    ctrl = sorted(glob.glob(os.path.join(ctrl_fs_dir, CTRL_FS_GLOB)))
    data = read_fs(ctrl[0]) if ctrl else None
    model = read_fs(os.path.join(demo_dir, MODEL_FS_NAME))
    metrics = None
    if data is not None and model is not None and len(data) == len(model):
        metrics = residual_metrics(model, data, rare_max_count=rare_max_count)

    status, reasons = assess_demo(metrics, top, cfg, na=na, n_fits=len(fits))

    row = {
        "item": base,
        "status": status,
        "reason": "; ".join(reasons) or "NA",
        "n_fits": len(fits),
        "n_converged": len(lls),
        "max_ll_key": best.get("max_ll_key", "NA") if best else "NA",
        "best_ll": lls[0] if lls else "NA",
        "ll_spread_top": round(max(top) - min(top), 4) if len(top) >= 2 else "NA",
        "Na": na if na is not None else "NA",
        "theta": best.get("theta", "NA") if best else "NA",
        "mutation_ctrl_diploid": best.get("mutation_ctrl_diploid", "NA") if best else "NA",
    }
    for key in ("max_abs_rare", "max_abs_rest", "mean_rest", "mean_rest_norm",
                "g2", "g2_per_snp", "n_snps", "n_bins_used",
                "n_masked_data_zero", "n_masked_model_zero"):
        value = metrics.get(key, "NA") if metrics else "NA"
        row[key] = round(value, 6) if isinstance(value, float) else value
    return row


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True, help="the demo output root, $OUT/demo")
    ap.add_argument("--run", required=True,
                    help="the run directory name, <FS_RUN>-<DEMO_MODEL>-<DEMO_RUN>")
    ap.add_argument("--sfs-out", required=True, help="the sfs output root, $OUT/sfs")
    ap.add_argument("--fs-run", required=True, help="FS_RUN")
    ap.add_argument("--rare-max-count", type=int,
                    default=conf.get_int("RESID_RARE_MAX_COUNT", 4))
    ap.add_argument("--max-rare-warn", type=float,
                    default=conf.get_float("RESID_MAX_RARE_WARN", 25.0))
    ap.add_argument("--max-rare-fail", type=float,
                    default=conf.get_float("RESID_MAX_RARE_FAIL", 40.0))
    ap.add_argument("--bulk-bias-warn", type=float,
                    default=conf.get_float("RESID_BULK_BIAS_WARN", 2.3))
    ap.add_argument("--rest-over-rare-min", type=float,
                    default=conf.get_float("RESID_REST_OVER_RARE_MIN", 5.0))
    ap.add_argument("--g2-warn", type=float,
                    default=conf.get_float("G2_PER_SNP_WARN", 0.008))
    ap.add_argument("--g2-fail", type=float,
                    default=conf.get_float("G2_PER_SNP_FAIL", 0.035))
    ap.add_argument("--num-convergence", type=int,
                    default=conf.get_int("NUM_CONVERGENCE", 5),
                    help="how many top fits the convergence spread covers")
    ap.add_argument("beds", nargs="+", help="the input beds the stage was run on")
    a = ap.parse_args(argv)

    cfg = {
        "max_rare_warn": a.max_rare_warn,
        "max_rare_fail": a.max_rare_fail,
        "bulk_bias_warn": a.bulk_bias_warn,
        "rest_over_rare_min": a.rest_over_rare_min,
        "g2_warn": a.g2_warn,
        "g2_fail": a.g2_fail,
    }

    rows = []
    for bed in a.beds:
        base = _bed_base(bed)
        rows.append(row_for(
            os.path.join(a.out, base, a.run),
            os.path.join(a.sfs_out, base, a.fs_run),
            base, cfg, a.rare_max_count, a.num_convergence))

    os.makedirs(a.out, exist_ok=True)
    dest = os.path.join(a.out, "check_demo.%s.tsv" % a.run)
    write_check(dest, rows)
    print("wrote %s" % dest)
    print("  for reference only: these flag fits worth a second look, not an "
          "accept/reject test")
    for r in rows:
        if r["status"] not in ("PASS", "SKIP"):
            print("  %s %s: %s" % (r["status"], r["item"], r["reason"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
