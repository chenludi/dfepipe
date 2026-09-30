#!/usr/bin/env python3
"""Checkpoint for the dfe stage: does the fitted DFE describe the test SFS?

FOR REFERENCE ONLY, as everywhere in this pipeline: a WARN means "unusual
compared with 137 previous fits", not "wrong". The figure remains the arbiter.

FAIL is reserved for structural problems, which do block a chain:
    no fits at all                     no *.InferDFE.bestfits in the run
                                       directory, so nothing was inferred
    large_ll.txt missing/unparseable   the collected fits the summary stage and
                                       the discretised DFE both read
    the spectra could not be read      the fit was never measured
    no usable bins, or the model gives no probability where variants are
                                       observed (from residuals.classify)

WARN is advisory: the calibrated residual and deviance thresholds, which never
FAIL however extreme, plus three remarks about the inference itself - too few
starting points converged, the top fits disagree, or the DFE does no better than
the neutral null. None of them blocks anything.

The residual metrics are not recomputed here: run_bestdfe.sh (and run_fitneu.sh
for the null) writes the best fit's folded model spectrum beside the fits as
`model_sfs_rank1.fs`, and this script compares it with the test spectrum through
scripts/checks/residuals.py.
"""
import argparse
import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checkfile import write_check  # noqa: E402
from residuals import classify, residual_metrics  # noqa: E402

# Top fits within this many log-likelihood units count as converged on one
# optimum. One ll unit is the conventional tolerance, as in check_demo.
CONVERGENCE_LL_SPREAD = 1.0

# Written by run_bestdfe.sh (rank 1 = the maximum-likelihood fit) and by
# run_fitneu.sh (the demography-only expectation), folded either way.
MODEL_FS_NAME = "model_sfs_rank1.fs"
LARGE_LL_NAME = "large_ll.txt"
# The test spectrum the DFE was fitted to. The tail of the name carries the
# bed's short name and CTRL_TAG, so it is globbed rather than rebuilt: the file
# the sfs stage wrote is the authority.
TEST_FS_GLOB = "has_qui.*.chrs_1t22.fs"
# run_dfe.sh's per-starting-point log, and the packs that define which starting
# points exist at all.
OUTCOMES_NAME = "fit_outcomes.tsv"
PACK_GLOB = "dadi-cli_dfe.*.sh"
BESTFITS_SUFFIX = ".InferDFE.bestfits"
OUTCOME_KINDS = ("converged", "no_likelihood", "timeout", "not_run")
# pipeline.conf's DFE_MAX_NOTRUN_FRAC / DFE_MIN_CONVERGED_FRAC /
# DFE_MIN_CONVERGED_N. Repeated here only as the CLI defaults; the caller passes
# the configured values through --max-notrun-frac and friends.
DFE_MAX_NOTRUN_FRAC = 0.20
DFE_MIN_CONVERGED_FRAC = 0.20
DFE_MIN_CONVERGED_N = 10
# run_fitneu.sh keeps the historical log name; its `ll: ` line is the null.
NEU_LOG_GLOB = "*.dadi_cal_ll_neu_0723.log"


# ------------------------------------------------------------------- verdict --

def assess_dfe(metrics, top_lls, cfg, ll_dfe=None, ll_neu=None, n_fits=None,
               n_converged_block=None, num_convergence=5,
               ll_spread_warn=CONVERGENCE_LL_SPREAD, outcomes=None,
               max_notrun_frac=DFE_MAX_NOTRUN_FRAC,
               min_converged_frac=DFE_MIN_CONVERGED_FRAC,
               min_converged_n=DFE_MIN_CONVERGED_N):
    """(status, reasons) for one bed and one pdf. See the module docstring.

    `metrics` is a residual_metrics() mapping, or None when the spectra could
    not be read. `top_lls` are the best log-likelihoods across starting points,
    `ll_dfe` the winning one (None when large_ll.txt is unusable), `ll_neu` the
    neutral null (None when run_fitneu.sh has not run) and `cfg` the threshold
    mapping residuals.classify takes.

    `n_fits` and `n_converged_block` are skipped when None, which is how the
    `neu` row - a single closed-form likelihood, not an optimisation - is
    assessed. `n_converged_block` counts the .bestfits files that carry a
    "# Converged results" block, which is stricter than count_outcomes()'s
    `n_converged` (a non-empty .bestfits): about half the files have no such
    block.

    `outcomes` is a count_outcomes() mapping, or None where there are no
    starting points to account for. no_likelihood and timeout are reported by
    row_for and never judged here: a bad starting point that fails fast is part
    of the design.
    """
    fails = []
    if n_fits is not None and n_fits < 1:
        fails.append(
            "no fits produced: no *.InferDFE.bestfits in the run directory. "
            "Either run_dfe.sh has not run, or every command in every pack "
            "failed - check the SGE logs and the pack files")
    if outcomes and outcomes["n_starting_points"] and not outcomes["n_converged"]:
        fails.append(
            "not one of the %d starting points converged: there is no best fit "
            "to take. %d produced no likelihood, %d timed out, %d never ran"
            % (outcomes["n_starting_points"], outcomes["n_no_likelihood"],
               outcomes["n_timeout"], outcomes["n_not_run"]))
    if ll_dfe is None:
        fails.append(
            "%s is missing or holds no fit: the summary stage and the "
            "discretised DFE both read it. Re-run scripts/dfe/run_bestdfe.sh"
            % LARGE_LL_NAME)
    if metrics is None:
        fails.append(
            "the model and data spectra could not be read, so the fit was not "
            "measured: check that %s sits beside the fits" % MODEL_FS_NAME)
    if fails:
        return "FAIL", fails

    # Residual and deviance thresholds only ever WARN; classify FAILs only on
    # structural problems. A structural FAIL is the whole verdict: the remarks
    # below are about an inference that at least produced something measurable.
    status, reasons = classify(metrics, cfg)
    if status == "FAIL":
        return status, reasons

    def warn(msg):
        reasons.append(msg)

    lls = [float(x) for x in top_lls]
    if len(lls) >= 2:
        spread = max(lls) - min(lls)
        if spread > ll_spread_warn:
            warn("the top %d fits span %.1f log-likelihood units (> %.1f): they "
                 "did not converge on one optimum, so the best fit may be local"
                 % (len(lls), spread, ll_spread_warn))
    if n_converged_block is not None and n_converged_block < num_convergence:
        warn("only %d starting point(s) converged (want %d): the optimum may be "
             "local" % (n_converged_block, num_convergence))
    if outcomes and outcomes["n_starting_points"]:
        n = outcomes["n_starting_points"]
        if outcomes["frac_not_run"] > max_notrun_frac:
            warn("%d of %d starting points never ran (%.0f%% > %.0f%%): the "
                 "array task hit its h_rt, so the parameter space was "
                 "under-explored. RESUME=1 finishes them"
                 % (outcomes["n_not_run"], n, 100 * outcomes["frac_not_run"],
                    100 * max_notrun_frac))
        if (outcomes["n_converged"] < min_converged_n
                or outcomes["frac_converged"] < min_converged_frac):
            warn("only %d of %d starting points converged (%.0f%%; want at "
                 "least %d and %.0f%%): the method looks unstable on this bed"
                 % (outcomes["n_converged"], n, 100 * outcomes["frac_converged"],
                    min_converged_n, 100 * min_converged_frac))
    if ll_neu is not None and ll_dfe <= float(ll_neu):
        warn("the DFE log-likelihood %.1f is no better than the neutral null "
             "%.1f: selection buys this model nothing"
             % (ll_dfe, float(ll_neu)))

    if reasons:
        return "WARN", reasons
    return status, reasons


# ------------------------------------------------------------------- readers --

def read_large_ll(path):
    """(header, rows) from large_ll.txt, or (None, None) if it is unusable.

    Rows come back in file order, which run_bestdfe.sh writes sorted by
    descending log-likelihood, so rows[0] is the winning fit. Numeric fields are
    converted; one that will not convert is left as the raw string.
    """
    try:
        with open(path) as fh:
            rd = csv.reader(fh, delimiter="\t")
            header = next(rd)
            rows = []
            for raw in rd:
                if not raw:
                    continue
                row = dict(zip(header, raw))
                for key, value in list(row.items()):
                    try:
                        row[key] = float(value)
                    except ValueError:
                        pass
                rows.append(row)
    except (OSError, StopIteration):
        return None, None
    return header, rows


def read_neu_ll(directory):
    """The neutral log-likelihood from run_fitneu.sh's log, or None.

    The log is the historical `<test fs base>.dadi_cal_ll_neu_0723.log`, and its
    `ll: ` line is the same one the summary stage greps.
    """
    hits = sorted(glob.glob(os.path.join(str(directory), NEU_LOG_GLOB)))
    if not hits:
        return None
    try:
        with open(hits[0]) as fh:
            for line in fh:
                if line.startswith("ll:"):
                    return float(line.split(":", 1)[1].strip())
    except (OSError, ValueError):
        return None
    return None


def read_fit_outcomes(path):
    """marker -> outcome from run_dfe.sh's fit_outcomes.tsv, or {} if absent.

    The log is APPEND-ONLY: a RESUME re-run adds rows rather than rewriting
    them, so a marker can appear more than once and the LAST row wins. Header
    lines, blank lines and unrecognised outcomes are skipped, which is what
    makes a header repeated by a second array task harmless.
    """
    out = {}
    try:
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                marker, outcome = parts[0].strip(), parts[1].strip()
                if not marker or marker == "marker":
                    continue
                if outcome in OUTCOME_KINDS:
                    out[marker] = outcome
    except OSError:
        return {}
    return out


def read_pack_markers(dfe_dir):
    """Every starting point the packs define, as markers, in pack order.

    The marker is the basename of --output-prefix, which run_dfe.sh writes as
    the starting point's index, so it names the .bestfits and the .opts files
    too. The packs are the authority on how many starting points exist: a
    marker in no pack was never meant to run.
    """
    markers = []
    for pack in sorted(glob.glob(os.path.join(str(dfe_dir), PACK_GLOB))):
        try:
            with open(pack) as fh:
                lines = fh.read().splitlines()
        except OSError:
            continue
        for line in lines:
            fields = line.split()
            for k in range(len(fields) - 1):
                if fields[k] == "--output-prefix":
                    markers.append(os.path.basename(fields[k + 1]))
    return list(dict.fromkeys(markers))


def count_outcomes(dfe_dir):
    """The four per-starting-point outcome counts for one run directory.

    A non-empty .bestfits is `converged` whatever the log says - the file
    system is the authority, and it is also the only evidence a run made before
    this log existed leaves behind. Otherwise the log's last row for that
    marker decides. A marker in the packs with neither is `not_run`: the array
    task never reached it.
    """
    dfe_dir = str(dfe_dir)
    markers = read_pack_markers(dfe_dir)
    logged = read_fit_outcomes(os.path.join(dfe_dir, OUTCOMES_NAME))
    converged = set()
    for path in glob.glob(os.path.join(dfe_dir, "*" + BESTFITS_SUFFIX)):
        try:
            if os.path.getsize(path) > 0:
                converged.add(os.path.basename(path)[:-len(BESTFITS_SUFFIX)])
        except OSError:
            continue
    if not markers:
        markers = sorted(set(logged) | converged)

    counts = dict.fromkeys(("n_" + k for k in OUTCOME_KINDS), 0)
    for marker in markers:
        if marker in converged:
            counts["n_converged"] += 1
        elif logged.get(marker) in ("no_likelihood", "timeout"):
            counts["n_" + logged[marker]] += 1
        elif marker in logged:
            # Logged converged, but the .bestfits is gone: it produced nothing
            # readable now, which is what no_likelihood means to every reader.
            counts["n_no_likelihood"] += 1
        else:
            counts["n_not_run"] += 1

    n = len(markers)
    counts["n_starting_points"] = n
    for kind in OUTCOME_KINDS:
        counts["frac_" + kind] = counts["n_" + kind] / n if n else 0.0
    return counts


def read_bestfits_lls(path):
    """Every converged log-likelihood in one dadi-cli .bestfits file.

    A file with no "# Converged results" block did not converge and contributes
    nothing, which is the rule plot_converge_dfe_0430.py:282-321 applies.
    """
    try:
        with open(path) as fh:
            lines = fh.read().splitlines()
    except OSError:
        return []
    out, found = [], False
    for line in lines:
        stripped = line.strip()
        if stripped == "# Converged results":
            found = True
            continue
        if not found:
            continue
        if stripped == "#":              # end of the converged block
            break
        if stripped.startswith("#") or not stripped:
            continue
        try:
            out.append(float(stripped.split("\t")[0]))
        except ValueError:
            break
    return out


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

def row_for(dfe_dir, test_fs_dir, neu_dir, base, pdf, cfg, rare_max_count,
            num_convergence, max_notrun_frac=DFE_MAX_NOTRUN_FRAC,
            min_converged_frac=DFE_MIN_CONVERGED_FRAC,
            min_converged_n=DFE_MIN_CONVERGED_N):
    """The single check row for one bed and one pdf."""
    is_neu = pdf == "neu"

    outcomes = None
    n_fits = n_converged_block = None
    top = []
    header, rows = (None, None)
    ll_dfe = None
    if is_neu:
        # The null has no optimisation and no large_ll.txt; its ll IS the fit.
        ll_dfe = read_neu_ll(dfe_dir)
        ll_neu = None
        best_file = "NA"
    else:
        fits = sorted(glob.glob(os.path.join(dfe_dir, "*.InferDFE.bestfits")))
        n_fits = len(fits)
        per_file = [read_bestfits_lls(f) for f in fits]
        n_converged_block = sum(1 for lls in per_file if lls)
        all_lls = sorted((max(lls) for lls in per_file if lls), reverse=True)
        top = all_lls[:num_convergence]
        header, rows = read_large_ll(os.path.join(dfe_dir, LARGE_LL_NAME))
        best = rows[0] if rows else None
        if best is not None and isinstance(best.get("Log(likelihood)"), float):
            ll_dfe = best["Log(likelihood)"]
        best_file = best.get("filename", "NA") if best else "NA"
        ll_neu = read_neu_ll(neu_dir)
        outcomes = count_outcomes(dfe_dir)

    test = sorted(glob.glob(os.path.join(test_fs_dir, TEST_FS_GLOB)))
    data = read_fs(test[0]) if test else None
    model = read_fs(os.path.join(dfe_dir, MODEL_FS_NAME))
    metrics = None
    if data is not None and model is not None and len(data) == len(model):
        metrics = residual_metrics(model, data, rare_max_count=rare_max_count)

    status, reasons = assess_dfe(
        metrics, top, cfg, ll_dfe=ll_dfe, ll_neu=ll_neu, n_fits=n_fits,
        n_converged_block=n_converged_block, num_convergence=num_convergence,
        outcomes=outcomes, max_notrun_frac=max_notrun_frac,
        min_converged_frac=min_converged_frac,
        min_converged_n=min_converged_n)

    row = {
        "item": base,
        "status": status,
        "pdf": pdf,
        "reason": "; ".join(reasons) or "NA",
        "n_fits": n_fits if n_fits is not None else "NA",
        "n_converged_block": (n_converged_block
                              if n_converged_block is not None else "NA"),
        "best_file": best_file,
        "ll": ll_dfe if ll_dfe is not None else "NA",
        "ll_neu": ll_neu if ll_neu is not None else "NA",
        "delta_ll_vs_neu": (round(ll_dfe - ll_neu, 4)
                            if ll_dfe is not None and ll_neu is not None else "NA"),
        "ll_spread_top": round(max(top) - min(top), 6) if len(top) >= 2 else "NA",
    }
    # Every starting point accounted for, in the four outcomes of
    # fit_outcomes.tsv. no_likelihood and timeout are reported, never judged.
    row["n_starting_points"] = outcomes["n_starting_points"] if outcomes else "NA"
    for kind in OUTCOME_KINDS:
        row["n_" + kind] = outcomes["n_" + kind] if outcomes else "NA"
        row["frac_" + kind] = (round(outcomes["frac_" + kind], 4)
                               if outcomes else "NA")
    # The fitted parameters, whichever this pdf has, straight from large_ll.txt.
    if header and rows:
        for key in header:
            if key in ("filename", "Log(likelihood)"):
                continue
            value = rows[0].get(key, "NA")
            row[key] = value
    for key in ("max_abs_rare", "max_abs_rest", "mean_rest", "mean_rest_norm",
                "g2", "g2_per_snp", "n_snps", "n_bins_used",
                "n_masked_data_zero", "n_masked_model_zero"):
        value = metrics.get(key, "NA") if metrics else "NA"
        row[key] = round(value, 6) if isinstance(value, float) else value
    return row


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True, help="the dfe output root, $OUT/dfe")
    ap.add_argument("--run", required=True,
                    help="the run directory name, <CACHE_RUN>-<pdf>-<DFE_RUN>")
    ap.add_argument("--pdf", required=True, help="the pdf this run fitted")
    ap.add_argument("--neu-run", default=None,
                    help="the neutral run directory name, for the null ll")
    ap.add_argument("--sfs-out", required=True, help="the sfs output root, $OUT/sfs")
    ap.add_argument("--fs-run", required=True, help="FS_RUN")
    ap.add_argument("--rare-max-count", type=int, default=4)
    ap.add_argument("--max-rare-warn", type=float, default=25.0)
    ap.add_argument("--max-rare-fail", type=float, default=40.0)
    ap.add_argument("--bulk-bias-warn", type=float, default=2.3)
    ap.add_argument("--rest-over-rare-min", type=float, default=5.0)
    ap.add_argument("--g2-warn", type=float, default=0.008)
    ap.add_argument("--g2-fail", type=float, default=0.035)
    ap.add_argument("--num-convergence", type=int, default=5,
                    help="how many starting points are wanted to converge")
    ap.add_argument("--max-notrun-frac", type=float,
                    default=DFE_MAX_NOTRUN_FRAC,
                    help="DFE_MAX_NOTRUN_FRAC: warn above this not_run share")
    ap.add_argument("--min-converged-frac", type=float,
                    default=DFE_MIN_CONVERGED_FRAC,
                    help="DFE_MIN_CONVERGED_FRAC: warn below this converged share")
    ap.add_argument("--min-converged-n", type=int, default=DFE_MIN_CONVERGED_N,
                    help="DFE_MIN_CONVERGED_N: warn below this many converged")
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
            os.path.join(a.out, base, a.neu_run) if a.neu_run else "",
            base, a.pdf, cfg, a.rare_max_count, a.num_convergence,
            max_notrun_frac=a.max_notrun_frac,
            min_converged_frac=a.min_converged_frac,
            min_converged_n=a.min_converged_n))

    os.makedirs(a.out, exist_ok=True)
    dest = os.path.join(a.out, "check_dfe.%s.tsv" % a.run)
    write_check(dest, rows)
    print("wrote %s" % dest)
    print("  for reference only: these flag fits worth a second look, not an "
          "accept/reject test")
    for r in rows:
        if r["status"] not in ("PASS", "SKIP"):
            print("  %s %s (%s): %s" % (r["status"], r["item"], r["pdf"],
                                        r["reason"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
