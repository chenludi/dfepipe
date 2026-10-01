#!/usr/bin/env python3
"""Checkpoint for the summary stage: is the table complete and self-consistent?

The earlier stages already judged their own output, so this does not re-measure
the fit. It asks whether the table can be read at face value:

    FAIL  sumdfe.tsv is missing, or a bed asked for has no row
    FAIL  a row's field count does not match the header
    WARN  ll_neu or a pdf's likelihood is NA - that column came from nothing
    WARN  a pdf's DFE bins do not sum to 1 within BIN_SUM_TOL
    WARN  a pdf's likelihood is no better than ll_neu: selection buys nothing
    WARN  u, l or ul is NA, so the mutational target is unknown
    PASS  one row per bed, every value present, bins summing to 1

Everything but a structurally unusable table is a WARN, for the same reason as
the other stages: a thin result is still a result, and the user decides.

    python3 check_summary.py --out <OUT>/summary --run <SUMMARY_RUN> \
        --pdfs "neu gamma neugamma" <bed> [bed...]
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checkfile import write_check  # noqa: E402
from conf import conf  # noqa: E402

TABLE_NAME = "sumdfe.tsv"
# Keyed by the run, as every other stage's check file is: two chains with
# different SUMMARY_RUNs must not overwrite each other.
CHECK_TEMPLATE = "check_summary.%s.tsv"
BIN_PREFIX = "bin_"
# The six bins are a partition of the DFE, so they sum to 1. The slack is for
# the CDF arithmetic, not for a real shortfall.
BIN_SUM_TOL = 1e-6


def _float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def assess_row(row, fit_pdfs, bin_sum_tol=BIN_SUM_TOL):
    """(status, reasons) for one bed's row."""
    reasons = []

    if _float(row.get("ll_neu")) is None:
        reasons.append(
            "ll_neu is NA: run_fitneu.sh has not run, or its log is unreadable. "
            "Every model comparison in the summary needs it")

    for pdf in fit_pdfs:
        ll = _float(row.get("ll_%s" % pdf))
        if ll is None:
            reasons.append(
                "ll_%s is NA: large_ll.txt is missing or holds no fit, so this "
                "pdf contributes nothing to the table" % pdf)
            continue
        ll_neu = _float(row.get("ll_neu"))
        if ll_neu is not None and ll <= ll_neu:
            reasons.append(
                "ll_%s %.1f is no better than ll_neu %.1f: selection buys this "
                "model nothing" % (pdf, ll, ll_neu))

        prefix = "%s_%s" % (pdf, BIN_PREFIX)
        bins = [_float(v) for k, v in row.items() if k.startswith(prefix)]
        if not bins:
            reasons.append(
                "%s has no DFE bins: dfe_*_%s.csv is missing" % (pdf, pdf))
        elif any(b is None for b in bins):
            reasons.append("%s has an unreadable DFE bin" % pdf)
        else:
            total = sum(bins)
            if abs(total - 1.0) > bin_sum_tol:
                reasons.append(
                    "%s DFE bins sum to %.6f, not 1 (tolerance %g): the "
                    "discretisation lost mass" % (pdf, total, bin_sum_tol))

    for name in ("u", "l", "ul"):
        if _float(row.get(name)) is None:
            reasons.append(
                "%s is NA: the fsmut table was unreadable, so the mutational "
                "target is unknown and the DFE cannot be weighted by it" % name)

    return ("WARN" if reasons else "PASS"), reasons


def read_table(path):
    """(columns, {outputbase: row}, ragged) from sumdfe.tsv."""
    try:
        with open(path) as fh:
            rd = csv.reader(fh, delimiter="\t")
            columns = next(rd)
            rows, ragged = {}, []
            for raw in rd:
                if not raw:
                    continue
                if len(raw) != len(columns):
                    ragged.append((raw[0] if raw else "?", len(raw)))
                    continue
                row = dict(zip(columns, raw))
                rows[row.get("outputbase", "?")] = row
    except (OSError, StopIteration):
        return None, None, None
    return columns, rows, ragged


def short_base(base, ctrl_tag):
    """run_summary.sh's outputbase, so a bed argument can find its row."""
    stem = "s.nc_exHET_noq.%s" % base
    stem = stem.replace("s.nc_", "", 1)
    stem = stem.replace("nc_", "", 1)
    return "%s_%s" % (stem, ctrl_tag)


def _bed_base(bed):
    base = os.path.basename(bed)
    return base[:-4] if base.endswith(".bed") else base


def rows_for(out_dir, run, beds, pdfs, ctrl_tag, bin_sum_tol=BIN_SUM_TOL):
    path = os.path.join(out_dir, run, TABLE_NAME)
    fit_pdfs = [p for p in pdfs if p != "neu"]

    columns, table, ragged = read_table(path)
    if columns is None:
        return [{"item": "sumdfe", "status": "FAIL", "n_rows": "NA",
                 "reason": "%s is missing or unreadable: run "
                           "scripts/summary/run_summary.sh" % path}]

    out = []
    for name, n in ragged or []:
        out.append({"item": name, "status": "FAIL", "n_rows": len(table),
                    "reason": "the row has %d fields but the header has %d, so "
                              "its columns do not mean what they say"
                              % (n, len(columns))})
    for bed in beds:
        key = short_base(_bed_base(bed), ctrl_tag)
        row = table.get(key)
        if row is None:
            out.append({"item": key, "status": "FAIL", "n_rows": len(table),
                        "reason": "no row in %s. The table has: %s"
                                  % (TABLE_NAME, ", ".join(sorted(table)) or "nothing")})
            continue
        status, reasons = assess_row(row, fit_pdfs, bin_sum_tol=bin_sum_tol)
        out.append({"item": key, "status": status, "n_rows": len(table),
                    "reason": "; ".join(reasons) or "NA"})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="<OUT>/summary")
    ap.add_argument("--run", required=True, help="SUMMARY_RUN")
    ap.add_argument("--pdfs", required=True, help="PDFS, space separated")
    ap.add_argument("--ctrl-tag", default=conf.get("CTRL_TAG", "qui2"))
    ap.add_argument("--bin-sum-tol", type=float, default=BIN_SUM_TOL)
    ap.add_argument("beds", nargs="+")
    a = ap.parse_args(argv)

    rows = rows_for(a.out, a.run, a.beds, a.pdfs.split(), a.ctrl_tag,
                    bin_sum_tol=a.bin_sum_tol)
    path = os.path.join(a.out, CHECK_TEMPLATE % a.run)
    write_check(path, rows)
    print("wrote %s" % path)
    print("  for reference only: these flag a table worth a second look, not "
          "an accept/reject test")
    for r in rows:
        if r["status"] != "PASS":
            print("  %s %s: %s" % (r["status"], r["item"], r["reason"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
