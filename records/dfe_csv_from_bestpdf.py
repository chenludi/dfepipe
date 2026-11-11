#!/usr/bin/env python3
"""2026-10-27  bestpdf_sumdfe.tsv -> one csv per group for plot_DFE_mutload2_fromcsv.py.

That plotter reads each csv with index_col=0: the index is the legend label and
the remaining columns are the six |s| bins. This writes the best-fit discrete
DFE, i.e. the delta_50_* columns, which summary_getbestpdf.py has already filled
from whichever pdf the delta-50 rule chose (gamma, neugamma, or neutral). The
chosen pdf is printed so the figure can be audited.

  python records/dfe_csv_from_bestpdf.py -i bestpdf_sumdfe.tsv -o <dir> \
      --group exHET_noq.GRCh38-ELS_x_Prom_qui2 --label "ELS_x_Prom"

The csvs this wrote for the 2026-10-27 figures are still on the cluster, in
$OUT/figs/261027/figS2_csv/ and oldnew_csv/. Only the figures drawn from them
were moved locally, to
DFE_non-coding/a_results_weekly/mv_261027_dfepipe_els_v4_phylop/figs/.
"""
import argparse
import csv
import os
import sys

BINS = ["0.0 <= |s| < 1e-05", "1e-05 <= |s| < 0.0001", "0.0001 <= |s| < 0.001",
        "0.001 <= |s| < 0.01", "0.01 <= |s| < 0.5", "0.5 <= |s| < inf"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", nargs="+", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--group", action="append", required=True,
                    help="outputbase value to extract; repeatable")
    ap.add_argument("--label", action="append", required=True,
                    help="legend label for the matching --group; repeatable")
    a = ap.parse_args()
    if len(a.group) != len(a.label):
        sys.exit("--group and --label must be given the same number of times")

    table = {}
    for p in a.input:
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                table[r["outputbase"]] = r

    os.makedirs(a.outdir, exist_ok=True)
    for g, lab in zip(a.group, a.label):
        if g not in table:
            sys.exit("group %r not found in any input (have %d groups)"
                     % (g, len(table)))
        r = table[g]
        vals = [r["delta_50_" + b] for b in BINS]
        out = os.path.join(a.outdir, "dfe_%s.csv" % lab)
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow([""] + BINS)
            w.writerow([lab] + vals)
        print("  %-26s <- %-52s best=%s" % (lab, g, r["pdf_plotting_delta_50"]))


if __name__ == "__main__":
    main()
