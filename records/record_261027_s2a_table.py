#!/usr/bin/env python3
"""2026-10-27  bestpdf_sumdfe.tsv -> the column layout of Data_S2_DFE.xlsx, S2a.

Reads one or more bestpdf_sumdfe.tsv files and writes a tsv whose columns are
the 47 of the "DFEs of ELS" section of S2a (sheet S2a.enhancers, header row 38),
in that order, ready to paste as new rows under that section.

Nothing is read from or written to the xlsx. S2a's own column names are
hard-coded below, taken from row 38 of the sheet on 2026-10-27.

Three notes on the mapping:

  * S2a's `pdf` and `pdf.1` hold the model name, not a file name, so they are
    written as the literal "gamma" and "neugamma". bestpdf's gamma_filename /
    neugamma_filename (e.g. "129.InferDFE.bestfits") have no S2a column.
  * Column names are not uniform across S2a's sections: the ELS section ends
    `ul`, `num_nonneutral`, while the binned sections use `num_non-neutral` and
    add `cumu_num-mutation(ul)`, `cumu_non-neutral`. --section picks which.
  * The cumu_* columns are cumulative over ordered rows of a whole binned
    series, so they are left empty for a handful of added rows; computing them
    would need the full series in its plotting order.

bestpdf columns with no S2a counterpart are dropped and listed with --verbose:
gamma_theta, gamma_sdadi, neugamma_theta, neugamma_sdadi, neugamma_Ne, File,
l_inputbed.

Output. This runs locally, not on the cluster. The table it wrote on
2026-10-27 was moved on 2026-11-10 to
DFE_non-coding/a_results_weekly/mv_261027_dfepipe_els_v4_phylop/S2a_rows_ELS_261027.tsv,
together with the bestpdf_sumdfe.*.tsv it was built from. It was in
dfepipe/exchange/ before that; that directory no longer exists.

  python records/record_261027_s2a_table.py -i <bestpdf_sumdfe.tsv> [...] -o out.tsv
  python records/record_261027_s2a_table.py -i a.tsv b.tsv -o out.tsv --only ELS --verbose
"""
import argparse
import csv
import sys

BINS = ["0.0 <= |s| < 1e-05", "1e-05 <= |s| < 0.0001", "0.0001 <= |s| < 0.001",
        "0.001 <= |s| < 0.01", "0.01 <= |s| < 0.5", "0.5 <= |s| < inf"]

# (S2a column name, bestpdf_sumdfe.tsv column name or a ("literal", value) pair)
MAP = [
    ("groups", "outputbase"),
    ("ll_neu", "ll_neu"),
    ("pdf", ("literal", "gamma")),
    ("ll_gamma", "ll_gamma"),
    ("shape", "gamma_shape"),
    ("scale", "gamma_scale"),
    ("scale_div2Na", "gamma_scale_div2Na"),
    ("Ne", "gamma_Ne"),
] + [(b, "gamma_bin_" + b) for b in BINS] + [
    ("pdf.1", ("literal", "neugamma")),
    ("ll_neugamma", "ll_neugamma"),
    ("shape.1", "neugamma_shape"),
    ("scale.1", "neugamma_scale"),
    ("pneu", "neugamma_pneu"),
    ("scale_div2Na.1", "neugamma_scale_div2Na"),
] + [(b + ".1", "neugamma_bin_" + b) for b in BINS] + [
    ("u", "u"),
    ("l", "l"),
    ("ll_gamma-ll_neu", "ll_gamma-ll_neu"),
    ("ll_neugamma-ll_gamma", "ll_neugamma-ll_gamma"),
    ("ll_neugamma-ll_neu", "ll_neugamma-ll_neu"),
    ("pdf_best_chi_5%", "pdf_best_chi_5%"),
    ("pdf_best_chi_0.10%", "pdf_best_chi_0.10%"),
    ("pdf_best_delta_10", "pdf_best_delta_10"),
    ("pdf_best_delta_30", "pdf_best_delta_30"),
    ("pdf_best_delta_50", "pdf_best_delta_50"),
    ("pdf_best_delta_100", "pdf_best_delta_100"),
    ("pdf_best_allconfidence", "pdf_best_allconfidence"),
    ("pdf_plotting_delta_50", "pdf_plotting_delta_50"),
] + [("delta_50_" + b, "delta_50_" + b) for b in BINS] + [
    ("ul", "ul"),
]

# the trailing columns, which differ between S2a's sections
TAIL = {
    "ELS":    [("num_nonneutral", "num_nonneutral")],
    "binned": [("num_non-neutral", "num_nonneutral"),
               ("cumu_num-mutation(ul)", None),
               ("cumu_non-neutral", None)],
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", nargs="+", required=True,
                    help="one or more bestpdf_sumdfe.tsv files")
    ap.add_argument("-o", "--output", required=True, help="output tsv")
    ap.add_argument("--section", choices=sorted(TAIL), default="ELS",
                    help="which S2a section's trailing columns to use (default ELS)")
    ap.add_argument("--only", default=None,
                    help="keep only rows whose groups value contains this substring")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()

    cols = MAP + TAIL[a.section]
    rows, seen, dropped = [], set(), None

    for path in a.input:
        with open(path, newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            if dropped is None:
                wanted = {src for _, src in MAP + TAIL[a.section]
                          if isinstance(src, str)}
                dropped = [c for c in rd.fieldnames if c not in wanted]
            for r in rd:
                if a.only and a.only not in r["outputbase"]:
                    continue
                if r["outputbase"] in seen:
                    print("duplicate groups value, keeping the first: %s"
                          % r["outputbase"], file=sys.stderr)
                    continue
                seen.add(r["outputbase"])
                out = []
                for name, src in cols:
                    if src is None:
                        out.append("")                       # cumu_*, see docstring
                    elif isinstance(src, tuple):
                        out.append(src[1])                   # literal pdf name
                    else:
                        if src not in r:
                            raise SystemExit(
                                "column %r missing from %s" % (src, path))
                        out.append(r[src])
                rows.append(out)

    with open(a.output, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow([name for name, _ in cols])
        w.writerows(rows)

    print("wrote %s: %d row(s), %d columns (S2a section %r)"
          % (a.output, len(rows), len(cols), a.section))
    if a.verbose:
        print("dropped, no S2a column: %s" % ", ".join(dropped))
        for r in rows:
            print("  %-54s best=%s" % (r[0], r[cols.index(("pdf_plotting_delta_50",
                                                           "pdf_plotting_delta_50"))]))


if __name__ == "__main__":
    main()
