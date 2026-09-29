#!/usr/bin/env python3
"""Checkpoint for the prep stage: how much sequence survived each filter.

prep is four bedtools calls in a fixed order, ported from
submit_exHETquiencmask.sh (subtract heterochromatin, then the quiescent
candidate-neutral set) and filt_nc-unmask.sh (subtract CDS, then intersect the
accessibility mask). run_prep.sh records base pairs after every one of them in
lengths.tsv, and this script turns that into a verdict plus a row that says
where the sequence went, so a surprising result localises to a step instead of
sending you back to the beds.

Verdicts, deliberately few:
    FAIL  nothing survived. The stage produced no test sites at all; anything
          downstream is meaningless, so this blocks a chain.
    WARN  fewer than `min_frac` of the input bp survived (default 5%). Heavy
          attrition is normal for some annotations - a bed sitting largely in
          heterochromatin or largely inside the quiescent control set loses most
          of itself legitimately - so this is advisory, never blocking.
    PASS  otherwise.

No threshold is invented beyond those two. In particular there is no minimum bp:
how much sequence is "enough" is a question about SNP counts, and the sfs stage
answers it against MIN_SNPS on the data itself.
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checkfile import write_check  # noqa: E402

# The filter chain, in the order run_prep.sh runs it. `minus_quiescent` is the
# candidate-neutral set: those sites become the controls, so a test site may
# never also be one.
ORDER = ["input", "minus_HET", "minus_quiescent", "minus_CDS", "after_mask"]

OUTPUT_TEMPLATE = "s.nc_exHET_noq.{base}.bed"


def assess(rows, min_frac=0.05):
    """(status, reasons) from a {stage: bp} mapping. See the module docstring."""
    missing = [s for s in ORDER if s not in rows]
    if missing:
        return "FAIL", ["lengths.tsv is missing step(s): %s" % ", ".join(missing)]

    total = rows["input"]
    final = rows["after_mask"]
    if total <= 0:
        return "FAIL", ["input bed has no sequence (%d bp)" % total]
    if final <= 0:
        return "FAIL", ["no sites left after filtering (0 bp after the mask)"]

    frac = final / total
    if frac < min_frac:
        return "WARN", ["only %.1f%% of input bp survived filtering (< %.1f%%)"
                        % (100 * frac, 100 * min_frac)]
    return "PASS", []


def _bed_base(bed):
    """The bed's base name, with or without a directory and with or without .bed."""
    base = os.path.basename(bed)
    return base[:-4] if base.endswith(".bed") else base


def _read_lengths(path):
    """{stage: bp} from lengths.tsv, or None if it is unreadable."""
    try:
        with open(path) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            return {r["stage"]: int(r["bp"]) for r in rd if r.get("stage")}
    except (OSError, KeyError, TypeError, ValueError):
        return None


def row_for(out_dir, run, bed, min_frac=0.05):
    """One check-file row for one bed."""
    base = _bed_base(bed)
    d = os.path.join(out_dir, base, run)
    lengths = _read_lengths(os.path.join(d, "lengths.tsv"))
    row = {"item": base}

    if lengths is None:
        row["status"] = "FAIL"
        row["reason"] = "no readable lengths.tsv in %s" % d
        return row

    status, reasons = assess(lengths, min_frac=min_frac)
    bed_out = os.path.join(d, OUTPUT_TEMPLATE.format(base=base))
    if not os.path.isfile(bed_out):
        status = "FAIL"
        reasons = reasons + ["output bed is missing: %s" % bed_out]

    row["status"] = status
    row["reason"] = "; ".join(reasons) or "NA"
    for step in ORDER:
        row[step] = lengths.get(step, "NA")
    total = lengths.get("input", 0)
    for step in ORDER[1:]:
        bp = lengths.get(step)
        name = "frac_retained" if step == "after_mask" else "frac_%s" % step
        row[name] = round(bp / total, 4) if (total and bp is not None) else "NA"
    return row


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True,
                    help="the prep output root, i.e. $OUT/prep")
    ap.add_argument("--run", required=True, help="PREP_RUN, the run key")
    ap.add_argument("--min-frac", type=float, default=0.05,
                    help="WARN below this fraction of input bp surviving (default 0.05)")
    ap.add_argument("beds", nargs="+", help="the input beds prep was run on")
    a = ap.parse_args(argv)

    rows = [row_for(a.out, a.run, bed, min_frac=a.min_frac) for bed in a.beds]
    dest = os.path.join(a.out, "check_prep.%s.tsv" % a.run)
    write_check(dest, rows)
    print("wrote %s" % dest)
    for r in rows:
        if r["status"] != "PASS":
            print("  %s %s: %s" % (r["status"], r["item"], r["reason"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
