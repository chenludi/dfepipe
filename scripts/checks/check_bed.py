#!/usr/bin/env python3
"""Preflight on an input bed, before any stage spends time on it.

This exists because the two most likely user errors both fail late and quietly:

  * a bed named `1 2 3` instead of `chr1 chr2 chr3` passes prep unharmed, then
    produces nothing in all 22 sfs tasks, and only surfaces at concat - an hour
    or more later, as "cannot concatenate".
  * coordinates from the wrong assembly produce plausible output that is simply
    wrong, with nothing anywhere to say so.

So this is a gate on the input, not a judgement of the result:

    FAIL  unreadable, empty, fewer than 3 columns, a non-integer or reversed
          coordinate, or no chromosome the run will analyse
    WARN  contigs in the bed that CHROMS excludes - those bases are silently
          dropped at concat and the reported l will not be the bed's l
    WARN  intervals longer than 2 x CTRL_FLANK: their control window falls
          entirely inside them, where controls cannot be, so they can never
          find one
    WARN  unsorted, or overlapping intervals
    WARN  a coordinate far beyond that chromosome's last masked base - a corrupt
          bed, or one from a different genome
    PASS   none of the above

The genome build is NOT checked, and the summary says so every time rather than
letting a PASS imply it was. Coordinates are compared against the mask only to
catch something grossly out of range - a corrupt bed, or one from a different
genome. Two builds of the same genome differ by far less than that, so no
coordinate test here can tell them apart.

    python3 check_bed.py --chroms "1 2 3 ..." [--ctrl-flank 25000]
                        [--mask-bed <path>] [--out <dir> --run <key>] <bed>...

"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checkfile import write_check  # noqa: E402
from conf import conf  # noqa: E402

CHECK_TEMPLATE = "check_bed.%s.tsv"
# A coordinate further than this past its chromosome's last masked base is out of
# range: a corrupt file, or a bed from a different genome altogether. The figure
# is a judgement, not a measurement - the mask does not reach the end of a
# chromosome, so some slack is needed, and 1 Mb is loose enough never to fire on
# a normal bed.
#
# This is NOT an assembly check and is not presented as one. Two builds of the
# same genome differ by far less than this, so it cannot tell them apart; that
# the build is unverifiable here is said once, in the summary line, rather than
# hedged on every bed.
OUT_OF_RANGE_SLACK = 1_000_000


def normalise(chrom):
    """CHROMS is bare numbers in pipeline.conf; every bed is chr-prefixed."""
    c = str(chrom).strip()
    return c if c.startswith("chr") else "chr" + c


def read_bed(path, window=0):
    """(per-chrom stats, problems). Streams, so a 100 MB bed is fine.

    `window` is 2 x CTRL_FLANK; intervals longer than it are counted as they go,
    because how many there are and how much sequence they hold is what decides
    whether the control window needs rethinking.
    """
    stats = {}
    problems = []
    n = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip() or line[0] == "#" or line.startswith("track"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                problems.append("line %d has %d column(s), needs at least 3 "
                                "(chrom, start, end), tab separated"
                                % (lineno, len(f)))
                if len(problems) > 5:
                    break
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                problems.append("line %d: start/end are not integers: %r, %r"
                                % (lineno, f[1], f[2]))
                if len(problems) > 5:
                    break
                continue
            if end <= start:
                problems.append("line %d: end %d is not after start %d"
                                % (lineno, end, start))
                if len(problems) > 5:
                    break
                continue
            s = stats.setdefault(f[0], {"n": 0, "bp": 0, "max_end": 0,
                                        "longest": 0, "n_long": 0, "bp_long": 0,
                                        "unsorted": 0, "overlaps": 0,
                                        "prev_end": -1, "prev_start": -1})
            s["n"] += 1
            s["bp"] += end - start
            s["max_end"] = max(s["max_end"], end)
            s["longest"] = max(s["longest"], end - start)
            if window and end - start > window:
                s["n_long"] += 1
                s["bp_long"] += end - start
            if start < s["prev_start"]:
                s["unsorted"] += 1
            elif start < s["prev_end"]:
                s["overlaps"] += 1
            s["prev_start"], s["prev_end"] = start, end
            n += 1
    if not n and not problems:
        problems.append("no usable interval found")
    return stats, problems


def mask_extents(path):
    """{chrom: last base present} from a bed, or None if unreadable."""
    if not path or not os.path.exists(path):
        return None
    out = {}
    try:
        with open(path) as fh:
            for line in fh:
                f = line.split("\t", 3)
                if len(f) < 3:
                    continue
                try:
                    end = int(f[2])
                except ValueError:
                    continue
                if end > out.get(f[0], 0):
                    out[f[0]] = end
    except OSError:
        return None
    return out or None


def assess(bed, stats, problems, chroms, ctrl_flank, mask=None):
    """(status, reasons, bp) for one bed."""
    if problems:
        return "FAIL", problems[:5], 0

    wanted = {normalise(c) for c in chroms}
    present = set(stats)
    usable = present & wanted
    bp = sum(s["bp"] for c, s in stats.items() if c in usable)

    if not usable:
        return "FAIL", [
            "none of the bed's chromosomes are in CHROMS. The bed has %s; the "
            "run analyses %s. A bed named without the 'chr' prefix is the usual "
            "cause, and it would otherwise fail hours later at concat."
            % (", ".join(sorted(present)[:5]) or "nothing",
               ", ".join(sorted(wanted)[:5]))], 0

    reasons = []
    extra = sorted(present - wanted)
    if extra:
        extra_bp = sum(stats[c]["bp"] for c in extra)
        reasons.append(
            "%d contig(s) are not in CHROMS and will be dropped without "
            "warning at concat: %s (%d bp, %.1f%% of the bed). The reported l "
            "will be smaller than the bed's own length."
            % (len(extra), ", ".join(extra[:6]), extra_bp,
               100.0 * extra_bp / (bp + extra_bp) if bp + extra_bp else 0.0))

    window = 2 * ctrl_flank
    n_long = sum(stats[c]["n_long"] for c in usable)
    if n_long:
        bp_long = sum(stats[c]["bp_long"] for c in usable)
        longest = max(stats[c]["longest"] for c in usable)
        reasons.append(
            "%d interval(s) are longer than 2 x CTRL_FLANK (%d bp), holding "
            "%d bp (%.1f%% of the analysed sequence); longest %d bp. The control "
            "window is the midpoint +/- CTRL_FLANK, so for these it lies "
            "entirely inside the interval, where controls cannot be - they can "
            "never find one and are dropped. Either raise CTRL_FLANK above half "
            "the longest interval, or the window needs to extend from the "
            "interval EDGES rather than its centre, which is a change to what "
            "\"nearby\" means. todo/tocheck.md item 1."
            % (n_long, window, bp_long,
               100.0 * bp_long / bp if bp else 0.0, longest))

    unsorted_n = sum(stats[c]["unsorted"] for c in usable)
    if unsorted_n:
        reasons.append("%d interval(s) are out of coordinate order; sort the "
                       "bed (sort -k1,1 -k2,2n)" % unsorted_n)
    overlaps = sum(stats[c]["overlaps"] for c in usable)
    if overlaps:
        reasons.append(
            "%d interval(s) overlap the one before. With MERGE_BEFORE_WINDOW=0 "
            "each is windowed separately, so overlapping intervals are counted "
            "once but windowed twice" % overlaps)

    # Coordinates grossly out of range: a corrupt file, or an entirely different
    # genome. Not an assembly check - see OUT_OF_RANGE_SLACK.
    if mask:
        beyond = []
        for c in sorted(usable):
            limit = mask.get(c)
            if limit and stats[c]["max_end"] > limit + OUT_OF_RANGE_SLACK:
                beyond.append("%s ends at %d, more than %d past the last masked "
                              "base at %d"
                              % (c, stats[c]["max_end"], OUT_OF_RANGE_SLACK,
                                 limit))
        if beyond:
            reasons.append(
                "coordinates are out of range for this reference data - a "
                "corrupt bed, or one from a different genome: %s"
                % "; ".join(beyond[:4]))

    return ("WARN" if reasons else "PASS"), reasons, bp


def rows_for(beds, chroms, ctrl_flank, mask_bed=None):
    mask = mask_extents(mask_bed)
    out = []
    for bed in beds:
        item = os.path.basename(bed)
        if not os.path.exists(bed) or os.path.getsize(bed) == 0:
            out.append({"item": item, "status": "FAIL", "n_intervals": 0,
                        "bp": 0, "reason": "missing or empty: %s" % bed})
            continue
        try:
            stats, problems = read_bed(bed, window=2 * ctrl_flank)
        except OSError as exc:
            out.append({"item": item, "status": "FAIL", "n_intervals": 0,
                        "bp": 0, "reason": "cannot read %s: %s" % (bed, exc)})
            continue
        status, reasons, bp = assess(bed, stats, problems, chroms, ctrl_flank,
                                     mask)
        wanted = {normalise(c) for c in chroms}
        out.append({
            "item": item, "status": status,
            "n_intervals": sum(s["n"] for c, s in stats.items() if c in wanted),
            "bp": bp, "reason": "; ".join(reasons) or "NA"})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=conf.get("CHROMS", " ".join(
        str(i) for i in range(1, 23))))
    ap.add_argument("--ctrl-flank", type=int,
                    default=conf.get_int("CTRL_FLANK", 25000))
    ap.add_argument("--mask-bed", default=os.environ.get("MASK_BED", ""))
    ap.add_argument("--out", default="", help="where to write the check file")
    ap.add_argument("--run", default="", help="run key, for the file name")
    ap.add_argument("beds", nargs="+")
    a = ap.parse_args(argv)

    chroms = a.chroms.split()
    if not chroms:
        # Empty means the caller passed an unset $CHROMS. Analysing nothing would
        # FAIL every bed and blame the chr prefix, so fall back to the config.
        chroms = conf.get("CHROMS", " ".join(str(i) for i in range(1, 23))).split()
        sys.stderr.write("check_bed: --chroms was empty, using pipeline.conf's "
                         "CHROMS (%d chromosomes)\n" % len(chroms))
    if not chroms:
        sys.stderr.write("check_bed: no chromosomes to analyse\n")
        return 2
    rows = rows_for(a.beds, chroms, a.ctrl_flank, a.mask_bed or None)
    if a.out and a.run:
        os.makedirs(a.out, exist_ok=True)
        path = os.path.join(a.out, CHECK_TEMPLATE % a.run)
        write_check(path, rows)
        print("wrote %s" % path)
    for r in rows:
        print("  %s %s: %d intervals, %d bp" % (r["status"], r["item"],
                                                r["n_intervals"], r["bp"]))
        if r["status"] != "PASS":
            print("      %s" % r["reason"])
    # Said once, not hedged on every bed: the build is the thing this cannot
    # check, and a reader should hear it exactly as often as it is true.
    print("  note: the genome build is NOT verified by this check. Confirm by "
          "hand that the bed matches the reference data.")
    # FAIL blocks; a WARN is information the user should have before waiting.
    return 1 if any(r["status"] == "FAIL" for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
