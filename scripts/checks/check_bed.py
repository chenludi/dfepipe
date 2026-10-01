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
    WARN  a coordinate beyond the end of that chromosome in MASK_BED, which
          usually means the bed is not the assembly the references are
    PASS   none of the above

    python3 check_bed.py --chroms "1 2 3 ..." [--ctrl-flank 25000]
                        [--mask-bed <path>] [--out <dir> --run <key>] <bed>...

The assembly is the one thing this cannot settle on its own: without MASK_BED it
says so rather than implying the bed was checked.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checkfile import write_check  # noqa: E402
from conf import conf  # noqa: E402

CHECK_TEMPLATE = "check_bed.%s.tsv"
# A coordinate this far past the last masked base on its chromosome is more
# likely a different assembly than a real interval.
BEYOND_MASK_SLACK = 1_000_000


def normalise(chrom):
    """CHROMS is bare numbers in pipeline.conf; every bed is chr-prefixed."""
    c = str(chrom).strip()
    return c if c.startswith("chr") else "chr" + c


def read_bed(path):
    """(per-chrom stats, problems). Streams, so a 100 MB bed is fine."""
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
                                        "longest": 0, "n_long": 0,
                                        "unsorted": 0, "overlaps": 0,
                                        "prev_end": -1, "prev_start": -1})
            s["n"] += 1
            s["bp"] += end - start
            s["max_end"] = max(s["max_end"], end)
            s["longest"] = max(s["longest"], end - start)
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
    n_long = sum(1 for c in usable for _ in [0] if stats[c]["longest"] > window)
    long_chroms = [c for c in sorted(usable) if stats[c]["longest"] > window]
    if long_chroms:
        longest = max(stats[c]["longest"] for c in long_chroms)
        reasons.append(
            "intervals longer than 2 x CTRL_FLANK (%d bp) on %s, longest %d bp: "
            "their control window lies entirely inside them, where controls "
            "cannot be, so they can never find one. Consider a larger "
            "CTRL_FLANK." % (window, ", ".join(long_chroms[:6]), longest))

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

    if mask:
        beyond = []
        for c in sorted(usable):
            limit = mask.get(c)
            if limit and stats[c]["max_end"] > limit + BEYOND_MASK_SLACK:
                beyond.append("%s ends at %d, past %d in the mask"
                              % (c, stats[c]["max_end"], limit))
        if beyond:
            reasons.append(
                "coordinates run past the end of the chromosome in MASK_BED, so "
                "the bed is probably not the assembly the reference data is: %s"
                % "; ".join(beyond[:4]))
    else:
        reasons.append(
            "the assembly was not checked: no readable MASK_BED. Confirm by "
            "hand that the bed matches the reference data's build")

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
            stats, problems = read_bed(bed)
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

    rows = rows_for(a.beds, a.chroms.split(), a.ctrl_flank, a.mask_bed or None)
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
    # FAIL blocks; a WARN is information the user should have before waiting.
    return 1 if any(r["status"] == "FAIL" for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
