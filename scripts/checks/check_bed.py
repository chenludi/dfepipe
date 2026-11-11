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

The genome build is not checked, and the summary says so every time rather than
letting a PASS imply it was. Coordinates are compared against the mask only to
catch something grossly out of range - a corrupt bed, or one from a different
genome. Two builds of the same genome differ by far less than that, so no
coordinate test here can tell them apart.

    python3 check_bed.py --chroms "1 2 3 ..." [--ctrl-flank 25000]
                        [--mask-bed <path>] [--out <dir> --run <key>]
                        [--reader auto|lines] <bed>...

"""
import argparse
import os
import sys

try:
    import numpy as np
except ImportError:        # the bulk reader is skipped, not the check
    np = None

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
# This is not an assembly check and is not presented as one. Two builds of the
# same genome differ by far less than this, so it cannot tell them apart; that
# the build is unverifiable here is said once, in the summary line, rather than
# hedged on every bed.
OUT_OF_RANGE_SLACK = 1_000_000


def normalise(chrom):
    """CHROMS is bare numbers in pipeline.conf; every bed is chr-prefixed."""
    c = str(chrom).strip()
    return c if c.startswith("chr") else "chr" + c


class Irregular(Exception):
    """The bulk reader met something only the line reader handles."""


# Measured, not guessed: 1 MB is both the fastest and the leanest of 1, 2, 4, 8,
# 16 and 32 MB on a per-base bed, because the arrays a chunk produces then stay
# in cache. Larger chunks cost time and memory together.
CHUNK = 1 << 20
MAX_CHROM_BYTES = 32
# Chromosome blocks in one chunk, not chromosomes in the bed: a chunk of a bed
# grouped by chromosome holds one or two, an interleaved one holds as many as it
# has lines. Past this the Python loop over blocks costs more than the bulk
# parse saves, so the reader hands over. Per chunk, so a bed does not trip it by
# being long.
MAX_BLOCKS = 4096


def blank_stats():
    return {"n": 0, "bp": 0, "max_end": 0, "longest": 0, "n_long": 0,
            "bp_long": 0, "unsorted": 0, "overlaps": 0,
            "prev_end": -1, "prev_start": -1}


def read_bed(path, window=0, reader="auto"):
    """(per-chrom stats, problems). Streams, so a 100 MB bed is fine.

    `window` is 2 x CTRL_FLANK; intervals longer than it are counted as they go,
    because how many there are and how much sequence they hold is what decides
    whether the control window needs rethinking.

    Two readers, one result. read_bed_lines is the reference: it owns every
    message a malformed bed produces. read_bed_numpy parses in bulk, measured
    13x faster on a bed of 5 million intervals, but handles only a well-formed
    bed and returns None for anything else, so the reference reader still
    decides what is wrong with it. reader="lines" forces the reference, for the
    tests and for a bed whose verdict is ever in doubt. tests/
    test_check_bed_readers.py holds them to the same answer.
    """
    if reader != "lines":
        bulk = read_bed_numpy(path, window)
        if bulk is not None:
            return bulk
    return read_bed_lines(path, window)


def read_bed_lines(path, window=0):
    """The reference reader: one Python pass, every message this check emits."""
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
            # get, not setdefault: setdefault builds the default dict on every
            # line whether or not the chromosome is new.
            s = stats.get(f[0])
            if s is None:
                s = stats[f[0]] = blank_stats()
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


def read_bed_numpy(path, window=0):
    """The same stats as read_bed_lines, or None if this reader cannot be sure.

    It hands over - rather than guessing - for numpy missing, a comment, track
    or blank line, a column count that is not the same on every line, a
    coordinate that is empty or not all digits, an interval whose end is not
    after its start, a line longer than CHUNK, or a bed with no usable interval.
    So a bed that reaches here is well formed, and anything else is reported by
    the reference reader in its own words.
    """
    if np is None:
        return None
    stats = {}
    try:
        with open(path, "rb") as fh:
            tail = b""
            while True:
                block = fh.read(CHUNK)
                if not block:
                    break
                if tail:
                    block = tail + block
                cut = block.rfind(b"\n") + 1
                if not cut:
                    raise Irregular("a line longer than CHUNK")
                tail = block[cut:]
                scan_chunk(block[:cut], window, stats)
            if tail.strip():
                scan_chunk(tail + b"\n", window, stats)
    except Irregular:
        return None
    if not sum(s["n"] for s in stats.values()):
        return None                # empty or all comments: the reference says so
    return stats, []


def parse_digits(a, lo, hi):
    """Decimal value of a[lo:hi] for every row, right to left.

    One pass per digit position. Up to the shortest field every row is still
    live, which is the common case because coordinates on a chromosome mostly
    have the same number of digits, so those passes skip the bounds mask
    entirely.
    """
    width = hi - lo
    wmin, wmax = int(width.min()), int(width.max())
    if wmin < 1:
        raise Irregular("an empty coordinate field")
    if wmax > 18:
        raise Irregular("a coordinate of more than 18 digits")
    # A coordinate of 9 digits reaches 10^9; the running value, and 10^j times a
    # digit, both stay inside int32 well past any real chromosome.
    out = np.zeros(lo.size, dtype=np.int32 if wmax <= 9 else np.int64)
    for j in range(wmax):
        at = hi - 1 - j
        live = None if j < wmin else at >= lo
        # uint8 subtraction wraps, so one comparison catches every non-digit.
        digit = a[at if live is None else np.maximum(at, lo)] - np.uint8(48)
        bad = digit > np.uint8(9)
        if live is not None:
            bad &= live
        if bad.any():
            raise Irregular("a coordinate that is not all digits")
        step = digit.astype(out.dtype)
        if live is not None:
            step *= live
        out += step * out.dtype.type(10 ** j)
    return out.astype(np.int64)


def scan_chunk(buf, window, stats):
    """Fold one whole-line chunk into `stats`."""
    a = np.frombuffer(buf, dtype=np.uint8)
    ends = np.flatnonzero(a == 10)
    if not ends.size:
        raise Irregular("no complete line")
    starts = np.empty(ends.size, dtype=np.int64)
    starts[0] = 0
    starts[1:] = ends[:-1] + 1
    # Any line opening with '#' or 't' goes to the reference reader, which is
    # wider than the comment and track lines it is aimed at but costs nothing: no
    # chromosome name starts with either. 10 and 13 are blank and bare-CR lines;
    # a CRLF bed is caught later, by the stray '\r' in its end coordinate.
    head = a[starts]
    if ((head == 35) | (head == 116) | (head == 10) | (head == 13)).any():
        raise Irregular("a comment, track, blank or bare-CR line")

    tabs = np.flatnonzero(a == 9)
    if tabs.size < 2 * ends.size or tabs.size % ends.size:
        raise Irregular("the column count is not the same on every line")
    ncol = tabs.size // ends.size
    tabs = tabs.reshape(ends.size, ncol)
    # These two make the reshape above exact, not merely plausible. tabs is
    # sorted, so if row k's first tab is past line k's start and its last is
    # before line k's end, all ncol of them are inside line k - and line k holds
    # at least ncol. Every tab is inside some line, and there are exactly
    # ncol x lines of them, so each line holds exactly ncol. A line with a
    # different count cannot satisfy both.
    if not ((tabs[:, 0] > starts).all() and (tabs[:, -1] < ends).all()):
        raise Irregular("the column count is not the same on every line")

    beg = parse_digits(a, tabs[:, 0] + 1, tabs[:, 1])
    fin = parse_digits(a, tabs[:, 1] + 1, tabs[:, 2] if ncol > 2 else ends)
    span = fin - beg
    if (span <= 0).any():
        raise Irregular("an interval whose end is not after its start")

    # Chromosome blocks: compare each name with the one above, so a bed grouped
    # by chromosome costs a few passes however long it is.
    width = tabs[:, 0] - starts
    wmax = int(width.max())
    if wmax > MAX_CHROM_BYTES:
        raise Irregular("a chromosome name over %d bytes" % MAX_CHROM_BYTES)
    same = width[1:] == width[:-1]
    for i in range(wmax):
        col = np.where(width > i, a[np.minimum(starts + i, a.size - 1)], 0)
        same &= col[1:] == col[:-1]
    edge = np.flatnonzero(~same) + 1
    blk_beg = np.concatenate(([0], edge))
    blk_end = np.concatenate((edge, [ends.size]))

    if blk_beg.size > MAX_BLOCKS:
        raise Irregular("more than %d chromosome blocks in one chunk"
                        % MAX_BLOCKS)

    for lo, hi in zip(blk_beg.tolist(), blk_end.tolist()):
        try:
            chrom = a[starts[lo]:tabs[lo, 0]].tobytes().decode("ascii")
        except UnicodeDecodeError:
            raise Irregular("a chromosome name that is not ascii")
        s = stats.get(chrom)
        if s is None:
            s = stats[chrom] = blank_stats()
        b, f, d = beg[lo:hi], fin[lo:hi], span[lo:hi]
        s["n"] += int(d.size)
        s["bp"] += int(d.sum())
        s["max_end"] = max(s["max_end"], int(f.max()))
        s["longest"] = max(s["longest"], int(d.max()))
        if window:
            long_one = d > window
            if long_one.any():
                s["n_long"] += int(long_one.sum())
                s["bp_long"] += int(d[long_one].sum())
        # The line before, per chromosome, carried across blocks and chunks.
        was_beg = np.empty(d.size, dtype=np.int64)
        was_beg[0] = s["prev_start"]
        was_beg[1:] = b[:-1]
        was_fin = np.empty(d.size, dtype=np.int64)
        was_fin[0] = s["prev_end"]
        was_fin[1:] = f[:-1]
        unsorted = b < was_beg
        s["unsorted"] += int(unsorted.sum())
        s["overlaps"] += int((~unsorted & (b < was_fin)).sum())
        s["prev_start"], s["prev_end"] = int(b[-1]), int(f[-1])


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
            "interval edges rather than its centre, which is a change to what "
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
            "%d interval(s) overlap the one before. prep merges its output "
            "(run_prep.sh:112), so each overlapping group becomes one interval "
            "and gets one control window, centred on the merged span - not one "
            "window per input interval. Merge them yourself first if you want "
            "control over where those midpoints fall" % overlaps)

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


def rows_for(beds, chroms, ctrl_flank, mask_bed=None, reader="auto"):
    mask = mask_extents(mask_bed)
    out = []
    for bed in beds:
        item = os.path.basename(bed)
        if not os.path.exists(bed) or os.path.getsize(bed) == 0:
            out.append({"item": item, "status": "FAIL", "n_intervals": 0,
                        "bp": 0, "reason": "missing or empty: %s" % bed})
            continue
        try:
            stats, problems = read_bed(bed, window=2 * ctrl_flank,
                                       reader=reader)
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
    ap.add_argument("--reader", default="auto", choices=("auto", "lines"),
                    help="auto parses in bulk and falls back; lines forces "
                         "the reference reader")
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
    rows = rows_for(a.beds, chroms, a.ctrl_flank, a.mask_bed or None,
                    reader=a.reader)
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
    print("  note: the genome build is not verified by this check. Confirm by "
          "hand that the bed matches the reference data.")
    # FAIL blocks; a WARN is information the user should have before waiting.
    return 1 if any(r["status"] == "FAIL" for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
