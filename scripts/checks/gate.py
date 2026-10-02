#!/usr/bin/env python3
"""Upstream gate: report problems in the previous stage's check file.

Usage: python3 gate.py <check_file> <upstream_stage_name>

FAIL blocks, and so does a missing check file - inside a chain the upstream
stage has just run, so no check file means its check crashed rather than that
checking was skipped. WARN does not block. A FAIL means the data is broken - a chromosome
missing in error, too few SNPs, a cache with too few usable spectra, a model
giving no probability where variants are observed - and carrying that into a
ten-hour job wastes a day, so the gate exits non-zero and the chain stops.
Set SKIP_GATE=1 to proceed anyway; the override is logged. FORCE=1 also works,
for compatibility, but it means six other things in the stage scripts and is
forwarded to every later job, so it rewrites finished work from here on.

A WARN is advisory and never blocks. Some of those thresholds are measured and
some are judgements; each one says which in config/pipeline.conf.

Running a stage on its own does not consult a gate at all, so re-running one
stage against older output - a normal workflow here - is unaffected.

Parsing goes through checkfile.read_status(), so the same library both writes
and reads the file and the quoting can never desync: a free-text `reason`
holding a tab or a newline is quoted on write and honoured on read, where a
line-and-tab parser would mis-column it.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checkfile import STATUSES, read_status  # noqa: E402

USAGE = ("usage: gate.py <check_file> <upstream_stage_name>\n"
         "       (a missing argument is a caller bug, not an upstream failure,\n"
         "        so this is the one case that does not exit 0)")

REFERENCE_ONLY = ("         for reference only: these flag runs worth a second "
                  "look, not an accept/reject test")


def _override():
    """SKIP_GATE=1 bypasses this gate and nothing else.

    FORCE=1 also works, because it always did and the docs said so - but FORCE
    means six other things in the stage scripts (overwrite a finished directory,
    redo every fit, rebuild a cache, re-attempt points logged `timeout`, submit a
    bed that already has a cache, overwrite a table), and it is forwarded to every
    job in every later phase. So using it to get past one gate also rewrites work
    from that stage onward. SKIP_GATE is the narrow one, and the messages below
    recommend it.
    """
    if os.environ.get("SKIP_GATE") == "1":
        return "SKIP_GATE=1"
    if os.environ.get("FORCE") == "1":
        return "FORCE=1"
    return None


def main(argv):
    if len(argv) != 3:
        sys.stderr.write(USAGE + "\n")
        return 2
    path, stage = argv[1], argv[2]

    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        # Blocking, like a FAIL. A gate is only consulted inside a chain, where
        # the upstream stage has just run and must have written one - so an
        # absent file means its check crashed, not that checking was skipped.
        # It used to proceed with a warning, which is how a TypeError in
        # check_summary.py stayed invisible: the stage warns and exits 0 when
        # its check dies, so nothing downstream could tell a crashed check from
        # a clean one.
        sys.stderr.write("ERROR: upstream stage '%s' wrote no check file "
                         "(%s)\n" % (stage, path))
        sys.stderr.write("         Its check did not run, so nothing verified "
                         "that stage. Re-run it and read the check output, or "
                         "set SKIP_GATE=1 to continue unverified.\n")
        ov = _override()
        if ov:
            sys.stderr.write("         %s: continuing unverified\n" % ov)
            return 0
        return 1

    status = read_status(path)
    failed = [i for i, s in status.items() if s == "FAIL"]
    warned = [i for i, s in status.items() if s == "WARN"]
    # Anything outside the vocabulary is reported rather than ignored: an
    # unreadable status must not be indistinguishable from "nothing wrong".
    unknown = [i for i, s in status.items() if s not in STATUSES]
    if not failed and not warned and not unknown:
        return 0

    if failed:
        sys.stderr.write("ERROR: upstream stage '%s' FAILED for: %s\n"
                         % (stage, " ".join(failed)))
    if warned:
        sys.stderr.write("WARNING: upstream stage '%s' WARNED for: %s\n"
                         % (stage, " ".join(warned)))
    if unknown:
        sys.stderr.write("WARNING: upstream stage '%s' wrote a status this gate does "
                         "not understand for: %s\n" % (stage, " ".join(unknown)))
    sys.stderr.write("         see %s\n" % path)
    if warned or unknown:
        sys.stderr.write(REFERENCE_ONLY + "\n")

    if failed:
        ov = _override()
        if ov:
            sys.stderr.write("         %s: continuing despite the failure(s) "
                             "above\n" % ov)
            return 0
        sys.stderr.write("         stopping. Fix the problem, or set SKIP_GATE=1 "
                         "to continue anyway. FORCE=1 also works but means six "
                         "other things in the stage scripts, so it would rewrite "
                         "finished work from here on.\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
