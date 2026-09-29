#!/usr/bin/env python3
"""Upstream gate: report problems in the previous stage's check file.

Usage: python3 gate.py <check_file> <upstream_stage_name>

Warns and ALWAYS exits 0. Re-running one stage against older output is a normal
workflow here, so an upstream problem is reported, never enforced. The only
non-zero exit is a missing command-line argument, which is a bug in the caller
rather than anything about the upstream stage.

Parsing goes through checkfile.read_status(), so the same library both writes
and reads the file and the quoting can never desync: a free-text `reason`
holding a tab or a newline is quoted on write and honoured on read, where a
line-and-tab parser would mis-column it.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from checkfile import read_status  # noqa: E402

USAGE = ("usage: gate.py <check_file> <upstream_stage_name>\n"
         "       (a missing argument is a caller bug, not an upstream failure,\n"
         "        so this is the one case that does not exit 0)")

REFERENCE_ONLY = ("         for reference only: these flag runs worth a second "
                  "look, not an accept/reject test")


def main(argv):
    if len(argv) != 3:
        sys.stderr.write(USAGE + "\n")
        return 2
    path, stage = argv[1], argv[2]

    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        sys.stderr.write("WARNING: no check file for upstream stage "
                         "'%s' (%s); proceeding unverified\n" % (stage, path))
        return 0

    status = read_status(path)
    failed = [i for i, s in status.items() if s == "FAIL"]
    warned = [i for i, s in status.items() if s == "WARN"]
    if not failed and not warned:
        return 0

    if failed:
        sys.stderr.write("WARNING: upstream stage '%s' FAILED for: %s\n"
                         % (stage, " ".join(failed)))
    if warned:
        sys.stderr.write("WARNING: upstream stage '%s' WARNED for: %s\n"
                         % (stage, " ".join(warned)))
    sys.stderr.write("         see %s; proceeding anyway\n" % path)
    sys.stderr.write(REFERENCE_ONLY + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
