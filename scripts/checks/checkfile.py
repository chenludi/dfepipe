#!/usr/bin/env python3
"""Read and write the per-stage check files.

A check file is a TSV whose first two columns are always `item` (usually the bed
name) and `status` (PASS / WARN / FAIL). Remaining columns are whatever the stage
measured; a row missing one writes NA rather than shifting the columns.
"""
import csv

REQUIRED = ["item", "status"]


def write_check(path, rows):
    cols = list(REQUIRED)
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, "NA") for c in cols])


def read_status(path):
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        return {r["item"]: r["status"] for r in rd}
