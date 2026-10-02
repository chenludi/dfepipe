#!/usr/bin/env python3
"""Read and write the per-stage check files.

A check file is a TSV whose first two columns are always `item` and `status`.
Remaining columns are whatever the stage measured; a row missing one writes NA
rather than shifting the columns.

`item` identifies the unit the row is about and must be unique within a file.
Where a stage works per chromosome, the unit is the pair, so the item is
"<bed>:<chrom>" and the file also carries one "<bed>" summary row. Duplicate
items would otherwise let a later row hide an earlier failure.

Status vocabulary:
    PASS  nothing to report
    WARN  worth a look, not an error; never blocks
    FAIL  almost certainly wrong
    SKIP  legitimately not applicable, e.g. a chromosome the bed has no intervals on

`status` is validated when the file is written, so a missing or mistyped status
is an error in the check script rather than a silent NA that later reads as
"nothing wrong".
"""
import csv

REQUIRED = ["item", "status"]
STATUSES = ("PASS", "WARN", "FAIL", "SKIP")
# worst-first, for aggregating duplicate items
_SEVERITY = {"FAIL": 3, "WARN": 2, "SKIP": 1, "PASS": 0}


def write_check(path, rows):
    """Write rows to `path`, validating item/status and sanitising values.

    Raises ValueError if a row lacks `item` or `status`, if `status` is not one
    of STATUSES, or if two rows share an item. Tabs and newlines inside values
    are replaced with spaces: a newline in a free-text `reason` would otherwise
    forge an extra row for anything reading the file line by line.
    """
    seen = set()
    for i, r in enumerate(rows):
        for k in REQUIRED:
            if not str(r.get(k, "")).strip():
                raise ValueError(f"row {i} has no {k!r}: {r!r}")
        if r["status"] not in STATUSES:
            raise ValueError(
                f"row {i} has status {r['status']!r}; expected one of {', '.join(STATUSES)}")
        if r["item"] in seen:
            raise ValueError(
                f"row {i} repeats item {r['item']!r}; items must be unique "
                "(use '<bed>:<chrom>' where a stage works per chromosome)")
        seen.add(r["item"])

    cols = list(REQUIRED)
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)

    def clean(v):
        return str(v).replace("\t", " ").replace("\r", " ").replace("\n", " ")

    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE,
                       escapechar="\\")
        w.writerow(cols)
        for r in rows:
            w.writerow([clean(r.get(c, "NA")) for c in cols])


def read_status(path):
    """item -> status. Duplicate items collapse to the WORST status, not the last.

    Nothing should write duplicates (write_check refuses), but a hand-edited or
    truncated file must not be able to hide a failure behind a later PASS.
    """
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            item, status = r.get("item"), r.get("status")
            if item is None:
                continue
            if item in out and _SEVERITY.get(out[item], -1) >= _SEVERITY.get(status, -1):
                continue
            out[item] = status
    return out
