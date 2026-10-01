"""Every check script must survive an empty run and write a usable check file.

This exists because of a bug the rest of the suite could not catch.
check_summary.py called `write_check(path, rows, fieldnames=[...])`, which that
function does not accept. The stage scripts run their check with

    if ! python3 .../check_x.py ...; then echo "WARNING: ..." ; fi

so a check that dies only prints a warning and the stage still exits 0 - the
failure is invisible unless someone reads the output. And the unit tests called
`rows_for()` and `assess_*()` directly, never `main()`, so the broken call was
never executed. It surfaced only on the cluster.

An empty run directory is the cheapest input that still exercises the whole
plumbing: argparse, the row builders, the write_check call and the output path.
Every check should report the absence rather than crash, so this doubles as the
"what happens when a stage has not run yet" test.
"""
import csv
import importlib
import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "checks"))

# module name -> the argument list, with {out} filled in per test
CHECKS = {
    "check_prep": ["--out", "{out}", "--run", "r", "A.bed"],
    "check_sfs": ["--out", "{out}", "--run", "r", "--prep-out", "{out}",
                  "--prep-run", "r", "--chroms", "1 2", "A.bed"],
    "check_demo": ["--out", "{out}", "--run", "r", "--sfs-out", "{out}",
                   "--fs-run", "r", "A.bed"],
    "check_cache": ["--out", "{out}", "--run", "r", "--gamma-pts", "300",
                    "--projection", "160", "A.bed"],
    "check_dfe": ["--out", "{out}", "--run", "r", "--pdf", "gamma",
                  "--sfs-out", "{out}", "--fs-run", "r", "A.bed"],
    "check_summary": ["--out", "{out}", "--run", "r", "--pdfs",
                      "neu gamma neugamma", "A.bed"],
    # check_bed is the one that runs BEFORE a stage, so its "empty run" is a
    # bed that does not exist. It must still write a check file and say FAIL.
    "check_bed": ["--chroms", "1 2", "--out", "{out}", "--run", "r",
                  "{out}/absent.bed"],
}
STATUSES = {"PASS", "WARN", "FAIL", "SKIP"}


@pytest.mark.parametrize("name", sorted(CHECKS))
def test_check_main_survives_an_empty_run(tmp_path, name):
    """main() must return an int and not raise. A TypeError in the write_check
    call, a missing argparse option or a bad import all fail here."""
    mod = importlib.import_module(name)
    argv = [a.format(out=str(tmp_path)) for a in CHECKS[name]]
    rc = mod.main(argv)
    assert isinstance(rc, int), "%s.main() returned %r" % (name, rc)


@pytest.mark.parametrize("name", sorted(CHECKS))
def test_check_writes_a_readable_check_file(tmp_path, name):
    """The file is the stage's record, so it must exist and parse, with every
    status drawn from the shared vocabulary."""
    mod = importlib.import_module(name)
    mod.main([a.format(out=str(tmp_path)) for a in CHECKS[name]])

    written = sorted(tmp_path.rglob("check_*.tsv"))
    assert written, "%s wrote no check file under %s" % (name, tmp_path)

    for path in written:
        with open(path) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert rows, "%s wrote an empty check file" % name
        for row in rows:
            assert row.get("item"), (name, row)
            assert row.get("status") in STATUSES, (name, row)


@pytest.mark.parametrize("name", sorted(CHECKS))
def test_an_empty_run_is_not_reported_as_PASS(tmp_path, name):
    """Nothing ran, so nothing can be fine. A check that says PASS on an empty
    directory would hide a stage that never started."""
    mod = importlib.import_module(name)
    mod.main([a.format(out=str(tmp_path)) for a in CHECKS[name]])
    statuses = set()
    for path in sorted(tmp_path.rglob("check_*.tsv")):
        with open(path) as fh:
            statuses |= {r["status"] for r in csv.DictReader(fh, delimiter="\t")}
    assert statuses and statuses != {"PASS"}, (name, statuses)
