"""Tests for the prep checkpoint.

prep's whole job is a four-filter chain, so its checkpoint is bp accounting: one
row per step, in the order the filters actually run, so a wrong result localises
to the step that ate the sequence.
"""
import csv
import pathlib
import subprocess
import sys

import pytest

from scripts.checks.check_prep import ORDER, assess, main
from scripts.checks.checkfile import read_status

ROOT = pathlib.Path(__file__).resolve().parent.parent

# A healthy run: nearly everything survives.
GOOD = {"input": 1000, "minus_HET": 980, "minus_quiescent": 950,
        "minus_CDS": 940, "after_mask": 900}


def _run_dir(out, base, run, rows, bed=True):
    """Fake one prep run directory: lengths.tsv plus the output bed."""
    d = out / base / run
    d.mkdir(parents=True, exist_ok=True)
    if rows is not None:
        with open(d / "lengths.tsv", "w") as fh:
            fh.write("stage\tbp\n")
            for k, v in rows.items():
                fh.write("%s\t%s\n" % (k, v))
    if bed:
        (d / ("s.nc_exHET_noq.%s.bed" % base)).write_text("chr1\t0\t1\n")
    return d


# --------------------------------------------------------------- filter order --

def test_order_is_the_four_filters_in_the_order_they_run():
    """submit_exHETquiencmask.sh: HET then quiescent; filt_nc-unmask.sh: CDS then mask."""
    assert ORDER == ["input", "minus_HET", "minus_quiescent", "minus_CDS", "after_mask"]


# ------------------------------------------------------------------- assess() --

def test_flags_total_loss():
    rows = dict(GOOD, after_mask=0)
    status, reasons = assess(rows, min_frac=0.05)
    assert status == "FAIL"
    assert any("no sites" in r for r in reasons)


def test_flags_heavy_loss():
    rows = dict(GOOD, after_mask=20)
    status, reasons = assess(rows, min_frac=0.05)
    assert status == "WARN"
    assert reasons


def test_normal_retention_passes():
    status, reasons = assess(GOOD, min_frac=0.05)
    assert status == "PASS"
    assert reasons == []


def test_exactly_min_frac_passes():
    """The WARN is for FEWER than min_frac, so the boundary itself is fine."""
    assert assess(dict(GOOD, after_mask=50), min_frac=0.05)[0] == "PASS"


def test_min_frac_is_a_knob_not_a_constant():
    rows = dict(GOOD, after_mask=40)          # 4% survives
    assert assess(rows, min_frac=0.05)[0] == "WARN"
    assert assess(rows, min_frac=0.01)[0] == "PASS"


def test_missing_lengths_rows_are_not_a_crash():
    """A truncated lengths.tsv must produce a verdict, not a traceback."""
    status, reasons = assess({"input": 1000}, min_frac=0.05)
    assert status == "FAIL"


def test_no_input_bp_is_a_fail_not_a_zero_division():
    status, _ = assess({"input": 0, "after_mask": 0}, min_frac=0.05)
    assert status == "FAIL"


# --------------------------------------------------------------------- main() --

def test_writes_check_file_with_every_step_and_frac(tmp_path):
    out = tmp_path / "prep"
    _run_dir(out, "EnhA", "prep_t", GOOD)
    assert main(["--out", str(out), "--run", "prep_t", "EnhA.bed"]) == 0

    dest = out / "check_prep.prep_t.tsv"
    assert dest.is_file()
    rows = list(csv.DictReader(dest.open(), delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["item"] == "EnhA"
    assert row["status"] == "PASS"
    for step in ORDER:                        # every step's bp is reported
        assert row[step] == str(GOOD[step])
    assert float(row["frac_retained"]) == pytest.approx(0.9)
    assert float(row["frac_minus_HET"]) == pytest.approx(0.98)
    assert float(row["frac_minus_quiescent"]) == pytest.approx(0.95)
    assert float(row["frac_minus_CDS"]) == pytest.approx(0.94)
    assert read_status(dest) == {"EnhA": "PASS"}


def test_check_file_covers_every_bed_and_keeps_the_worst_verdict(tmp_path):
    out = tmp_path / "prep"
    _run_dir(out, "EnhA", "prep_t", GOOD)
    _run_dir(out, "ELS", "prep_t", dict(GOOD, after_mask=0))
    assert main(["--out", str(out), "--run", "prep_t", "EnhA.bed", "ELS.bed"]) == 0
    assert read_status(out / "check_prep.prep_t.tsv") == {"EnhA": "PASS", "ELS": "FAIL"}


def test_missing_lengths_file_is_a_fail_row_not_a_crash(tmp_path):
    """The check file must still be written, or the gate has nothing to read."""
    out = tmp_path / "prep"
    _run_dir(out, "EnhA", "prep_t", GOOD)
    _run_dir(out, "Gone", "prep_t", None, bed=False)
    assert main(["--out", str(out), "--run", "prep_t", "EnhA.bed", "Gone.bed"]) == 0
    dest = out / "check_prep.prep_t.tsv"
    assert read_status(dest) == {"EnhA": "PASS", "Gone": "FAIL"}
    text = dest.read_text()
    assert "lengths.tsv" in text


def test_missing_output_bed_is_a_fail(tmp_path):
    out = tmp_path / "prep"
    _run_dir(out, "EnhA", "prep_t", GOOD, bed=False)
    assert main(["--out", str(out), "--run", "prep_t", "EnhA.bed"]) == 0
    assert read_status(out / "check_prep.prep_t.tsv") == {"EnhA": "FAIL"}


def test_bed_argument_may_be_a_path_and_may_lack_the_suffix(tmp_path):
    out = tmp_path / "prep"
    _run_dir(out, "EnhA", "prep_t", GOOD)
    assert main(["--out", str(out), "--run", "prep_t",
                 str(tmp_path / "beds" / "EnhA")]) == 0
    assert read_status(out / "check_prep.prep_t.tsv") == {"EnhA": "PASS"}


def test_runs_as_a_script(tmp_path):
    """run_prep.sh invokes it with `python3 <path>`, so the CLI must stand alone."""
    out = tmp_path / "prep"
    _run_dir(out, "EnhA", "prep_t", GOOD)
    r = subprocess.run([sys.executable, str(ROOT / "scripts/checks/check_prep.py"),
                        "--out", str(out), "--run", "prep_t", "EnhA.bed"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (out / "check_prep.prep_t.tsv").is_file()
