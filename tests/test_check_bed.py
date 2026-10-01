"""The input-bed preflight.

Both failures it exists to catch are slow and quiet without it: a bed named `1`
instead of `chr1` wastes the whole sfs array before concat notices, and a bed
from the wrong assembly produces plausible numbers that are simply wrong.
"""
import os
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "checks"))
from check_bed import main, rows_for, normalise  # noqa: E402

CHAIN = str(ROOT / "scripts" / "chains" / "run_all.sh")
CHROMS = [str(i) for i in range(1, 23)]


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def one(tmp_path, text, chroms=CHROMS, flank=25000, mask=None):
    bed = write(tmp_path, "t.bed", text)
    return rows_for([bed], chroms, flank, mask)[0]


# ------------------------------------------------------------------ structure --

def test_a_good_bed_passes(tmp_path):
    mask = write(tmp_path, "mask.bed", "chr1\t0\t1000000\nchr2\t0\t1000000\n")
    r = one(tmp_path, "chr1\t100\t200\nchr2\t300\t400\n", mask=mask)
    assert r["status"] == "PASS", r["reason"]
    assert r["n_intervals"] == 2 and r["bp"] == 200


@pytest.mark.parametrize("text,fragment", [
    ("chr1\t100\n", "needs at least 3"),
    ("chr1\tabc\t200\n", "not integers"),
    ("chr1\t200\t100\n", "not after start"),
])
def test_a_malformed_line_fails(tmp_path, text, fragment):
    r = one(tmp_path, text)
    assert r["status"] == "FAIL"
    assert fragment in r["reason"]


def test_an_empty_bed_fails(tmp_path):
    r = one(tmp_path, "")
    assert r["status"] == "FAIL"


def test_a_missing_bed_fails(tmp_path):
    r = rows_for([str(tmp_path / "absent.bed")], CHROMS, 25000)[0]
    assert r["status"] == "FAIL" and "missing or empty" in r["reason"]


def test_comments_and_track_lines_are_skipped(tmp_path):
    mask = write(tmp_path, "mask.bed", "chr1\t0\t1000000\n")
    r = one(tmp_path, "# a comment\ntrack name=x\nchr1\t100\t200\n", mask=mask)
    assert r["status"] == "PASS" and r["n_intervals"] == 1


# ------------------------------------------------------- the chromosome names --

def test_a_bed_without_the_chr_prefix_fails(tmp_path):
    """The failure this check exists for: it would otherwise cost the whole sfs
    array and surface as 'cannot concatenate' an hour later."""
    r = one(tmp_path, "1\t100\t200\n2\t300\t400\n")
    assert r["status"] == "FAIL"
    assert "'chr' prefix" in r["reason"] and "concat" in r["reason"]


def test_contigs_outside_CHROMS_warn_with_the_bp_they_cost(tmp_path):
    """chrX survives prep and is dropped at concat without a word, so the
    reported l is not the bed's l."""
    r = one(tmp_path, "chr1\t0\t100\nchrX\t0\t900\n")
    assert r["status"] == "WARN"
    assert "chrX" in r["reason"] and "dropped without" in r["reason"]
    assert r["bp"] == 100                      # only the analysed bases


def test_normalise_accepts_both_conventions():
    assert normalise("1") == "chr1" and normalise("chr1") == "chr1"


# --------------------------------------------------------- the science warnings --

def test_an_interval_longer_than_the_window_warns(tmp_path):
    """Its control window falls inside itself, where controls cannot be."""
    r = one(tmp_path, "chr1\t0\t60000\n", flank=25000)
    assert r["status"] == "WARN"
    assert "can never find one" in r["reason"]


def test_the_long_interval_warning_quantifies_the_loss(tmp_path):
    """How many and how much bp is what decides whether CTRL_FLANK needs raising
    or the window needs to extend from the edges instead of the centre."""
    r = one(tmp_path, "chr1\t0\t60000\nchr1\t100000\t101000\n"
                      "chr2\t0\t200000\n", flank=25000)
    assert "2 interval(s) are longer" in r["reason"]
    assert "260000 bp" in r["reason"]
    assert "longest 200000 bp" in r["reason"]
    assert "EDGES" in r["reason"]          # names the alternative


def test_an_interval_just_under_the_window_does_not(tmp_path):
    mask = write(tmp_path, "mask.bed", "chr1\t0\t1000000\n")
    r = one(tmp_path, "chr1\t0\t49999\n", flank=25000, mask=mask)
    assert r["status"] == "PASS", r["reason"]


def test_unsorted_and_overlapping_intervals_warn(tmp_path):
    assert "out of coordinate order" in one(
        tmp_path, "chr1\t500\t600\nchr1\t100\t200\n")["reason"]
    r = one(tmp_path, "chr1\t100\t200\nchr1\t150\t400\n")["reason"]
    assert "overlap the one before" in r
    # prep merges, so an overlapping group gets ONE window, not one each. The
    # warning used to claim the opposite.
    assert "ONE control window" in r and "run_prep.sh:112" in r


def test_a_coordinate_far_past_the_mask_warns_as_out_of_range(tmp_path):
    """Catches a corrupt bed or one from a different genome. NOT an assembly
    check: two builds of the same genome differ by far less than the slack."""
    mask = write(tmp_path, "mask.bed", "chr1\t0\t1000\n")
    r = one(tmp_path, "chr1\t100\t9000000\n", mask=mask)
    assert "out of range" in r["reason"]


def test_the_summary_says_the_build_is_not_verified(tmp_path, capsys):
    """Said once, where it is true, rather than hedged on every bed."""
    good = write(tmp_path, "good.bed", "chr1\t100\t200\n")
    main(["--chroms", "1", good])
    out = capsys.readouterr().out
    assert "genome build is NOT verified" in out
    assert "hg19" not in out              # claims nothing unsourced


# -------------------------------------------------------------------- the CLI --

def test_main_returns_1_on_FAIL_and_0_otherwise(tmp_path):
    bad = write(tmp_path, "bad.bed", "1\t100\t200\n")
    good = write(tmp_path, "good.bed", "chr1\t100\t200\n")
    assert main(["--chroms", "1 2", bad]) == 1
    assert main(["--chroms", "1 2", good]) == 0


def test_main_writes_a_check_file_when_asked(tmp_path):
    good = write(tmp_path, "good.bed", "chr1\t100\t200\n")
    out = tmp_path / "prep"
    assert main(["--chroms", "1", "--out", str(out), "--run", "r", good]) == 0
    written = out / "check_bed.r.tsv"
    assert written.exists() and "good.bed" in written.read_text()


# -------------------------------------------------- the chain refuses up front --

def run_chain(tmp_path, bed, **extra):
    env = dict(os.environ)
    env.update({"OUT": str(tmp_path / "out"), "REFS_DIR": str(tmp_path),
                "PREP_RUN": "p", "FS_RUN": "f", "DEMO_RUN": "d",
                "CACHE_RUN": "c", "DFE_RUN": "x", "SUMMARY_RUN": "s",
                "DRY_RUN": "1"})
    env.update(extra)
    return subprocess.run(["bash", CHAIN, bed], capture_output=True, text=True,
                          env=env)


def test_the_chain_submits_nothing_for_a_misnamed_bed(tmp_path):
    bed = write(tmp_path, "bad.bed", "1\t100\t200\n")
    r = run_chain(tmp_path, bed)
    assert r.returncode == 1
    assert "nothing was submitted" in r.stderr
    assert "prep" not in r.stdout          # no stage line was printed


def test_FORCE_submits_anyway(tmp_path):
    bed = write(tmp_path, "bad.bed", "1\t100\t200\n")
    r = run_chain(tmp_path, bed, FORCE="1")
    assert r.returncode == 0
    assert "submitting despite the preflight failure" in r.stderr


def test_a_warning_does_not_stop_the_chain(tmp_path):
    """Every warning it raises is a judgement for the user, not a blocker."""
    bed = write(tmp_path, "ok.bed", "chr1\t0\t100\nchrX\t0\t900\n")
    r = run_chain(tmp_path, bed)
    assert r.returncode == 0
    assert "phase 1 submitted" in r.stdout


# ------------------------------------------- the assembly is always stated --




def test_an_empty_chroms_falls_back_to_the_config(tmp_path, capsys):
    """A profile does not export CHROMS, so `--chroms "$CHROMS"` in a shell that
    sourced only a profile passed an empty string - and an empty chromosome list
    FAILED every bed while blaming the chr prefix, which is the most misleading
    message this package can emit."""
    good = write(tmp_path, "good.bed", "chr1\t100\t200\n")
    assert main(["--chroms", "", good]) == 0
    assert "--chroms was empty" in capsys.readouterr().err
