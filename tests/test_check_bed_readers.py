"""The two bed readers must agree.

read_bed_numpy parses in bulk and is the one that runs; read_bed_lines is the
reference and owns every message. A disagreement between them is a wrong verdict
on a real bed, so this compares them field by field on generated beds, and
checks that the bulk reader hands over - rather than guessing - on everything it
does not handle.
"""
import pathlib
import random
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "checks"))
import check_bed                                           # noqa: E402
from check_bed import read_bed_lines, read_bed_numpy        # noqa: E402

FIELDS = ("n", "bp", "max_end", "longest", "n_long", "bp_long",
          "unsorted", "overlaps", "prev_start", "prev_end")
WINDOW = 50000


def agree(tmp_path, text, window=WINDOW, name="t.bed"):
    """Both readers on the same bed; returns the bulk result, or None if it
    handed over."""
    bed = tmp_path / name
    bed.write_text(text)
    ref_stats, ref_problems = read_bed_lines(str(bed), window)
    bulk = read_bed_numpy(str(bed), window)
    if bulk is None:
        return None, ref_stats, ref_problems
    stats, problems = bulk
    assert problems == [], "the bulk reader never reports problems"
    assert not ref_problems, "it should have handed over: %s" % ref_problems
    assert set(stats) == set(ref_stats)
    for c in ref_stats:
        for k in FIELDS:
            assert stats[c][k] == ref_stats[c][k], (
                "%s %s: bulk %r, reference %r"
                % (c, k, stats[c][k], ref_stats[c][k]))
    return stats, ref_stats, ref_problems


# ----------------------------------------------------------- the shapes seen --

def test_a_plain_sorted_bed(tmp_path):
    s, _, _ = agree(tmp_path, "chr1\t100\t200\nchr1\t300\t400\nchr2\t0\t50\n")
    assert s is not None and s["chr1"]["n"] == 2


def test_per_base_intervals(tmp_path):
    """The phyloP beds: about one bp per interval, which is why this path
    exists."""
    rows = "".join("chr1\t%d\t%d\n" % (i * 2, i * 2 + 1) for i in range(5000))
    s, _, _ = agree(tmp_path, rows)
    assert s is not None and s["chr1"]["n"] == 5000 and s["chr1"]["bp"] == 5000


def test_unsorted_overlapping_and_long(tmp_path):
    s, _, _ = agree(tmp_path, "chr1\t0\t60000\n"          # longer than window
                              "chr1\t100000\t200000\n"
                              "chr1\t150000\t160000\n"    # overlaps
                              "chr1\t1000\t2000\n")       # out of order
    assert s is not None
    assert s["chr1"]["unsorted"] == 1 and s["chr1"]["overlaps"] == 1
    assert s["chr1"]["n_long"] == 2


def test_chromosomes_interleaved(tmp_path):
    """prev_start/prev_end are per chromosome, so interleaving must not leak."""
    rows = "".join("chr%d\t%d\t%d\n" % (1 + i % 3, i * 10, i * 10 + 5)
                   for i in range(300))
    s, _, _ = agree(tmp_path, rows)
    assert s is not None and len(s) == 3


def test_more_than_three_columns(tmp_path):
    s, _, _ = agree(tmp_path, "chr1\t100\t200\tname\t0\t+\n"
                              "chr1\t300\t400\tname\t0\t-\n")
    assert s is not None and s["chr1"]["bp"] == 200


def test_no_trailing_newline(tmp_path):
    s, _, _ = agree(tmp_path, "chr1\t100\t200\nchr1\t300\t400")
    assert s is not None and s["chr1"]["n"] == 2


def test_contigs_outside_the_main_set(tmp_path):
    s, _, _ = agree(tmp_path, "chr1\t0\t100\nchrX\t0\t900\n"
                              "chrUn_GL000220v1\t0\t10\n")
    assert s is not None and len(s) == 3


# ------------------------------------------------- where it must hand over --

@pytest.mark.parametrize("text,why", [
    ("chr1\t100\n", "too few columns"),
    ("chr1\tabc\t200\n", "not an integer"),
    ("chr1\t200\t100\n", "reversed"),
    ("chr1\t-5\t200\n", "negative"),
    ("chr1\t\t200\n", "empty field"),
    ("", "empty file"),
    ("# comment\nchr1\t100\t200\n", "comment line"),
    ("track name=x\nchr1\t100\t200\n", "track line"),
    ("chr1\t100\t200\n\nchr1\t300\t400\n", "blank line"),
    ("chr1\t100\t200\r\n", "CRLF"),
    ("chr1\t100\t200\nchr1\t300\t400\tx\n", "column count varies"),
])
def test_the_bulk_reader_hands_over(tmp_path, text, why):
    bulk, _, _ = agree(tmp_path, text)
    assert bulk is None, "should have handed over: %s" % why


def test_handing_over_still_produces_the_reference_message(tmp_path):
    """What the user sees for a bad bed must come from the reference reader."""
    bed = tmp_path / "bad.bed"
    bed.write_text("chr1\t200\t100\n")
    stats, problems = check_bed.read_bed(str(bed), WINDOW)
    assert problems and "not after start" in problems[0]


# ------------------------------------------------------ the chunk boundary --

def test_agreement_across_chunk_boundaries(tmp_path, monkeypatch):
    """A chunk ends mid-bed, so the carried-over state is what is being tested.
    Every chunk size from tiny to larger than the bed."""
    rows = "".join("chr%d\t%d\t%d\n" % (1 + i % 2, i * 7, i * 7 + 3)
                   for i in range(400))
    bed = tmp_path / "c.bed"
    bed.write_text(rows)
    ref, _ = read_bed_lines(str(bed), WINDOW)
    for chunk in (16, 37, 64, 127, 256, 1024, 4096, 1 << 20):
        monkeypatch.setattr(check_bed, "CHUNK", chunk)
        got = read_bed_numpy(str(bed), WINDOW)
        assert got is not None, "chunk=%d" % chunk
        for c in ref:
            for k in FIELDS:
                assert got[0][c][k] == ref[c][k], "chunk=%d %s %s" % (chunk, c, k)


def test_a_line_longer_than_the_chunk_hands_over(tmp_path, monkeypatch):
    bed = tmp_path / "long.bed"
    bed.write_text("chr1\t100\t200\n")
    monkeypatch.setattr(check_bed, "CHUNK", 4)
    assert read_bed_numpy(str(bed), WINDOW) is None


# -------------------------------------------- the guards, at their own scope --

def test_a_long_bed_still_uses_the_bulk_reader(tmp_path, monkeypatch):
    """The block guard is per chunk. It was once a budget spent across them, so
    a bed grouped by chromosome lost the bulk reader purely for being long."""
    rows = "".join("chr1\t%d\t%d\n" % (i * 4, i * 4 + 2) for i in range(600))
    bed = tmp_path / "long.bed"
    bed.write_text(rows)
    monkeypatch.setattr(check_bed, "CHUNK", 64)       # ~80 chunks
    monkeypatch.setattr(check_bed, "MAX_BLOCKS", 3)   # 1 block per chunk, so fine
    got = read_bed_numpy(str(bed), WINDOW)
    assert got is not None and got[0]["chr1"]["n"] == 600


def test_an_interleaved_bed_hands_over(tmp_path, monkeypatch):
    """Many blocks in one chunk is the case the guard is for: the Python loop
    over blocks then costs more than the bulk parse saves."""
    rows = "".join("chr%d\t%d\t%d\n" % (1 + i % 8, i * 10, i * 10 + 5)
                   for i in range(200))
    bed = tmp_path / "mixed.bed"
    bed.write_text(rows)
    monkeypatch.setattr(check_bed, "MAX_BLOCKS", 3)
    assert read_bed_numpy(str(bed), WINDOW) is None


@pytest.mark.parametrize("text", [
    "chr1\t100\t200\tx\nchr1\t300\n",        # 3 tabs then 1: total divides
    "chr1\t100\nchr1\t300\t400\tx\n",        # 1 tab then 3
    "chr1\t100\t200\nchr1\t300\t400\tx\ty\n",
])
def test_a_varying_column_count_hands_over(tmp_path, text):
    """Each line's tabs are counted exactly, so a total that happens to divide
    by the line count cannot be mistaken for a uniform bed."""
    bulk, _, _ = agree(tmp_path, text)
    assert bulk is None


# ------------------------------------------------------------------ random --

@pytest.mark.parametrize("seed", range(40))
def test_random_beds_agree(tmp_path, seed):
    """Coordinates, widths and chromosome order all varying, so sortedness,
    overlaps and the per-chromosome carry are all exercised together."""
    rnd = random.Random(seed)
    rows = []
    for _ in range(rnd.randint(1, 400)):
        c = "chr%s" % rnd.choice([1, 2, 7, 22, "X", "Un_GL000220v1"])
        start = rnd.randint(0, 3_000_000)
        rows.append("%s\t%d\t%d\n" % (c, start,
                                      start + rnd.choice([1, 1, 1, 2, 10, 999,
                                                          60000, 200000])))
    if rnd.random() < 0.5:
        rows.sort()
    text = "".join(rows)
    if rnd.random() < 0.3:
        text = text.rstrip("\n")
    bulk, _, problems = agree(tmp_path, text)
    assert bulk is not None and not problems


# --------------------------------------------------- the fallback is audible --

def test_a_missing_numpy_is_said_out_loud(tmp_path, monkeypatch, capsys):
    """The bulk reader going missing is an 11x slowdown. Silent, it is how an
    optimisation comes to not run at all: the hoffman2 login node's python has
    no numpy."""
    monkeypatch.setattr(check_bed, "np", None)
    bed = tmp_path / "ok.bed"
    bed.write_text("chr1\t100\t200\n")
    assert check_bed.main(["--chroms", "1", str(bed)]) == 0
    assert "numpy is not available" in capsys.readouterr().out


def test_nothing_is_said_when_numpy_is_there(tmp_path, capsys):
    bed = tmp_path / "ok.bed"
    bed.write_text("chr1\t100\t200\n")
    assert check_bed.main(["--chroms", "1", str(bed)]) == 0
    assert "numpy is not available" not in capsys.readouterr().out
