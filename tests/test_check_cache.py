"""Tests for the cache checkpoint.

The plan's Task 8 sketched `assess_cache(nbytes, usable_frac, gamma_lo, gamma_hi,
na, min_bytes, min_usable)` and a `load_cache_stats` guessed from memory. The
shipped module keeps that positional signature (the plan's four tests are below,
unchanged) and adds what inspecting a real cache showed was needed:

  * `s_high`, so the WARN compares against 2*GAMMA_S_HIGH*Na from config rather
    than a hardcoded 0.5;
  * `cache_stats(gammas, spectra)`, pure and numpy-free, so the usable-spectrum
    count is testable without dadi installed. dadi stores `spectra` as a plain
    (gamma_pts x ns+1) float array, never a list of Spectrum objects, and a
    failed extrapolation shows up as a row whose sum is <= 0, not as None.

`assess_cache` still takes an absolute `min_bytes`, but nothing supplies one
directly any more: `row_for` derives it as CACHE_MIN_SIZE_FRAC of
`expected_bytes(gamma_pts, projection)`, because a cache's size follows from the
config and not from the data. The plan's fixed 600000 was larger than a correct
marker-1 cache.
"""
import pickle

import pytest

from scripts.checks.check_cache import (
    PICKLE_OVERHEAD,
    RANGE_TOLERANCE,
    assess_cache,
    cache_stats,
    expected_bytes,
    load_cache_stats,
    read_spec,
    row_for,
)

# The reference run: PROJECTION=160, marker 1's gamma_pts=300. The published
# phastCons17way_0-5 cache for these is 559,739 bytes.
REAL_GAMMA_PTS = 300
REAL_PROJECTION = 160
REAL_CACHE_BYTES = 559739


# ----------------------------------------------------- the plan's four tests --

def test_small_file_fails():
    status, reasons = assess_cache(nbytes=1000, usable_frac=1.0, gamma_lo=0.01,
                                   gamma_hi=1e4, na=10000, min_bytes=600000,
                                   min_usable=0.9)
    assert status == "FAIL" and any("size" in r for r in reasons)


def test_many_failed_spectra_fails():
    status, reasons = assess_cache(2_000_000, 0.4, 0.01, 1e4, 10000, 600000, 0.9)
    assert status == "FAIL" and any("extrapolation" in r for r in reasons)


def test_narrow_gamma_range_warns():
    """Upper bound well short of 2*0.5*Na leaves strong selection unexplorable."""
    status, reasons = assess_cache(2_000_000, 1.0, 0.01, 50.0, 10000, 600000, 0.9)
    assert status == "WARN" and any("range" in r for r in reasons)


def test_good_cache_passes():
    assert assess_cache(2_000_000, 1.0, 0.02, 10000.0, 10000, 600000, 0.9)[0] == "PASS"


# ------------------------------------------------------------ verdict, extra --

def test_missing_file_fails():
    status, reasons = assess_cache(None, None, None, None, 10000, 600000, 0.9)
    assert status == "FAIL" and any("missing" in r for r in reasons)


def test_usable_frac_exactly_at_the_threshold_passes():
    assert assess_cache(2_000_000, 0.9, 0.02, 10000.0, 10000, 600000, 0.9)[0] == "PASS"


def test_the_range_warning_uses_the_configured_s_high():
    """GAMMA_S_HIGH=0.05 makes 2*0.05*Na=1000 the target, so 500 is fine."""
    assert assess_cache(2_000_000, 1.0, 0.02, 500.0, 10000, 600000, 0.9,
                        s_high=0.05)[0] == "PASS"
    assert assess_cache(2_000_000, 1.0, 0.02, 500.0, 10000, 600000, 0.9,
                        s_high=0.5)[0] == "WARN"


def test_range_warning_boundary_is_the_tolerance():
    """RANGE_TOLERANCE of 2*s_high*Na is the cutoff, not an epsilon either side."""
    target = 2 * 0.5 * 10000 * RANGE_TOLERANCE
    assert assess_cache(2_000_000, 1.0, 0.02, target, 10000, 600000, 0.9)[0] == "PASS"
    assert assess_cache(2_000_000, 1.0, 0.02, target * 0.99, 10000, 600000,
                        0.9)[0] == "WARN"


def test_an_unknown_na_does_not_warn_about_the_range():
    """Na comes from the demo stage; without it there is nothing to compare to."""
    assert assess_cache(2_000_000, 1.0, 0.02, 50.0, None, 600000, 0.9)[0] == "PASS"


def test_a_failing_size_reports_only_the_size():
    """The other numbers are meaningless once the file is truncated."""
    status, reasons = assess_cache(10, 0.1, 0.01, 1.0, 10000, 600000, 0.9)
    assert status == "FAIL" and len(reasons) == 1


# --------------------------------------------------------- the expected size --

def test_expected_bytes_is_the_array_plus_overhead():
    payload = REAL_GAMMA_PTS * (REAL_PROJECTION + 1) * 8
    assert payload == 386400
    assert expected_bytes(REAL_GAMMA_PTS, REAL_PROJECTION) == int(
        payload * PICKLE_OVERHEAD)


def test_expected_bytes_brackets_the_real_cache():
    """The estimate must land near the measured 559,739 bytes, not below half."""
    expected = expected_bytes(REAL_GAMMA_PTS, REAL_PROJECTION)
    assert 0.8 * REAL_CACHE_BYTES < expected < 1.2 * REAL_CACHE_BYTES


def test_expected_bytes_scales_with_both_axes():
    one = expected_bytes(300, 160)
    assert expected_bytes(600, 160) == pytest.approx(2 * one, rel=1e-6)
    assert expected_bytes(1000, 160) > expected_bytes(1000, 100)


def test_the_real_cache_clears_the_derived_floor():
    """The published cache is the validation target: it must not FAIL on size."""
    floor = int(0.5 * expected_bytes(REAL_GAMMA_PTS, REAL_PROJECTION))
    assert REAL_CACHE_BYTES > floor
    status, _ = assess_cache(REAL_CACHE_BYTES, 289 / 300.0, 0.0369201141,
                             18460.057042534212, 18460.057042534212, floor, 0.9)
    assert status == "PASS"


def test_a_cache_built_with_too_few_gamma_points_fails_on_size():
    """Asking for 1000 points and finding a 300-point cache is caught."""
    floor = int(0.5 * expected_bytes(1000, REAL_PROJECTION))
    status, reasons = assess_cache(REAL_CACHE_BYTES, 1.0, 0.02, 1e4, 10000,
                                   floor, 0.9)
    assert status == "FAIL" and any("size" in r for r in reasons)


# ---------------------------------------------------------------- statistics --

def test_cache_stats_counts_usable_rows():
    """A row summing to <= 0 is where the extrapolation broke down."""
    gammas = [-100.0, -10.0, -1.0, -0.1]
    spectra = [[1.0, 2.0], [-3.0, 1.0], [0.0, 0.0], [5.0, 5.0]]
    n, usable, lo, hi = cache_stats(gammas, spectra)
    assert (n, usable) == (4, 2)
    # the range covers the USABLE rows: |−100| and |−0.1|
    assert (lo, hi) == (0.1, 100.0)


def test_cache_stats_skips_none_and_nonfinite_rows():
    gammas = [-4.0, -3.0, -2.0, -1.0]
    spectra = [[1.0, 1.0], None, [float("nan"), 5.0], [float("inf"), 1.0]]
    n, usable, lo, hi = cache_stats(gammas, spectra)
    assert (n, usable, lo, hi) == (4, 1, 4.0, 4.0)


def test_cache_stats_with_nothing_usable():
    n, usable, lo, hi = cache_stats([-1.0, -2.0], [[0.0], [-1.0]])
    assert (n, usable) == (2, 0)
    assert lo is None and hi is None


def test_cache_stats_empty():
    assert cache_stats([], []) == (0, 0, None, None)


def test_cache_stats_accepts_numpy():
    numpy = pytest.importorskip("numpy")
    gammas = -numpy.logspace(4, -2, 6)
    spectra = numpy.ones((6, 3))
    spectra[2] = -1.0
    n, usable, lo, hi = cache_stats(gammas, spectra)
    assert (n, usable) == (6, 5)
    assert hi == pytest.approx(1e4)


# ------------------------------------------------------------------- readers --

class FakeCache(object):
    """Stands in for dadi.DFE.Cache1D_mod.Cache1D, which needs dadi to unpickle.

    Only the two attributes the real class defines and this check reads are here:
    `gammas` (negative, log-spaced) and `spectra` (one row per gamma).
    """

    def __init__(self, gammas, spectra):
        self.gammas = gammas
        self.spectra = spectra


def test_load_cache_stats_reads_a_pickle(tmp_path):
    f = tmp_path / "cache_1.spectra.bpkl"
    with open(str(f), "wb") as fh:
        pickle.dump(FakeCache([-10.0, -1.0, -0.1], [[1.0], [2.0], [-1.0]]), fh)
    assert load_cache_stats(str(f)) == (3, 2, 1.0, 10.0)


def test_load_cache_stats_on_a_corrupt_pickle_is_none(tmp_path):
    f = tmp_path / "cache_1.spectra.bpkl"
    f.write_bytes(b"not a pickle at all")
    assert load_cache_stats(str(f)) is None


def test_read_spec_returns_na(tmp_path):
    (tmp_path / "spec.txt").write_text(
        "CACHE_RUN=1\nNa=18460.057042534212\ngamma_hi=18460.0570425342\n")
    spec = read_spec(str(tmp_path / "spec.txt"))
    assert spec["Na"] == pytest.approx(18460.057042534212)
    assert spec["CACHE_RUN"] == "1"


def test_read_spec_missing_file_is_empty(tmp_path):
    assert read_spec(str(tmp_path / "absent.txt")) == {}


def test_read_spec_ignores_an_unparseable_na(tmp_path):
    (tmp_path / "spec.txt").write_text("Na=NA\n")
    assert "Na" not in read_spec(str(tmp_path / "spec.txt"))


# ---------------------------------------------------------------------- rows --

def _cache_dir(tmp_path, base, run):
    d = tmp_path / base / run
    d.mkdir(parents=True)
    return d


def _stub(d, name, gammas, spectra):
    with open(str(d / name), "wb") as fh:
        pickle.dump(FakeCache(gammas, spectra), fh)


def _row(d, marker="1", gamma_pts=1, projection=1, min_size_frac=0.5,
         min_usable=0.9, s_high=0.5):
    """row_for() with the size floor scaled to the stub caches below.

    The stubs hold a handful of one-element rows, so gamma_pts/projection are
    given as 1: the point of these tests is the rest of the row, and the derived
    floor gets its own tests above and below.
    """
    return row_for(str(d), "ELS", marker, min_size_frac, min_usable, s_high,
                   gamma_pts, projection)


def test_row_for_a_missing_cache_fails(tmp_path):
    d = _cache_dir(tmp_path, "ELS", "1-1-1")
    row = _row(d, gamma_pts=REAL_GAMMA_PTS, projection=REAL_PROJECTION)
    assert row["status"] == "FAIL"
    assert "missing" in row["reason"]
    assert row["item"] == "ELS"
    assert row["file_bytes"] == "NA"
    # The expected size is knowable without the file, and is reported anyway.
    assert row["expected_bytes"] == expected_bytes(REAL_GAMMA_PTS, REAL_PROJECTION)
    assert row["min_bytes"] == int(0.5 * row["expected_bytes"])


def test_row_for_a_good_cache_reports_the_numbers(tmp_path):
    d = _cache_dir(tmp_path, "ELS", "1-1-1")
    _stub(d, "cache_1.spectra.bpkl", [-18460.0, -1000.0, -0.0369],
          [[1.0], [1.0], [1.0]])
    (d / "spec.txt").write_text("Na=18460.057042534212\n")
    row = _row(d)
    assert row["status"] == "PASS"
    assert row["n_gamma"] == 3 and row["n_usable"] == 3
    assert row["gamma_lo"] == pytest.approx(0.0369)
    assert row["gamma_hi"] == pytest.approx(18460.0)
    assert row["Na"] == pytest.approx(18460.057042534212)
    assert row["file_bytes"] > 0


def test_row_for_derives_the_size_floor_from_gamma_pts(tmp_path):
    """The same small cache passes as a 1-point cache and fails as a 300-point one."""
    d = _cache_dir(tmp_path, "ELS", "1-1-1")
    _stub(d, "cache_1.spectra.bpkl", [-1.0], [[1.0]])
    assert _row(d)["status"] == "PASS"
    row = _row(d, gamma_pts=REAL_GAMMA_PTS, projection=REAL_PROJECTION)
    assert row["status"] == "FAIL" and "size" in row["reason"]


def test_row_for_honours_min_size_frac(tmp_path):
    d = _cache_dir(tmp_path, "ELS", "1-1-1")
    _stub(d, "cache_1.spectra.bpkl", [-1.0], [[1.0]])
    nbytes = _row(d)["file_bytes"]
    # A gamma_pts whose expected size is just above the file: the file clears
    # half of it but not all of it.
    gamma_pts = int(2 * nbytes / (2 * 8 * PICKLE_OVERHEAD))
    assert _row(d, gamma_pts=gamma_pts, min_size_frac=0.5)["status"] == "PASS"
    assert _row(d, gamma_pts=gamma_pts, min_size_frac=1.0)["status"] == "FAIL"


def test_row_for_an_unreadable_cache_fails(tmp_path):
    d = _cache_dir(tmp_path, "ELS", "1-1-1")
    (d / "cache_1.spectra.bpkl").write_bytes(b"x" * 100)
    row = _row(d)
    assert row["status"] == "FAIL" and "unpickle" in row["reason"]


def test_row_for_uses_the_marker_in_the_file_name(tmp_path):
    d = _cache_dir(tmp_path, "ELS", "1-large-1")
    _stub(d, "cache_large.spectra.bpkl", [-1.0], [[1.0]])
    row = _row(d, marker="large")
    assert row["status"] == "PASS" and row["n_gamma"] == 1


def test_main_writes_the_check_file_beside_the_cache(tmp_path):
    from scripts.checks.check_cache import main
    d = _cache_dir(tmp_path, "ELS", "1-1-1")
    _stub(d, "cache_1.spectra.bpkl", [-18460.0, -0.04], [[1.0], [1.0]])
    (d / "spec.txt").write_text("Na=18460.0\n")
    rc = main(["--out", str(tmp_path), "--run", "1-1-1", "--marker", "1",
               "--min-size-frac", "0.5", "--gamma-pts", "1", "--projection", "1",
               "--min-usable", "0.9", "--gamma-s-high", "0.5", "ELS.bed"])
    assert rc == 0
    text = (d / "check_cache.tsv").read_text()
    assert text.splitlines()[0].startswith("item\tstatus")
    assert "\tPASS\t" in text
    assert "ELS" in text


def test_main_requires_gamma_pts_and_projection(tmp_path):
    """The size floor is derived from them, so there is no sane default."""
    from scripts.checks.check_cache import main
    _cache_dir(tmp_path, "ELS", "1-1-1")
    with pytest.raises(SystemExit):
        main(["--out", str(tmp_path), "--run", "1-1-1", "ELS.bed"])


def test_main_returns_zero_even_when_a_bed_fails(tmp_path):
    """The checkpoint is a diagnostic; the gate, not the exit status, blocks."""
    from scripts.checks.check_cache import main
    _cache_dir(tmp_path, "ELS", "1-1-1")
    assert main(["--out", str(tmp_path), "--run", "1-1-1", "--marker", "1",
                 "--gamma-pts", str(REAL_GAMMA_PTS),
                 "--projection", str(REAL_PROJECTION), "ELS.bed"]) == 0
    assert "FAIL" in (tmp_path / "ELS" / "1-1-1" / "check_cache.tsv").read_text()
