"""run_bestpdf.sh end to end: the model-selection half of the summary stage.

Every branch of the comparison is driven by writing a sumdfe.tsv with chosen
likelihoods, so the test says what the decision should be rather than echoing
what the code does.
"""
import csv
import os
import pathlib
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
SCRIPT = str(ROOT / "scripts" / "summary" / "run_bestpdf.sh")
CUTS = str(ROOT / "config" / "cutvalue_modelsel.csv")

BINS = ["0.0 <= |s| < 1e-05", "1e-05 <= |s| < 0.0001", "0.0001 <= |s| < 0.001",
        "0.001 <= |s| < 0.01", "0.01 <= |s| < 0.5", "0.5 <= |s| < inf"]
GAMMA_BINS = [0.40, 0.25, 0.25, 0.09, 0.01, 0.0]
NEUGAMMA_BINS = [0.70, 0.15, 0.10, 0.04, 0.01, 0.0]
UL = 2.0


def write_sumdfe(path, ll_neu, ll_gamma, ll_neugamma, l=5_000_000.0, ul=UL):
    cols = ["outputbase", "ll_neu", "ll_gamma"]
    cols += ["gamma_bin_%s" % b for b in BINS]
    cols += ["ll_neugamma"]
    cols += ["neugamma_bin_%s" % b for b in BINS]
    cols += ["u", "l", "l_inputbed", "ul"]
    vals = ["theBed", ll_neu, ll_gamma] + GAMMA_BINS + [ll_neugamma] \
        + NEUGAMMA_BINS + [1e-8, l, l, ul]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerow(["NA" if v is None else v for v in vals])


def run(tmp_path, primary="delta_50", small_primary="chi_5%", small_l="1500000",
        force=False):
    env = dict(os.environ)
    env.update({
        "OUT": str(tmp_path / "out"), "REFS_DIR": str(tmp_path),
        "PREP_RUN": "p", "FS_RUN": "f", "DEMO_RUN": "d", "CACHE_RUN": "c",
        "DFE_RUN": "x", "SUMMARY_RUN": "s", "MODELSEL_TABLE": CUTS,
        "MODELSEL_PRIMARY": primary, "MODELSEL_SMALL_PRIMARY": small_primary,
        "MODELSEL_SMALL_L": small_l,
    })
    if force:
        env["FORCE"] = "1"
    return subprocess.run(["bash", SCRIPT], capture_output=True, text=True,
                          env=env)


def result(tmp_path):
    path = tmp_path / "out" / "summary" / "s" / "bestpdf_sumdfe.tsv"
    with open(path) as fh:
        return next(csv.DictReader(fh, delimiter="\t"))


def setup(tmp_path, **kw):
    write_sumdfe(tmp_path / "out" / "summary" / "s" / "sumdfe.tsv", **kw)


# -------------------------------------------------------------- the branches --

def test_gamma_wins_when_it_clears_the_margin_and_neugamma_adds_nothing(tmp_path):
    # gamma - neu = 100 > 50; neugamma - gamma = 0 < 50
    setup(tmp_path, ll_neu=-1100.0, ll_gamma=-1000.0, ll_neugamma=-1000.0)
    assert run(tmp_path).returncode == 0
    r = result(tmp_path)
    assert r["pdf_best_delta_50"] == "gamma"
    assert r["pdf_plotting_delta_50"] == "gamma"


def test_neutral_wins_when_gamma_does_not_clear_the_margin(tmp_path):
    # gamma - neu = 10, under delta_50's 50; neugamma adds nothing
    setup(tmp_path, ll_neu=-1010.0, ll_gamma=-1000.0, ll_neugamma=-1000.0)
    assert run(tmp_path).returncode == 0
    r = result(tmp_path)
    assert r["pdf_best_delta_50"] == "neutral"
    # ... while the looser chi_5% cutoff of 3 does accept gamma
    assert r["pdf_best_chi_5%"] == "gamma"


def test_neugamma_wins_from_gamma(tmp_path):
    # gamma - neu = 1000 > 50, and neugamma - gamma = 100 > 50
    setup(tmp_path, ll_neu=-2000.0, ll_gamma=-1000.0, ll_neugamma=-900.0)
    assert run(tmp_path).returncode == 0
    assert result(tmp_path)["pdf_best_delta_50"] == "neugamma"


def test_neugamma_wins_from_neutral(tmp_path):
    """gamma fails against neutral, but neugamma beats neutral directly - the
    branch summary_getbestpdf.py:52-56 exists for."""
    # gamma - neu = 10 (< 50) so neutral; neugamma - neu = 100 (> 50)
    setup(tmp_path, ll_neu=-1010.0, ll_gamma=-1000.0, ll_neugamma=-910.0)
    assert run(tmp_path).returncode == 0
    assert result(tmp_path)["pdf_best_delta_50"] == "neugamma"


def test_the_three_deltas_are_recorded(tmp_path):
    setup(tmp_path, ll_neu=-2000.0, ll_gamma=-1000.0, ll_neugamma=-900.0)
    run(tmp_path)
    r = result(tmp_path)
    assert float(r["ll_gamma-ll_neu"]) == pytest.approx(1000.0)
    assert float(r["ll_neugamma-ll_gamma"]) == pytest.approx(100.0)
    assert float(r["ll_neugamma-ll_neu"]) == pytest.approx(1100.0)


# ------------------------------------------------------- agreement across levels --

def test_allconfidence_is_filled_when_every_level_agrees(tmp_path):
    setup(tmp_path, ll_neu=-3000.0, ll_gamma=-1000.0, ll_neugamma=-1000.0)
    run(tmp_path)
    assert result(tmp_path)["pdf_best_allconfidence"] == "gamma"


def test_allconfidence_is_blank_when_levels_disagree(tmp_path):
    """gamma - neu = 10 clears chi_5% (3) and delta_10 (10 is not > 10, so no)
    but not delta_50, so the levels cannot agree."""
    setup(tmp_path, ll_neu=-1010.0, ll_gamma=-1000.0, ll_neugamma=-1000.0)
    run(tmp_path)
    r = result(tmp_path)
    assert r["pdf_best_allconfidence"] == ""
    assert r["pdf_best_chi_5%"] == "gamma" and r["pdf_best_delta_50"] == "neutral"


# ------------------------------------------------------------ the reported DFE --

def test_the_winning_pdfs_bins_are_copied(tmp_path):
    setup(tmp_path, ll_neu=-2000.0, ll_gamma=-1000.0, ll_neugamma=-900.0)
    run(tmp_path)
    r = result(tmp_path)
    assert r["pdf_plotting_delta_50"] == "neugamma"
    got = [float(r["delta_50_%s" % b]) for b in BINS]
    assert got == pytest.approx(NEUGAMMA_BINS)


def test_a_neutral_verdict_gives_all_mass_to_the_first_bin(tmp_path):
    """summary_getbestpdf.py:124-126 - not the gamma bins, which is the bug its
    2025/04/02 header records fixing."""
    setup(tmp_path, ll_neu=-1010.0, ll_gamma=-1000.0, ll_neugamma=-1000.0)
    run(tmp_path)
    r = result(tmp_path)
    assert r["pdf_plotting_delta_50"] == "neutral"
    got = [float(r["delta_50_%s" % b]) for b in BINS]
    assert got == pytest.approx([1.0, 0, 0, 0, 0, 0])


def test_num_nonneutral_is_one_minus_the_neutral_bin_times_ul(tmp_path):
    setup(tmp_path, ll_neu=-1100.0, ll_gamma=-1000.0, ll_neugamma=-1000.0)
    run(tmp_path)
    r = result(tmp_path)
    assert float(r["num_nonneutral"]) == pytest.approx((1 - GAMMA_BINS[0]) * UL)


def test_num_nonneutral_is_zero_for_a_neutral_verdict(tmp_path):
    setup(tmp_path, ll_neu=-1010.0, ll_gamma=-1000.0, ll_neugamma=-1000.0)
    run(tmp_path)
    assert float(result(tmp_path)["num_nonneutral"]) == pytest.approx(0.0)


# --------------------------------------------------------- the small-target rule --

def test_a_small_target_uses_the_small_cutoff(tmp_path):
    """Below MODELSEL_SMALL_L a flat 50-unit margin is too strict, so the
    nominal chi-square cutoff decides instead."""
    setup(tmp_path, ll_neu=-1010.0, ll_gamma=-1000.0, ll_neugamma=-1000.0,
          l=100_000.0)
    run(tmp_path)
    r = result(tmp_path)
    assert r["modelsel_level"] == "chi_5%"
    assert r["pdf_plotting_delta_50"] == "gamma"          # delta_50 would say neutral
    assert r["pdf_best_delta_50"] == "neutral"   # still recorded


def test_a_large_target_uses_the_primary_cutoff(tmp_path):
    setup(tmp_path, ll_neu=-1010.0, ll_gamma=-1000.0, ll_neugamma=-1000.0,
          l=5_000_000.0)
    run(tmp_path)
    r = result(tmp_path)
    assert r["modelsel_level"] == "delta_50"
    assert r["pdf_plotting_delta_50"] == "neutral"


def test_small_l_zero_restores_the_original_single_cutoff(tmp_path):
    setup(tmp_path, ll_neu=-1010.0, ll_gamma=-1000.0, ll_neugamma=-1000.0,
          l=1.0)
    run(tmp_path, small_l="0")
    assert result(tmp_path)["modelsel_level"] == "delta_50"


# ------------------------------------------------------------------ guard rails --

def test_a_missing_likelihood_gives_NA_not_a_verdict(tmp_path):
    """The original would raise on the subtraction; a wrong answer would be
    worse than no answer."""
    setup(tmp_path, ll_neu=None, ll_gamma=-1000.0, ll_neugamma=-1000.0)
    assert run(tmp_path).returncode == 0
    r = result(tmp_path)
    assert r["ll_gamma-ll_neu"] == "NA"
    assert r["pdf_best_delta_50"] == "NA" and r["pdf_plotting_delta_50"] == "NA"
    assert r["pdf_best_allconfidence"] == ""
    assert r["num_nonneutral"] == "NA"


def test_an_unknown_primary_level_is_refused(tmp_path):
    setup(tmp_path, ll_neu=-1100.0, ll_gamma=-1000.0, ll_neugamma=-1000.0)
    r = run(tmp_path, primary="delta_999")
    assert r.returncode != 0
    assert "not a confidence level" in r.stderr


def test_no_summary_table_is_refused(tmp_path):
    (tmp_path / "out" / "summary" / "s").mkdir(parents=True)
    r = run(tmp_path)
    assert r.returncode == 1
    assert "no summary table" in r.stderr


def test_an_existing_output_is_not_overwritten(tmp_path):
    setup(tmp_path, ll_neu=-1100.0, ll_gamma=-1000.0, ll_neugamma=-1000.0)
    assert run(tmp_path).returncode == 0
    again = run(tmp_path)
    assert again.returncode == 1 and "already exists" in again.stderr
    assert run(tmp_path, force=True).returncode == 0


def test_every_input_column_survives(tmp_path):
    """The table is appended to, never rebuilt, so nothing upstream is lost."""
    setup(tmp_path, ll_neu=-1100.0, ll_gamma=-1000.0, ll_neugamma=-1000.0)
    run(tmp_path)
    src = tmp_path / "out" / "summary" / "s" / "sumdfe.tsv"
    with open(src) as fh:
        original = next(csv.reader(fh, delimiter="\t"))
    r = result(tmp_path)
    for col in original:
        assert col in r, col


def test_the_small_target_rule_is_off_by_default(tmp_path):
    """It changes which pdf is reported for a short bed, and the original applied
    one cutoff to every bed. A result-affecting invention must be opt-in."""
    conf = (ROOT / "config" / "pipeline.conf").read_text()
    assert 'MODELSEL_SMALL_L="${MODELSEL_SMALL_L:-0}"' in conf


def test_with_the_rule_off_a_short_bed_uses_the_primary_cutoff(tmp_path):
    setup(tmp_path, ll_neu=-1010.0, ll_gamma=-1000.0, ll_neugamma=-1000.0,
          l=1000.0)
    run(tmp_path, small_l="0")
    r = result(tmp_path)
    assert r["modelsel_level"] == "delta_50"
    assert r["pdf_plotting_delta_50"] == "neutral"
