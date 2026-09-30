"""Tests for the demo checkpoint.

The plan's Task 7 was written against an older residuals API (`classify(metrics,
good, warn, bias)` and a `max_abs_resid` metric). The calibrated module that
actually shipped takes `classify(metrics, cfg)` over `max_abs_rare` /
`g2_per_snp` / `mean_rest_norm` / `max_abs_rest`, and its residual thresholds
never FAIL. These tests follow the shipped module.
"""
import textwrap

import pytest

from scripts.checks.check_demo import (
    CONVERGENCE_LL_SPREAD,
    assess_demo,
    parse_best_inference,
    read_bestfits_ll,
    read_fs,
)

# pipeline.conf's RESID_*/G2_* defaults, the values run_bestfit.sh passes through.
CFG = {
    "max_rare_warn": 25.0,
    "max_rare_fail": 40.0,
    "bulk_bias_warn": 2.3,
    "rest_over_rare_min": 5.0,
    "g2_warn": 0.008,
    "g2_fail": 0.035,
}

GOOD = {
    "max_abs_rare": 1.8, "max_abs_rest": 1.1, "mean_rest": 0.05,
    "mean_rest_norm": 0.2, "g2": 4.0, "g2_per_snp": 1e-5,
    "n_snps": 4.0e6, "n_bins_used": 80,
    "n_masked_data_zero": 79, "n_masked_model_zero": 0,
}


def metrics(**over):
    m = dict(GOOD)
    m.update(over)
    return m


# ------------------------------------------------------------------ verdicts --

def test_good_fit_passes():
    status, reasons = assess_demo(metrics(), [-1000.0, -1000.4, -1000.9], CFG,
                                 na=18460.0, n_fits=256)
    assert status == "PASS" and reasons == []


def test_extreme_residual_only_warns():
    """Residual thresholds never FAIL, however extreme (Global Constraints)."""
    status, reasons = assess_demo(metrics(max_abs_rare=99.0), [-1000.0], CFG,
                                  na=18460.0, n_fits=256)
    assert status == "WARN"
    assert any("max_abs_rare" in r for r in reasons)


def test_extreme_deviance_only_warns():
    status, reasons = assess_demo(metrics(g2_per_snp=0.4), [-1000.0], CFG,
                                  na=18460.0, n_fits=256)
    assert status == "WARN"
    assert any("g2_per_snp" in r for r in reasons)


def test_model_zero_where_variants_observed_fails():
    status, reasons = assess_demo(metrics(n_masked_model_zero=3), [-1000.0], CFG,
                                  na=18460.0, n_fits=256)
    assert status == "FAIL"
    assert any("no usable probability" in r for r in reasons)


def test_no_usable_bins_fails():
    status, reasons = assess_demo(metrics(n_bins_used=0), [-1000.0], CFG,
                                  na=18460.0, n_fits=256)
    assert status == "FAIL"
    assert any("no usable bins" in r for r in reasons)


def test_unreadable_spectra_fails():
    """No metrics at all means nothing was measured: that blocks, like no bins."""
    status, reasons = assess_demo(None, [-1000.0], CFG, na=18460.0, n_fits=256)
    assert status == "FAIL"
    assert any("spectr" in r for r in reasons)


def test_no_fits_fails():
    status, reasons = assess_demo(None, [], CFG, na=None, n_fits=0)
    assert status == "FAIL"
    assert any("no fits" in r for r in reasons)


def test_missing_best_inference_fails():
    status, reasons = assess_demo(metrics(), [-1000.0], CFG, na=None, n_fits=256)
    assert status == "FAIL"
    assert any("best_inference.txt" in r for r in reasons)


@pytest.mark.parametrize("na", [0.0, -1.0])
def test_non_positive_na_fails(na):
    status, reasons = assess_demo(metrics(), [-1000.0], CFG, na=na, n_fits=256)
    assert status == "FAIL"
    assert any("Na" in r for r in reasons)


def test_poor_convergence_only_warns():
    """Convergence is desirable, not required: it warns, never fails."""
    status, reasons = assess_demo(metrics(), [-1000.0, -1200.0, -1500.0], CFG,
                                  na=18460.0, n_fits=256)
    assert status == "WARN"
    assert any("converge" in r for r in reasons)


def test_convergence_spread_just_inside_threshold_passes():
    top = [-1000.0, -1000.0 - CONVERGENCE_LL_SPREAD + 0.01]
    assert assess_demo(metrics(), top, CFG, na=1.0, n_fits=2)[0] == "PASS"


def test_single_fit_cannot_be_judged_on_convergence():
    status, reasons = assess_demo(metrics(), [-1000.0], CFG, na=1.0, n_fits=1)
    assert status == "PASS" and reasons == []


def test_bulk_bias_is_reported_in_words():
    status, reasons = assess_demo(metrics(mean_rest=1.4, mean_rest_norm=3.0),
                                  [-1.0, -1.1], CFG, na=1.0, n_fits=2)
    assert status == "WARN"
    assert any("overestimates" in r for r in reasons)


def test_structural_fail_wins_over_convergence_warning():
    status, reasons = assess_demo(metrics(n_bins_used=0), [-1000.0, -2000.0], CFG,
                                  na=1.0, n_fits=2)
    assert status == "FAIL"
    assert not any("converge" in r for r in reasons)


# ------------------------------------------------------------------- readers --

HISTORICAL_BEST_INFERENCE = textwrap.dedent("""\
    dir_outdemog:/somewhere/quiecent_ctrl.x_qui2.chrs_1t22-demo_three_epoch-1
    model:three_epoch
    likelyhood_cutoff:500.0
    mutfile:/somewhere/fsmut.quiecent_ctrl.x_qui2.chrs_1t22.tab
    mutation_ctrl_diploid:13.629336445542616
    max_ll_key:29
    ll:-519.8079138965113
    nuB:1.716150083036583
    nuF:1000.0
    TB:0.7241850572875557
    TF:0.005606701763295142
    theta:503196.6564732152
    Na:18460.05704253424
    """)


def test_parse_best_inference_reads_the_historical_format(tmp_path):
    """The two keys later stages grep for must come back, exactly."""
    f = tmp_path / "best_inference.txt"
    f.write_text(HISTORICAL_BEST_INFERENCE)
    d = parse_best_inference(f)
    assert d["max_ll_key"] == "29"
    assert d["Na"] == pytest.approx(18460.05704253424)
    assert d["theta"] == pytest.approx(503196.6564732152)
    assert d["ll"] == pytest.approx(-519.8079138965113)


def test_parse_best_inference_missing_file_is_none(tmp_path):
    assert parse_best_inference(tmp_path / "absent.txt") is None


def test_parse_best_inference_without_the_grepped_keys_is_none(tmp_path):
    f = tmp_path / "best_inference.txt"
    f.write_text("model:three_epoch\ntheta:1.0\n")
    assert parse_best_inference(f) is None


BESTFITS = textwrap.dedent("""\
    #
    # Converged results
    # Log(likelihood)\tnuB\tnuF\tTB\tTF\ttheta
    -519.8079138965113\t1.716\t1000.0\t0.724\t0.0056\t503196.65
    -519.8141750932264\t1.716\t996.28\t0.724\t0.0056\t503134.15
    #
    # Top 100 results
    # Log(likelihood)\tnuB\tnuF\tTB\tTF\ttheta
    -519.8079138965113\t1.716\t1000.0\t0.724\t0.0056\t503196.65
    -900.0\t1.0\t1.0\t1.0\t1.0\t1.0
    """)


def test_read_bestfits_ll_takes_the_top_converged_row(tmp_path):
    f = tmp_path / "three_epoch.ctrl.29.InferDM.bestfits"
    f.write_text(BESTFITS)
    assert read_bestfits_ll(f) == pytest.approx(-519.8079138965113)


def test_read_bestfits_ll_ignores_an_unconverged_file(tmp_path):
    """No "# Converged results" block means the run did not converge at all."""
    f = tmp_path / "three_epoch.ctrl.7.InferDM.bestfits"
    f.write_text("#\n# Top 100 results\n# Log(likelihood)\tnuB\n-3.0\t1.0\n")
    assert read_bestfits_ll(f) is None


def test_read_fs_zeroes_masked_bins(tmp_path):
    f = tmp_path / "x.fs"
    f.write_text('5 folded "YRI"\n10 20 30 40 50\n1 0 0 1 1\n')
    vals = read_fs(f)
    assert list(vals) == [0.0, 20.0, 30.0, 0.0, 0.0]


def test_read_fs_without_a_mask_row(tmp_path):
    f = tmp_path / "x.fs"
    f.write_text('4 unfolded "YRI"\n1 2 3 4\n')
    assert list(read_fs(f)) == [1.0, 2.0, 3.0, 4.0]


def test_read_fs_missing_file_is_none(tmp_path):
    assert read_fs(tmp_path / "absent.fs") is None
