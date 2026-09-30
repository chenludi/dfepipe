"""Tests for the dfe checkpoint.

The plan's Task 9 was written against an older residuals API
(`classify(metrics, good, warn, bias)` over a `max_abs_resid` metric, with
residuals able to FAIL). The calibrated module that shipped takes
`classify(metrics, cfg)` over `max_abs_rare` / `g2_per_snp` / `mean_rest_norm` /
`max_abs_rest`, and its residual thresholds never FAIL. These tests follow the
shipped module, as tests/test_check_demo.py does.
"""
import textwrap

import pytest

from scripts.checks.check_dfe import (
    CONVERGENCE_LL_SPREAD,
    assess_dfe,
    read_large_ll,
    read_neu_ll,
)

# pipeline.conf's RESID_*/G2_* defaults, the values run_bestdfe.sh passes on.
CFG = {
    "max_rare_warn": 25.0,
    "max_rare_fail": 40.0,
    "bulk_bias_warn": 2.3,
    "rest_over_rare_min": 5.0,
    "g2_warn": 0.008,
    "g2_fail": 0.035,
}

GOOD = {
    "max_abs_rare": 1.5, "max_abs_rest": 1.1, "mean_rest": 0.1,
    "mean_rest_norm": 0.2, "g2": 4.0, "g2_per_snp": 1e-5,
    "n_snps": 4.0e6, "n_bins_used": 80,
    "n_masked_data_zero": 79, "n_masked_model_zero": 0,
}


def metrics(**over):
    m = dict(GOOD)
    m.update(over)
    return m


def assess(m=None, top=(-900.0, -900.1, -900.2, -900.3, -900.4), **kw):
    kw.setdefault("ll_dfe", -900.0)
    kw.setdefault("ll_neu", -1000.0)
    kw.setdefault("n_fits", 144)
    kw.setdefault("n_converged", 5)
    return assess_dfe(metrics() if m is None else m, list(top), CFG, **kw)


# ------------------------------------------------------------------ verdicts --

def test_good_fit_passes():
    status, reasons = assess()
    assert status == "PASS" and reasons == []


def test_no_better_than_neutral_warns():
    status, reasons = assess(ll_dfe=-1000.5, top=[-1000.5])
    assert status == "WARN"
    assert any("neutral" in r for r in reasons)


def test_equal_to_neutral_warns():
    """`<=`, not `<`: a DFE that only ties the null has bought nothing."""
    status, reasons = assess(ll_dfe=-1000.0, top=[-1000.0])
    assert status == "WARN"
    assert any("neutral" in r for r in reasons)


def test_few_converged_warns_only():
    status, reasons = assess(n_converged=1, top=[-900.0])
    assert status == "WARN"
    assert any("converge" in r for r in reasons)


def test_scattered_top_fits_warn_only():
    status, reasons = assess(top=[-900.0, -905.0, -910.0])
    assert status == "WARN"
    assert any("optimum" in r for r in reasons)


def test_extreme_residual_only_warns():
    """Residual thresholds never FAIL, however extreme (Global Constraints)."""
    status, reasons = assess(metrics(max_abs_rare=99.0))
    assert status == "WARN"
    assert any("max_abs_rare" in r for r in reasons)


def test_extreme_deviance_only_warns():
    status, reasons = assess(metrics(g2_per_snp=0.4))
    assert status == "WARN"
    assert any("g2_per_snp" in r for r in reasons)


def test_model_zero_where_variants_observed_fails():
    status, reasons = assess(metrics(n_masked_model_zero=2))
    assert status == "FAIL"
    assert any("no usable probability" in r for r in reasons)


def test_no_usable_bins_fails():
    status, reasons = assess(metrics(n_bins_used=0))
    assert status == "FAIL"
    assert any("no usable bins" in r for r in reasons)


def test_no_fits_fails():
    status, reasons = assess(n_fits=0, top=[])
    assert status == "FAIL"
    assert any("no fits" in r for r in reasons)


def test_unparseable_large_ll_fails():
    status, reasons = assess(ll_dfe=None, top=[])
    assert status == "FAIL"
    assert any("large_ll" in r for r in reasons)


def test_unreadable_spectra_fails():
    status, reasons = assess_dfe(None, [-900.0], CFG, ll_dfe=-900.0,
                                 ll_neu=-1000.0, n_fits=144, n_converged=5)
    assert status == "FAIL"
    assert any("spectra" in r for r in reasons)


def test_missing_neutral_ll_is_not_a_failure():
    """run_fitneu.sh may not have run yet; that is not a broken DFE fit."""
    status, reasons = assess(ll_neu=None)
    assert status == "PASS" and reasons == []


def test_structural_failure_is_the_whole_verdict():
    status, reasons = assess(metrics(n_masked_model_zero=1), n_converged=1,
                             top=[-900.0])
    assert status == "FAIL"
    assert not any("converge" in r for r in reasons)


# ------------------------------------------------------------------- readers --

LARGE_LL = textwrap.dedent("""\
    filename\tLog(likelihood)\tshape\tscale\ttheta\tscale_div2Na\tsdadi\tNe
    178.InferDFE.bestfits\t-610.1733440139265\t0.16899270309552797\t65.58966927719874\t54890.81236797567\t0.001776529431249104\t0.00030022051071554697\t18460.05704253424
    158.InferDFE.bestfits\t-610.1733440800572\t0.16899134410792874\t65.59132448373833\t54890.81236797567\t0.0017765742633570377\t0.0003002256726722592\t18460.05704253424
    """)


def test_read_large_ll_returns_the_top_row_first(tmp_path):
    f = tmp_path / "large_ll.txt"
    f.write_text(LARGE_LL)
    header, rows = read_large_ll(f)
    assert header[:2] == ["filename", "Log(likelihood)"]
    assert len(rows) == 2
    assert rows[0]["filename"] == "178.InferDFE.bestfits"
    assert rows[0]["Log(likelihood)"] == pytest.approx(-610.1733440139265)
    assert rows[0]["shape"] == pytest.approx(0.16899270309552797)


def test_read_large_ll_missing_file_is_none(tmp_path):
    assert read_large_ll(tmp_path / "absent.txt") == (None, None)


def test_read_large_ll_header_only_is_no_rows(tmp_path):
    f = tmp_path / "large_ll.txt"
    f.write_text(LARGE_LL.splitlines()[0] + "\n")
    header, rows = read_large_ll(f)
    assert header is not None and rows == []


def test_read_neu_ll_greps_the_historical_log(tmp_path):
    (tmp_path / "has_qui.foo.chrs_1t22.dadi_cal_ll_neu_0723.log").write_text(
        "run_fitneu.sh ...\nll: -28858.706648022526\nmutratio: 0.109\n")
    assert read_neu_ll(tmp_path) == pytest.approx(-28858.706648022526)


def test_read_neu_ll_absent_is_none(tmp_path):
    assert read_neu_ll(tmp_path) is None


def test_convergence_spread_default_is_one_ll_unit():
    assert CONVERGENCE_LL_SPREAD == 1.0
