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
    count_outcomes,
    read_fit_outcomes,
    read_large_ll,
    read_neu_ll,
    read_pack_markers,
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
    kw.setdefault("n_converged_block", 5)
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
    status, reasons = assess(n_converged_block=1, top=[-900.0])
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
                                 ll_neu=-1000.0, n_fits=144, n_converged_block=5)
    assert status == "FAIL"
    assert any("spectra" in r for r in reasons)


def test_missing_neutral_ll_is_not_a_failure():
    """run_fitneu.sh may not have run yet; that is not a broken DFE fit."""
    status, reasons = assess(ll_neu=None)
    assert status == "PASS" and reasons == []


def test_structural_failure_is_the_whole_verdict():
    status, reasons = assess(metrics(n_masked_model_zero=1), n_converged_block=1,
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


# --------------------------------------------------- per-starting-point log --

def write_pack(d, name, markers):
    """A pack file holding one InferDFE command line per marker."""
    lines = ["dadi-cli InferDFE --fs x.fs --output-prefix %s --optimizations 100"
             % (d / m) for m in markers]
    (d / name).write_text("\n".join(lines) + "\n")


def write_log(d, rows, header=True):
    text = "marker\toutcome\tseconds\tnote\n" if header else ""
    text += "".join("%s\t%s\t%s\t%s\n" % r for r in rows)
    (d / "fit_outcomes.tsv").write_text(text)


def test_read_fit_outcomes_skips_the_header(tmp_path):
    write_log(tmp_path, [("0", "converged", "12", "rc=0"),
                         ("1", "timeout", "3600", "killed after 60m")])
    assert read_fit_outcomes(tmp_path / "fit_outcomes.tsv") == {
        "0": "converged", "1": "timeout"}


def test_read_fit_outcomes_takes_the_last_row_per_marker(tmp_path):
    """RESUME appends; the reader resolves a repeated marker to its last row."""
    write_log(tmp_path, [("7", "timeout", "3600", "killed after 60m"),
                         ("7", "converged", "44", "rc=0")])
    assert read_fit_outcomes(tmp_path / "fit_outcomes.tsv")["7"] == "converged"


def test_read_fit_outcomes_ignores_junk_and_repeated_headers(tmp_path):
    (tmp_path / "fit_outcomes.tsv").write_text(
        "marker\toutcome\tseconds\tnote\n"
        "0\tconverged\t1\trc=0\n"
        "\n"
        "marker\toutcome\tseconds\tnote\n"
        "1\tnonsense\t1\t\n"
        "2\tno_likelihood\t9\trc=0, no .bestfits\n")
    assert read_fit_outcomes(tmp_path / "fit_outcomes.tsv") == {
        "0": "converged", "2": "no_likelihood"}


def test_read_fit_outcomes_absent_file_is_empty(tmp_path):
    assert read_fit_outcomes(tmp_path / "fit_outcomes.tsv") == {}


def test_read_pack_markers_reads_every_pack_in_order(tmp_path):
    write_pack(tmp_path, "dadi-cli_dfe.0.sh", ["0", "1"])
    write_pack(tmp_path, "dadi-cli_dfe.50.sh", ["50", "51"])
    assert read_pack_markers(tmp_path) == ["0", "1", "50", "51"]


def test_read_pack_markers_without_packs_is_empty(tmp_path):
    assert read_pack_markers(tmp_path) == []


def bestfits(d, marker, text="#\n# Converged results\n-1.0\t0.1\n"):
    (d / ("%s.InferDFE.bestfits" % marker)).write_text(text)


def test_count_outcomes_classifies_all_four(tmp_path):
    write_pack(tmp_path, "dadi-cli_dfe.0.sh", ["0", "1", "2", "3", "4"])
    bestfits(tmp_path, "0")
    bestfits(tmp_path, "1")
    write_log(tmp_path, [("0", "converged", "12", "rc=0"),
                         ("1", "converged", "13", "rc=0"),
                         ("2", "no_likelihood", "40", "rc=0, no .bestfits"),
                         ("3", "timeout", "3600", "killed after 60m")])
    c = count_outcomes(tmp_path)
    assert c["n_starting_points"] == 5
    assert (c["n_converged"], c["n_no_likelihood"], c["n_timeout"],
            c["n_not_run"]) == (2, 1, 1, 1)
    assert c["frac_converged"] == pytest.approx(0.4)
    assert c["frac_not_run"] == pytest.approx(0.2)


def test_count_outcomes_an_empty_bestfits_is_not_converged(tmp_path):
    """dadi-cli can leave a zero-byte file behind; run_dfe.sh tests -s too."""
    write_pack(tmp_path, "dadi-cli_dfe.0.sh", ["0"])
    (tmp_path / "0.InferDFE.bestfits").write_text("")
    write_log(tmp_path, [("0", "no_likelihood", "40", "rc=0, no .bestfits")])
    c = count_outcomes(tmp_path)
    assert (c["n_converged"], c["n_no_likelihood"]) == (0, 1)


def test_count_outcomes_bestfits_beats_the_log(tmp_path):
    """A .bestfits present now is converged, whatever an older row claimed."""
    write_pack(tmp_path, "dadi-cli_dfe.0.sh", ["0"])
    bestfits(tmp_path, "0")
    write_log(tmp_path, [("0", "timeout", "3600", "killed after 60m")])
    assert count_outcomes(tmp_path)["n_converged"] == 1


def test_count_outcomes_without_a_log_is_converged_or_not_run(tmp_path):
    """A run from before the log existed still classifies, just more coarsely."""
    write_pack(tmp_path, "dadi-cli_dfe.0.sh", ["0", "1", "2"])
    bestfits(tmp_path, "0")
    c = count_outcomes(tmp_path)
    assert (c["n_converged"], c["n_not_run"]) == (1, 2)


def test_count_outcomes_without_packs_falls_back_to_log_and_bestfits(tmp_path):
    write_log(tmp_path, [("0", "timeout", "3600", "killed after 60m")])
    bestfits(tmp_path, "1")
    c = count_outcomes(tmp_path)
    assert c["n_starting_points"] == 2
    assert (c["n_converged"], c["n_timeout"], c["n_not_run"]) == (1, 1, 0)


def test_count_outcomes_on_an_empty_directory(tmp_path):
    c = count_outcomes(tmp_path)
    assert c["n_starting_points"] == 0
    assert c["frac_converged"] == 0.0


# --------------------------------------------------------- outcome verdicts --

def outcomes(n=100, converged=65, no_likelihood=30, timeout=5, not_run=0):
    c = {"n_starting_points": n, "n_converged": converged,
         "n_no_likelihood": no_likelihood, "n_timeout": timeout,
         "n_not_run": not_run}
    for key in ("converged", "no_likelihood", "timeout", "not_run"):
        c["frac_" + key] = c["n_" + key] / n if n else 0.0
    return c


def test_typical_outcome_mix_passes():
    status, reasons = assess(outcomes=outcomes())
    assert status == "PASS" and reasons == []


def test_too_many_never_ran_warns():
    """4 of 190 is fine; 46 of 190 is a job that ran out of time."""
    status, reasons = assess(outcomes=outcomes(190, 100, 44, 0, 46))
    assert status == "WARN"
    assert any("never ran" in r for r in reasons)


def test_a_few_never_ran_does_not_warn():
    status, reasons = assess(outcomes=outcomes(190, 120, 66, 0, 4))
    assert status == "PASS" and reasons == []


def test_too_few_converged_by_fraction_warns():
    """dELS: 3 of 158 converged."""
    status, reasons = assess(outcomes=outcomes(158, 3, 155, 0, 0))
    assert status == "WARN"
    assert any("converged" in r for r in reasons)


def test_too_few_converged_by_count_warns():
    """Above the fraction but below the floor: 9 of 20."""
    status, reasons = assess(outcomes=outcomes(20, 9, 11, 0, 0))
    assert status == "WARN"
    assert any("converged" in r for r in reasons)


def test_no_converged_fits_fails():
    status, reasons = assess(outcomes=outcomes(190, 0, 190, 0, 0), n_fits=0,
                             top=[])
    assert status == "FAIL"
    assert any("converged" in r for r in reasons)


def test_no_likelihood_is_reported_never_judged():
    """A bad starting point that fails fast is part of the design."""
    status, reasons = assess(outcomes=outcomes(190, 100, 90, 0, 0))
    assert status == "PASS" and reasons == []


def test_timeout_is_reported_never_judged():
    status, reasons = assess(outcomes=outcomes(190, 100, 0, 90, 0))
    assert status == "PASS" and reasons == []


def test_outcome_thresholds_are_arguments():
    status, _ = assess(outcomes=outcomes(190, 100, 44, 0, 46),
                       max_notrun_frac=0.5)
    assert status == "PASS"


def test_outcomes_absent_is_not_judged():
    """The `neu` row has no starting points at all."""
    status, reasons = assess(outcomes=None)
    assert status == "PASS" and reasons == []
