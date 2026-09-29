"""Tests for the trimmed, calibrated residual check.

Only the metrics that were shown to discriminate good from bad fits are
computed, so these tests pin both the metric definitions and every branch of
the classifier.
"""
import math

import numpy as np
import pytest

from scripts.checks import residuals as R
from scripts.checks.residuals import residual_metrics, classify


EXPECTED_KEYS = {
    "max_abs_rare",
    "max_abs_rest",
    "mean_rest",
    "mean_rest_norm",
    "g2",
    "g2_per_snp",
    "n_snps",
    "n_bins_used",
    "n_masked_data_zero",
    "n_masked_model_zero",
}

# Thresholds as calibrated on 137 previous fits (plan, Global Constraints).
CFG = {
    "max_rare_warn": 25.0,
    "max_rare_fail": 40.0,
    "bulk_bias_warn": 2.3,
    "rest_over_rare_min": 5.0,
    "g2_warn": 0.008,
    "g2_fail": 0.035,
}


def base_metrics(**over):
    """A metrics dict that classifies PASS, with overrides applied."""
    m = {
        "max_abs_rare": 1.0,
        "max_abs_rest": 0.5,
        "mean_rest": 0.1,
        "mean_rest_norm": 0.2,
        "g2": 1.0,
        "g2_per_snp": 0.0001,
        "n_snps": 10000.0,
        "n_bins_used": 20,
        "n_masked_data_zero": 0,
        "n_masked_model_zero": 0,
    }
    m.update(over)
    return m


# --------------------------------------------------------------- metrics --


def test_returns_exactly_the_expected_keys():
    data = np.array([0, 100.0, 80.0, 60.0, 40.0, 20.0, 10.0, 0])
    m = residual_metrics(data.copy(), data)
    assert set(m) == EXPECTED_KEYS


def test_using_dadi_flag_is_a_bool():
    assert isinstance(R.USING_DADI, bool)


def test_positive_residual_means_model_overestimates():
    """dadi returns -1.5*(datatrans-modeltrans)/model**(1/6): model > data is positive."""
    data = np.array([0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 0])
    model = data * 2.0                      # model above data everywhere
    model[0] = model[-1] = 0.0
    m = residual_metrics(model, data, rare_max_count=4)
    assert m["mean_rest"] > 0
    assert m["max_abs_rare"] > 0


def test_perfect_fit_is_all_zero_and_passes():
    data = np.array([0, 500.0, 400.0, 300.0, 200.0, 100.0, 50.0, 0])
    m = residual_metrics(data.copy(), data)
    assert m["max_abs_rare"] == pytest.approx(0.0, abs=1e-9)
    assert m["max_abs_rest"] == pytest.approx(0.0, abs=1e-9)
    assert m["mean_rest"] == pytest.approx(0.0, abs=1e-9)
    assert m["mean_rest_norm"] == pytest.approx(0.0, abs=1e-9)
    assert m["g2"] == pytest.approx(0.0, abs=1e-9)
    assert m["g2_per_snp"] == pytest.approx(0.0, abs=1e-12)
    assert m["n_snps"] == pytest.approx(1550.0)
    assert m["n_bins_used"] == 6
    assert classify(m, CFG)[0] == "PASS"


def test_fixed_categories_are_never_used():
    """Index 0 and the last index are the fixed categories: excluded from everything."""
    data = np.array([1e9, 10.0, 10.0, 10.0, 1e9])
    model = np.array([1.0, 10.0, 10.0, 10.0, 1.0])   # wildly off in the fixed bins
    m = residual_metrics(model, data)
    assert m["n_bins_used"] == 3
    assert m["n_snps"] == pytest.approx(30.0)
    assert m["max_abs_rare"] == pytest.approx(0.0, abs=1e-9)


def test_rare_and_rest_are_split_at_the_configured_count():
    data = np.array([0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 0])
    model = data.copy()
    model[1:5] = 40.0                      # only allele counts 1-4 are off
    m = residual_metrics(model, data, rare_max_count=4)
    assert m["max_abs_rare"] > 1.0
    assert m["max_abs_rest"] == pytest.approx(0.0, abs=1e-9)

    # Same arrays, split moved down to 2: counts 3 and 4 now land in "rest".
    m2 = residual_metrics(model, data, rare_max_count=2)
    assert m2["max_abs_rest"] == pytest.approx(m2["max_abs_rare"], rel=1e-12)


def test_data_zero_bins_are_counted_and_excluded():
    """data == 0 is legitimate (no variants observed): masked, counted, not a failure."""
    data = np.array([0, 10.0, 0.0, 0.0, 10.0, 0])
    model = np.array([0, 10.0, 99.0, 99.0, 10.0, 0])
    m = residual_metrics(model, data)
    assert m["n_masked_data_zero"] == 2
    assert m["n_masked_model_zero"] == 0
    assert m["n_bins_used"] == 2
    assert m["max_abs_rare"] == pytest.approx(0.0, abs=1e-9)
    assert m["g2"] == pytest.approx(0.0, abs=1e-9)
    assert classify(m, CFG)[0] == "PASS"


def test_model_zero_bin_is_counted_and_fails():
    """data > 0 with a non-finite residual: the model gives ~0 probability there."""
    data = np.array([0, 10.0, 10.0, 10.0, 0])
    model = np.array([0, 10.0, 0.0, 10.0, 0])
    m = residual_metrics(model, data)
    assert m["n_masked_model_zero"] == 1
    assert m["n_masked_data_zero"] == 0
    assert m["n_bins_used"] == 2
    status, reasons = classify(m, CFG)
    assert status == "FAIL"
    assert any("variants are observed" in r for r in reasons)


def test_g2_matches_a_hand_computation():
    # One interior bin used: data 10, model 5.
    data = np.array([0, 10.0, 0])
    model = np.array([0, 5.0, 0])
    m = residual_metrics(model, data)
    expected = 2.0 * (10.0 * math.log(10.0 / 5.0) - (10.0 - 5.0))
    assert m["g2"] == pytest.approx(expected)
    assert m["n_snps"] == pytest.approx(10.0)
    assert m["g2_per_snp"] == pytest.approx(expected / 10.0)

    # Two interior bins, one over- and one under-predicted.
    data2 = np.array([0, 10.0, 20.0, 0])
    model2 = np.array([0, 5.0, 40.0, 0])
    expected2 = (2.0 * (10.0 * math.log(10.0 / 5.0) - 5.0)
                 + 2.0 * (20.0 * math.log(20.0 / 40.0) + 20.0))
    m2 = residual_metrics(model2, data2)
    assert m2["g2"] == pytest.approx(expected2)


def test_mean_rest_norm_is_normalised_by_sqrt_n_snps():
    # Rest bins (counts 5,6,7) identical in both cases, so mean_rest is identical.
    # A rare bin fitted perfectly carries the extra SNPs: 1000 -> 4000 total.
    def build(rare_count):
        data = np.array([0, float(rare_count), 1.0, 1.0, 1.0,
                         100.0, 100.0, 100.0, 0])
        model = data.copy()
        model[5:8] = 150.0
        return residual_metrics(model, data, rare_max_count=4)

    small = build(700 - 3)     # n_snps = 697 + 3 + 300 = 1000
    big = build(3700 - 3)      # n_snps = 3697 + 3 + 300 = 4000
    assert small["n_snps"] == pytest.approx(1000.0)
    assert big["n_snps"] == pytest.approx(4000.0)
    assert big["mean_rest"] == pytest.approx(small["mean_rest"])
    assert big["mean_rest_norm"] == pytest.approx(small["mean_rest_norm"] / 2.0)
    # And the definition itself.
    assert small["mean_rest_norm"] == pytest.approx(
        abs(small["mean_rest"]) / math.sqrt(1000.0 / 1e5))


def test_mean_rest_norm_uses_absolute_value():
    data = np.array([0, 10.0, 10.0, 10.0, 10.0, 100.0, 100.0, 0])
    model = data.copy()
    model[5:7] = 50.0                      # model under-predicts the rest bins
    m = residual_metrics(model, data, rare_max_count=4)
    assert m["mean_rest"] < 0
    assert m["mean_rest_norm"] > 0


# -------------------------------------------------------------- classify --


def test_classify_pass():
    status, reasons = classify(base_metrics(), CFG)
    assert status == "PASS"
    assert reasons == []


def test_classify_fail_model_zero_beats_everything():
    status, reasons = classify(base_metrics(n_masked_model_zero=1), CFG)
    assert status == "FAIL"
    assert any("variants are observed" in r for r in reasons)


def test_classify_warns_not_fails_on_extreme_max_abs_rare():
    status, reasons = classify(base_metrics(max_abs_rare=47.2), CFG)
    assert status == "WARN"
    assert any("max_abs_rare" in r for r in reasons)


def test_classify_warns_not_fails_on_extreme_g2_per_snp():
    status, reasons = classify(base_metrics(g2_per_snp=0.0395), CFG)
    assert status == "WARN"
    assert any("g2_per_snp" in r for r in reasons)


def test_classify_warn_on_max_abs_rare():
    status, reasons = classify(base_metrics(max_abs_rare=30.0), CFG)
    assert status == "WARN"
    assert any("max_abs_rare" in r for r in reasons)


def test_classify_warn_on_g2_per_snp():
    status, reasons = classify(base_metrics(g2_per_snp=0.02), CFG)
    assert status == "WARN"
    assert any("g2_per_snp" in r for r in reasons)


def test_classify_warn_on_bulk_bias():
    status, reasons = classify(base_metrics(mean_rest=-3.0, mean_rest_norm=3.0), CFG)
    assert status == "WARN"
    assert any("bias" in r.lower() for r in reasons)
    assert any("underestimate" in r for r in reasons)


def test_classify_warn_on_rest_worse_than_rare():
    status, reasons = classify(base_metrics(max_abs_rare=3.0, max_abs_rest=8.0), CFG)
    assert status == "WARN"
    assert any("max_abs_rest" in r for r in reasons)


def test_classify_rest_over_rare_needs_both_conditions():
    # Bigger than rare but below rest_over_rare_min: not a warning.
    assert classify(base_metrics(max_abs_rare=1.0, max_abs_rest=4.0), CFG)[0] == "PASS"
    # Above rest_over_rare_min but smaller than rare: not a warning either.
    assert classify(base_metrics(max_abs_rare=20.0, max_abs_rest=6.0), CFG)[0] == "PASS"


def test_classify_thresholds_are_strict_inequalities():
    assert classify(base_metrics(max_abs_rare=CFG["max_rare_warn"]), CFG)[0] == "PASS"
    assert classify(base_metrics(max_abs_rare=CFG["max_rare_fail"]), CFG)[0] == "WARN"
    assert classify(base_metrics(g2_per_snp=CFG["g2_warn"]), CFG)[0] == "PASS"
    assert classify(base_metrics(mean_rest_norm=CFG["bulk_bias_warn"]), CFG)[0] == "PASS"


def test_classify_reports_several_reasons_at_once():
    status, reasons = classify(
        base_metrics(max_abs_rare=30.0, g2_per_snp=0.02, mean_rest_norm=5.0), CFG)
    assert status == "WARN"
    assert len(reasons) == 3


def test_everything_masked_fails_rather_than_passing():
    """A fit where no bin could be measured must not come out clean."""
    cfg = dict(max_rare_warn=25, max_rare_fail=40, bulk_bias_warn=2.3,
               rest_over_rare_min=5, g2_warn=0.008, g2_fail=0.035)
    m = {"max_abs_rare": 0.0, "max_abs_rest": 0.0, "mean_rest": 0.0,
         "mean_rest_norm": 0.0, "g2_per_snp": 0.0, "n_masked_model_zero": 0,
         "n_bins_used": 0}
    status, reasons = classify(m, cfg)
    assert status == "FAIL"
    assert any("no usable bins" in r for r in reasons)


def test_structural_problems_still_fail_while_residuals_only_warn():
    """FAIL is reserved for broken data; residual thresholds are advisory."""
    cfg = dict(max_rare_warn=25, max_rare_fail=40, bulk_bias_warn=2.3,
               rest_over_rare_min=5, g2_warn=0.008, g2_fail=0.035)
    wild = {"max_abs_rare": 999.0, "max_abs_rest": 0.0, "mean_rest": 0.0,
            "mean_rest_norm": 0.0, "g2_per_snp": 99.0, "n_masked_model_zero": 0,
            "n_bins_used": 70}
    assert classify(wild, cfg)[0] == "WARN"
    broken = dict(wild, n_masked_model_zero=3)
    assert classify(broken, cfg)[0] == "FAIL"
