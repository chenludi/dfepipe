import numpy as np
import pytest
from scripts.checks.residuals import residual_metrics, classify


def test_positive_residual_means_model_overestimates():
    """dadi returns -1.5*(datatrans-modeltrans)/model**(1/6), so model>data is positive."""
    data = np.array([0, 100.0, 100.0, 0])
    model = np.array([0, 200.0, 200.0, 0])   # model above data everywhere
    m = residual_metrics(model, data)
    assert m["mean_resid"] > 0


def test_perfect_fit_is_all_zero():
    data = np.array([0, 50.0, 40.0, 30.0, 20.0, 0])
    m = residual_metrics(data.copy(), data)
    assert m["max_abs_resid"] == pytest.approx(0.0, abs=1e-9)
    assert m["rms_resid"] == pytest.approx(0.0, abs=1e-9)


def test_rare_and_rest_are_split_at_count_four():
    data = np.array([0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 0])
    model = data.copy()
    model[1:5] = 40.0            # counts 1-4 badly off
    m = residual_metrics(model, data, rare_max_count=4)
    assert m["max_abs_resid_rare"] > m["max_abs_resid_rest"]
    assert m["max_abs_resid_rest"] == pytest.approx(0.0, abs=1e-9)


def test_longest_same_sign_run():
    data = np.array([0, 10.0, 10.0, 10.0, 10.0, 10.0, 0])
    model = np.array([0, 20.0, 20.0, 20.0, 5.0, 20.0, 0])   # +++-+
    m = residual_metrics(model, data)
    assert m["longest_same_sign_run"] == 3


def test_masked_bins_are_counted_and_excluded():
    """data == 0 is masked by dadi's Anscombe residual, so it must not enter the metrics."""
    data = np.array([0, 10.0, 0.0, 10.0, 0])
    model = np.array([0, 10.0, 99.0, 10.0, 0])
    m = residual_metrics(model, data)
    assert m["n_masked_bins"] == 1
    assert m["max_abs_resid"] == pytest.approx(0.0, abs=1e-9)


def test_classify_thresholds():
    assert classify({"max_abs_resid": 1.0, "mean_resid": 0.1}, 2, 10, 0.5)[0] == "PASS"
    assert classify({"max_abs_resid": 5.0, "mean_resid": 0.1}, 2, 10, 0.5)[0] == "WARN"
    assert classify({"max_abs_resid": 20.0, "mean_resid": 0.1}, 2, 10, 0.5)[0] == "FAIL"
    status, reasons = classify({"max_abs_resid": 1.0, "mean_resid": 3.0}, 2, 10, 0.5)
    assert status == "WARN" and any("overestimate" in r for r in reasons)
