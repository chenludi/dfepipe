"""Names and messages must not claim more than the code does.

Two real defects motivated this file, both found by the user reading output
rather than by any test:

  * `RESID_MAX_RARE_FAIL` and `G2_PER_SNP_FAIL` were named as though crossing
    them stopped a run. Nothing stopped: residuals.classify only ever appends to
    `warns` for these. They are now `..._LARGE` and `..._BADFIT`.
  * the warning read "beyond anything seen in past runs (>40.0)". That threshold
    is a percentile of gamma and neugamma fits only - the calibration globbed
    `*_model_sfs.fs`, and the `_neu` run directories hold none - so for the
    neutral null the sentence asserted a comparison that was never made.

Both were wording, not arithmetic, which is exactly why no test caught them.
"""
import pathlib
import re
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts" / "checks"))
from residuals import classify  # noqa: E402

CONF = (ROOT / "config" / "pipeline.conf").read_text()
RESID = (ROOT / "scripts" / "checks" / "residuals.py").read_text()

CFG = dict(max_rare_large=25.0, max_rare_badfit=40.0, g2_large=0.008,
           g2_badfit=0.035, bulk_bias_warn=2.3, rest_over_rare_min=5.0)


def metrics(**kw):
    base = dict(max_abs_rare=1.0, max_abs_rest=1.0, mean_rest=0.0,
                mean_rest_norm=0.0, g2_per_snp=0.0, n_masked_model_zero=0,
                n_bins_used=80)
    base.update(kw)
    return base


# ------------------------------------------------- a _FAIL key must be able to fail --

def test_no_residual_threshold_is_named_FAIL():
    """A key named FAIL that only warns is the defect this file exists for."""
    offenders = [k for k in re.findall(r'^(RESID_\w+|G2_\w+)=', CONF, re.M)
                 if k.endswith("_FAIL")]
    assert offenders == [], (
        "these only ever warn, so the name must not say FAIL: %s" % offenders)


@pytest.mark.parametrize("key,metric,value", [
    ("max_rare_badfit", "max_abs_rare", 1e6),
    ("g2_badfit", "g2_per_snp", 1e6),
])
def test_crossing_the_badfit_threshold_warns_and_never_fails(key, metric, value):
    """However extreme the residual, it stays advisory. FAIL is reserved for
    structural problems, which is the invariant the old name contradicted."""
    status, reasons = classify(metrics(**{metric: value}), CFG)
    assert status == "WARN", reasons
    assert reasons


def test_a_structural_problem_still_fails():
    """The reserved FAIL must still work, or the distinction is meaningless."""
    status, reasons = classify(metrics(n_masked_model_zero=3), CFG)
    assert status == "FAIL"
    assert any("no usable probability" in r for r in reasons)


# ------------------------------------------------------- no unsupported claims --

def test_no_message_claims_a_comparison_with_past_runs():
    """The thresholds have no calibration set for the neutral null, and the
    package does not ship the calibration data for any model, so no message may
    assert where past runs landed."""
    for phrase in ("beyond anything seen in past runs",
                   "unusual compared with"):
        assert phrase not in RESID, (
            "%r claims a comparison the code cannot support" % phrase)


def test_the_badfit_message_says_what_it_means():
    """The user asked for wording a reader understands without the source."""
    _, reasons = classify(metrics(max_abs_rare=49.6), CFG)
    joined = " ".join(reasons)
    assert "residual is large" in joined
    assert "not a good fit" in joined


@pytest.mark.parametrize("metric,value,expect", [
    ("max_abs_rare", 49.6, "rare-bin residual is large"),
    ("max_abs_rare", 30.0, "rare-bin residual is on the large side"),
    ("g2_per_snp", 0.1681, "deviance per SNP is large"),
    ("g2_per_snp", 0.01, "deviance per SNP is on the large side"),
])
def test_both_tiers_are_distinguishable_in_words(metric, value, expect):
    """Two tiers that read identically would be worse than one."""
    _, reasons = classify(metrics(**{metric: value}), CFG)
    assert any(expect in r for r in reasons), reasons


def test_no_message_leads_with_an_internal_variable_name():
    """`max_abs_rest 8.0 is worse than max_abs_rare 3.0` told the reader nothing
    unless they had the source open."""
    cases = [metrics(max_abs_rare=49.6), metrics(g2_per_snp=0.1681),
             metrics(mean_rest_norm=13.76, mean_rest=1.0),
             metrics(max_abs_rare=3.0, max_abs_rest=8.0)]
    for m in cases:
        _, reasons = classify(m, CFG)
        for r in reasons:
            first = r.split()[0].rstrip(":,")
            assert "_" not in first, "message opens with a variable name: %r" % r


# ------------------------------------------- no check may overclaim its calibration --

CHECK_FILES = sorted((ROOT / "scripts" / "checks").glob("*.py"))


def test_the_check_files_were_found():
    assert len(CHECK_FILES) >= 10, [p.name for p in CHECK_FILES]


@pytest.mark.parametrize("path", CHECK_FILES, ids=lambda p: p.name)
def test_no_check_claims_a_past_run_comparison(path):
    """Four docstrings asserted that a WARN meant "unusual compared with past
    runs". The calibration set is 137 gamma and neugamma DFE fits: it does not
    cover the demographic fit, the cache, the neutral null, or any of the
    judgement thresholds, so the claim was wrong in each of those files."""
    text = path.read_text()
    for phrase in ("unusual compared with", "beyond anything seen in past runs",
                   "percentiles of past runs"):
        assert phrase not in text, "%s: %r" % (path.name, phrase)


def test_every_conf_fallback_matches_the_configured_value():
    """conf.py promises a threshold "exists in exactly one place". A fallback that
    disagrees with the config means a hand-run check and a pipeline-run check use
    different numbers - HASQUI_MIN_FRAC fell back to 0.90 after the configured
    value moved to 0.75."""
    conf_text = (ROOT / "config" / "pipeline.conf").read_text()
    mismatches = []
    for path in CHECK_FILES:
        for key, fallback in re.findall(
                r'conf\.get(?:_float|_int)\("([A-Z0-9_]+)",\s*([0-9.eE+-]+)\)',
                path.read_text()):
            m = re.search(r'^%s="\$\{%s:-([^}]*)\}"' % (key, key), conf_text, re.M)
            if not m:
                continue
            if float(fallback) != float(m.group(1)):
                mismatches.append("%s: %s fallback %s but config has %s"
                                  % (path.name, key, fallback, m.group(1)))
    assert mismatches == [], "\n".join(mismatches)
