#!/usr/bin/env python3
"""Goal tests for a built environment.

These do not ask "is dadi 2.3.3 installed". They ask "does this environment do
what we need it to do": the custom pdfs are registered and numerically correct,
and the external tools the pipeline shells out to are callable.

That framing is deliberate. When a future dadi release appears, run this against
it: if everything passes, the vendored patch has been made redundant by upstream
and the pin can move. If it fails, keep the vendored dadi 2.3.3. Reproducing the
published results always uses the vendored copy either way.

    python3 test_env.py [--prefix ENV_PREFIX]

Exit 0 if every check passes, 1 otherwise, with one line per check.
"""
import argparse
import os
import shutil
import subprocess
import sys

# A fixed input with values computed from the patched neugamma in dadi0116.
# gamma(shape=0.2, scale=1000) mixed with pneu=0.3 of neutral mass below 1e-4.
NEUGAMMA_X = [1e-5, 1e-3, 1.0, 100.0]
NEUGAMMA_PARAMS = (0.2, 1000.0, 0.3)

RESULTS = []


def check(name, fn):
    try:
        detail = fn()
        RESULTS.append((True, name, detail or ""))
    except Exception as exc:                       # noqa: BLE001 - report everything
        RESULTS.append((False, name, f"{type(exc).__name__}: {exc}"))


def _import_dadi():
    import dadi
    return dadi


def check_versions():
    import dadi
    import dadi_cli
    return f"dadi {dadi.__version__}, dadi_cli {getattr(dadi_cli, '__version__', 'unknown')}"


def check_pdfs_registered():
    """The patch's whole purpose: the custom pdfs are dispatchable by name."""
    from dadi_cli.Pdfs import get_dadi_pdf
    missing = [n for n in ("neugamma", "neugammazero", "gammaneud1")
               if get_dadi_pdf(n) is None]
    if missing:
        raise AssertionError(f"not registered in dadi-cli: {', '.join(missing)}")
    return "neugamma, neugammazero, gammaneud1 all dispatchable"


def check_neugamma_numeric():
    """neugamma must be the additive mixture, not something upstream renamed."""
    import numpy as np
    import scipy.stats.distributions as ssd
    from dadi.DFE import PDFs

    xx = np.array(NEUGAMMA_X)
    got = np.asarray(PDFs.neugamma(xx, NEUGAMMA_PARAMS), dtype=float)

    alpha, beta, pneu = NEUGAMMA_PARAMS
    want = (1 - pneu) * ssd.gamma.pdf(xx, alpha, scale=beta)
    want = want + np.where(xx < 1e-4, pneu / 1e-4, 0.0)

    if not np.allclose(got, want, rtol=1e-10, atol=0):
        raise AssertionError(f"neugamma returned {got}, expected {want}")
    return "additive mixture, neutral spike below 1e-4"


def check_neugammazero_is_disjoint():
    """neugammazero renormalises instead of adding: the two must differ."""
    import numpy as np
    from dadi.DFE import PDFs

    xx = np.array(NEUGAMMA_X)
    a = np.asarray(PDFs.neugamma(xx, NEUGAMMA_PARAMS), dtype=float)
    b = np.asarray(PDFs.neugammazero(xx, NEUGAMMA_PARAMS), dtype=float)
    if np.allclose(a, b):
        raise AssertionError("neugammazero is identical to neugamma; the patch is wrong")
    return "distinct from neugamma, as intended"


def check_residual_uses_dadi():
    """The pipeline's own metrics must be running on dadi, not the fallback."""
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.dirname(here))
    from checks.residuals import USING_DADI
    if not USING_DADI:
        raise AssertionError("residuals.py fell back to its local formula; dadi not importable")
    return "residuals.py is using dadi's own implementation"


def check_tool(name):
    def fn():
        path = shutil.which(name)
        if not path:
            raise AssertionError(f"{name} not on PATH (module not loaded?)")
        return path
    return fn


def check_cli_runs(prefix):
    def fn():
        exe = os.path.join(prefix, "bin", "dadi-cli") if prefix else "dadi-cli"
        out = subprocess.run([exe, "Pdf", "--names"], capture_output=True, text=True)
        if "neugamma" not in (out.stdout + out.stderr):
            raise AssertionError("`dadi-cli Pdf --names` does not list neugamma")
        return "dadi-cli lists neugamma"
    return fn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", default="", help="conda env prefix, for the dadi-cli check")
    a = ap.parse_args()

    check("dadi importable", _import_dadi)
    check("versions", check_versions)
    check("custom pdfs registered", check_pdfs_registered)
    check("neugamma numerically correct", check_neugamma_numeric)
    check("neugammazero distinct", check_neugammazero_is_disjoint)
    check("dadi-cli lists neugamma", check_cli_runs(a.prefix))
    check("residuals.py uses dadi", check_residual_uses_dadi)
    check("bedtools on PATH", check_tool("bedtools"))
    check("Rscript on PATH", check_tool("Rscript"))

    width = max(len(n) for _, n, _ in RESULTS)
    failed = 0
    for ok, name, detail in RESULTS:
        print(f"{'PASS' if ok else 'FAIL'}  {name:<{width}}  {detail}")
        failed += 0 if ok else 1
    if failed:
        print(f"\n{failed} check(s) failed", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
