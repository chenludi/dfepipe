#!/bin/bash
# demo, the pick: choose the maximum-likelihood fit, write best_inference.txt,
# and draw the diagnostics.
#
#   bash scripts/demo/run_bestdemo.sh <bed> [bed...]
#
# <bed> is the ORIGINAL bed, the same argument prep, sfs and run_demo.sh took.
# Run it once every array task of run_demo.sh has finished.
#
# Ported from piperef_stddfe_4_plotdemo.sh and plot_bestll_demo_args_0421.py.
# The selection rule, unchanged:
#
#   1. every *.InferDM.bestfits in the run directory contributes the FIRST row
#      of its "# Converged results" block, keyed by the p0 marker in its name
#      (plot_bestll_demo_args_0421.py:139-180). A file with no converged block
#      contributes nothing.
#   2. Na = theta / (2 * mutation_ctrl_diploid), where mutation_ctrl_diploid =
#      2 * u * L comes from the control fsmut table's genome-wide row:
#      u = ROULETTE_SCALE * avepersite_diploi_Roulette_mutrate / 2 and
#      L = length_adj_inputbed (:96-114).
#   3. fits with ll <= the cutoff are dropped. The cutoff is DEMO_LL_CUTOFF,
#      given positive and negated (500 in piperef_stddfe_4_plotdemo.sh:43); if
#      it is stricter than the best fit it falls back to max(ll) - 100 (:223-224).
#   4. the winner is the largest ll among what is left; best_inference.txt
#      records it as `key:value` lines, in the historical order, and the later
#      stages grep it for `max_ll_key:` and `Na:`.
#
# Differences from the original, all deliberate:
#   * the parameter names come from the .bestfits header line dadi-cli writes,
#     not from a per-model table in the script. plot_bestll_demo_args_0421.py:119-122
#     lists bottlegrowth_1d's parameters as (nu, T), which is wrong - it has
#     three - and would abort on that model at :163.
#   * figure names drop the cutoff value. The original wrote
#     `ll_nuB_619.8079138965113.pdf`; the cutoff is recorded in
#     best_inference.txt, so the file name does not have to carry it.
#   * the model-vs-data figure is named sfs_bestfit.pdf rather than
#     `<the whole bestfits file name>.pdf`, and is written into the run
#     directory rather than copied to a dated ~/mv_<date> viewing directory.
#     The original's Plot call did run and its figures are on disk (667
#     *.InferDM.bestfits.pdf, including the Tx run's); only the name differs.
#   * the best-fit model spectrum is written out as bestfit_model.fs, which is
#     what check_demo.py measures the residuals on.
#   * PDF only; PNG=1 additionally writes PNG.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# shellcheck source=../../config/env.sh
source "${REPO_DIR}/config/env.sh"

: "${OUT:?OUT not set}"
: "${FS_RUN:?FS_RUN not set}"
: "${DEMO_RUN:?DEMO_RUN not set}"
: "${DEMO_MODEL:?DEMO_MODEL not set}"
: "${CTRL_TAG:?CTRL_TAG not set}"
: "${NUM_CONVERGENCE:?NUM_CONVERGENCE not set}"
: "${RESID_RARE_MAX_COUNT:?RESID_RARE_MAX_COUNT not set}"
: "${RESID_MAX_RARE_WARN:?RESID_MAX_RARE_WARN not set}"
: "${RESID_MAX_RARE_FAIL:?RESID_MAX_RARE_FAIL not set}"
: "${RESID_BULK_BIAS_WARN:?RESID_BULK_BIAS_WARN not set}"
: "${RESID_REST_OVER_RARE_MIN:?RESID_REST_OVER_RARE_MIN not set}"
: "${G2_PER_SNP_WARN:?G2_PER_SNP_WARN not set}"
: "${G2_PER_SNP_FAIL:?G2_PER_SNP_FAIL not set}"

# Positive here, as in piperef_stddfe_4_plotdemo.sh:43; negated internally.
: "${DEMO_LL_CUTOFF:?DEMO_LL_CUTOFF not set}"
: "${ROULETTE_SCALE:?ROULETTE_SCALE not set}"

if [[ $# -eq 0 ]]; then
  echo "usage: run_bestdemo.sh <bed> [bed...]" >&2
  exit 2
fi

case "$DEMO_MODEL" in
  three_epoch)  dadi_model="three_epoch" ;;
  two_epoch)    dadi_model="two_epoch" ;;
  growth)       dadi_model="growth" ;;
  bottlegrowth) dadi_model="bottlegrowth_1d" ;;
  *) echo "run_bestdemo: DEMO_MODEL='${DEMO_MODEL}' is not supported" >&2; exit 2 ;;
esac

activate_conda

failed=0
for bed in "$@"; do
  base="$(basename "$bed" .bed)"
  # The spectrum and mutation table names, as run_concat.sh builds them.
  shortbase="s.nc_exHET_noq.${base}"
  shortbase="${shortbase/s.nc_/}"
  shortbase="${shortbase/nc_/}"
  fsbase="quiecent_ctrl.${shortbase}_${CTRL_TAG}.chrs_1t22"
  sfs_dir="${OUT}/sfs/${base}/${FS_RUN}"
  data_fs="${sfs_dir}/${fsbase}.fs"
  mutfile="${sfs_dir}/fsmut.${fsbase}.tab"
  d="${OUT}/demo/${base}/${FS_RUN}-${DEMO_MODEL}-${DEMO_RUN}"

  for f in "$data_fs" "$mutfile"; do
    [[ -s "$f" ]] || {
      echo "run_bestdemo: missing sfs output for '${base}': ${f}" >&2
      echo "      Run scripts/sfs/run_concat.sh first, or check FS_RUN." >&2
      failed=1; continue 2; }
  done
  if ! compgen -G "${d}/*.InferDM.bestfits" > /dev/null; then
    echo "run_bestdemo: no fits for '${base}' in ${d}" >&2
    echo "      Run scripts/demo/run_demo.sh (or scripts/sge/submit_demo.sh)" >&2
    echo "      first, or check DEMO_RUN." >&2
    failed=1; continue
  fi

  echo "run_bestdemo: ${base}: picking the maximum-likelihood fit in ${d}"
  if ! python3 - "$d" "$dadi_model" "$DEMO_MODEL" "$data_fs" "$mutfile" \
                "$DEMO_LL_CUTOFF" "${PNG:-0}" "$ROULETTE_SCALE" <<'PY'
"""Pick the maximum-likelihood demographic fit, write best_inference.txt, plot.

A faithful port of plot_bestll_demo_args_0421.py; the differences are listed in
run_bestdemo.sh's header.
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")                     # no display on a compute node
import matplotlib.pyplot as plt           # noqa: E402
import dadi                               # noqa: E402

demo_dir, dadi_model, model, data_fs, mutfile, cutoff_arg, png = sys.argv[1:8]
want_png = png not in ("", "0", "no", "false")

# Roulette stores per-site rates scaled; ROULETTE_SCALE (pipeline.conf) converts
# a table value to a diploid per-site rate (plot_bestll_demo_args_0421.py:97).
# Passed in rather than hardcoded, so every stage uses one value.
ROULETTE_SCALE = float(sys.argv[8])


def die(msg):
    sys.stderr.write("run_bestdemo: %s\n" % msg)
    sys.exit(1)


def save(stem):
    """PDF always, PNG as well under PNG=1."""
    plt.savefig(os.path.join(demo_dir, stem + ".pdf"))
    if want_png:
        plt.savefig(os.path.join(demo_dir, stem + ".png"), dpi=200)
    plt.close()


def plain_ticks():
    """Na runs into the tens of thousands; scientific notation hides the value."""
    try:
        plt.ticklabel_format(style="plain")
    except AttributeError:
        pass


# ------------------------------------------------- the control mutation count --
# mutation_ctrl_diploid = 2*u*L over the CONTROL sequence; theta = 4*Na*u*L, so
# Na = theta / (2 * mutation_ctrl_diploid). The genome-wide row of the
# concatenated fsmut table is the one whose fsfile column ends in ".etal";
# the per-chromosome rows follow it.
with open(mutfile) as fh:
    header = fh.readline().rstrip("\n").split("\t")
    rows = [dict(zip(header, line.rstrip("\n").split("\t")))
            for line in fh if line.strip()]
etal = [r for r in rows if "etal" in r.get("fsfile", "")]
if not etal:
    die("no genome-wide ('.etal') row in %s: it holds only per-chromosome rows, "
        "so the control mutation count cannot be formed" % mutfile)
try:
    length = float(etal[0]["length_adj_inputbed"])
    rate = float(etal[0]["avepersite_diploi_Roulette_mutrate"])
except (KeyError, ValueError) as exc:
    die("cannot read length_adj_inputbed / avepersite_diploi_Roulette_mutrate "
        "from %s (%s)" % (mutfile, exc))
u = ROULETTE_SCALE * rate / 2
mutation_ctrl_diploid = 2 * (u * length)
if mutation_ctrl_diploid <= 0:
    die("the control mutation count came out as %r; check %s"
        % (mutation_ctrl_diploid, mutfile))

# --------------------------------------------------------- collect the fits --
# One entry per .bestfits file, keyed by the p0 marker in its name (the third
# field from the end: <model>.<fs base>.<marker>.InferDM.bestfits).
files = sorted(glob.glob(os.path.join(demo_dir, "*.InferDM.bestfits")))
data_dict = {}
names = None
for path in files:
    with open(path) as fh:
        lines = fh.read().splitlines()
    found = False
    for line in lines:
        stripped = line.strip()
        if stripped == "# Converged results":
            found = True
            continue
        if not found:
            continue
        if stripped.startswith("# Log(likelihood)"):
            # dadi-cli names the columns; that is the authority on this model's
            # parameters, so no per-model table is needed here.
            names = ["ll"] + stripped.split("\t")[1:]
            continue
        if stripped.startswith("#") or not stripped:
            if stripped == "#":
                break
            continue
        try:
            values = [float(v) for v in stripped.split("\t")]
        except ValueError:
            break
        if names is None or len(values) != len(names):
            sys.stderr.write(
                "run_bestdemo: %s has %d values but %s columns; skipping\n"
                % (os.path.basename(path), len(values),
                   len(names) if names else "no"))
            break
        key = os.path.basename(path).split(".")[-3]
        data_dict[key] = dict(zip(names, values))
        break

if not data_dict:
    die("none of the %d .bestfits file(s) in %s holds a converged result"
        % (len(files), demo_dir))
if "theta" not in names:
    die("the .bestfits header has no theta column: %s" % names)

for entry in data_dict.values():
    entry["Na"] = entry["theta"] / (2 * mutation_ctrl_diploid)

lls = [e["ll"] for e in data_dict.values()]
nas = [e["Na"] for e in data_dict.values()]

# ------------------------------------------------------------- ll ~ Na, all --
plt.scatter(nas, lls)
plain_ticks()
plt.xlabel("Na")
plt.ylabel("ll")
plt.title("%s: ll ~ Na (%d converged fits)" % (model, len(lls)))
save("ll_Na")

# ------------------------------------------------------------------- cutoff --
# Given positive, used negative. If it is stricter than the best fit nothing
# would survive, so fall back to max(ll) - 100, as the original did.
cutoff_input = float(cutoff_arg)
cutoff = -cutoff_input if cutoff_input > 0 else cutoff_input
if cutoff > max(lls):
    cutoff = max(lls) - 100
    sys.stderr.write(
        "run_bestdemo: the cutoff %g is stricter than the best fit %g; using "
        "max(ll) - 100 = %g instead\n" % (-cutoff_input, max(lls), cutoff))

kept = {k: v for k, v in data_dict.items() if v["ll"] > cutoff}
lls2 = sorted((e["ll"] for e in kept.values()), reverse=True)
best_key = max(kept, key=lambda k: kept[k]["ll"])
best = kept[best_key]

# -------------------------------------------------------- best_inference.txt --
# `key:value` per line, in the historical order. The order is cosmetic; the key
# names are not - the cache and DFE stages grep this file for `max_ll_key:` and
# `Na:`, and `dict(x.split(":") for x in f)` is how it is meant to be read back.
out = {
    "dir_outdemog": demo_dir,
    "model": model,
    "likelyhood_cutoff": cutoff_input,
    "mutfile": mutfile,
    "mutation_ctrl_diploid": mutation_ctrl_diploid,
    "max_ll_key": best_key,
}
out.update(best)
with open(os.path.join(demo_dir, "best_inference.txt"), "w") as fh:
    for key, value in out.items():
        fh.write("%s:%s\n" % (key, value))

# Kept for parity with the original; nothing downstream is required to read it.
with open(os.path.join(demo_dir, "maxll_demoparam_file.txt"), "w") as fh:
    fh.write("dir_outdemog model bestfit_demogfile mutation_ctrl_diploid theta\n")
    fh.write("%s %s %s%s.InferDM.bestfits %s %s\n"
             % (demo_dir, model, model, best_key, mutation_ctrl_diploid,
                best["theta"]))
print("run_bestdemo: max_ll_key=%s ll=%s Na=%s" % (best_key, best["ll"], best["Na"]))

# ------------------------------------------------------------ ll ~ parameter --
second = kept[sorted(kept, key=lambda k: kept[k]["ll"], reverse=True)[1]] \
    if len(lls2) > 1 else None
params = [n for n in names if n not in ("ll", "theta")]


def circle(entry, x, y, label):
    plt.scatter(entry[x], entry[y], edgecolors="red", marker="o",
                linewidths=1, facecolors="none", label=label)


for param in ["Na"] + params:
    xs = [e[param] for e in kept.values()]
    ys = [e["ll"] for e in kept.values()]
    plt.scatter(xs, ys)
    plain_ticks()
    if second is not None:
        circle(second, param, "ll", "2nd best ll")
    circle(best, param, "ll", "best ll")
    plt.xlabel(param)
    plt.ylabel("ll")
    plt.title("%s: ll ~ %s (ll > %.1f)" % (model, param, cutoff))
    save("ll_%s" % param)

for param in params:
    xs = [e[param] for e in kept.values()]
    ys = [e["Na"] for e in kept.values()]
    sc = plt.scatter(xs, ys, c=[e["ll"] for e in kept.values()],
                     cmap="viridis", alpha=0.5)
    if second is not None:
        circle(second, param, "Na", "2nd best ll")
    circle(best, param, "Na", "best ll")
    plt.colorbar(sc).set_label("log likelihood")
    plt.xlabel(param)
    plt.ylabel("Na")
    plt.title("%s: Na ~ %s" % (model, param))
    save("Na_%s" % param)

# ------------------------------------------------------- the model spectrum --
# The best fit's .bestfits, globbed by marker exactly as
# piperef_stddfe_4_plotdemo.sh:51 does.
hits = glob.glob(os.path.join(demo_dir, "*.%s.InferDM.bestfits" % best_key))
if not hits:
    die("no .bestfits file for the winning marker %s in %s" % (best_key, demo_dir))
bestfits = hits[0]

# The model spectrum as a file, for check_demo.py to measure residuals on. The
# parameters are read and the grids chosen by dadi-cli's own helpers, and the
# result is folded and optimally scaled to the data. That matches what the figure
# shows: dadi-cli's Plot.py never folds, but the two dadi functions it calls both
# do - Inference.optimal_sfs_scaling:664-665 folds before computing the scaling,
# and Plotting.plot_1d_comp_Poisson folds before plotting - so the residuals here
# and the residuals in sfs_bestfit.pdf are measured on the same array.
from dadi_cli.Models import get_model                                  # noqa: E402
from dadi_cli.utilities import get_opts_and_theta, pts_l_func          # noqa: E402

func, _ = get_model(dadi_model, None)
popt, _theta = get_opts_and_theta(bestfits)
observed = dadi.Spectrum.from_file(data_fs)
func_ex = dadi.Numerics.make_extrap_func(func)
expected = func_ex(popt, observed.sample_sizes, pts_l_func(observed.sample_sizes))
if observed.folded and not expected.folded:
    expected = expected.fold()
dadi.Inference.optimally_scaled_sfs(expected, observed).to_file(
    os.path.join(demo_dir, "bestfit_model.fs"))
PY
  then
    echo "run_bestdemo: ${base}: best-fit selection failed" >&2
    failed=1
    continue
  fi

  # The model-vs-data spectrum, from dadi-cli itself. The winning .bestfits is
  # recovered through best_inference.txt's max_ll_key, exactly as
  # piperef_stddfe_4_plotdemo.sh:50-52 does.
  key="$(awk -F: '/^max_ll_key:/ {print $2}' "${d}/best_inference.txt")"
  bestfits="$(ls -1 "${d}"/*."${key}".InferDM.bestfits 2>/dev/null | head -1 || true)"
  if [[ -z "$bestfits" ]]; then
    echo "run_bestdemo: ${base}: no .bestfits for max_ll_key=${key}" >&2
    failed=1
    continue
  fi
  dadi-cli Plot --fs "$data_fs" --demo-popt "$bestfits" \
    --output "${d}/sfs_bestfit.pdf" --model "$dadi_model" --nomisid
  if [[ "${PNG:-0}" != "0" ]]; then
    dadi-cli Plot --fs "$data_fs" --demo-popt "$bestfits" \
      --output "${d}/sfs_bestfit.png" --model "$dadi_model" --nomisid
  fi
  echo "run_bestdemo: ${base}: wrote ${d}/sfs_bestfit.pdf"
done

# ---------------------------------------------------------------- checkpoint --
# Diagnostic only: a bug in the check must never discard a finished stage, so its
# failure warns and the stage's own exit status is left alone.
if ! python3 "${REPO_DIR}/scripts/checks/check_demo.py" \
       --out "${OUT}/demo" --run "${FS_RUN}-${DEMO_MODEL}-${DEMO_RUN}" \
       --sfs-out "${OUT}/sfs" --fs-run "${FS_RUN}" \
       --rare-max-count "${RESID_RARE_MAX_COUNT}" \
       --max-rare-warn "${RESID_MAX_RARE_WARN}" \
       --max-rare-fail "${RESID_MAX_RARE_FAIL}" \
       --bulk-bias-warn "${RESID_BULK_BIAS_WARN}" \
       --rest-over-rare-min "${RESID_REST_OVER_RARE_MIN}" \
       --g2-warn "${G2_PER_SNP_WARN}" --g2-fail "${G2_PER_SNP_FAIL}" \
       --num-convergence "${NUM_CONVERGENCE}" "$@"; then
  echo "WARNING: check_demo.py failed; the demo output above is unaffected." >&2
  echo "         Inspect best_inference.txt and sfs_bestfit.pdf by hand before" >&2
  echo "         trusting this run." >&2
fi

exit "$failed"
