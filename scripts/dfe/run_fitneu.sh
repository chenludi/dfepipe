#!/bin/bash
# dfe, the neutral null: the log-likelihood of the TEST spectrum under the
# demographic model alone, with no selection.
#
#   bash scripts/dfe/run_fitneu.sh <bed> [bed...]
#
# <bed> is the ORIGINAL bed, the same argument prep, sfs and demo took. This is
# the `neu` entry of PDFS and the null the summary stage's model comparison
# subtracts from the gamma / neugamma likelihoods, so it must run for every bed
# whose DFE is fitted.
#
# Ported from piperef_stddfe_5_fitneu.sh -> dadi_cal_ll_neu_0723.py. The
# calculation, unchanged:
#
#   1. read the winning demographic fit's parameters and theta from the
#      .bestfits file named by best_inference.txt's max_ll_key
#      (dadi_cal_ll_neu_0723.py:96-135, piperef_stddfe_5_fitneu.sh:56-57)
#   2. mutratio = u_test*L_test / u_ctrl*L_ctrl from the two fsmut tables'
#      genome-wide (".etal") rows, u = 1.015e-7 * rate / 2         (:154-171)
#   3. the expected spectrum from the demographic model on dadi-cli's CACHE
#      grid sizes (2.2n+2, 2.4n+4, 2.6n+6), scaled by
#      theta_test = theta_demog * mutratio                         (:138-210)
#   4. ll = dadi.Inference.ll(expected, test spectrum)             (:214)
#
# It writes `<test fs base>.dadi_cal_ll_neu_0723.log`, whose `ll: ` line is what
# the summary stage greps; the name is the historical one verbatim.
#
# The expected spectrum is handed to dadi unfolded while the test spectrum is
# folded (dadi builds it with polarized=False). That is safe, and was checked
# rather than assumed: dadi.Inference.ll_per_bin begins
# `if data.folded and not model.folded: model = model.fold()`, so the comparison
# is folded either way. Verified on this bed - unfolded and pre-folded inputs
# give the same ll to 1e-11 - so the null IS comparable with the DFE
# likelihoods, which dadi-cli folds the same way inside InferDFE.
#
# Differences from the original, all deliberate:
#   * no `rm -r $dirout` (piperef_stddfe_5_fitneu.sh:51). Results are keyed by
#     DFE_RUN and are never deleted implicitly; a non-empty run directory stops
#     the stage unless FORCE=1, as in every other stage.
#   * nothing is copied to a dated ~/mv_<date> viewing directory.
#   * the figure is `sfs_neu.pdf`, not `ll_<round(ll)>-<fs base>.pdf`; PDF only,
#     PNG=1 additionally writes PNG.
#   * the folded expected spectrum is written as `model_sfs_rank1.fs` so
#     check_dfe.py can measure residuals on it, the same name run_bestdfe.sh
#     uses for a fitted DFE. Folded, so the residuals are bin-for-bin against
#     the folded data - the same folding dadi applies inside ll().
#   * the original's figure title ran `';'.join("neu")` -> "n;e;u"; the title
#     here just says neu.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# shellcheck source=../../config/env.sh
source "${REPO_DIR}/config/env.sh"

: "${OUT:?OUT not set}"
: "${FS_RUN:?FS_RUN not set: the dfe stage reads the sfs output}"
: "${DEMO_RUN:?DEMO_RUN not set}"
: "${DEMO_MODEL:?DEMO_MODEL not set}"
: "${CACHE_RUN:?CACHE_RUN not set: it keys the dfe run directory}"
: "${DFE_RUN:?DFE_RUN not set}"
: "${CTRL_TAG:?CTRL_TAG not set}"
: "${ROULETTE_SCALE:?ROULETTE_SCALE not set}"
: "${NUM_CONVERGENCE:?NUM_CONVERGENCE not set}"
: "${RESID_RARE_MAX_COUNT:?RESID_RARE_MAX_COUNT not set}"
: "${RESID_MAX_RARE_WARN:?RESID_MAX_RARE_WARN not set}"
: "${RESID_MAX_RARE_FAIL:?RESID_MAX_RARE_FAIL not set}"
: "${RESID_BULK_BIAS_WARN:?RESID_BULK_BIAS_WARN not set}"
: "${RESID_REST_OVER_RARE_MIN:?RESID_REST_OVER_RARE_MIN not set}"
: "${G2_PER_SNP_WARN:?G2_PER_SNP_WARN not set}"
: "${G2_PER_SNP_FAIL:?G2_PER_SNP_FAIL not set}"

if [[ $# -eq 0 ]]; then
  echo "usage: run_fitneu.sh <bed> [bed...]" >&2
  exit 2
fi

# dadi_cal_ll_neu_0723.py:186-197 accepts these; pipeline.conf's DEMO_MODEL
# spells bottlegrowth without the dadi suffix, so the two names meet here.
case "$DEMO_MODEL" in
  three_epoch|two_epoch|growth) demo_func="$DEMO_MODEL" ;;
  bottlegrowth)                 demo_func="bottlegrowth_1d" ;;
  *) echo "run_fitneu: DEMO_MODEL='${DEMO_MODEL}' is not supported" >&2; exit 2 ;;
esac

# --------------------------------------------------------------- preflight ----
# Every input for every bed is located before any work starts.
beds=(); test_fs=(); ctrl_fs=(); mut_t=(); mut_c=(); popts=(); dirs=()
for bed in "$@"; do
  base="$(basename "$bed" .bed)"
  # The genome-wide spectrum names, as run_concat.sh builds them: prep's
  # "s.nc_exHET_noq.<base>" with the "s.nc_" and then "nc_" prefixes dropped.
  shortbase="s.nc_exHET_noq.${base}"
  shortbase="${shortbase/s.nc_/}"
  shortbase="${shortbase/nc_/}"
  fsbase2="${shortbase}_${CTRL_TAG}.chrs_1t22"
  sfs_dir="${OUT}/sfs/${base}/${FS_RUN}"
  demo_dir="${OUT}/demo/${base}/${FS_RUN}-${DEMO_MODEL}-${DEMO_RUN}"
  best="${demo_dir}/best_inference.txt"

  t_fs="${sfs_dir}/has_qui.${fsbase2}.fs"
  c_fs="${sfs_dir}/quiecent_ctrl.${fsbase2}.fs"
  t_mut="${sfs_dir}/fsmut.has_qui.${fsbase2}.tab"
  c_mut="${sfs_dir}/fsmut.quiecent_ctrl.${fsbase2}.tab"
  for f in "$t_fs" "$c_fs" "$t_mut" "$c_mut"; do
    [[ -s "$f" ]] || {
      echo "run_fitneu: missing sfs output for '${base}': ${f}" >&2
      echo "      Run scripts/sfs/run_concat.sh first, or check FS_RUN." >&2
      exit 1; }
  done
  [[ -s "$best" ]] || {
    echo "run_fitneu: no demographic best fit for '${base}': ${best}" >&2
    echo "      Run scripts/demo/run_bestdemo.sh first, or check DEMO_RUN." >&2
    exit 1; }

  # The winning fit, recovered through max_ll_key exactly as
  # piperef_stddfe_5_fitneu.sh:56-57 does.
  key="$(awk -F: '/^max_ll_key:/ {print $2}' "$best")"
  popt="$(ls -1 "${demo_dir}"/*."${key}".InferDM.bestfits 2>/dev/null | head -1 || true)"
  [[ -n "$popt" ]] || {
    echo "run_fitneu: no .bestfits for max_ll_key=${key} in ${demo_dir}" >&2
    exit 1; }

  d="${OUT}/dfe/${base}/${CACHE_RUN}-${CACHE_MARKER}-neu-${DFE_RUN}"
  if [[ -d "$d" && -n "$(ls -A "$d" 2>/dev/null)" ]]; then
    if [[ "${FORCE:-0}" == "1" ]]; then
      echo "run_fitneu: FORCE=1 - overwriting the existing run in ${d}" >&2
    else
      echo "run_fitneu: run directory exists and is not empty: ${d}" >&2
      echo "      Results are keyed by DFE_RUN, so an existing run is never" >&2
      echo "      overwritten implicitly. Pick a new DFE_RUN, or re-run with" >&2
      echo "      FORCE=1 to replace this one." >&2
      exit 1
    fi
  fi

  beds+=("$base"); test_fs+=("$t_fs"); ctrl_fs+=("$c_fs")
  mut_t+=("$t_mut"); mut_c+=("$c_mut"); popts+=("$popt"); dirs+=("$d")
done

activate_conda

# ------------------------------------------------------------------ the null ----
failed=0
for i in "${!beds[@]}"; do
  base="${beds[$i]}"
  d="${dirs[$i]}"
  mkdir -p "$d"

  { echo "DFE_RUN=${DFE_RUN}"
    echo "CACHE_RUN=${CACHE_RUN}"
    echo "pdf=neu"
    echo "FS_RUN=${FS_RUN}"
    echo "DEMO_RUN=${DEMO_RUN}"
    echo "DEMO_MODEL=${DEMO_MODEL}"
    echo "dadi_model=${demo_func}"
    echo "test_fs=${test_fs[$i]}"
    echo "ctrl_fs=${ctrl_fs[$i]}"
    echo "fsmut_test=${mut_t[$i]}"
    echo "fsmut_ctrl=${mut_c[$i]}"
    echo "demopopt=${popts[$i]}"
    echo "date=$(date -Iseconds 2>/dev/null || date)"
  } > "${d}/spec.txt"

  echo "run_fitneu: ${base}: neutral likelihood of $(basename "${test_fs[$i]}")"
  if python3 - "$d" "${test_fs[$i]}" "${ctrl_fs[$i]}" "${popts[$i]}" \
               "$demo_func" "${mut_t[$i]}" "${mut_c[$i]}" "${PNG:-0}" \
               "$ROULETTE_SCALE" <<'PY'
"""The neutral log-likelihood, a faithful port of dadi_cal_ll_neu_0723.py.

Differences are listed in run_fitneu.sh's header. The expected spectrum is
handed to dadi.Inference.ll unfolded, as the original had it; ll_per_bin folds
it against the folded data itself, and only the copy written out for the
residual diagnostic is folded here.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")                     # no display on a compute node
import matplotlib.pyplot as plt           # noqa: E402
import dadi                               # noqa: E402
import dadi.Demographics1D                # noqa: E402
import dadi.Numerics                      # noqa: E402

dirout, testfs, ctrlfs, demopopt, model, testfsmut, ctrlfsmut, png = sys.argv[1:9]
ROULETTE_SCALE = float(sys.argv[9])
want_png = png not in ("", "0", "no", "false")

def die(msg):
    sys.stderr.write("run_fitneu: %s\n" % msg)
    sys.exit(1)


def get_opts_and_theta(filename, gen_cache=False):
    """dadi-cli's reader, copied from dadi_cal_ll_neu_0723.py:96-135.

    The FIRST data row of the file is taken, which after the "# Converged
    results" header is the best converged fit.
    """
    opts, param_names, is_converged = [], [], False
    with open(filename) as fid:
        for line in fid.readlines():
            if line.startswith("# Converged"):
                is_converged = True
            elif line.startswith("# Log(likelihood)"):
                param_names = line.rstrip().split("\t")
            elif line.startswith("#"):
                continue
            else:
                try:
                    opts = [float(v) for v in line.rstrip().split("\t")]
                    break
                except ValueError:
                    pass
    if not opts:
        die("no numeric fit row in %s" % filename)
    theta = opts[-1]
    if gen_cache and "misid" in param_names:
        opts = opts[1:-2]
    else:
        opts = opts[1:-1]
    if not is_converged:
        sys.stderr.write("run_fitneu: no converged results in %s\n" % filename)
    return opts, theta


def cache_pts_l_func(sample_sizes):
    """dadi-cli's CACHE grid sizes (dadi_cal_ll_neu_0723.py:138-151).

    Not dadi-cli's pts_l_func: the neutral model is evaluated on the same grid
    the selection cache was built on, so the null and the DFE fits share it.
    """
    n = max(sample_sizes)
    return (int(n * 2.2) + 2, int(n * 2.4) + 4, int(n * 2.6) + 6)


def getul(fsmutfile):
    """u*L from the genome-wide row of an fsmut table (:154-171).

    The row is the FIRST after the header and its fsfile column must contain
    "etal" (chr1 "et al", the concatenation marker run_concat.sh writes).
    """
    import csv
    with open(fsmutfile) as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    if "etal" not in row.get("fsfile", ""):
        die("%s has no genome-wide ('etal') row first: it holds only "
            "per-chromosome rows, so u*L cannot be formed" % fsmutfile)
    length = float(row["length_adj_inputbed"])
    u = ROULETTE_SCALE * float(row["avepersite_diploi_Roulette_mutrate"]) / 2
    return u * length


MODELS = {
    "three_epoch": dadi.Demographics1D.three_epoch,
    "two_epoch": dadi.Demographics1D.two_epoch,
    "growth": dadi.Demographics1D.growth,
    "bottlegrowth_1d": dadi.Demographics1D.bottlegrowth_1d,
    "snm": getattr(dadi.Demographics1D, "snm", None),
    "snm_1d": getattr(dadi.Demographics1D, "snm_1d", None),
}
func = MODELS.get(model)
if func is None:
    die("model '%s' is not supported" % model)

demoopts, theta_demog = get_opts_and_theta(demopopt, gen_cache=True)

ul_test, ul_ctrl = getul(testfsmut), getul(ctrlfsmut)
if ul_ctrl <= 0:
    die("the control u*L came out as %r; check %s" % (ul_ctrl, ctrlfsmut))
mutratio = ul_test / ul_ctrl

fs_ctrl = dadi.Spectrum.from_file(ctrlfs)
fs_test = dadi.Spectrum.from_file(testfs)

func_ex = dadi.Numerics.make_extrap_func(func)
pts_l = cache_pts_l_func(fs_ctrl.sample_sizes)
sfs_demog = func_ex(demoopts, fs_ctrl.sample_sizes, pts_l)

theta_test = theta_demog * mutratio
sfs_test_expected = sfs_demog * theta_test

ll = dadi.Inference.ll(sfs_test_expected, fs_test)

# ------------------------------------------------------------------- outputs --
# The log name is the historical one: the summary stage greps its `ll: ` line.
stem = os.path.basename(testfs)
if stem.endswith(".fs"):
    stem = stem[:-3]
logpath = os.path.join(dirout, stem + ".dadi_cal_ll_neu_0723.log")
with open(logpath, "w") as fh:
    fh.write("run_fitneu.sh %s\n" % " ".join(sys.argv[1:]))
    fh.write("ll: {}\n".format(ll))
    fh.write("mutratio: {}\n".format(mutratio))
    fh.write("theta_demog: {}\n".format(theta_demog))
    fh.write("theta_test: {}\n".format(theta_test))
    fh.write("demopopt: {}\n".format(demopopt))

# Folded, for the residual diagnostic only - see the header. The name matches
# the one run_bestdfe.sh gives a fitted DFE's best model, so check_dfe.py reads
# one file whatever the pdf.
sfs_test_expected.fold().to_file(os.path.join(dirout, "model_sfs_rank1.fs"))

fig = plt.figure()
fig.clear()
dadi.Plotting.plot_1d_comp_multinom(sfs_test_expected, fs_test)
plt.title("neu  ll:%s" % ll)
fig.savefig(os.path.join(dirout, "sfs_neu.pdf"))
if want_png:
    fig.savefig(os.path.join(dirout, "sfs_neu.png"), dpi=200)
plt.close(fig)

print("run_fitneu: ll=%s  mutratio=%s" % (ll, mutratio))
print("run_fitneu: wrote %s" % logpath)
PY
  then
    :
  else
    echo "run_fitneu: ${base}: the neutral likelihood failed" >&2
    failed=1
  fi
done

# ---------------------------------------------------------------- checkpoint --
# Diagnostic only: a bug in the check must never discard a finished stage, so its
# failure warns and the stage's own exit status is left alone.
if ! python3 "${REPO_DIR}/scripts/checks/check_dfe.py" \
       --out "${OUT}/dfe" --run "${CACHE_RUN}-${CACHE_MARKER}-neu-${DFE_RUN}" --pdf neu \
       --sfs-out "${OUT}/sfs" --fs-run "${FS_RUN}" \
       --rare-max-count "${RESID_RARE_MAX_COUNT}" \
       --max-rare-warn "${RESID_MAX_RARE_WARN}" \
       --max-rare-fail "${RESID_MAX_RARE_FAIL}" \
       --bulk-bias-warn "${RESID_BULK_BIAS_WARN}" \
       --rest-over-rare-min "${RESID_REST_OVER_RARE_MIN}" \
       --g2-warn "${G2_PER_SNP_WARN}" --g2-fail "${G2_PER_SNP_FAIL}" \
       --num-convergence "${NUM_CONVERGENCE}" "$@"; then
  echo "WARNING: check_dfe.py failed; the neutral likelihood above is" >&2
  echo "         unaffected. Inspect the .dadi_cal_ll_neu_0723.log by hand." >&2
fi

exit "$failed"
