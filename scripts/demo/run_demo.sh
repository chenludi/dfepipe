#!/bin/bash
# demo, the fits: infer the demographic history from the CONTROL spectrum, once
# per line of one initial-parameter file.
#
#   bash scripts/demo/run_demo.sh <p0file> <bed> [bed...]
#
# <p0file> is a file in DEMO_P0_DIR, by name or by full path; <bed> is the
# ORIGINAL bed, the same argument prep and the sfs stage took. submit_demo.sh
# submits one array task per p0 file, so this script is the body of one task.
#
# Ported from piperef_stddfe_3_inferdemo.sh (which p0 files, which spectrum) and
# dadi-cli_demog_fsfile.sh (the InferDM call). Each p0 file holds one starting
# point per line:
#
#   <p0marker> <p0 x nparam> <lbounds x nparam> <ubounds x nparam>
#
# and each line becomes its own InferDM run, wrapped in `timeout`, exactly as
# submit_dadi-cli_demog_fsfile.sh:60 did. The output prefix keeps the historical
# name, <DEMO_MODEL>.<fs basename>.<p0marker>, because the p0marker is the
# `max_ll_key` that best_inference.txt reports and run_bestdemo.sh recovers the
# winning file by globbing for it.
#
# The spectrum fitted is the quiescent control one, not the test one: the
# demography is estimated from putatively neutral sites and then held fixed while
# the DFE is fitted to the test sites.
#
# Differences from the original, all deliberate:
#   * `--output-prefix` and `--cpus` are spelled in full. The original wrote
#     `--output` and `--cpu`, which worked only as argparse abbreviations.
#   * a fit whose .bestfits already exists is skipped, so an array task that was
#     killed or timed out can simply be re-run; FORCE=1 redoes them.
#   * dadi-cli_demog_fsfile.sh:32 computes its output prefix before reading
#     $model_demog and $p0marker, so every file it writes is named ".<fs>."  -
#     a hidden file with no marker in it. The published results predate that
#     edit and carry the full name, which is what is reproduced here.
#   * a failed or timed-out fit is reported and the rest of the file continues
#     (as the original's while-loop did), but a task in which NOTHING succeeded
#     exits non-zero instead of pretending to have worked.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# shellcheck source=../../config/env.sh
source "${REPO_DIR}/config/env.sh"

: "${OUT:?OUT not set}"
: "${FS_RUN:?FS_RUN not set: the demo stage reads the sfs output}"
: "${DEMO_RUN:?DEMO_RUN not set}"
: "${DEMO_MODEL:?DEMO_MODEL not set}"
: "${DEMO_P0_DIR:?DEMO_P0_DIR not set: the directory of initial-parameter files}"
: "${CTRL_TAG:?CTRL_TAG not set}"
: "${MAXEVAL:?MAXEVAL not set}"
: "${NUM_CONVERGENCE:?NUM_CONVERGENCE not set}"
: "${SGE_DEMO_PE:?SGE_DEMO_PE not set: the core count InferDM is given}"

# Demography's own optcount and per-fit clock: 40 and 30m, from
# dadi-cli_demog_fsfile.sh:97 and submit_dadi-cli_demog_fsfile.sh:60.
# OPTCOUNT (100) is the DFE stage's, not this one's.
: "${DEMO_OPTCOUNT:?DEMO_OPTCOUNT not set}"
: "${DEMO_TIMEOUT:?DEMO_TIMEOUT not set}"

if [[ $# -lt 2 ]]; then
  echo "usage: run_demo.sh <p0file> <bed> [bed...]" >&2
  exit 2
fi

p0arg="$1"; shift
beds=("$@")

# dadi-cli's model name and its parameter count. pipeline.conf accepts
# `bottlegrowth`; dadi-cli calls the 1-population version `bottlegrowth_1d`.
# pipeline.conf accepts only three_epoch and two_epoch - the two with a p0 grid -
# so the branches below for growth and bottlegrowth are unreachable today. They
# are kept because adding a p0 grid is all it takes to need them again.
# (dadi-cli_demog_fsfile.sh:74-75), so the two names meet here.
case "$DEMO_MODEL" in
  three_epoch)  dadi_model="three_epoch";     nparam=4 ;;   # nuB nuF TB TF
  two_epoch)    dadi_model="two_epoch";       nparam=2 ;;   # nu T
  growth)       dadi_model="growth";          nparam=2 ;;   # nu T
  bottlegrowth) dadi_model="bottlegrowth_1d"; nparam=3 ;;   # nuB nuF T
  *) echo "run_demo: DEMO_MODEL='${DEMO_MODEL}' is not supported" >&2; exit 2 ;;
esac

case "$p0arg" in
  */*) p0file="$p0arg" ;;
  *)   p0file="${DEMO_P0_DIR}/${p0arg}" ;;
esac
[[ -s "$p0file" ]] || {
  echo "run_demo: initial-parameter file missing or empty: ${p0file}" >&2
  echo "      DEMO_P0_DIR is ${DEMO_P0_DIR}" >&2
  exit 1; }

# --------------------------------------------------------------- preflight ----
# Every spectrum is located before any fitting starts: a task that runs for an
# hour and then finds a typo in the second bed has wasted the hour.
fs_files=(); demo_dirs=(); fs_bases=()
for bed in "${beds[@]}"; do
  base="$(basename "$bed" .bed)"
  # The genome-wide spectrum's name, as run_concat.sh builds it: prep's
  # "s.nc_exHET_noq.<base>" with the "s.nc_" and then "nc_" prefixes dropped.
  shortbase="s.nc_exHET_noq.${base}"
  shortbase="${shortbase/s.nc_/}"
  shortbase="${shortbase/nc_/}"
  fsbase="quiecent_ctrl.${shortbase}_${CTRL_TAG}.chrs_1t22"
  fs="${OUT}/sfs/${base}/${FS_RUN}/${fsbase}.fs"
  [[ -s "$fs" ]] || {
    echo "run_demo: no control spectrum for '${base}': ${fs}" >&2
    echo "      Run scripts/sfs/run_concat.sh first, or check FS_RUN." >&2
    exit 1; }
  fs_files+=("$fs")
  fs_bases+=("$fsbase")
  demo_dirs+=("${OUT}/demo/${base}/${FS_RUN}-${DEMO_MODEL}-${DEMO_RUN}")
done

activate_conda

# --------------------------------------------------------------------- fits ----
n_done=0; n_skipped=0; n_failed=0

for i in "${!beds[@]}"; do
  fs="${fs_files[$i]}"
  fsbase="${fs_bases[$i]}"
  d="${demo_dirs[$i]}"
  mkdir -p "$d"

  { echo "FS_RUN=${FS_RUN}"
    echo "DEMO_RUN=${DEMO_RUN}"
    echo "DEMO_MODEL=${DEMO_MODEL}"
    echo "dadi_cli_model=${dadi_model}"
    echo "DEMO_P0_DIR=${DEMO_P0_DIR}"
    echo "control_fs=${fs}"
    echo "MAXEVAL=${MAXEVAL}"
    echo "DEMO_OPTCOUNT=${DEMO_OPTCOUNT}"
    echo "NUM_CONVERGENCE=${NUM_CONVERGENCE}"
    echo "DEMO_TIMEOUT=${DEMO_TIMEOUT}"
    echo "cpus=${SGE_DEMO_PE}"
    echo "date=$(date -Iseconds 2>/dev/null || date)"
  } > "${d}/spec.txt"

  lineno=0
  while IFS= read -r line || [[ -n "$line" ]]; do
    lineno=$((lineno + 1))
    [[ -z "${line// /}" ]] && continue
    case "$line" in \#*) continue ;; esac

    read -r -a f <<< "$line"
    want=$((1 + 3 * nparam))
    if [[ "${#f[@]}" -ne "$want" ]]; then
      echo "run_demo: $(basename "$p0file") line ${lineno} has ${#f[@]} fields," \
           "expected ${want} for ${DEMO_MODEL} (marker + ${nparam} p0 +" \
           "${nparam} lbounds + ${nparam} ubounds)" >&2
      n_failed=$((n_failed + 1))
      continue
    fi

    marker="${f[0]}"
    p0=("${f[@]:1:nparam}")
    lb=("${f[@]:1+nparam:nparam}")
    ub=("${f[@]:1+2*nparam:nparam}")
    # -1 means "this parameter is not fixed", one per parameter.
    consts=(); for ((k = 0; k < nparam; k++)); do consts+=("-1"); done

    prefix="${d}/${dadi_model}.${fsbase}.${marker}"
    if [[ -s "${prefix}.InferDM.bestfits" && "${FORCE:-0}" != "1" ]]; then
      n_skipped=$((n_skipped + 1))
      continue
    fi

    cmd=(dadi-cli InferDM
         --fs "$fs"
         --nomisid
         --model "$dadi_model"
         --lbounds "${lb[@]}"
         --ubounds "${ub[@]}"
         --p0 "${p0[@]}"
         --constants "${consts[@]}"
         --output-prefix "$prefix"
         --maxeval "$MAXEVAL"
         --optimizations "$DEMO_OPTCOUNT"
         --cpus "$SGE_DEMO_PE"
         --check-convergence "$NUM_CONVERGENCE")

    echo "run_demo: ${DEMO_MODEL} p0marker=${marker} (line ${lineno})"
    printf '%s\n' "${cmd[*]}" > "${prefix}.log"
    if timeout "$DEMO_TIMEOUT" "${cmd[@]}" >> "${prefix}.log" 2>&1; then
      echo "Finished running dadi-cli InferDM" >> "${prefix}.log"
      n_done=$((n_done + 1))
    else
      rc=$?
      echo "dadi-cli InferDM exited ${rc}" >> "${prefix}.log"
      if [[ "$rc" -eq 124 ]]; then
        echo "run_demo: p0marker=${marker} hit the ${DEMO_TIMEOUT} timeout" >&2
      else
        echo "run_demo: p0marker=${marker} failed (exit ${rc}); see ${prefix}.log" >&2
      fi
      n_failed=$((n_failed + 1))
    fi
  done < "$p0file"
done

echo "run_demo: $(basename "$p0file"): ${n_done} fitted, ${n_skipped} already" \
     "present, ${n_failed} failed"

if [[ "$n_done" -eq 0 && "$n_skipped" -eq 0 ]]; then
  echo "run_demo: not one fit succeeded for $(basename "$p0file")" >&2
  exit 1
fi
exit 0
