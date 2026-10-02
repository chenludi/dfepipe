#!/bin/bash
# submit_dfe.sh <pdf> <bed> [bed...]
#
# One submission per pdf, as piperef_stddfe_7_fitgamneug.sh was invoked (once per
# entry of PDFS; record_251006_repeat_phyloP243_0-5.sh:47-49 shows the loop).
#
#   pdf = neu   a single job running scripts/dfe/run_fitneu.sh for every bed.
#   otherwise   the command packs are generated HERE, on the submitting host,
#               and then one array task runs each (bed, pack) pair. Generating
#               them up front is what makes the array size knowable, and it
#               fails fast on a missing cache or demography instead of leaving
#               every task to discover it.
#
#   TERSE=1   print only the job id, for chaining
#   HOLD=<id> hold until that job finishes (e.g. the cache job)
#   FORCE=1   forwarded: replaces an existing run directory, and redoes fits
#   RESUME=1  forwarded: keeps finished fits, runs only the missing starting
#             points - what to use after an h_rt kill
#             whose .InferDFE.bestfits already exist
#
# The best-fit pick is not submitted here: it can only run once every pack has
# finished. Run it after this job, holding on its id:
#
#   jid=$(TERSE=1 bash scripts/sge/submit_dfe.sh gamma mybed.bed)
#   HOLD="${jid%%.*}" ... scripts/dfe/run_bestdfe.sh gamma mybed.bed
#
# The original submitted one job per pack with no dependency tracking at all
# (piperef_stddfe_7_fitgamneug.sh:96-103), and each job ran up to PACKNUMBER
# commands of up to DFE_TIMEOUT each - 50 x 60m against an h_rt of 23:50:00, so
# a full pack could not finish inside its own job. See the report.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# qsub is not on PATH in a non-interactive shell; env.sh knows where it lives.
# shellcheck source=../../config/env.sh
source "${REPO_DIR}/config/env.sh"
ensure_scheduler

: "${OUT:?OUT not set}"
: "${FS_RUN:?FS_RUN not set}"
: "${DEMO_RUN:?DEMO_RUN not set}"
: "${DEMO_MODEL:?DEMO_MODEL not set}"
: "${CACHE_RUN:?CACHE_RUN not set}"
: "${CACHE_MARKER:?CACHE_MARKER not set}"   # part of the run key: two markers must not share a directory
: "${DFE_RUN:?DFE_RUN not set}"
: "${SGE_DFE_RT:?SGE_DFE_RT not set}"

# PACKNUMBER fits of up to DFE_TIMEOUT each can outlast h_rt - 50 x 60m against
# 23:50 in the shipped config. The worker's deadline guard means the overflow is
# recorded as not_run and finished by RESUME=1 rather than killed part-way, but
# say so up front so the extra pass is expected, not a surprise.
# How many fits go in one task. DFE_TARGET_TASK_HOURS, if set, decides it from
# DFE_SECONDS_PER_FIT; otherwise PACKNUMBER stands as configured. The estimate
# is deliberately crude: getting it wrong costs one extra RESUME pass, because
# run_dfe.sh stops before h_rt and logs the remainder as not_run.
if [[ -n "${DFE_TARGET_TASK_HOURS:-}" ]]; then
  _want=$(( DFE_TARGET_TASK_HOURS * 3600 / DFE_SECONDS_PER_FIT ))
  (( _want < 1 )) && _want=1
  echo "submit_dfe: DFE_TARGET_TASK_HOURS=${DFE_TARGET_TASK_HOURS} at" \
       "${DFE_SECONDS_PER_FIT}s per fit -> PACKNUMBER=${_want}" \
       "(was ${PACKNUMBER})" >&2
  PACKNUMBER=$_want
  export PACKNUMBER
fi

# h_rt is the hard ceiling whatever the preference says.
_rt_h=${SGE_DFE_RT%%:*}
_est_h=$(( PACKNUMBER * DFE_SECONDS_PER_FIT / 3600 ))
case "$DFE_TIMEOUT" in
  *h) _to_min=$(( ${DFE_TIMEOUT%h} * 60 )) ;;
  *m) _to_min=${DFE_TIMEOUT%m} ;;
  *)  _to_min=60 ;;
esac
echo "submit_dfe: ${PACKNUMBER} fits per task, estimated ${_est_h}h" \
     "(worst case $(( PACKNUMBER * _to_min / 60 ))h), h_rt ${SGE_DFE_RT}" >&2
if (( _est_h * 60 > 10#${_rt_h} * 60 )); then
  echo "submit_dfe: the ESTIMATE already exceeds h_rt. Each task will stop" \
       "early and log the rest as not_run; RESUME=1 finishes them. For one" \
       "pass use DFE_TARGET_TASK_HOURS=$(( 10#${_rt_h} * 4 / 5 ))." >&2
fi
: "${SGE_DFE_MEM:?SGE_DFE_MEM not set}"
: "${SGE_DFE_PE:?SGE_DFE_PE not set}"

if [[ $# -lt 2 ]]; then
  echo "usage: submit_dfe.sh <pdf> <bed> [bed...]" >&2
  echo "       pdf is one of: ${PDFS}" >&2
  exit 2
fi
pdf="$1"; shift

# The job's working directory is not the submitting one, so bed arguments are
# resolved while the caller's $PWD still means something.
beds=()
for bed in "$@"; do
  case "$bed" in
    /*) beds+=("$bed") ;;
    *)  beds+=("${PWD}/${bed}") ;;
  esac
done

# Logs live under $OUT, not in the package: the package is synced to the cluster
# and must stay clean, and the logs belong with the results they describe.
log_dir="${OUT}/dfe/log"
mkdir -p "$log_dir"

# The config the launcher validated is forwarded explicitly, so a task cannot
# silently pick up a different run key from the worker node's environment. -V
# carries the rest of the exported environment (the site profile).
job_env="REPO_DIR=${REPO_DIR},OUT=${OUT},FS_RUN=${FS_RUN},DEMO_RUN=${DEMO_RUN}"
job_env="${job_env},DEMO_MODEL=${DEMO_MODEL},CACHE_RUN=${CACHE_RUN}"
job_env="${job_env},DFE_RUN=${DFE_RUN},DFE_PDF=${pdf}"
# The worker turns h_rt into a deadline so run_dfe.sh can stop before a fit
# it cannot finish, instead of having one killed part-way.
job_env="${job_env},SGE_DFE_RT=${SGE_DFE_RT}"
# The chain sets these; a stage submitted by hand leaves them unset and
# no gate runs. scripts/sge/gate_upstream.sh
[[ -n "${GATE_FILES:-}" ]] && job_env="${job_env},GATE_FILES=${GATE_FILES}"
[[ -n "${GATE_STAGE:-}" ]] && job_env="${job_env},GATE_STAGE=${GATE_STAGE}"
[[ -n "${FORCE:-}" ]] && job_env="${job_env},FORCE=${FORCE}"
[[ -n "${SKIP_GATE:-}" ]] && job_env="${job_env},SKIP_GATE=${SKIP_GATE}"
[[ -n "${RESUME:-}" ]] && job_env="${job_env},RESUME=${RESUME}"
[[ -n "${DFE_CACHE_FILE:-}" ]] && job_env="${job_env},DFE_CACHE_FILE=${DFE_CACHE_FILE}"
[[ -n "${DFE_TIMEOUT:-}" ]] && job_env="${job_env},DFE_TIMEOUT=${DFE_TIMEOUT}"

# ------------------------------------------------------------- the null, once --
if [[ "$pdf" == "neu" ]]; then
  qsub ${TERSE:+-terse} \
    -N "dfe_neu_${DFE_RUN}" \
    -cwd -j y -V \
    -o "${log_dir}/neu.${DFE_RUN}.\$JOB_ID.log" \
    -l h_rt="${SGE_DFE_RT}",h_data="${SGE_DFE_MEM}" \
    -pe shared "${SGE_DFE_PE}" \
    ${HOLD:+-hold_jid "${HOLD}"} \
    -v "${job_env}" \
    "${here}/worker_dfe.sh" "${beds[@]}"
  exit 0
fi

# ----------------------------------------------------------------- the packs --
# Generated here, in this process, with the config this script just validated.
echo "submit_dfe: generating the ${pdf} command packs" >&2
DFE_GENERATE_ONLY=1 bash "${REPO_DIR}/scripts/dfe/run_dfe.sh" "$pdf" "${beds[@]}" >&2

# One task per (bed, pack), in one fixed order both sides agree on. LC_ALL=C so
# the launcher and the worker cannot disagree about collation.
task_beds=(); task_packs=(); n_points=0
for bed in "${beds[@]}"; do
  base="$(basename "$bed" .bed)"
  d="${OUT}/dfe/${base}/${CACHE_RUN}-${CACHE_MARKER}-${pdf}-${DFE_RUN}"
  found=0
  while IFS= read -r p; do
    task_beds+=("$bed"); task_packs+=("$(basename "$p")"); found=1
    # `|| true`: grep -c exits 1 when it counts zero lines, and under `set -e`
    # that aborted this script with no message at all - the "no command packs"
    # branch below was never reached. A pack with no fit lines is possible when a
    # previous generation was killed between creating the file and appending.
    n_points=$(( n_points + $(grep -cv '^[[:space:]]*\(#\|$\)' "$p" || true) ))
  done < <(LC_ALL=C ls -1 "${d}"/dadi-cli_dfe.*.sh 2>/dev/null)
  if [[ "$found" -eq 0 ]]; then
    echo "submit_dfe: no command packs were generated for '${base}' in ${d}" >&2
    exit 1
  fi
done

if [[ -n "${DFE_MAX_TASKS:-}" ]] && (( ${#task_packs[@]} > DFE_MAX_TASKS )); then
  _need=$(( (n_points + DFE_MAX_TASKS - 1) / DFE_MAX_TASKS ))
  echo "submit_dfe: this would be ${#task_packs[@]} array tasks but" \
       "DFE_MAX_TASKS=${DFE_MAX_TASKS}. Nothing was submitted." >&2
  echo "      ${n_points} starting points over ${DFE_MAX_TASKS} tasks needs" \
       "PACKNUMBER=${_need} (about $(( _need * DFE_SECONDS_PER_FIT / 3600 ))h" \
       "per task). Re-run with FORCE=1 PACKNUMBER=${_need}." >&2
  exit 1
fi

job_env="${job_env},DFE_TASK_BEDS=${task_beds[*]},DFE_PACKS=${task_packs[*]}"

qsub ${TERSE:+-terse} \
  -N "dfe_${pdf}_${DFE_RUN}" \
  -t "1-${#task_packs[@]}" \
  -cwd -j y -V \
  -o "${log_dir}/${pdf}.${DFE_RUN}.\$JOB_ID.\$TASK_ID.log" \
  -l h_rt="${SGE_DFE_RT}",h_data="${SGE_DFE_MEM}" \
  -pe shared "${SGE_DFE_PE}" \
  ${HOLD:+-hold_jid "${HOLD}"} \
  -v "${job_env}" \
  "${here}/worker_dfe.sh"
