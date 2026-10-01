#!/bin/bash
# submit_cache.sh <bed> [bed...]
#
# One array job over the BEDS: one cache is ~10 hours on 10 cores, so two beds in
# one task would not fit inside h_rt. Every task builds one bed's cache.
#
#   TERSE=1   print only the job id, for chaining
#   HOLD=<id> hold until that job finishes (e.g. the demo best-fit job)
#   FORCE=1   forwarded to run_cache.sh, which then rebuilds an existing cache
#
# The cache is the most expensive step in the pipeline, so the launcher checks
# what it can before queueing: every bed must already have a best_inference.txt,
# and a bed whose cache exists is not submitted at all unless FORCE=1.
#
# Cache cost is data-independent: gamma_pts (linear) x grid size (~quadratic),
# divided by cpus. Measured 20 min at 300 pts / 1540 / 10 cpus. We estimate from
# that and raise SGE_CACHE_RT automatically when a marker needs more, so a big
# marker cannot silently inherit a request sized for the standard one.
# on highp, -pe shared 10. SGE_CACHE_EXTRA carries the highp; empty asks for no
# extra resource.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# qsub is not on PATH in a non-interactive shell; env.sh knows where it lives.
source "${REPO_DIR}/config/env.sh"
ensure_scheduler

: "${OUT:?OUT not set}"
: "${FS_RUN:?FS_RUN not set}"
: "${DEMO_RUN:?DEMO_RUN not set}"
: "${DEMO_MODEL:?DEMO_MODEL not set}"
: "${CACHE_RUN:?CACHE_RUN not set}"
: "${CACHE_MARKER:?CACHE_MARKER not set}"
: "${SGE_CACHE_RT:?SGE_CACHE_RT not set}"
: "${SGE_CACHE_MEM:?SGE_CACHE_MEM not set}"
: "${SGE_CACHE_PE:?SGE_CACHE_PE not set}"

if [[ $# -eq 0 ]]; then
  echo "usage: submit_cache.sh <bed> [bed...]" >&2
  exit 2
fi

run_dir_name="${DEMO_RUN}-${CACHE_MARKER}-${CACHE_RUN}"

# The job's working directory is not the submitting one for relative lookups we
# can rely on, so bed arguments are resolved here, while the caller's $PWD still
# means something. A bed with nothing to do is dropped from the array.
beds=()
for bed in "$@"; do
  case "$bed" in
    /*) abs="$bed" ;;
    *)  abs="${PWD}/${bed}" ;;
  esac
  base="$(basename "$abs" .bed)"
  best="${OUT}/demo/${base}/${FS_RUN}-${DEMO_MODEL}-${DEMO_RUN}/best_inference.txt"
  if [[ ! -s "$best" ]]; then
    echo "submit_cache: '${base}' has no best_inference.txt: ${best}" >&2
    echo "      The gamma bounds are 2*s*Na and Na comes from there. Run" >&2
    echo "      scripts/demo/run_bestfit.sh first, or check FS_RUN /" >&2
    echo "      DEMO_MODEL / DEMO_RUN." >&2
    exit 1
  fi
  out_bpkl="${OUT}/cache/${base}/${run_dir_name}/cache_${CACHE_MARKER}.spectra.bpkl"
  if [[ -s "$out_bpkl" && "${FORCE:-0}" != "1" ]]; then
    echo "submit_cache: '${base}' already has a cache, not submitting it:" >&2
    echo "      ${out_bpkl}  (FORCE=1 rebuilds it)" >&2
    continue
  fi
  beds+=("$abs")
done

if [[ "${#beds[@]}" -eq 0 ]]; then
  echo "submit_cache: nothing to do; every bed already has a cache for" \
       "${run_dir_name}" >&2
  exit 0
fi

# Logs live under $OUT, not in the package: the package is synced to the cluster
# and must stay clean, and the logs belong with the results they describe.
log_dir="${OUT}/cache/log"
mkdir -p "$log_dir"

# The config the launcher validated is forwarded explicitly, so a task cannot
# silently pick up a different run key - or a different bed list, which would
# shift every task's index - from the worker node's environment. -V carries the
# rest of the exported environment (the site profile).
job_env="REPO_DIR=${REPO_DIR},OUT=${OUT},FS_RUN=${FS_RUN},DEMO_RUN=${DEMO_RUN}"
job_env="${job_env},DEMO_MODEL=${DEMO_MODEL},CACHE_RUN=${CACHE_RUN}"
job_env="${job_env},CACHE_MARKER=${CACHE_MARKER}"
job_env="${job_env},CACHE_BEDS=${beds[*]}"
# The chain sets these; a stage submitted by hand leaves them unset and
# no gate runs. scripts/sge/gate_upstream.sh
[[ -n "${GATE_FILES:-}" ]] && job_env="${job_env},GATE_FILES=${GATE_FILES}"
[[ -n "${GATE_STAGE:-}" ]] && job_env="${job_env},GATE_STAGE=${GATE_STAGE}"
[[ -n "${FORCE:-}" ]] && job_env="${job_env},FORCE=${FORCE}"

# Estimate, then request 3x it (a killed cache leaves nothing to resume from,
# but 6x an extrapolated estimate exceeds the queue limits for big markers).
max_grid="$(printf '%s\n' ${GRIDS} | sort -n | tail -1)"
est_s="$(awk -v g="${GAMMA_PTS}" -v m="${max_grid}" -v c="${SGE_CACHE_PE}" \
    'BEGIN {printf "%d", 40 * g * (m/1540)^2 / c}')"
req_s=$(( est_s * 3 ))
# Never go below the configured request; raise it if the estimate says so.
conf_s="$(awk -F: '{print ($1*3600)+($2*60)+$3}' <<< "${SGE_CACHE_RT}")"
if (( req_s > conf_s )); then
  SGE_CACHE_RT="$(printf '%d:%02d:00' $(( req_s / 3600 )) $(( (req_s % 3600) / 60 )))"
  echo "submit_cache: ${GAMMA_PTS} gamma points on grid ${max_grid} with ${SGE_CACHE_PE} cpus" >&2
  echo "              estimates ~$(( est_s / 60 )) min; requesting h_rt=${SGE_CACHE_RT}" >&2
else
  echo "submit_cache: estimated ~$(( est_s / 60 )) min; requesting h_rt=${SGE_CACHE_RT}" >&2
fi

qsub ${TERSE:+-terse} \
  -N "dfe_cache_${CACHE_RUN}" \
  -t "1-${#beds[@]}" \
  -cwd -j y -V \
  -o "${log_dir}/cache.${run_dir_name}.\$JOB_ID.\$TASK_ID.log" \
  -l h_rt="${SGE_CACHE_RT}",h_data="${SGE_CACHE_MEM}"${SGE_CACHE_EXTRA:+,${SGE_CACHE_EXTRA}} \
  -pe shared "${SGE_CACHE_PE}" \
  ${HOLD:+-hold_jid "${HOLD}"} \
  -v "${job_env}" \
  "${here}/worker_cache.sh"
