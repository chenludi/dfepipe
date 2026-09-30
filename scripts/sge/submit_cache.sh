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
# The queue request is submit_dadi-cli_cache_simp.sh's: h_rt=10:10:00,h_data=15G
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
[[ -n "${FORCE:-}" ]] && job_env="${job_env},FORCE=${FORCE}"

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
