#!/bin/bash
# submit_demo.sh <bed> [bed...]
#
# One array job over the initial-parameter files in DEMO_P0_DIR: every starting
# point is an independent optimisation, so they run in parallel. Every task
# handles every bed for its own p0 file, as piperef_stddfe_3_inferdemo.sh's inner
# loop did (one qsub per p0 file, all beds inside).
#
#   TERSE=1   print only the job id, for chaining
#   HOLD=<id> hold until that job finishes (e.g. the sfs concatenation)
#   FORCE=1   forwarded to run_demo.sh, which then redoes fits that already exist
#
# The best-fit pick is not submitted here: it can only run once every task has
# finished. Submit it after this job, holding on its id:
#
#   jid=$(TERSE=1 bash scripts/sge/submit_demo.sh mybed.bed)
#   HOLD="${jid%%.*}" ... scripts/demo/run_bestdemo.sh mybed.bed
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
: "${DEMO_P0_DIR:?DEMO_P0_DIR not set}"
: "${SGE_DEMO_RT:?SGE_DEMO_RT not set}"
: "${SGE_DEMO_MEM:?SGE_DEMO_MEM not set}"
: "${SGE_DEMO_PE:?SGE_DEMO_PE not set}"

if [[ $# -eq 0 ]]; then
  echo "usage: submit_demo.sh <bed> [bed...]" >&2
  exit 2
fi

# piperef_stddfe_3_inferdemo.sh:36-40 refuses more than 500 spectra in one go and
# tells the user to write a job file instead. One bed is one spectrum here.
if [[ $# -gt 500 ]]; then
  echo "submit_demo: $# beds is more than 500 spectra for one submission." >&2
  echo "      Split the list and submit it in batches." >&2
  exit 2
fi

# The p0 file list, in one fixed order the worker reindexes from. LC_ALL=C so
# the launcher and the worker cannot disagree about collation.
[[ -d "$DEMO_P0_DIR" ]] || {
  echo "submit_demo: DEMO_P0_DIR does not exist: ${DEMO_P0_DIR}" >&2
  echo "      It holds one file per batch of starting points; see" >&2
  echo "      config/pipeline.conf's DEMO_MODEL / DEMO_P0_SET." >&2
  exit 1; }
p0files=()
while IFS= read -r f; do
  [[ -s "${DEMO_P0_DIR}/${f}" ]] && p0files+=("$f")
done < <(LC_ALL=C ls -1 "$DEMO_P0_DIR")
if [[ "${#p0files[@]}" -eq 0 ]]; then
  echo "submit_demo: no initial-parameter files in ${DEMO_P0_DIR}" >&2
  exit 1
fi

# The job's working directory is not the submitting one for relative lookups we
# can rely on, so bed arguments are resolved here, while the caller's $PWD still
# means something.
beds=()
for bed in "$@"; do
  case "$bed" in
    /*) beds+=("$bed") ;;
    *)  beds+=("${PWD}/${bed}") ;;
  esac
done

# Logs live under $OUT, not in the package: the package is synced to the cluster
# and must stay clean, and the logs belong with the results they describe.
log_dir="${OUT}/demo/log"
mkdir -p "$log_dir"

# The config the launcher validated is forwarded explicitly, so a task cannot
# silently pick up a different run key - or a different p0 file list, which would
# shift every task's index - from the worker node's environment. -V carries the
# rest of the exported environment (the site profile).
job_env="REPO_DIR=${REPO_DIR},OUT=${OUT},FS_RUN=${FS_RUN},DEMO_RUN=${DEMO_RUN}"
job_env="${job_env},DEMO_MODEL=${DEMO_MODEL},DEMO_P0_DIR=${DEMO_P0_DIR}"
job_env="${job_env},DEMO_P0_FILES=${p0files[*]}"
# The chain sets these; a stage submitted by hand leaves them unset and
# no gate runs. scripts/sge/gate_upstream.sh
[[ -n "${GATE_FILES:-}" ]] && job_env="${job_env},GATE_FILES=${GATE_FILES}"
[[ -n "${GATE_STAGE:-}" ]] && job_env="${job_env},GATE_STAGE=${GATE_STAGE}"
[[ -n "${FORCE:-}" ]] && job_env="${job_env},FORCE=${FORCE}"

[[ -n "${SKIP_GATE:-}" ]] && job_env="${job_env},SKIP_GATE=${SKIP_GATE}"
qsub ${TERSE:+-terse} \
  -N "dfe_demo_${DEMO_RUN}" \
  -t "1-${#p0files[@]}" \
  -cwd -j y -V \
  -o "${log_dir}/demo.${DEMO_RUN}.\$JOB_ID.\$TASK_ID.log" \
  -l h_rt="${SGE_DEMO_RT}",h_data="${SGE_DEMO_MEM}" \
  -pe shared "${SGE_DEMO_PE}" \
  ${HOLD:+-hold_jid "${HOLD}"} \
  -v "${job_env}" \
  "${here}/worker_demo.sh" "${beds[@]}"
