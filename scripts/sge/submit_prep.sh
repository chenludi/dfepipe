#!/bin/bash
# submit_prep.sh <bed> [bed...]
#
# One SGE job for the whole list: prep is a handful of bedtools calls, so the
# queue wait dominates and splitting per bed buys nothing.
#
#   TERSE=1   print only the job id, for chaining
#   HOLD=<id> hold until that job finishes
#   FORCE=1   forwarded to run_prep.sh, which then replaces an existing run
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# qsub is not on PATH in a non-interactive shell; env.sh knows where it lives.
source "${REPO_DIR}/config/env.sh"
ensure_scheduler

: "${OUT:?OUT not set}"
: "${PREP_RUN:?PREP_RUN not set}"
: "${SGE_PREP_RT:?SGE_PREP_RT not set}"
: "${SGE_PREP_MEM:?SGE_PREP_MEM not set}"

if [[ $# -eq 0 ]]; then
  echo "usage: submit_prep.sh <bed> [bed...]" >&2
  exit 2
fi

# The job's working directory is not the submitting one for relative lookups we
# can rely on, so bed arguments are resolved here, while the caller's $PWD still
# means something.
# Deliberately not an array over beds, unlike sfs and demo. run_prep.sh runs
# check_bed.py and check_prep.py itself (run_prep.sh:153,157) and both write one
# file for the whole stage, so parallel tasks would race on it. Prep is also the
# cheap stage - 64 s a bed, about 3 min for three - so splitting it would buy
# minutes and cost a check file. Move the checks to their own step first.
beds=()
for bed in "$@"; do
  case "$bed" in
    /*) beds+=("$bed") ;;
    *)  beds+=("${PWD}/${bed}") ;;
  esac
done

# Logs live under $OUT, not in the package: the package is synced to the cluster
# and must stay clean, and the logs belong with the results they describe.
log_dir="${OUT}/prep/log"
mkdir -p "$log_dir"

# The config the launcher validated is forwarded explicitly, so the job cannot
# silently pick up a different run key from the worker node's environment. -V
# carries the rest of the exported environment (the site profile).
job_env="REPO_DIR=${REPO_DIR},OUT=${OUT},PREP_RUN=${PREP_RUN}"
# The chain sets these; a stage submitted by hand leaves them unset and
# no gate runs. scripts/sge/gate_upstream.sh
[[ -n "${GATE_FILES:-}" ]] && job_env="${job_env},GATE_FILES=${GATE_FILES}"
[[ -n "${GATE_STAGE:-}" ]] && job_env="${job_env},GATE_STAGE=${GATE_STAGE}"
[[ -n "${FORCE:-}" ]] && job_env="${job_env},FORCE=${FORCE}"

[[ -n "${SKIP_GATE:-}" ]] && job_env="${job_env},SKIP_GATE=${SKIP_GATE}"
qsub ${TERSE:+-terse} \
  -N "dfe_prep_${PREP_RUN}" \
  -cwd -j y -V \
  -o "${log_dir}/prep.${PREP_RUN}.\$JOB_ID.log" \
  -l h_rt="${SGE_PREP_RT}",h_data="${SGE_PREP_MEM}" \
  -pe shared "${SGE_PREP_PE}" \
  ${HOLD:+-hold_jid "${HOLD}"} \
  -v "${job_env}" \
  "${here}/worker_prep.sh" "${beds[@]}"
