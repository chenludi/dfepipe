#!/bin/bash
# submit_sfs.sh <bed> [bed...]
#
# One array job over CHROMS: control pairing and the VCF intersect are the slow
# part of this pipeline and are independent per chromosome, so they run in
# parallel. Every task handles every bed for its own chromosome.
#
#   TERSE=1   print only the job id, for chaining
#   HOLD=<id> hold until that job finishes (e.g. the prep job)
#   FORCE=1   forwarded to run_sfs_chrom.sh, which then redoes finished chromosomes
#
# Concatenation is not submitted here: it can only run once every task has
# finished. Submit it after this job, holding on its id:
#
#   jid=$(TERSE=1 bash scripts/sge/submit_sfs.sh mybed.bed)
#   HOLD="${jid%%.*}" ... scripts/sfs/run_concat.sh mybed.bed
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
: "${PREP_RUN:?PREP_RUN not set}"
: "${CHROMS:?CHROMS not set}"
: "${SGE_SFS_RT:?SGE_SFS_RT not set}"
: "${SGE_SFS_MEM:?SGE_SFS_MEM not set}"

if [[ $# -eq 0 ]]; then
  echo "usage: submit_sfs.sh <bed> [bed...]" >&2
  exit 2
fi

n_chroms=$(echo "$CHROMS" | wc -w | tr -d ' ')
if [[ "$n_chroms" -lt 1 ]]; then
  echo "submit_sfs: CHROMS is empty" >&2
  exit 2
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
# check_sfs.py names this directory when a chromosome produced nothing.
log_dir="${OUT}/sfs/log"
mkdir -p "$log_dir"

# The config the launcher validated is forwarded explicitly, so a task cannot
# silently pick up a different run key - or a different chromosome list, which
# would shift every task's index - from the worker node's environment. -V
# carries the rest of the exported environment (the site profile).
job_env="REPO_DIR=${REPO_DIR},OUT=${OUT},FS_RUN=${FS_RUN},PREP_RUN=${PREP_RUN},CHROMS=${CHROMS}"
# The chain sets these; a stage submitted by hand leaves them unset and
# no gate runs. scripts/sge/gate_upstream.sh
[[ -n "${GATE_FILES:-}" ]] && job_env="${job_env},GATE_FILES=${GATE_FILES}"
[[ -n "${GATE_STAGE:-}" ]] && job_env="${job_env},GATE_STAGE=${GATE_STAGE}"
[[ -n "${FORCE:-}" ]] && job_env="${job_env},FORCE=${FORCE}"

[[ -n "${SKIP_GATE:-}" ]] && job_env="${job_env},SKIP_GATE=${SKIP_GATE}"
qsub ${TERSE:+-terse} \
  -N "dfe_sfs_${FS_RUN}" \
  -t "1-${n_chroms}" \
  -cwd -j y -V \
  -o "${log_dir}/sfs.${FS_RUN}.\$JOB_ID.\$TASK_ID.log" \
  -l h_rt="${SGE_SFS_RT}",h_data="${SGE_SFS_MEM}" \
  -pe shared "${SGE_SFS_PE}" \
  ${HOLD:+-hold_jid "${HOLD}"} \
  -v "${job_env}" \
  "${here}/worker_sfs.sh" "${beds[@]}"
