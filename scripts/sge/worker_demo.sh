#!/bin/bash
# SGE worker for the demo stage. Not run by hand: submit_demo.sh submits it as an
# array job, one task per initial-parameter file in DEMO_P0_DIR.
#
# UGE spools a COPY of the script into its own spool directory, so $BASH_SOURCE
# here points at that copy and cannot locate the package. REPO_DIR is therefore
# passed in by the launcher and is the only way back to the repo.
#
# SGE_TASK_ID indexes into the p0 file list. The launcher forwards DEMO_P0_FILES,
# the list it counted the array size from, so the task's idea of which file it
# owns cannot drift from the one that was checked - a plain `ls` here would shift
# every index if a file were added while the job sat in the queue.
set -euo pipefail
: "${REPO_DIR:?REPO_DIR must be exported by the launcher (submit_demo.sh)}"
: "${DEMO_P0_FILES:?DEMO_P0_FILES must be exported by the launcher (submit_demo.sh)}"
: "${SGE_TASK_ID:?SGE_TASK_ID not set: worker_demo.sh runs as an array task}"

read -r -a _p0files <<< "$DEMO_P0_FILES"
idx=$(( SGE_TASK_ID - 1 ))
if (( idx < 0 || idx >= ${#_p0files[@]} )); then
  echo "worker_demo: SGE_TASK_ID=${SGE_TASK_ID} is outside DEMO_P0_FILES" \
       "(${#_p0files[@]} entries)" >&2
  exit 1
fi
p0file="${_p0files[$idx]}"

# shellcheck source=gate_upstream.sh
source "${REPO_DIR}/scripts/sge/gate_upstream.sh"
dfepipe_gate_or_exit "worker_demo"

echo "worker_demo: host=$(hostname) job=${JOB_ID:-none} task=${SGE_TASK_ID}" \
     "p0file=${p0file} started $(date)"
bash "${REPO_DIR}/scripts/demo/run_demo.sh" "$p0file" "$@"
echo "worker_demo: p0file=${p0file} finished $(date)"
