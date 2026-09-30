#!/bin/bash
# SGE worker for the cache stage. Not run by hand: submit_cache.sh submits it as
# an array job, one task per bed.
#
# UGE spools a COPY of the script into its own spool directory, so $BASH_SOURCE
# here points at that copy and cannot locate the package. REPO_DIR is therefore
# passed in by the launcher and is the only way back to the repo.
#
# SGE_TASK_ID indexes into the bed list. The launcher forwards CACHE_BEDS, the
# list it counted the array size from, so the task's idea of which bed it owns
# cannot drift from the one that was checked - re-deriving it here would shift
# every index if a bed's cache appeared while the job sat in the queue.
set -euo pipefail
: "${REPO_DIR:?REPO_DIR must be exported by the launcher (submit_cache.sh)}"
: "${CACHE_BEDS:?CACHE_BEDS must be exported by the launcher (submit_cache.sh)}"
: "${SGE_TASK_ID:?SGE_TASK_ID not set: worker_cache.sh runs as an array task}"

read -r -a _beds <<< "$CACHE_BEDS"
idx=$(( SGE_TASK_ID - 1 ))
if (( idx < 0 || idx >= ${#_beds[@]} )); then
  echo "worker_cache: SGE_TASK_ID=${SGE_TASK_ID} is outside CACHE_BEDS" \
       "(${#_beds[@]} entries)" >&2
  exit 1
fi
bed="${_beds[$idx]}"

echo "worker_cache: host=$(hostname) job=${JOB_ID:-none} task=${SGE_TASK_ID}" \
     "bed=$(basename "$bed") started $(date)"
bash "${REPO_DIR}/scripts/cache/run_cache.sh" "$bed"
echo "worker_cache: bed=$(basename "$bed") finished $(date)"
