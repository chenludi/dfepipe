#!/bin/bash
# SGE worker for the dfe stage. Not run by hand: submit_dfe.sh submits it as an
# array job, one task per (bed, command pack) pair - or as a single job for the
# `neu` null, which is one closed-form likelihood per bed.
#
# UGE spools a COPY of the script into its own spool directory, so $BASH_SOURCE
# here points at that copy and cannot locate the package. REPO_DIR is therefore
# passed in by the launcher and is the only way back to the repo.
#
# SGE_TASK_ID indexes into two aligned lists the launcher forwards, DFE_TASK_BEDS
# and DFE_PACKS. They are the lists it counted the array size from, so a task's
# idea of what it owns cannot drift - a plain `ls` here would shift every index
# if the packs were regenerated while the job sat in the queue. The number of
# packs depends on how many starting points survive the bed's own gamma range,
# so it can differ between beds; pairing them explicitly is what keeps that safe.
set -euo pipefail
: "${REPO_DIR:?REPO_DIR must be exported by the launcher (submit_dfe.sh)}"
: "${DFE_PDF:?DFE_PDF must be exported by the launcher (submit_dfe.sh)}"

if [[ "$DFE_PDF" == "neu" ]]; then
  echo "worker_dfe: host=$(hostname) job=${JOB_ID:-none} pdf=neu started $(date)"
  bash "${REPO_DIR}/scripts/dfe/run_fitneu.sh" "$@"
  echo "worker_dfe: pdf=neu finished $(date)"
  exit 0
fi

: "${DFE_PACKS:?DFE_PACKS must be exported by the launcher (submit_dfe.sh)}"
: "${DFE_TASK_BEDS:?DFE_TASK_BEDS must be exported by the launcher (submit_dfe.sh)}"
: "${SGE_TASK_ID:?SGE_TASK_ID not set: worker_dfe.sh runs as an array task}"

read -r -a _packs <<< "$DFE_PACKS"
read -r -a _beds <<< "$DFE_TASK_BEDS"
if (( ${#_packs[@]} != ${#_beds[@]} )); then
  echo "worker_dfe: DFE_PACKS has ${#_packs[@]} entries but DFE_TASK_BEDS has" \
       "${#_beds[@]}; the launcher built them inconsistently" >&2
  exit 1
fi
idx=$(( SGE_TASK_ID - 1 ))
if (( idx < 0 || idx >= ${#_packs[@]} )); then
  echo "worker_dfe: SGE_TASK_ID=${SGE_TASK_ID} is outside DFE_PACKS" \
       "(${#_packs[@]} entries)" >&2
  exit 1
fi
pack="${_packs[$idx]}"
bed="${_beds[$idx]}"

echo "worker_dfe: host=$(hostname) job=${JOB_ID:-none} task=${SGE_TASK_ID}" \
     "pdf=${DFE_PDF} bed=$(basename "$bed") pack=${pack} started $(date)"
DFE_PACK="$pack" bash "${REPO_DIR}/scripts/dfe/run_dfe.sh" "$DFE_PDF" "$bed"
echo "worker_dfe: ${pack} finished $(date)"
