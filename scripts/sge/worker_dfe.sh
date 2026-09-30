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

# When this task's h_rt runs out, as an epoch second. run_dfe.sh uses it to stop
# before a fit it cannot finish, rather than have the scheduler kill one
# part-way: PACKNUMBER x DFE_TIMEOUT can exceed h_rt, so that is not a rare
# case. Best effort - without SGE_DFE_RT the guard is simply off.
if [[ -n "${SGE_DFE_RT:-}" ]]; then
  IFS=: read -r _h _m _s <<< "$SGE_DFE_RT"
  _rt=$(( 10#${_h:-0} * 3600 + 10#${_m:-0} * 60 + 10#${_s:-0} ))
  if (( _rt > 0 )); then
    export DFE_TASK_DEADLINE=$(( $(date +%s) + _rt - 300 ))   # 5 min of margin
    echo "worker_dfe: h_rt=${SGE_DFE_RT}; will start no fit after" \
         "$(date -d "@${DFE_TASK_DEADLINE}" 2>/dev/null || echo "+${_rt}s")"
  fi
fi

DFE_PACK="$pack" bash "${REPO_DIR}/scripts/dfe/run_dfe.sh" "$DFE_PDF" "$bed"
echo "worker_dfe: ${pack} finished $(date)"
