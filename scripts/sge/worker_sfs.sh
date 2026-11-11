#!/bin/bash
# SGE worker for the sfs stage. Not run by hand: submit_sfs.sh submits it as an
# array job, one task per (chromosome, bed) pair.
#
# UGE spools a COPY of the script into its own spool directory, so $BASH_SOURCE
# here points at that copy and cannot locate the package. REPO_DIR is therefore
# passed in by the launcher and is the only way back to the repo.
#
# SGE_TASK_ID indexes into CHROMS x the bed list, beds fastest. The launcher
# forwards CHROMS verbatim and passes the same bed list it sized the array
# from, so the task's idea of either cannot drift from what was checked.
set -euo pipefail
: "${REPO_DIR:?REPO_DIR must be exported by the launcher (submit_sfs.sh)}"
: "${CHROMS:?CHROMS must be exported by the launcher (submit_sfs.sh)}"
: "${SGE_TASK_ID:?SGE_TASK_ID not set: worker_sfs.sh runs as an array task}"

read -r -a _chroms <<< "$CHROMS"
_beds=("$@")
if (( ${#_beds[@]} == 0 )); then
  echo "worker_sfs: no bed arguments" >&2
  exit 1
fi
_total=$(( ${#_chroms[@]} * ${#_beds[@]} ))
idx=$(( SGE_TASK_ID - 1 ))
if (( idx < 0 || idx >= _total )); then
  echo "worker_sfs: SGE_TASK_ID=${SGE_TASK_ID} is outside CHROMS x beds" \
       "(${#_chroms[@]} x ${#_beds[@]} = ${_total} tasks)" >&2
  exit 1
fi
chrom="${_chroms[$(( idx / ${#_beds[@]} ))]}"
bed="${_beds[$(( idx % ${#_beds[@]} ))]}"

# shellcheck source=gate_upstream.sh
source "${REPO_DIR}/scripts/sge/gate_upstream.sh"
dfepipe_gate_or_exit "worker_sfs"

echo "worker_sfs: host=$(hostname) job=${JOB_ID:-none} task=${SGE_TASK_ID}" \
     "chrom=${chrom} bed=$(basename "$bed") started $(date)"
bash "${REPO_DIR}/scripts/sfs/run_sfs_chrom.sh" "$chrom" "$bed"
echo "worker_sfs: chrom=${chrom} bed=$(basename "$bed") finished $(date)"
