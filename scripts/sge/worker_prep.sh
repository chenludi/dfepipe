#!/bin/bash
# SGE worker for the prep stage. Not run by hand: submit_prep.sh submits it.
#
# UGE spools a COPY of the script into its own spool directory, so $BASH_SOURCE
# here points at that copy and cannot locate the package. REPO_DIR is therefore
# passed in by the launcher and is the only way back to the repo.
set -euo pipefail
: "${REPO_DIR:?REPO_DIR must be exported by the launcher (submit_prep.sh)}"

# shellcheck source=gate_upstream.sh
source "${REPO_DIR}/scripts/sge/gate_upstream.sh"
dfepipe_gate_or_exit "worker_prep"

echo "worker_prep: host=$(hostname) job=${JOB_ID:-none} started $(date)"
bash "${REPO_DIR}/scripts/prep/run_prep.sh" "$@"
echo "worker_prep: finished $(date)"
