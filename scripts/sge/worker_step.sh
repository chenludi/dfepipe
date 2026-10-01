#!/bin/bash
# Runs one pipeline step as its own SGE job, after gating on the upstream check.
#
#   qsub ... worker_step.sh <script> [args...]
#
# For the steps that are not arrays and so have no submit_*.sh of their own:
# run_concat.sh, run_bestfit.sh, run_fitneu.sh, run_bestdfe.sh, run_summary.sh
# and run_bestpdf.sh. Each can only start once the array before it has finished,
# which is what `qsub -hold_jid` arranges.
#
# Exported by the chain, all optional:
#   GATE_FILES  space-separated check files to gate on before running the step
#   GATE_STAGE  the upstream stage's name, for the gate's message
#
# The gate runs here rather than at submission time because the check file does
# not exist yet when the chain is submitted - every job goes in at once, and SGE
# orders them. A failed gate stops this step; it does not stop the jobs held on
# it, because -hold_jid releases on completion whatever the exit status. That is
# deliberate: each stage gates on its own upstream, so a break stops every step
# that depends on the broken data, each saying so in its own log.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"

if [[ $# -lt 1 ]]; then
  echo "usage: worker_step.sh <script> [args...]" >&2
  exit 2
fi
script="$1"; shift
[[ "$script" = /* ]] || script="${REPO_DIR}/${script}"
[[ -f "$script" ]] || { echo "worker_step: no such script: ${script}" >&2; exit 2; }

echo "worker_step: host=$(hostname) job=${JOB_ID:-none} started $(date)"
echo "worker_step: $(basename "$script") $*"

# shellcheck source=gate_upstream.sh
source "${REPO_DIR}/scripts/sge/gate_upstream.sh"
dfepipe_gate_or_exit "worker_step"

bash "$script" "$@"
rc=$?
echo "worker_step: $(basename "$script") finished $(date) rc=${rc}"
exit "$rc"
