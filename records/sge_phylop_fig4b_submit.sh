#!/bin/bash
# Launcher job for record_261026_phylop_fig4b_run.sh. It runs that record script
# with SUBMIT=1 on a compute node and does nothing else, so the record script's
# own qsubs come from here instead of from the login node.
#
# This exists because of the preflight, not the submission. run_all.sh:98 reads
# every interval of every input bed before submitting anything, and the record
# script invokes the chain twice - once for the dry run at its step 2, once for
# the submit at step 3 - so the preflight is paid twice. The phyloP cross-bin
# beds are near-per-base: 39 of them hold 1,390,879,772 intervals at 1.03 bp
# each, which is about 6 minutes a pass with numpy and about 71 without. The
# login node's /usr/bin/python3 is 3.6.8 and has no numpy, and its governor
# killed a long pass at ~13 minutes on 2026-10-27. That is why the 39-bed run
# could not be launched from there at all. activate_conda here gives python3
# 3.12.1 with numpy 1.26.4, measured at 6 s on one 21.9 million interval bed
# against 67 s for the login node's python.
#
# Not run by hand. UGE spools a copy of this script, so $BASH_SOURCE points at
# the spool copy; everything it needs arrives through qsub -v and -V.
set -euo pipefail
: "${REPO_DIR:?REPO_DIR must come from the driver}"
: "${SUMMARY_RUN:?SUMMARY_RUN must come from the driver}"

cd "$REPO_DIR"

# shellcheck source=../config/env.sh
source "${REPO_DIR}/config/env.sh"
# Only the preflight needs this. Every stage script activates its own env
# (run_sfs_chrom.sh, run_demo.sh, run_dfe.sh and the rest), so inheriting this
# one through run_all's -V changes nothing downstream.
activate_conda

echo "sge_phylop_fig4b_submit: host=$(hostname) job=${JOB_ID:-none} started $(date)"
echo "  python3   $(command -v python3)"
if python3 -c 'import numpy' 2>/dev/null; then
  echo "  numpy     $(python3 -c 'import numpy; print(numpy.__version__)')"
else
  # Said here rather than discovered later: without numpy the two preflights
  # take about 71 minutes each instead of 6.
  echo "  numpy     absent - each preflight will fall back to the line reader"
fi
echo "  selection ONLY='${ONLY:-}' EXCEPT='${EXCEPT:-}' PROBE='${PROBE:-0}'"
echo "  run keys  SUMMARY_RUN=${SUMMARY_RUN} PHYLOP_SFS_RT=${PHYLOP_SFS_RT:-unset}"

SUBMIT=1 VIA_JOB=0 bash records/record_261026_phylop_fig4b_run.sh
