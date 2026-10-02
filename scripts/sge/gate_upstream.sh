#!/bin/bash
# Sourced by every worker. Runs the upstream gate, if the chain asked for one.
#
#   GATE_FILES  space-separated check files to inspect before doing any work
#   GATE_STAGE  the upstream stage's name, for the message
#
# The gate runs inside the job, not at submission time, because the check file
# does not exist yet when the chain is submitted: every job goes in at once and
# SGE orders them with -hold_jid.
#
# An array task gates like any other, so all 22 tasks of a stage each consult
# the same file. That is a few milliseconds each and it is what makes the gate
# work at all: -hold_jid releases when the held job FINISHES, whatever its exit
# status, so SGE cannot cancel what follows a failure. Each task has to stop
# itself.
#
# Nothing happens when GATE_FILES is unset, which is the case for a stage run by
# hand - re-running one stage against older output is a normal workflow here.

dfepipe_gate_upstream() {
  : "${REPO_DIR:?REPO_DIR not set}"
  local f
  for f in ${GATE_FILES:-}; do
    echo "gate: ${GATE_STAGE:-upstream} -> ${f}"
    if ! python3 "${REPO_DIR}/scripts/checks/gate.py" \
           "$f" "${GATE_STAGE:-upstream}"; then
      return 1
    fi
  done
  return 0
}

# The one line every worker needs. Takes the worker's own name for the message.
dfepipe_gate_or_exit() {
  if ! dfepipe_gate_upstream; then
    echo "${1:-worker}: the upstream gate stopped this task; nothing was run." >&2
    echo "      Fix the stage it names, or set SKIP_GATE=1 to carry on." >&2
    echo "      Not FORCE=1: that also means overwrite a finished directory, redo" >&2
    echo "      every fit, rebuild a cache and overwrite a table, and it is" >&2
    echo "      forwarded to every later job, so it rewrites work from here on." >&2
    exit 1
  fi
}
