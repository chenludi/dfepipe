#!/bin/bash
# The whole pipeline, submitted in one go: bed files in, DFE table out.
#
#   bash scripts/chains/run_all.sh <bed> [bed...]
#
#   DRY_RUN=1   print what would be submitted and exit, submitting nothing
#   FROM=<stage>  start at that stage, holding on nothing. One of
#                 prep sfs demo cache dfe summary. The stages before it must
#                 already have finished for the same run keys.
#
# Submitted in three PHASES, because two stages cannot be submitted before the
# stage above them has actually run:
#
#   submit_cache.sh reads best_inference.txt (:58) to find the demography
#   submit_dfe.sh   generates the command packs (:144), which needs the test
#                   spectrum and Na from the demographic fit
#
# An array's size must be known at qsub time, so neither can be deferred into
# the job. Instead each phase ends by submitting a LAUNCHER job, held on the
# phase's last job, which re-runs this script with the next FROM. The chain
# submits itself forward:
#
#   phase 1  FROM=prep    prep sfs concat demo bestdemo  -> launcher FROM=cache
#   phase 2  FROM=cache   cache                          -> launcher FROM=dfe
#   phase 3  FROM=dfe     fitneu dfe bestdfe summary bestpdf
#
# So `qstat` shows only the current phase. That is the cost of a dynamically
# sized array; the alternative is submitting a dfe array of the wrong size.
#
# Within a phase every job goes in at once, ordered by `qsub -hold_jid`. The
# order is forced by the data:
#
#   prep      array over beds          filters each bed to the analysable sites
#   sfs       array over chromosomes   per-chromosome spectra and mutation tables
#   concat    one job                  the genome-wide spectrum
#   demo      array over p0 files      the demographic fits
#   bestdemo  one job                  picks the winning demography
#   cache     array over beds          the selection-coefficient cache
#   fitneu    one job                  the neutral null likelihood
#   dfe       array per pdf            the DFE fits, one task per pack
#   bestdfe   one job per pdf          picks the winning fit, writes the DFE
#   summary   one job                  one row per bed
#   bestpdf   one job                  which pdf each spectrum supports
#
# Every job gates on its upstream check file before doing any work, array tasks
# included (scripts/sge/gate_upstream.sh). A FAIL upstream stops that job rather
# than carrying broken data into a long one; FORCE=1 overrides, as everywhere
# else. Two stages gate one further back than you might expect, because a
# stage's check file is written by the step AFTER its array - noted at each.
#
# Note on -hold_jid: it releases when the held job FINISHES, whatever its exit
# status. So a failed stage does not cancel what follows - each later step gates
# on its own upstream and stops itself, each saying so in its own log. That is
# why the gates matter more than the holds.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# qsub is not on PATH in a non-interactive shell - which is how this is normally
# run, over ssh or from a job. The submit_*.sh scripts each call this; the chain
# qsubs its own step jobs, so it needs it too.
# shellcheck source=../../config/env.sh
source "${REPO_DIR}/config/env.sh"
[[ "${DRY_RUN:-0}" == "1" ]] || ensure_scheduler

: "${OUT:?OUT not set}"
: "${PREP_RUN:?PREP_RUN not set}"
: "${FS_RUN:?FS_RUN not set}"
: "${DEMO_RUN:?DEMO_RUN not set}"
: "${DEMO_MODEL:?DEMO_MODEL not set}"
: "${CACHE_RUN:?CACHE_RUN not set}"
: "${CACHE_MARKER:?CACHE_MARKER not set}"
: "${DFE_RUN:?DFE_RUN not set}"
: "${SUMMARY_RUN:?SUMMARY_RUN not set}"
: "${PDFS:?PDFS not set}"
: "${SGE_STEP_RT:?SGE_STEP_RT not set}"
: "${SGE_STEP_MEM:?SGE_STEP_MEM not set}"
: "${SGE_STEP_PE:?SGE_STEP_PE not set}"

if [[ $# -eq 0 ]]; then
  echo "usage: run_all.sh <bed> [bed...]" >&2
  exit 2
fi

# Resolved while the caller's $PWD still means something: the jobs run elsewhere.
beds=()
for bed in "$@"; do
  [[ -f "$bed" ]] || { echo "run_all: no such bed: ${bed}" >&2; exit 1; }
  beds+=("$(cd "$(dirname "$bed")" && pwd)/$(basename "$bed")")
done

# Preflight the input before a single job is submitted. A bed named `1` instead
# of `chr1` would otherwise pass prep, produce nothing in all 22 sfs tasks, and
# surface as "cannot concatenate" an hour later. FAIL blocks; WARN is printed and
# the chain goes on, because every warning it raises is a judgement for the user.
if ! python3 "${REPO_DIR}/scripts/checks/check_bed.py" \
       --chroms "$CHROMS" --ctrl-flank "$CTRL_FLANK" \
       --mask-bed "${MASK_BED:-}" "${beds[@]}"; then
  echo "run_all: the input bed did not pass preflight; nothing was submitted." >&2
  echo "      Fix the bed, or set FORCE=1 to submit anyway." >&2
  [[ "${FORCE:-0}" == "1" ]] || exit 1
  echo "run_all: FORCE=1 - submitting despite the preflight failure" >&2
fi

STAGES=(prep sfs demo cache dfe summary)
from="${FROM:-prep}"
printf '%s\n' "${STAGES[@]}" | grep -qx "$from" || {
  echo "run_all: FROM=${from} is not one of: ${STAGES[*]}" >&2; exit 2; }

# Is this stage at or after the one we are starting from?
want() {
  local target="$1" seen=0
  for s in "${STAGES[@]}"; do
    [[ "$s" == "$from" ]] && seen=1
    [[ "$s" == "$target" ]] && { [[ "$seen" == 1 ]] && return 0 || return 1; }
  done
  return 1
}

log_dir="${OUT}/chains/${SUMMARY_RUN}"
# A dry run must leave the filesystem alone.
[[ "${DRY_RUN:-0}" == "1" ]] || mkdir -p "$log_dir"

# The check files each stage writes, which the next one gates on.
chk_prep="${OUT}/prep/check_prep.${PREP_RUN}.tsv"
chk_sfs="${OUT}/sfs/check_sfs.${FS_RUN}.tsv"
chk_demo="${OUT}/demo/check_demo.${FS_RUN}-${DEMO_MODEL}-${DEMO_RUN}.tsv"
# cache writes one per bed, inside that bed's cache directory.
chk_cache=()
for bed in "${beds[@]}"; do
  base="$(basename "$bed" .bed)"
  chk_cache+=("${OUT}/cache/${base}/${DEMO_RUN}-${CACHE_MARKER}-${CACHE_RUN}/check_cache.tsv")
done

dry="${DRY_RUN:-0}"
say()  { printf '  %-16s %s\n' "$1" "$2"; }
# Inside $( ) stdout is the job id, so anything explanatory goes to stderr.
note() { printf '  %-16s %s\n' "$1" "$2" >&2; }

# Submit one non-array step, printing nothing but its job id.
# step <name> <hold> <gate files> <gate stage> <script> [args...]
step() {
  local name="$1" hold="$2" gates="$3" gstage="$4"; shift 4
  local env="REPO_DIR=${REPO_DIR},OUT=${OUT},PREP_RUN=${PREP_RUN}"
  env="${env},FS_RUN=${FS_RUN},DEMO_RUN=${DEMO_RUN},DEMO_MODEL=${DEMO_MODEL}"
  env="${env},CACHE_RUN=${CACHE_RUN},CACHE_MARKER=${CACHE_MARKER}"
  env="${env},DFE_RUN=${DFE_RUN},SUMMARY_RUN=${SUMMARY_RUN}"
  [[ -n "$gates" ]] && env="${env},GATE_FILES=${gates},GATE_STAGE=${gstage}"
  [[ -n "${FORCE:-}" ]] && env="${env},FORCE=${FORCE}"
  [[ -n "${RESUME:-}" ]] && env="${env},RESUME=${RESUME}"
  [[ -n "${DFE_CACHE_FILE:-}" ]] && env="${env},DFE_CACHE_FILE=${DFE_CACHE_FILE}"

  if [[ "$dry" == "1" ]]; then
    note "$name" "qsub -N step_${name}_${SUMMARY_RUN}${hold:+ -hold_jid $hold} -- $*"
    echo "DRY-${name}"
    return 0
  fi
  qsub -terse \
    -N "step_${name}_${SUMMARY_RUN}" \
    -cwd -j y -V \
    -o "${log_dir}/${name}.\$JOB_ID.log" \
    -l h_rt="${SGE_STEP_RT}",h_data="${SGE_STEP_MEM}" \
    -pe shared "${SGE_STEP_PE}" \
    ${hold:+-hold_jid "${hold}"} \
    -v "$env" \
    "${REPO_DIR}/scripts/sge/worker_step.sh" "$@" | tail -1
}

# Submit one array stage through its own submit_*.sh, which knows its resources.
# array <name> <hold> <gate files> <gate stage> [args...]
array() {
  local name="$1" hold="$2" gates="$3" gstage="$4"; shift 4
  if [[ "$dry" == "1" ]]; then
    note "$name" "TERSE=1 HOLD=${hold:-none} GATE=${gstage:-none} submit_${name}.sh $*"
    echo "DRY-${name}"
    return 0
  fi
  TERSE=1 HOLD="${hold:-}" GATE_FILES="$gates" GATE_STAGE="$gstage" \
    bash "${REPO_DIR}/scripts/sge/submit_${name}.sh" "$@"
}

# The phase each stage belongs to, and which FROM starts the next one. A phase
# boundary sits wherever a submit script needs the previous stage's output.
next_phase_from() {
  case "$from" in
    prep|sfs|demo) echo cache ;;
    cache)         echo dfe ;;
    *)             echo "" ;;
  esac
}

# Submit the launcher for the next phase: this same script, FROM the next phase,
# as a job held on the last job of this one. Nothing is passed but FROM - the
# launcher re-reads the config, so there is one config read per phase.
launch_next() {
  local hold="$1" nextfrom="$2"
  local env="REPO_DIR=${REPO_DIR},OUT=${OUT},PREP_RUN=${PREP_RUN}"
  env="${env},FS_RUN=${FS_RUN},DEMO_RUN=${DEMO_RUN},DEMO_MODEL=${DEMO_MODEL}"
  env="${env},CACHE_RUN=${CACHE_RUN},CACHE_MARKER=${CACHE_MARKER}"
  env="${env},DFE_RUN=${DFE_RUN},SUMMARY_RUN=${SUMMARY_RUN},FROM=${nextfrom}"
  env="${env},PDFS=${PDFS},PACKNUMBER=${PACKNUMBER},DFE_TIMEOUT=${DFE_TIMEOUT}"
  env="${env},SGE_DFE_RT=${SGE_DFE_RT}"
  [[ -n "${FORCE:-}" ]] && env="${env},FORCE=${FORCE}"
  [[ -n "${RESUME:-}" ]] && env="${env},RESUME=${RESUME}"
  [[ -n "${DFE_CACHE_FILE:-}" ]] && env="${env},DFE_CACHE_FILE=${DFE_CACHE_FILE}"

  if [[ "$dry" == "1" ]]; then
    note "launch:${nextfrom}" "qsub -N chain_${nextfrom}_${SUMMARY_RUN} -hold_jid ${hold} -- run_all.sh FROM=${nextfrom}"
    echo "DRY-launch-${nextfrom}"
    return 0
  fi
  qsub -terse \
    -N "chain_${nextfrom}_${SUMMARY_RUN}" \
    -cwd -j y -V \
    -o "${log_dir}/launch_${nextfrom}.\$JOB_ID.log" \
    -l h_rt="00:30:00",h_data="4G" \
    -pe shared 1 \
    ${hold:+-hold_jid "${hold}"} \
    -v "$env" \
    "${REPO_DIR}/scripts/chains/run_all.sh" "${beds[@]}" | tail -1
}

# A job id from `qsub -terse`: an array prints `12345.1-22:1`, a plain job
# `12345`. -hold_jid wants the bare id.
bare() { local id="${1%%.*}"; echo "${id%% *}"; }

echo "run_all: ${#beds[@]} bed(s), from '${from}', logs in ${log_dir}"
[[ "$dry" == "1" ]] && echo "run_all: DRY_RUN=1, nothing will be submitted"

jid=""        # what the next stage holds on
jid_prep=""; jid_concat=""; jid_bestdemo=""; jid_cache=""

if want prep; then
  jid_prep="$(bare "$(array prep "" "" "" "${beds[@]}")")"
  say prep "$jid_prep"
  jid="$jid_prep"
fi

if want sfs; then
  jid_sfs="$(bare "$(array sfs "$jid" "$chk_prep" prep "${beds[@]}")")"
  say sfs "$jid_sfs"
  # concat gates on prep, not sfs: check_sfs.tsv is written by run_concat.sh
  # itself, so there is nothing from sfs to inspect yet. The sfs array has
  # already gated on prep, and each array task reports its own chromosome.
  jid_concat="$(bare "$(step concat "$jid_sfs" "$chk_prep" prep \
                        scripts/sfs/run_concat.sh "${beds[@]}")")"
  say concat "$jid_concat"
  jid="$jid_concat"
fi

if want demo; then
  jid_demo="$(bare "$(array demo "$jid" "$chk_sfs" sfs "${beds[@]}")")"
  say demo "$jid_demo"
  jid_bestdemo="$(bare "$(step bestdemo "$jid_demo" "$chk_sfs" sfs \
                         scripts/demo/run_bestdemo.sh "${beds[@]}")")"
  say bestdemo "$jid_bestdemo"
  jid="$jid_bestdemo"
fi

# Phase 1 ends here: submit_cache.sh reads best_inference.txt, which bestdemo
# has not written yet.
if [[ "$from" != cache && "$from" != dfe && "$from" != summary ]]; then
  jid_launch="$(bare "$(launch_next "$jid" cache)")"
  say "launch:cache" "$jid_launch"
  echo "run_all: phase 1 submitted. The rest follows automatically; watch"
  echo "         ${log_dir}/launch_cache.*.log"
  exit 0
fi

if want cache; then
  jid_cache="$(bare "$(array cache "$jid" "$chk_demo" demo "${beds[@]}")")"
  say cache "$jid_cache"
  jid="$jid_cache"
  # Phase 2 ends here: submit_dfe.sh generates the packs, which needs the
  # spectrum and Na, and the array size depends on how many packs result.
  if [[ "$from" == cache ]]; then
    jid_launch="$(bare "$(launch_next "$jid_cache" dfe)")"
    say "launch:dfe" "$jid_launch"
    echo "run_all: phase 2 submitted. Watch ${log_dir}/launch_dfe.*.log"
    exit 0
  fi
fi

dfe_last=()
if want dfe; then
  jid_neu="$(bare "$(step fitneu "$jid" "${chk_cache[*]}" cache \
                     scripts/dfe/run_fitneu.sh "${beds[@]}")")"
  say fitneu "$jid_neu"
  dfe_last+=("$jid_neu")
  for pdf in $PDFS; do
    [[ "$pdf" == "neu" ]] && continue
    jid_fit="$(bare "$(array dfe "$jid" "${chk_cache[*]}" cache "$pdf" "${beds[@]}")")"
    say "dfe:${pdf}" "$jid_fit"
      # bestdfe gates on the cache, the nearest checked upstream: the dfe array
    # writes no check file of its own - check_dfe.py runs inside run_bestdfe.sh
    # and run_fitneu.sh, after the fits exist.
    jid_best="$(bare "$(step "bestdfe_${pdf}" "$jid_fit" "${chk_cache[*]}" cache \
                        scripts/dfe/run_bestdfe.sh "$pdf" "${beds[@]}")")"
    say "bestdfe:${pdf}" "$jid_best"
    dfe_last+=("$jid_best")
  done
fi

if want summary; then
  hold="$(IFS=,; echo "${dfe_last[*]:-}")"
  chk_dfe=""
  for pdf in $PDFS; do
    [[ "$pdf" == "neu" ]] && continue
    chk_dfe="${chk_dfe}${chk_dfe:+ }${OUT}/dfe/check_dfe.${CACHE_RUN}-${CACHE_MARKER}-${pdf}-${DFE_RUN}.tsv"
  done
  jid_sum="$(bare "$(step summary "$hold" "$chk_dfe" dfe \
                     scripts/summary/run_summary.sh "${beds[@]}")")"
  say summary "$jid_sum"
  jid_bp="$(bare "$(step bestpdf "$jid_sum" "${OUT}/summary/check_summary.${SUMMARY_RUN}.tsv" summary \
                    scripts/summary/run_bestpdf.sh)")"
  say bestpdf "$jid_bp"
  echo "run_all: the table will be ${OUT}/summary/${SUMMARY_RUN}/bestpdf_sumdfe.tsv"
fi

exit 0
