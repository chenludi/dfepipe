#!/bin/bash
# 2026-10-26  DFEs for the PhyloP version of Fig. 4B.
#
# Takes the 40 beds built by record_261026_phylop_crossbins.sh and runs them
# through dfepipe from the sfs stage. Two settings differ from the profile and
# both are deliberate.
#
# 1. The chain starts at sfs, not prep (FROM=sfs, SKIP_GATE=1).
#    These beds are already s.nc_exHET_noq.*: intersections of per-bin beds that
#    filt_nc-unmask.sh had already put through HET, quiescent, CDS and the
#    accessibility mask. An intersection of two such sets is still free of all
#    four, so prep's filters are no-ops on content. Running them anyway would
#    cost four bedtools passes over 40 files of 9-54 million intervals inside
#    prep's single job, which is hours of work for no change. Each bed is instead
#    hard-linked into place as prep's output, where run_sfs_chrom.sh:146 looks
#    for it, and SKIP_GATE=1 says there is no prep check file to gate on.
#
# 2. MERGE_BEFORE_WINDOW=1, against the profile's 0.
#    This is NOT cosmetic. The control window is centred on the midpoint of each
#    test interval, so whether the bed is merged first decides what the window is
#    centred on - and these beds are near-per-base, so the two give very
#    different control sets. The profile's 0 is safe only because prep merges
#    (run_prep.sh:112); skipping prep removes that, so sfs has to do it.
#    windctrl_splitchr.sh:41-42 merges, added 2024-11-04 16:59, and the published
#    PhastCons conditional run wrote its per-chromosome files at 17:11-17:13 the
#    same day, so the figure this one is being compared against was merged too.
#
#   bash records/record_261026_phylop_fig4b_run.sh              stage, check, dry run
#   PROBE=1 SUBMIT=1 bash records/...                           one bed only, to time it
#   SUBMIT=1 bash records/record_261026_phylop_fig4b_run.sh     all 40
#   EXCEPT="<bed base>" SUMMARY_RUN=<new> SUBMIT=1 bash records/...   a later batch
#   VIA_JOB=1 SUBMIT=1 bash records/...    launch from a compute node, not here
#
# Run PROBE=1 first. SGE_SFS_RT is a guess here: an sfs task's cost scales with
# interval count, and these beds carry 20-30 times what GRCh38-ELS did, where a
# task took about 14 minutes. 12h is generous rather than measured, and the ELS
# run on 2026-10-25 had to have its sfs runtime raised by hand for less reason
# than this. One bed through sfs gives the real number for the other 39.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
source config/profiles/hoffman_replicate.conf
source config/pipeline.conf
source config/env.sh

PRIMARY=phyloP243way
MAMMAL=phyloP447way
NCON_LABEL=40-100
BINS=(0-5 5-10 10-15 15-20 20-25 25-30 30-35 35-40 40-45 45-50 50-55 55-60 \
      60-65 65-70 70-75 75-80 80-85 85-90 90-95 95-100)

export PREP_RUN=prep_pp261026 FS_RUN=fs_pp261026 DEMO_RUN=demo_pp261026 \
       CACHE_RUN=cache_pp261026 DFE_RUN=dfe_pp261026
# run_summary.sh:77 and run_bestpdf.sh:75 refuse a table that already exists, so
# a batch run after an earlier one needs its own SUMMARY_RUN. The other five
# keys stay fixed, so every batch's outputs land together.
export SUMMARY_RUN="${SUMMARY_RUN:-sum_pp261026}"
export MERGE_BEFORE_WINDOW=1
export SKIP_GATE=1
# Set unconditionally: pipeline.conf has already run by this point and put
# 3:40:00 in SGE_SFS_RT, so "${SGE_SFS_RT:-12:00:00}" would silently keep 3:40.
# Override with PHYLOP_SFS_RT, not SGE_SFS_RT, for the same reason.
export SGE_SFS_RT="${PHYLOP_SFS_RT:-12:00:00}"

src="$REFS_DIR/data_bed/conserve_phylop_ncon_261026"
lens="${src}/lengths_261026.tsv"      # bed -> bp, written by the crossbins run
stage="$REFS_DIR/data_bed/phylop_ncon_pp261026"
log="$REFS_DIR/data_bed/conserve_phylop_ncon_261026/run_pp261026.log"
logdir="$REFS_DIR/data_bed/conserve_phylop_ncon_261026/log"

mkdir -p "$stage"
say () { echo "$@" | tee -a "$log"; }

say "=============================================================="
say " phyloP Fig. 4B DFEs   $(date '+%Y-%m-%d %H:%M:%S')   host=$(hostname)"
say " run keys: ${PREP_RUN} ${FS_RUN} ${DEMO_RUN} ${CACHE_RUN} ${DFE_RUN} ${SUMMARY_RUN}"
say " FROM=sfs  SKIP_GATE=1  MERGE_BEFORE_WINDOW=1  SGE_SFS_RT=${SGE_SFS_RT}"
say "=============================================================="

# ---------------------------------------------------- 1. name and link them --
# The union beds are named s.nc_exHET_noq.<X>.bed. dfepipe derives a bed's base
# from its filename and then looks for prep output at
# <base dir>/s.nc_exHET_noq.<base>.bed, so passing them under their own name
# would give every downstream file a doubled s.nc_exHET_noq. prefix and group
# names that no longer line up with Data S2. Dropping the prefix for the bed
# argument makes the union bed itself exactly what prep would have written.
# Hard links, so 48 GB is not copied and both names are the same bytes.
say
say "1. staging bed arguments in ${stage}"
beds=(); n_link=0; n_miss=0
link_one () {   # link_one <base>
  local base="$1" s="${src}/s.nc_exHET_noq.${1}.bed" a="${stage}/${1}.bed"
  if [[ ! -s "$s" ]]; then say "   missing: $s"; n_miss=$((n_miss + 1)); return 0; fi
  [[ -e "$a" ]] || ln "$s" "$a"
  local p="${OUT}/prep/${base}/${PREP_RUN}"
  mkdir -p "$p"
  [[ -e "${p}/s.nc_exHET_noq.${base}.bed" ]] || ln "$s" "${p}/s.nc_exHET_noq.${base}.bed"
  {
    echo "stage          prep (not run; staged by record_261026_phylop_fig4b_run.sh)"
    echo "date           $(date '+%Y-%m-%d %H:%M:%S')"
    echo "source         $s"
    echo "why            the bed is already s.nc_exHET_noq.*: HET, quiescent, CDS"
    echo "               and mask were applied to both per-bin beds before they"
    echo "               were intersected, so prep's filters are no-ops on it."
    echo "not merged     prep would also merge; sfs does it instead, with"
    echo "               MERGE_BEFORE_WINDOW=1. See the header of this record."
  } > "${p}/spec.txt"
  # check_sfs.py:218 reads after_mask out of prep's lengths.tsv for the totals
  # row; without it bp_prep, n_intervals_prep and hasqui_ratio come out NA for
  # the bed as a whole, which is the one number worth having per bed. The bp is
  # taken from the crossbins length table rather than rescanned, since these
  # beds are 40 files of over a gigabyte. input and after_mask are the same
  # number here because prep's four filters are no-ops on an already filtered
  # bed; spec.txt beside this file says prep did not run.
  local bp
  bp=$(awk -F'\t' -v n="s.nc_exHET_noq.${base}.bed" '$1 == n {print $6; exit}' "$lens")
  if [[ -z "$bp" ]]; then
    say "   no bp for ${base} in ${lens}"; n_miss=$((n_miss + 1)); return 0
  fi
  printf 'stage\tbp\ninput\t%s\nafter_mask\t%s\n' "$bp" "$bp" > "${p}/lengths.tsv"
  beds+=("$a"); n_link=$((n_link + 1))
}
for b in "${BINS[@]}"; do link_one "${PRIMARY}_${b}_${MAMMAL}_${NCON_LABEL}"; done
for b in "${BINS[@]}"; do link_one "${PRIMARY}_${NCON_LABEL}_${MAMMAL}_${b}"; done
if (( n_miss )); then
  say "   ${n_miss} bed(s) missing. Run record_261026_phylop_crossbins.sh first."
  exit 1
fi
say "   ${n_link} beds linked, and the same bytes linked as prep output"
say "   link count on one of them: $(stat -c %h "${beds[0]}")  (3 = src + argument + prep output)"

# ONLY and EXCEPT name beds explicitly, by base name, space or comma separated.
# Prefer them to SKIP_DONE: they state which beds run instead of inferring it
# from what is on disk, and the chosen list is written to the log. An unknown
# name is an error, so a typo cannot silently shrink the run.
#
#   ONLY="phyloP243way_0-5_phyloP447way_40-100"        run just that one
#   EXCEPT="phyloP243way_0-5_phyloP447way_40-100"      run the other 39
select_beds() {
  local mode="$1" list="$2"
  local a base w hit
  local keep=() want=() unknown=()
  read -r -a want <<< "${list//,/ }"
  for w in "${want[@]}"; do
    hit=0
    for a in "${beds[@]}"; do
      [[ "$(basename "$a" .bed)" == "$w" ]] && { hit=1; break; }
    done
    (( hit )) || unknown+=("$w")
  done
  if (( ${#unknown[@]} )); then
    say "   ${mode}: no such bed: ${unknown[*]}"
    exit 1
  fi
  for a in "${beds[@]}"; do
    base="$(basename "$a" .bed)"; hit=0
    for w in "${want[@]}"; do [[ "$base" == "$w" ]] && { hit=1; break; }; done
    if [[ "$mode" == ONLY ]]; then (( hit )) && keep+=("$a"); else (( hit )) || keep+=("$a"); fi
  done
  if (( ${#keep[@]} == 0 )); then say "   ${mode}: nothing left to run."; exit 1; fi
  beds=("${keep[@]}")
  say "   ${mode}=${list}"
  say "   ${#beds[@]} bed(s) selected:"
  for a in "${beds[@]}"; do say "     $(basename "$a" .bed)"; done
}
[[ -n "${ONLY:-}"   ]] && select_beds ONLY   "$ONLY"
[[ -n "${EXCEPT:-}" ]] && select_beds EXCEPT "$EXCEPT"

# SKIP_DONE=1 drops any bed that already has a concatenated test spectrum for
# this FS_RUN. That is what makes it safe to add beds to a run already in
# flight: run_sfs_chrom.sh:161 and run_concat.sh:146 refuse to overwrite a
# finished output rather than clobbering it. Without it the already-done bed's
# sfs tasks fail and concat refuses that one bed; concat still concatenates
# every other bed and only exits non-zero (run_concat.sh sets failed=1 and
# continues, exiting at line 247), so the cost is noise and a failed-looking
# job, not lost work.
# It keys on the concatenated spectrum, so a bed still part-way through sfs is
# NOT skipped and would still collide; only add beds between stages, not during
# one.
if [[ "${SKIP_DONE:-0}" == "1" ]]; then
  keep=(); n_skip=0
  for a in "${beds[@]}"; do
    base="$(basename "$a" .bed)"
    if compgen -G "${OUT}/sfs/${base}/${FS_RUN}/has_qui.*.chrs_1t22.fs" > /dev/null; then
      say "   skipping ${base}: already concatenated for ${FS_RUN}"
      n_skip=$((n_skip + 1))
    else
      keep+=("$a")
    fi
  done
  beds=("${keep[@]}")
  say "   SKIP_DONE=1: ${n_skip} bed(s) already done, ${#beds[@]} to run"
  if (( ${#beds[@]} == 0 )); then
    say "   nothing left to do."
    exit 0
  fi
fi

if [[ "${PROBE:-0}" == "1" ]]; then
  beds=("${beds[0]}")
  say
  say "   PROBE=1: only $(basename "${beds[0]}")"
  say "   This measures one sfs task so the other 39 can be sized from a number"
  say "   rather than a guess. It writes into the same run keys, so the rest can"
  say "   be added later with the same command and no PROBE."
fi

# ------------------------------------- 1b. hand the launch to a compute node --
# run_all.sh:98 preflights every interval of every bed, and the chain below is
# invoked twice, so the preflight is paid twice. On these near-per-base beds
# that is ~6 minutes a pass with numpy and ~71 without, and the login node's
# /usr/bin/python3 has neither numpy nor a governor that tolerates the wait.
# VIA_JOB=1 moves both passes to a compute node, where activate_conda has one
# that does. Same qsub shape as record_261026_phylop_crossbins.sh uses.
# The bed selection above has already been validated here, so a typo in ONLY or
# EXCEPT still fails now rather than inside the job.
if [[ "${VIA_JOB:-0}" == "1" ]]; then
  [[ "${SUBMIT:-0}" == "1" ]] || {
    say "VIA_JOB=1 submits the real run, so it needs SUBMIT=1."; exit 2; }
  mkdir -p "$logdir"
  jid=$(qsub -terse -N pxsubmit \
    -cwd -j y -V \
    -o "${logdir}/fig4b_submit.\$JOB_ID.log" \
    -l h_rt=2:00:00,h_data=8G \
    -v "REPO_DIR=${PWD},SUMMARY_RUN=${SUMMARY_RUN},PHYLOP_SFS_RT=${SGE_SFS_RT},ONLY=${ONLY:-},EXCEPT=${EXCEPT:-},PROBE=${PROBE:-0}" \
    records/sge_phylop_fig4b_submit.sh | cut -d. -f1)
  say
  say "1b. launch   ${jid}   runs this script with SUBMIT=1 on a compute node"
  say "    ${#beds[@]} bed(s): sfs array ${#beds[@]} x 22 = $(( ${#beds[@]} * 22 )) tasks"
  say "    log  ${logdir}/fig4b_submit.${jid}.log"
  say "    the sfs jobs appear in qstat once its two preflights are done, ~15 min"
  say
  qstat -u "$USER" 2>&1 | tee -a "$log"
  say "logged to ${log}"
  exit 0
fi

# ------------------------------------------------------------- 2. dry run --
say
say "2. dry run   ($(( ${#beds[@]} )) bed(s): sfs array ${#beds[@]} x 22 = $(( ${#beds[@]} * 22 )) tasks)"
DRY_RUN=1 FROM=sfs bash scripts/chains/run_all.sh "${beds[@]}" 2>&1 | tail -16 | tee -a "$log"

# -------------------------------------------------------------- 3. submit --
if [[ "${SUBMIT:-0}" != "1" ]]; then
  say
  say "Nothing was submitted."
  say "   PROBE=1 SUBMIT=1 bash records/record_261026_phylop_fig4b_run.sh   one bed"
  say "   SUBMIT=1 bash records/record_261026_phylop_fig4b_run.sh           all 40"
  exit 0
fi

say
say "3. submitting"
FROM=sfs bash scripts/chains/run_all.sh "${beds[@]}" 2>&1 | tee -a "$log"
say
qstat -u "$USER" 2>&1 | tee -a "$log"
say
say "logged to ${log}"
say "the table will be \$OUT/summary/${SUMMARY_RUN}/bestpdf_sumdfe.tsv"
