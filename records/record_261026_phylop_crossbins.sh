#!/bin/bash
# 2026-10-26  Cross-binned phyloP beds, for the PhyloP version of Fig. 4B.
#
# Fig. 4B has four series; two are plain constraint bins and two are bins of one
# score restricted to sites not conserved by the other. This builds the beds for
# the two restricted series, the same way
# ../../record_1103_pc17nonconserve_pc470.sh built pc17_40-100_pc470_<bin>:
# intersect the two scores' bins pairwise, then cat the cells.
#
#   primates, not conserved in mammals   <primate>_<bin>_<mammal>_40-100.bed
#   mammals,  not conserved in primates  <primate>_40-100_<mammal>_<bin>.bed
#
# "not conserved" is the bottom 60% of the other score, bins 40-45 through
# 95-100. That matches Fig. 4B, where the restriction is the PhastCons 470way
# zero-score block (its bottom 60%) and pc17_40-100. It does not match the
# phyloP243way_0-5_only sets behind Fig. 4D, which use "not in the top 5%".
#
#   bash records/record_261026_phylop_crossbins.sh              inventory only
#   SUBMIT=1 bash records/record_261026_phylop_crossbins.sh     inventory, then submit
#
# Four jobs, each held on the one before:
#   cells     array, one task per cell; all 336 are built fresh, into a staging
#             directory. Nothing is reused.
#   compare   the 144 cells that already exist in non_conserve/ are compared
#             against the fresh build. On a match the 336 are consolidated into
#             one directory, one copy each; on any mismatch it stops and the
#             union step refuses to run.
#   union     array, one task per output bed, 40 of them.
#   collect   gathers the per-bed length table.
#
# Why all 336 rather than reusing the 144: those cells are all dated
# 2025-10-21, no script in script/ or scripts/ names them, and nothing
# downstream ever consumed them, so their only provenance was their filenames.
# Rebuilding costs 144 extra array tasks and turns the question into a check.
#
# Overridable: PRIMATE, MAMMAL (score names), TC (array concurrency).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
source config/profiles/hoffman_replicate.conf
source config/pipeline.conf
source config/env.sh

PRIMATE="${PRIMATE:-phyloP243way}"
MAMMAL="${MAMMAL:-phyloP447way}"
TC="${TC:-40}"

BINS=(0-5 5-10 10-15 15-20 20-25 25-30 30-35 35-40 40-45 45-50 50-55 55-60 \
      60-65 65-70 70-75 75-80 80-85 85-90 90-95 95-100)
NCON=(40-45 45-50 50-55 55-60 60-65 65-70 70-75 75-80 80-85 85-90 90-95 95-100)
NCON_LABEL=40-100

bins_dir="$REFS_DIR/data_bed/conserve_bins_wg_nc"
bindir_p="${bins_dir}/${PRIMATE}"
bindir_m="${bins_dir}/${MAMMAL}"
cell_old="${bins_dir}/phyloP243-phyloP447/non_conserve"
cells="${bins_dir}/phyloP243-phyloP447/cells"
stage="$REFS_DIR/data_bed/conserve_bins_wg_ovlp/${PRIMATE}_x_${MAMMAL}_261026"
outdir="$REFS_DIR/data_bed/conserve_phylop_ncon_261026"
logdir="${outdir}/log"
log="${outdir}/build_261026.log"
check="${outdir}/cells_check.tsv"

mkdir -p "$stage" "$outdir" "$logdir"
say () { echo "$@" | tee -a "$log"; }
cellname () { echo "s.nc_exHET_noq.${PRIMATE}_${1}_${MAMMAL}_${2}.bed"; }

say "=============================================================="
say " phyloP cross-binned beds   $(date '+%Y-%m-%d %H:%M:%S')   host=$(hostname)"
say " primate score ${PRIMATE}   mammal score ${MAMMAL}"
say " not conserved = bins ${NCON[*]}  (${NCON_LABEL}, the bottom 60%)"
say "=============================================================="

# ------------------------------------------------------------- 1. the bins --
say
say "1. per-bin input beds"
missing_bins=0
for b in "${BINS[@]}"; do
  for d in "${bindir_p}:${PRIMATE}" "${bindir_m}:${MAMMAL}"; do
    f="${d%%:*}/s.nc_exHET_noq.${d##*:}_${b}.bed"
    [[ -s "$f" ]] || { say "   missing: $f"; missing_bins=$((missing_bins + 1)); }
  done
done
if (( missing_bins )); then
  say "   ${missing_bins} per-bin bed(s) missing. Stopping: build those first."
  exit 1
fi
say "   all $(( 2 * ${#BINS[@]} )) present"

# ------------------------------------------------------------ 2. the cells --
# Series A needs every (bin, ncon) pair, series B every (ncon, bin) pair. What
# the two requirements share is ncon x ncon, the 144 that already exist, so the
# two lists are concatenated and deduplicated rather than counted twice.
say
say "2. pairwise cells"
cells_tasks="${logdir}/tasks_cells.tsv"
compare_list="${logdir}/cells_to_compare.tsv"
{
  for a in "${BINS[@]}"; do for b in "${NCON[@]}"; do printf '%s\t%s\n' "$a" "$b"; done; done
  for b in "${BINS[@]}"; do for a in "${NCON[@]}"; do printf '%s\t%s\n' "$a" "$b"; done; done
} | sort -u > "$cells_tasks"

: > "$compare_list"
n_old=0
while IFS=$'\t' read -r a b; do
  [[ -s "${cell_old}/$(cellname "$a" "$b")" ]] || continue
  printf '%s\t%s\n' "$a" "$b" >> "$compare_list"
  n_old=$((n_old + 1))
done < "$cells_tasks"

n_cell_tasks=$(wc -l < "$cells_tasks" | tr -d " ")
# A finished run leaves all 336 in cells/ and removes the staging directory, so
# a second run would otherwise rebuild every cell and compare against a CELL_OLD
# that is now empty, which would pass for the wrong reason.
# find exits non-zero on a directory that does not exist, which under
# `set -e` with pipefail would abort here on the first run, so it is only
# called once the directory is known to be there.
n_done=0
if [[ -d "$cells" ]]; then
  n_done=$(find "$cells" -maxdepth 1 -name '*.bed' | wc -l | tr -d ' ')
fi
if (( n_done >= n_cell_tasks )) && [[ "${FORCE:-0}" != "1" ]]; then
  say "   ${cells} already holds ${n_done} cells, which is the full set."
  say "   Nothing to build. Use FORCE=1 to rebuild them anyway."
  say "   The union step alone can be rerun from those cells."
  exit 0
fi
say "   to build      ${n_cell_tasks}   all of them, into ${stage}"
say "   of those      ${n_old}   also exist in ${cell_old}, and will be compared"
say "   not reused    nothing; every cell is rebuilt from the per-bin beds"

# ------------------------------------------------------ 3. the output beds --
say
say "3. output beds"
union_tasks="${logdir}/tasks_union.tsv"
: > "$union_tasks"
for a in "${BINS[@]}"; do printf 'A\t%s\n' "$a" >> "$union_tasks"; done
for b in "${BINS[@]}"; do printf 'B\t%s\n' "$b" >> "$union_tasks"; done
n_union=$(wc -l < "$union_tasks" | tr -d " ")
say "   series A  ${PRIMATE}_<bin>_${MAMMAL}_${NCON_LABEL}   ${#BINS[@]} beds"
say "   series B  ${PRIMATE}_${NCON_LABEL}_${MAMMAL}_<bin>   ${#BINS[@]} beds"
say "   from      ${cells}   (after compare consolidates into it)"
say "   into      ${outdir}"
say "   each is a cat of ${#NCON[@]} cells, unsorted and unmerged, as the"
say "   PhastCons recipe left them; dfepipe prep sorts and merges."

# -------------------------------------------------------------- 4. submit --
if [[ "${SUBMIT:-0}" != "1" ]]; then
  say
  say "Nothing was submitted. To submit:"
  say "   SUBMIT=1 bash records/record_261026_phylop_crossbins.sh"
  exit 0
fi

ensure_scheduler
# Explicit name=value pairs, the convention the dfepipe submitters use, so a
# task's inputs cannot depend on what happened to be exported. NCON_LIST holds
# spaces; qsub -v preserves them (checked on 2026-10-25 with DEMO_P0_FILES).
job_env="REPO_DIR=${PWD},PRIMATE=${PRIMATE},MAMMAL=${MAMMAL}"
job_env="${job_env},NCON_LABEL=${NCON_LABEL},NCON_LIST=${NCON[*]}"
job_env="${job_env},CELL_OLD=${cell_old},CELLS=${cells},STAGE=${stage}"
job_env="${job_env},OUTDIR=${outdir},CHECK=${check}"
job_env="${job_env},BINDIR_P=${bindir_p},BINDIR_M=${bindir_m}"
job_env="${job_env},CELLS_TASKS=${cells_tasks},UNION_TASKS=${union_tasks}"
job_env="${job_env},COMPARE_LIST=${compare_list},BEDTOOLS_MOD=${BEDTOOLS_MOD}"

rm -f "$check"      # a stale PASS from an earlier attempt must not gate this one
say
say "4. submitting"

# FROM=compare resubmits the comparison onwards without rebuilding the cells,
# for the case where the cells are already staged and only the comparison's
# verdict needs revisiting (see ACCEPT_ORDER in sge_phylop_compare.sh).
if [[ "${FROM:-}" == "compare" ]]; then
  n_staged=0
  if [[ -d "$stage" ]]; then
    n_staged=$(find "$stage" -maxdepth 1 -name '*.bed' | wc -l | tr -d ' ')
  fi
  if (( n_staged != n_cell_tasks )); then
    say "   FROM=compare needs all ${n_cell_tasks} cells staged; ${stage} holds ${n_staged}."
    exit 1
  fi
  say "   cells     skipped, ${n_staged} already staged"
  jid_cmp=$(qsub -terse -N pxcompare \
    -cwd -j y -V \
    -o "${logdir}/compare.\$JOB_ID.log" \
    -l h_rt=4:00:00,h_data=4G \
    -v "$job_env,ACCEPT_ORDER=${ACCEPT_ORDER:-0},PURGE_OLD=${PURGE_OLD:-1}" \
    records/sge_phylop_compare.sh | cut -d. -f1)
  say "   compare   ${jid_cmp}      ACCEPT_ORDER=${ACCEPT_ORDER:-0} PURGE_OLD=${PURGE_OLD:-1}"
else
jid_cells=$(qsub -terse -N pxcell \
  -t "1-${n_cell_tasks}" -tc "$TC" \
  -cwd -j y -V \
  -o "${logdir}/cells.\$JOB_ID.\$TASK_ID.log" \
  -l h_rt=3:00:00,h_data=4G \
  -v "$job_env" records/sge_phylop_cell.sh | cut -d. -f1)
say "   cells     ${jid_cells}   1-${n_cell_tasks}   at most ${TC} at once"

jid_cmp=$(qsub -terse -N pxcompare \
  -cwd -j y -V \
  -o "${logdir}/compare.\$JOB_ID.log" \
  -l h_rt=4:00:00,h_data=4G \
  -hold_jid "${jid_cells}" \
  -v "$job_env,ACCEPT_ORDER=${ACCEPT_ORDER:-0},PURGE_OLD=${PURGE_OLD:-1}" \
  records/sge_phylop_compare.sh | cut -d. -f1)
say "   compare   ${jid_cmp}      ${n_old} cells, then consolidates"
fi

jid_union=$(qsub -terse -N pxunion \
  -t "1-${n_union}" -tc "$TC" \
  -cwd -j y -V \
  -o "${logdir}/union.\$JOB_ID.\$TASK_ID.log" \
  -l h_rt=1:00:00,h_data=2G \
  -hold_jid "${jid_cmp}" -v "$job_env" records/sge_phylop_union.sh | cut -d. -f1)
say "   union     ${jid_union}   1-${n_union}   gated on ${check}"

jid_collect=$(qsub -terse -N pxcollect \
  -cwd -j y -V \
  -o "${logdir}/collect.\$JOB_ID.log" \
  -l h_rt=0:30:00,h_data=2G \
  -hold_jid "${jid_union}" -v "$job_env" records/sge_phylop_collect.sh | cut -d. -f1)
say "   collect   ${jid_collect}"
say
qstat -u "$USER" 2>&1 | tee -a "$log"
say
say "logged to ${log}"
say "the comparison lands in ${outdir}/compare_261026.tsv"
say "the per-bed table lands in ${outdir}/lengths_261026.tsv"
