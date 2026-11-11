#!/bin/bash
# Array worker for record_261026_phylop_crossbins.sh, step "cells".
# One task per pairwise intersect of a primate constraint bin with a mammal one.
# All 336 cells are built here, into STAGE; nothing existing is reused.
#
# Not run by hand. UGE spools a copy of this script, so $BASH_SOURCE points at
# the spool copy; everything it needs arrives through qsub -v.
#
# -sorted is not optional: the per-bin beds are 0.8-1.7 GB of mostly
# single-base intervals, and plain `bedtools intersect` holds -b in memory.
# They come out of filt_nc-unmask.sh already sorted with `sort -k1,1 -k2,2n`,
# and bedtools fails loudly if that is not true, so a wrong assumption cannot
# pass silently. The records and their order are the same either way, which
# matters because the compare step hashes these files.
#
# The bed is written under a temporary name and moved into place only once
# complete, so a task killed partway cannot leave a truncated file that a later
# step would treat as finished.
set -euo pipefail
: "${REPO_DIR:?REPO_DIR must come from the driver}"
: "${CELLS_TASKS:?CELLS_TASKS must come from the driver}"
: "${STAGE:?STAGE must come from the driver}"
: "${SGE_TASK_ID:?sge_phylop_cell.sh runs as an array task}"

# shellcheck source=../config/env.sh
source "${REPO_DIR}/config/env.sh"
load_modules "$BEDTOOLS_MOD"

IFS=$'\t' read -r a b < <(sed -n "${SGE_TASK_ID}p" "$CELLS_TASKS")
[[ -n "${a:-}" && -n "${b:-}" ]] || {
  echo "no task on line ${SGE_TASK_ID} of ${CELLS_TASKS}" >&2; exit 1; }

name="s.nc_exHET_noq.${PRIMATE}_${a}_${MAMMAL}_${b}.bed"
file_a="${BINDIR_P}/s.nc_exHET_noq.${PRIMATE}_${a}.bed"
file_b="${BINDIR_M}/s.nc_exHET_noq.${MAMMAL}_${b}.bed"
for f in "$file_a" "$file_b"; do
  [[ -s "$f" ]] || { echo "missing input: $f" >&2; exit 1; }
done

echo "sge_phylop_cell: host=$(hostname) job=${JOB_ID:-none} task=${SGE_TASK_ID}" \
     "a=${a} b=${b} started $(date)"

out="${STAGE}/${name}"
tmp="${out}.part.${JOB_ID:-0}.${SGE_TASK_ID}"
bedtools intersect -sorted -a "$file_a" -b "$file_b" > "$tmp"
mv -f "$tmp" "$out"

n_int=$(wc -l < "$out")
n_bp=$(awk '{s += ($3 > $2) ? $3 - $2 : $2 - $3} END {printf "%d", s+0}' "$out")
printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$a" "$b" "$n_int" "$n_bp" > "${out}.track"
echo "sge_phylop_cell: wrote ${out}  ${n_int} intervals  ${n_bp} bp  $(date)"
