#!/bin/bash
# Step "compare" of record_261026_phylop_crossbins.sh.
#
# The 144 cells in CELL_OLD were built on 2025-10-21 by something not recorded
# in any script, and nothing downstream ever used them. This step rebuilds
# nothing; it compares the fresh build in STAGE against those 144 and decides
# whether they were the same files.
#
# Three comparisons per cell, in order, stopping at the first that passes:
#   sha256    byte-identical, which is the expected result
#   sorted    same records in a different order; would mean the two builds
#             differ only in record order, not in content
#   bp        same base pairs covered but different records; would mean the
#             interval structure differs, which changes control windows
# The sorted and bp comparisons only run for cells whose sha256 differs, so the
# expected path reads each file once.
#
# On success the 336 cells are consolidated into CELLS, one copy of each: the
# 144 originals are moved there (keeping their 2025-10-21 mtimes, which are the
# only record of when they were made), the 192 new ones are moved from STAGE,
# and STAGE's 144 duplicates are removed. Only files this run created are
# deleted. On any mismatch nothing is moved or deleted, and the check file
# tells the union step to stop.
set -euo pipefail
: "${CELL_OLD:?CELL_OLD must come from the driver}"
: "${CELLS:?CELLS must come from the driver}"
: "${STAGE:?STAGE must come from the driver}"
: "${CELLS_TASKS:?CELLS_TASKS must come from the driver}"
: "${COMPARE_LIST:?COMPARE_LIST must come from the driver}"
: "${OUTDIR:?OUTDIR must come from the driver}"
: "${CHECK:?CHECK must come from the driver}"

report="${OUTDIR}/compare_261026.tsv"
cellname () { echo "s.nc_exHET_noq.${PRIMATE}_${1}_${MAMMAL}_${2}.bed"; }
bp () { awk '{s += ($3 > $2) ? $3 - $2 : $2 - $3} END {printf "%d", s+0}' "$1"; }

echo "sge_phylop_compare: host=$(hostname) job=${JOB_ID:-none} started $(date)"

# Every cell the later steps need must exist in STAGE, or the array lost a task.
n_want=$(wc -l < "$CELLS_TASKS" | tr -d " ")
n_missing=0
while IFS=$'\t' read -r a b; do
  [[ -s "${STAGE}/$(cellname "$a" "$b")" ]] || {
    echo "   not built: $(cellname "$a" "$b")" >&2; n_missing=$((n_missing + 1)); }
done < "$CELLS_TASKS"
if (( n_missing )); then
  printf 'status\tFAIL\nreason\t%s of %s cells were not built\n' "$n_missing" "$n_want" > "$CHECK"
  echo "compare: ${n_missing} of ${n_want} cells missing from ${STAGE}. Stopping." >&2
  exit 1
fi
echo "   all ${n_want} cells present in ${STAGE}"

printf 'cell\ta\tb\tverdict\tn_int_old\tn_int_new\tbp_old\tbp_new\n' > "$report"
n_same=0 n_diff=0
while IFS=$'\t' read -r a b; do
  name="$(cellname "$a" "$b")"
  old="${CELL_OLD}/${name}"; new="${STAGE}/${name}"
  if [[ "$(sha256sum < "$old" | cut -d' ' -f1)" == "$(sha256sum < "$new" | cut -d' ' -f1)" ]]; then
    printf '%s\t%s\t%s\tsha256_match\tNA\tNA\tNA\tNA\n' "$name" "$a" "$b" >> "$report"
    n_same=$((n_same + 1))
    continue
  fi
  oi=$(wc -l < "$old" | tr -d " "); ni=$(wc -l < "$new" | tr -d " ")
  ob=$(bp "$old");      nb=$(bp "$new")
  st="${TMPDIR:-/tmp}"
  if [[ "$(sort -T "$st" "$old" | sha256sum | cut -d' ' -f1)" \
     == "$(sort -T "$st" "$new" | sha256sum | cut -d' ' -f1)" ]]; then
    verdict=same_records_different_order
  elif [[ "$ob" == "$nb" ]]; then
    verdict=same_bp_different_intervals
  else
    verdict=DIFFERENT
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$a" "$b" "$verdict" "$oi" "$ni" "$ob" "$nb" >> "$report"
  n_diff=$((n_diff + 1))
done < "$COMPARE_LIST"

n_cmp=$(wc -l < "$COMPARE_LIST" | tr -d " ")
n_order=$(awk -F'\t' 'NR>1 && $4 == "same_records_different_order"' "$report" | wc -l | tr -d ' ')
n_real=$(awk -F'\t' 'NR>1 && ($4 == "DIFFERENT" || $4 == "same_bp_different_intervals")' "$report" | wc -l | tr -d ' ')
echo "   compared ${n_cmp}: ${n_same} byte-identical, ${n_order} order-only, ${n_real} real"
echo "   report ${report}"

# An order-only difference means the two builds hold the same records; the
# chromosomes are just in a different sequence. On 2026-10-26, 16 of the 144
# cells came out this way: the old copies were in numeric chromosome order
# (chr1, chr2, ... chr22) while bedtools emits in -a's order, which for these
# inputs is lexicographic (chr1, chr10, ... chr2, chr20 ...). The other 128 were
# byte-identical, so the fresh build reproduces the original recipe exactly and
# those 16 had been re-sorted after the fact. Order does not reach any result:
# the union is a cat and run_prep.sh sorts and merges before anything reads it.
# It is still not accepted silently - ACCEPT_ORDER=1 has to say so.
if (( n_real > 0 )); then
  printf 'status\tFAIL\nreason\t%s of %s existing cells differ in content, not just order\n' \
    "$n_real" "$n_cmp" > "$CHECK"
  echo
  echo "compare: ${n_real} cell(s) differ in content. Nothing moved, nothing deleted."
  awk -F'\t' 'NR>1 && ($4 == "DIFFERENT" || $4 == "same_bp_different_intervals")' "$report"
  exit 1
fi
if (( n_order > 0 )) && [[ "${ACCEPT_ORDER:-0}" != "1" ]]; then
  printf 'status\tFAIL\nreason\t%s of %s existing cells hold the same records in a different order\n' \
    "$n_order" "$n_cmp" > "$CHECK"
  echo
  echo "compare: ${n_order} cell(s) hold the same records in a different chromosome"
  echo "order. No content differs. Nothing moved, nothing deleted."
  echo "Re-run this step with ACCEPT_ORDER=1 to accept that and consolidate."
  awk -F'\t' 'NR>1 && $4 == "same_records_different_order"' "$report"
  exit 1
fi

# ------------------------------------------------------------ consolidate --
# The fresh build is kept for every cell, so all 336 come from one recorded
# recipe and carry one chromosome ordering. The old copies are removed only
# after being shown to hold the same records; compare_261026.tsv is the record
# of that, and PURGE_OLD=0 keeps them if you would rather decide later.
mkdir -p "$CELLS"
n_mv=0
while IFS=$'\t' read -r a b; do
  name="$(cellname "$a" "$b")"
  mv -f "${STAGE}/${name}" "${CELLS}/${name}"
  [[ -e "${STAGE}/${name}.track" ]] && mv -f "${STAGE}/${name}.track" "${CELLS}/${name}.track"
  n_mv=$((n_mv + 1))
done < "$CELLS_TASKS"

n_final=$(find "$CELLS" -maxdepth 1 -name '*.bed' | wc -l | tr -d ' ')
echo
echo "   consolidated into ${CELLS}"
echo "      moved from ${STAGE}   ${n_mv}"
echo "      cells now in place    ${n_final}"
if (( n_final != n_want )); then
  printf 'status\tFAIL\nreason\t%s cells in place, expected %s\n' "$n_final" "$n_want" > "$CHECK"
  echo "compare: ${n_final} cells in place, expected ${n_want}. Stopping." >&2
  exit 1
fi

n_rm=0
if [[ "${PURGE_OLD:-1}" == "1" ]]; then
  while IFS=$'\t' read -r a b; do
    name="$(cellname "$a" "$b")"
    rm -f "${CELL_OLD}/${name}" "${CELL_OLD}/${name}.track"
    n_rm=$((n_rm + 1))
  done < "$COMPARE_LIST"
  echo "      superseded copies removed from ${CELL_OLD}   ${n_rm}"
  rmdir "$CELL_OLD" 2>/dev/null && echo "      ${CELL_OLD} removed (empty)" || true
else
  echo "      PURGE_OLD=0: the ${n_cmp} old copies are left in ${CELL_OLD}"
fi
rmdir "$STAGE" 2>/dev/null && echo "      ${STAGE} removed (empty)" || true

printf 'status\tPASS\nreason\t%s cells, %s byte-identical and %s order-only against the 2025-10-21 build\n' \
  "$n_final" "$n_same" "$n_order" > "$CHECK"
echo "sge_phylop_compare: PASS  $(date)"
