#!/bin/bash
# Array worker for record_261026_phylop_crossbins.sh, step "union".
# One task per output bed: cat the 12 not-conserved cells of one constraint bin.
#
#   series A   <primate>_<bin>_<mammal>_40-100     bin fixed on the primate side
#   series B   <primate>_40-100_<mammal>_<bin>     bin fixed on the mammal side
#
# Gated on the compare step's check file, the same way dfepipe's stage workers
# gate on theirs: -hold_jid releases a job whatever the held job's exit status
# was, so without this a failed comparison would still be followed by 40 unions
# built from cells whose provenance is in doubt.
#
# A missing cell is fatal rather than skipped: a quietly short union would be a
# wrong DFE input, not a smaller one.
#
# cat only, no sort and no merge, which is what
# record_1103_pc17nonconserve_pc470.sh did for the PhastCons series. The cells
# are individually sorted and disjoint, so the union holds the right set of
# bases; dfepipe prep sorts and merges it like any other bed. Merging here
# instead would join bases that are adjacent across two different bins, moving
# interval midpoints and so moving the +/-25 kb control windows.
set -euo pipefail
: "${CELLS:?CELLS must come from the driver}"
: "${UNION_TASKS:?UNION_TASKS must come from the driver}"
: "${CHECK:?CHECK must come from the driver}"
: "${OUTDIR:?OUTDIR must come from the driver}"
: "${SGE_TASK_ID:?sge_phylop_union.sh runs as an array task}"

if [[ ! -s "$CHECK" ]]; then
  echo "union: no check file at ${CHECK}; the compare step did not finish. Stopping." >&2
  exit 1
fi
status=$(awk -F'\t' '$1 == "status" {print $2; exit}' "$CHECK")
if [[ "$status" != "PASS" ]]; then
  reason=$(awk -F'\t' '$1 == "reason" {print $2; exit}' "$CHECK")
  echo "union: compare says ${status:-unknown}: ${reason:-no reason given}. Stopping." >&2
  exit 1
fi

IFS=$'\t' read -r series bin < <(sed -n "${SGE_TASK_ID}p" "$UNION_TASKS")
[[ -n "${series:-}" && -n "${bin:-}" ]] || {
  echo "no task on line ${SGE_TASK_ID} of ${UNION_TASKS}" >&2; exit 1; }
read -r -a ncon <<< "$NCON_LIST"

case "$series" in
  A) name="s.nc_exHET_noq.${PRIMATE}_${bin}_${MAMMAL}_${NCON_LABEL}.bed" ;;
  B) name="s.nc_exHET_noq.${PRIMATE}_${NCON_LABEL}_${MAMMAL}_${bin}.bed" ;;
  *) echo "unknown series '${series}'" >&2; exit 1 ;;
esac

cells=()
for n in "${ncon[@]}"; do
  if [[ "$series" == "A" ]]; then
    cn="s.nc_exHET_noq.${PRIMATE}_${bin}_${MAMMAL}_${n}.bed"
  else
    cn="s.nc_exHET_noq.${PRIMATE}_${n}_${MAMMAL}_${bin}.bed"
  fi
  [[ -s "${CELLS}/${cn}" ]] || { echo "missing cell ${CELLS}/${cn}" >&2; exit 1; }
  cells+=("${CELLS}/${cn}")
done

echo "sge_phylop_union: host=$(hostname) job=${JOB_ID:-none} task=${SGE_TASK_ID}" \
     "series=${series} bin=${bin} cells=${#cells[@]} started $(date)"

out="${OUTDIR}/${name}"
tmp="${out}.part.${JOB_ID:-0}.${SGE_TASK_ID}"
cat "${cells[@]}" > "$tmp"
mv -f "$tmp" "$out"

n_int=$(wc -l < "$out" | tr -d " ")
n_bp=$(awk '{s += ($3 > $2) ? $3 - $2 : $2 - $3} END {printf "%d", s+0}' "$out")
mkdir -p "${OUTDIR}/lengths"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$series" "$bin" "${#cells[@]}" "$n_int" "$n_bp" \
  > "${OUTDIR}/lengths/${name}.tsv"
echo "sge_phylop_union: wrote ${out}  ${n_int} intervals  ${n_bp} bp  $(date)"
