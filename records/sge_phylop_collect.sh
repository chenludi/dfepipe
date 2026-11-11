#!/bin/bash
# Final step of record_261026_phylop_crossbins.sh: gather the per-bed length
# files the union array wrote into one table, and flag anything suspicious.
#
# Each union task writes its own one-line file rather than appending to a shared
# one, because 40 concurrent appends over NFS is not safe. This job runs after
# the array and concatenates them.
#
# SMALL_BP is an attention threshold, not a rule. The PhastCons version of this
# analysis needed its sparsest cells pooled across three bins before a DFE could
# be fitted (record_1103_pc17nonconserve_pc470.sh:8), so a bin that comes out far
# smaller than its neighbours is worth looking at before 40 DFE chains are run.
set -euo pipefail
: "${OUTDIR:?OUTDIR must come from the driver}"
: "${UNION_TASKS:?UNION_TASKS must come from the driver}"
SMALL_BP="${SMALL_BP:-2000000}"

out="${OUTDIR}/lengths_261026.tsv"
printf 'bed\tseries\tbin\tn_cells\tn_intervals\tbp\n' > "$out"
cat "${OUTDIR}"/lengths/*.tsv 2>/dev/null | sort -k2,2 -k3,3V >> "$out"

n=$(( $(wc -l < "$out" | tr -d " ") - 1 ))
want=$(wc -l < "$UNION_TASKS" | tr -d " ")
echo "collected ${n} of ${want} beds into ${out}"
echo
column -t -s$'\t' "$out"
echo
if (( n != want )); then
  echo "warning: ${want} beds were expected; $(( want - n )) union task(s) did not report."
  echo "Check ${OUTDIR}/log/union.*.log before using this table."
fi
awk -F'\t' -v lim="$SMALL_BP" 'NR>1 && $6+0 < lim {print "   " $1 "   " $6 " bp"}' "$out" \
  > "${OUTDIR}/lengths/small.txt" || true
if [[ -s "${OUTDIR}/lengths/small.txt" ]]; then
  echo "beds under ${SMALL_BP} bp, worth a look before running the DFEs:"
  cat "${OUTDIR}/lengths/small.txt"
else
  echo "no bed is under ${SMALL_BP} bp"
fi
