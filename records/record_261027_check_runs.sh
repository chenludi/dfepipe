#!/bin/bash
# 2026-10-27  Did the three finished runs finish correctly?
#
# Read-only. Nothing here writes into a run directory; it only reads the
# check_*.tsv files each stage leaves behind and the two summary tables.
#
#   els_old  three old-ENCODE ELS beds        *_261025
#   els_v4   four ENCODE v4 beds              *_261026
#   phylop   one conserved x non-conserved    *_pp261026   (the Fig. 4B probe)
#
# Per stage it checks three things: the check file exists, every expected item
# has exactly one row, and no row says FAIL. It then adds the assertions that a
# status column cannot make on its own - no errored, broken or unrun dfe fits,
# and one summary row per bed in both final tables.
#
# Two WARN families are expected and are not failures:
#   check_bed   chrX/chrY are not in CHROMS and get dropped at concat, and some
#               beds are not coordinate-sorted (sfs sorts them anyway).
#   check_dfe neu  goodness-of-fit flags on the neutral fit. The neu pdf has no
#               starting-point grid, so its count columns are NA by design.
#
#   bash records/record_261027_check_runs.sh             all three
#   RUNS="phylop" bash records/record_261027_check_runs.sh
#   VERBOSE=1 bash records/record_261027_check_runs.sh   list PASS rows too
#
# Exit status is 0 only if all requested runs are clean.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
source config/profiles/hoffman_replicate.conf
source config/pipeline.conf          # PDFS, DEMO_MODEL, CHROMS
# env.sh is deliberately not sourced: no conda is needed to read tsv files.

RUNS="${RUNS:-els_old els_v4 phylop}"
VERBOSE="${VERBOSE:-0}"

# run key sets and bed lists. prep="-" means the run deliberately had no prep
# stage: the phylop beds entered at sfs, already filtered, via FROM=sfs.
run_spec() {
  case "$1" in
    els_old) echo "prep_261025|fs_261025|demo_261025|cache_261025|dfe_261025|sum_261025|GRCh38-ELS GRCh38-ELS_x_Prom GRCh38-ELS_not_Prom" ;;
    els_v4)  echo "prep_261026|fs_261026|demo_261026|cache_261026|dfe_261026|sum_261026|GRCh38-PLSv4 GRCh38-ELSv4 GRCh38-ELSv4_x_Prom GRCh38-ELSv4_not_Prom" ;;
    phylop)  echo "-|fs_pp261026|demo_pp261026|cache_pp261026|dfe_pp261026|sum_pp261026|phyloP243way_0-5_phyloP447way_40-100" ;;
    *) echo "" ;;
  esac
}

run_fail=0; run_warn=0        # reset per run
tot_fail=0

say()  { printf '%s\n' "$*"; }
bad()  { printf '  FAIL  %-28s %s\n' "$1" "$2"; run_fail=$((run_fail+1)); }
warn() { printf '  warn  %-28s %s\n' "$1" "$2"; run_warn=$((run_warn+1)); }
ok()   { printf '  ok    %-28s %s\n' "$1" "$2"; }

# scan_table <file> <stage label> <expected items, space separated>
# Prints one line per non-PASS / missing / duplicate row, then a #COUNT line.
scan_table() {
  local f="$1" label="$2" expected="$3" out line st item rsn
  if [ ! -f "$f" ]; then
    bad "$label" "check file missing: ${f#"$OUT"/}"
    return
  fi
  out=$(awk -F'\t' -v expected="$expected" '
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i
            rc = ("reason" in h) ? h["reason"] : 0
            next }
    {
      item=$1; st=$2; rsn = (rc>0 && rc<=NF) ? $rc : "NA"
      seen[item]++; n[st]++
      if (st == "PASS" || st == "SKIP") { if (v=="1") printf "%s\t%s\t%s\n", st, item, rsn }
      else printf "%s\t%s\t%s\n", st, item, rsn
    }
    END {
      ne = split(expected, e, " "); extra=0
      for (i=1;i<=ne;i++) if (!(e[i] in seen)) { printf "MISSING\t%s\t-\n", e[i]; extra++ }
      for (k in seen) if (seen[k] > 1) { printf "DUPLICATE\t%s\t%d rows\n", k, seen[k]; extra++ }
      printf "#COUNT\t%d\t%d\t%d\t%d\t%d\n", n["PASS"]+0, n["WARN"]+0, n["FAIL"]+0, extra, n["SKIP"]+0
    }' v="$VERBOSE" "$f")

  local np nw nf nx ns cnt skipnote
  cnt=$(printf '%s\n' "$out" | grep '^#COUNT' | head -1 || true)
  read -r _ np nw nf nx ns <<<"$cnt"
  np=${np:-0}; nw=${nw:-0}; nf=${nf:-0}; nx=${nx:-0}; ns=${ns:-0}
  skipnote=""; [ "$ns" -gt 0 ] && skipnote=", $ns skipped (no intervals there)"

  if [ "$nf" -gt 0 ] || [ "$nx" -gt 0 ]; then
    bad "$label" "$np pass, $nw warn, $nf fail, $nx missing or duplicated$skipnote"
  elif [ "$nw" -gt 0 ]; then
    warn "$label" "$np pass, $nw warn$skipnote"
  else
    ok "$label" "$np rows, all pass$skipnote"
  fi

  while IFS=$'\t' read -r st item rsn; do
    [ -n "${st:-}" ] || continue
    case "$st" in PASS|SKIP) printf '          %-50s %s\n' "$item" "$(printf '%s' "$st" | tr 'A-Z' 'a-z')"; continue ;; esac
    printf '        %-7s %-50s %s\n' "$st" "$item" "$(printf '%s' "$rsn" | cut -c1-150)"
  done <<<"$(printf '%s\n' "$out" | grep -v '^#COUNT' || true)"
}

# numeric column assertions the status flag does not cover
dfe_counts() {                       # <file> <label>
  local f="$1" label="$2" out
  [ -f "$f" ] || return 0
  out=$(awk -F'\t' '
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
    {
      e=$h["n_error"]+0; b=$h["n_broken"]+0; nr=$h["n_not_run"]+0; t=$h["n_timeout"]+0
      sp=$h["n_starting_points"]; bf=$h["n_bestfit"]
      if (e>0||b>0||nr>0) printf "BAD\t%s\terror=%d broken=%d not_run=%d\n", $1, e, b, nr
      if (t>0)            printf "TO\t%s\t%d fit(s) hit the time limit\n", $1, t
      if (sp ~ /^[0-9]+$/) printf "SP\t%s\t%s of %s starting points gave a bestfit\n", $1, bf, sp
    }' "$f")
  while IFS=$'\t' read -r k item msg; do
    case "${k:-}" in
      BAD) bad "$label" "$item: $msg" ;;
      TO)  warn "$label" "$item: $msg" ;;
      SP)  [ "$VERBOSE" = "1" ] && printf '          %-50s %s\n' "$item" "$msg" || true ;;
    esac
  done <<<"$out"
}

for rk in $RUNS; do
  spec=$(run_spec "$rk")
  if [ -z "$spec" ]; then say "unknown run key: $rk"; tot_fail=$((tot_fail+1)); continue; fi
  IFS='|' read -r PREP_RUN FS_RUN DEMO_RUN CACHE_RUN DFE_RUN SUMMARY_RUN BEDS <<<"$spec"
  run_fail=0; run_warn=0
  nbed=$(printf '%s\n' $BEDS | wc -l | tr -d ' ')

  say ""
  say "=============================================================================="
  say "$rk   $nbed bed(s)   keys: ${PREP_RUN} ${FS_RUN} ${DEMO_RUN} ${CACHE_RUN} ${DFE_RUN} ${SUMMARY_RUN}"
  say "=============================================================================="

  # prep
  if [ "$PREP_RUN" = "-" ]; then
    ok "prep" "skipped by design, this run entered at sfs"
  else
    exp=""; for b in $BEDS; do exp="$exp ${b}.bed"; done
    scan_table "$OUT/prep/check_bed.${PREP_RUN}.tsv"  "prep: check_bed"  "$exp"
    scan_table "$OUT/prep/check_prep.${PREP_RUN}.tsv" "prep: check_prep" "$BEDS"
  fi

  # sfs: one row per bed per chromosome, plus a per-bed total row
  exp=""
  for b in $BEDS; do
    exp="$exp $b"
    for c in $CHROMS; do
      case "$c" in chr*) cn="$c" ;; *) cn="chr$c" ;; esac   # check_sfs.py normalise_chrom
      exp="$exp ${b}:${cn}"
    done
  done
  scan_table "$OUT/sfs/check_sfs.${FS_RUN}.tsv" "sfs" "$exp"

  # demo
  scan_table "$OUT/demo/check_demo.${FS_RUN}-${DEMO_MODEL}-${DEMO_RUN}.tsv" "demo" "$BEDS"

  # cache: one check file per bed
  miss=0
  for b in $BEDS; do
    cf="$OUT/cache/${b}/${DEMO_RUN}-1-${CACHE_RUN}/check_cache.tsv"
    [ -f "$cf" ] || { bad "cache" "missing for $b"; miss=$((miss+1)); continue; }
    scan_table "$cf" "cache: $b" "$b"
  done

  # dfe, one file per pdf
  for pdf in $PDFS; do
    df="$OUT/dfe/check_dfe.${CACHE_RUN}-1-${pdf}-${DFE_RUN}.tsv"
    scan_table "$df" "dfe: $pdf" "$BEDS"
    [ "$pdf" = "neu" ] || dfe_counts "$df" "dfe: $pdf"
  done

  # summary check file plus the two tables
  exp=""; for b in $BEDS; do exp="$exp exHET_noq.${b}_qui2"; done
  scan_table "$OUT/summary/check_summary.${SUMMARY_RUN}.tsv" "summary" "$exp"

  for t in sumdfe.tsv bestpdf_sumdfe.tsv; do
    tf="$OUT/summary/${SUMMARY_RUN}/$t"
    if [ ! -f "$tf" ]; then bad "summary: $t" "missing"; continue; fi
    rows=$(( $(wc -l < "$tf") - 1 ))
    gone=""
    for b in $BEDS; do
      grep -q "^exHET_noq.${b}_qui2	" "$tf" || gone="$gone $b"
    done
    if [ -n "$gone" ]; then bad "summary: $t" "no row for:$gone"
    elif [ "$rows" -ne "$nbed" ]; then bad "summary: $t" "$rows rows, expected $nbed"
    else ok "summary: $t" "$rows rows, one per bed"; fi
  done

  say ""
  if [ "$run_fail" -gt 0 ]; then
    say ">> $rk: FAIL  ($run_fail problem(s), $run_warn warning(s))"
    tot_fail=$((tot_fail+1))
  else
    say ">> $rk: finished correctly  ($run_warn expected warning(s))"
  fi
done

say ""
say "=============================================================================="
if [ "$tot_fail" -eq 0 ]; then
  say "all requested runs finished correctly"
else
  say "$tot_fail run(s) have problems, see FAIL lines above"
fi
exit "$tot_fail"
