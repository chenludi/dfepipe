#!/bin/bash
# Delete the dfepipe test-run output on hoffman2, and record what was deleted.
#
# Runs on the cluster login node. Nothing is removed unless DELETE=1 is set, so
# the default is a listing you can read before committing to anything.
#
#   bash record_261025_cleanup_testruns.sh              list the dev tree
#   bash record_261025_cleanup_testruns.sh all          list dev + runs 1 and 2
#   DELETE=1 bash record_261025_cleanup_testruns.sh     delete the dev tree
#   DELETE=1 bash record_261025_cleanup_testruns.sh all delete both
#
# Scopes:
#   dev   output_dfepipe.dev_pre261001, the development tree. 8.0 G, of which
#         7.75 G is scratch left by one killed ELS sfs run (fs_260929nm): the
#         EXIT trap in run_sfs_chrom.sh does not fire on SIGKILL. No result in
#         it was ever used; dfe_260929 holds fits made against two different
#         caches, which is why CACHE_MARKER is now part of the directory name.
#   all   also output_dfepipe, the two validated acceptance runs. 543 M. Their
#         tables, DFE csv and figures are already under
#         a_results_weekly/mv_261002_dfepipe_run{1,2}/ locally, and the script
#         will not delete it unless --i-have-the-tables is passed as well.
#
# Never touched: data_bed, result_dadi, output_fs, dfepipe itself.
set -euo pipefail

ROOT="${ROOT:-$HOME/project-klohmuel/noncdoingdfe/test_0118}"
scope="${1:-dev}"
confirm_tables="${2:-}"
log="${ROOT}/cleanup_261025.log"

[[ -d "$ROOT" ]] || { echo "not a directory: $ROOT" >&2; exit 1; }

# Only these may ever be removed. A target not on this list is a bug, not a
# scope the caller forgot to add.
case "$scope" in
  dev) targets=(output_dfepipe.dev_pre261001) ;;
  all) targets=(output_dfepipe.dev_pre261001 output_dfepipe) ;;
  *)   echo "scope must be dev or all, not '${scope}'" >&2; exit 2 ;;
esac

if [[ "$scope" == "all" && "$confirm_tables" != "--i-have-the-tables" ]]; then
  echo "scope 'all' also deletes output_dfepipe, the two acceptance runs." >&2
  echo "Check that a_results_weekly/mv_261002_dfepipe_run{1,2}/ hold the" >&2
  echo "tables and figures you want, then pass --i-have-the-tables." >&2
  exit 2
fi

protected=(data_bed result_dadi output_fs dfepipe)
for p in "${protected[@]}"; do
  for t in "${targets[@]}"; do
    [[ "$t" == "$p" ]] && { echo "refusing to delete ${p}" >&2; exit 1; }
  done
done

mode="listing only"
[[ "${DELETE:-0}" == "1" ]] && mode="deleting"
echo "=============================================================="
echo " dfepipe test-run cleanup   scope=${scope}   ${mode}"
echo " root: ${ROOT}"
echo " date: $(date '+%Y-%m-%d %H:%M:%S')"
echo "=============================================================="
echo

total_kb=0
present=()
for t in "${targets[@]}"; do
  path="${ROOT}/${t}"
  if [[ ! -d "$path" ]]; then
    echo "  absent, nothing to do:  ${path}"
    continue
  fi
  present+=("$path")
  kb=$(du -sk "$path" | cut -f1)
  total_kb=$((total_kb + kb))
  echo "  ${path}"
  echo "      size:  $(du -sh "$path" | cut -f1)"
  echo "      files: $(find "$path" -type f | wc -l | tr -d ' ')"
  echo "      subdirectories to be removed:"
  # one line per second-level directory, which is the level a reader recognises
  find "$path" -mindepth 1 -maxdepth 2 -type d \
    | sort \
    | while read -r d; do printf '        %-62s %s\n' "${d#$ROOT/}" "$(du -sh "$d" | cut -f1)"; done
  echo
done

if [[ ${#present[@]} -eq 0 ]]; then
  echo "Nothing to delete."
  exit 0
fi

printf 'Total to remove: %s MB in %d top-level director%s\n\n' \
  "$((total_kb / 1024))" "${#present[@]}" \
  "$([[ ${#present[@]} -eq 1 ]] && echo y || echo ies)"

if [[ "${DELETE:-0}" != "1" ]]; then
  echo "Listing only. Nothing was removed."
  echo "To delete:  DELETE=1 bash record_261025_cleanup_testruns.sh ${scope} ${confirm_tables}"
  exit 0
fi

{
  echo "## $(date '+%Y-%m-%d %H:%M:%S')  scope=${scope}  host=$(hostname)"
  for path in "${present[@]}"; do
    echo "deleted ${path}  ($(du -sh "$path" | cut -f1), $(find "$path" -type f | wc -l | tr -d ' ') files)"
  done
} >> "$log"

for path in "${present[@]}"; do
  echo "deleting ${path} ..."
  rm -rf "$path"
  echo "  done; still present? $([[ -d "$path" ]] && echo yes || echo no)"
done

echo
echo "Recorded in ${log}:"
tail -n $(( ${#present[@]} + 1 )) "$log"
