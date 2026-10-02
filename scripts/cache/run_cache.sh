#!/bin/bash
# cache: build the selection-coefficient spectra cache from the demographic fit,
# one cache per bed.
#
#   bash scripts/cache/run_cache.sh <bed> [bed...]
#
# <bed> is the ORIGINAL bed, the same argument prep, sfs and the demo stage took.
# Run it after scripts/demo/run_bestdemo.sh: the gamma bounds are 2*s*Na and Na
# comes from that stage's best_inference.txt.
#
# Ported from submit_dadi-cli_cache_simp.sh (which fit, which Na, the output
# name) and dadi-cli_cache_simp.sh (the bounds at :79-81, the grid at :55, the
# gamma_pts default at :59, and the GenerateCache call at :117-118). The two
# non-default marker rows in config/cache_markers.tsv come from a third original,
# pipefs_exinitials_6_cache.sh:128-137:
#
#   max_ll_key   from best_inference.txt, then *.<max_ll_key>.InferDM.bestfits
#                beside it is the winning fit passed as --demo-popt
#   Na           from the same file
#   gamma bounds 2 * GAMMA_S_LOW * Na  to  2 * GAMMA_S_HIGH * Na  (1e-6, 0.5),
#                formatted with awk "%.10f" exactly as the original
#   --model      <dadi model>_sel, e.g. three_epoch_sel
#
# Measured 17 min on 10 cores for phastCons17way 0-5 at the default marker
# (grids 1500 1520 1540, 300 gamma points). The original asked for h_rt=10:10:00
# (submit_dadi-cli_cache_simp.sh:8), which was a queue request, not a measurement;
# an earlier version of this header repeated it as "a ~10 hour job". An existing
# cache is skipped rather than rebuilt; FORCE=1 rebuilds it.
#
# CACHE_MARKER names the parameter set and appears both in the run directory and
# in the output file name. It is resolved ONCE, before any bed, against
# CACHE_MARKER_TABLE (config/cache_markers.tsv), whose columns are
#
#   marker  gamma_s_low  gamma_s_high  gamma_lo  gamma_hi  gamma_pts  grids  note
#
# A marker sets EITHER gamma_s_low/gamma_s_high, and the bounds are 2*s*Na per
# bed (marker 1 does this), OR gamma_lo/gamma_hi, and they are absolute and the
# same for every bed. The unused pair is NA; setting neither is an error. The
# row's gamma_pts and grids override the config values when they are not NA.
#
# Differences from the original, all deliberate:
#   * `--sample-sizes` is spelled in full. The original wrote `--sample-size`,
#     which worked only as an argparse abbreviation (dadi-cli __main__.py's
#     add_sample_sizes_argument), the same slip run_demo.sh records for
#     `--output-prefix` / `--cpus`.
#   * the dead positive-gamma branch is gone. In dadi-cli_cache_simp.sh:113-114
#     its dadi-cli call is commented out, so only the else branch ever ran.
#   * command_cache.sh holds the command that ACTUALLY RAN, and is a runnable
#     script. The original always wrote the --additional-gammas variant
#     (:124), i.e. not the command it had just executed - and with the array
#     empty that line is not even valid.
#   * the cache lands in its own $OUT/cache/<bed>/... directory rather than on
#     top of the demo run's, so a second marker or run cannot overwrite a first.
#   * a failed GenerateCache has its output removed, so a half-written pickle
#     cannot be mistaken for a finished cache on the next run's skip check. A
#     job killed outright (h_rt) can still leave one; check_cache.py catches it.
#   * DRYRUN=1 prints and records the command without running it, for checking
#     what a submission would do before spending the ten hours.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# shellcheck source=../../config/env.sh
source "${REPO_DIR}/config/env.sh"

: "${OUT:?OUT not set}"
: "${FS_RUN:?FS_RUN not set: the cache stage reads the demo output}"
: "${DEMO_RUN:?DEMO_RUN not set}"
: "${DEMO_MODEL:?DEMO_MODEL not set}"
: "${CACHE_RUN:?CACHE_RUN not set}"
: "${CACHE_MARKER:?CACHE_MARKER not set}"
: "${CACHE_MARKER_TABLE:?CACHE_MARKER_TABLE not set}"
: "${GAMMA_PTS:?GAMMA_PTS not set}"
: "${GRIDS:?GRIDS not set}"
: "${GAMMA_S_LOW:?GAMMA_S_LOW not set}"
: "${GAMMA_S_HIGH:?GAMMA_S_HIGH not set}"
: "${PROJECTION:?PROJECTION not set}"
: "${SGE_CACHE_PE:?SGE_CACHE_PE not set: the core count GenerateCache is given}"
: "${CACHE_MIN_SIZE_FRAC:?CACHE_MIN_SIZE_FRAC not set}"
: "${CACHE_MIN_USABLE_FRAC:?CACHE_MIN_USABLE_FRAC not set}"

if [[ $# -eq 0 ]]; then
  echo "usage: run_cache.sh <bed> [bed...]" >&2
  exit 2
fi

# dadi-cli's model name, as run_demo.sh maps it; the cache needs the selection
# version, which dadi.DFE.DemogSelModels spells <name>_sel.
case "$DEMO_MODEL" in
  three_epoch)  dadi_model="three_epoch" ;;
  two_epoch)    dadi_model="two_epoch" ;;
  growth)       dadi_model="growth" ;;
  bottlegrowth) dadi_model="bottlegrowth_1d" ;;
  *) echo "run_cache: DEMO_MODEL='${DEMO_MODEL}' is not supported" >&2; exit 2 ;;
esac
model_sel="${dadi_model}_sel"

# ----------------------------------------------------------------- the marker --
[[ -s "$CACHE_MARKER_TABLE" ]] || {
  echo "run_cache: CACHE_MARKER_TABLE missing or empty: ${CACHE_MARKER_TABLE}" >&2
  exit 1; }
row="$(awk -F'\t' -v m="$CACHE_MARKER" 'NR > 1 && $1 == m {print; exit}' \
       "$CACHE_MARKER_TABLE")"
if [[ -z "$row" ]]; then
  echo "run_cache: CACHE_MARKER='${CACHE_MARKER}' is not in ${CACHE_MARKER_TABLE}" >&2
  echo "      Known markers:" \
       "$(awk -F'\t' 'NR > 1 && NF {printf "%s ", $1}' "$CACHE_MARKER_TABLE")" >&2
  exit 1
fi
IFS=$'\t' read -r _m row_s_low row_s_high row_gamma_lo row_gamma_hi \
        row_gamma_pts row_grids row_note <<< "$row"

# A marker either scales the bounds by Na (gamma_s_low/high) or gives them
# absolutely (gamma_lo/hi). Exactly one pair must be set.
if [[ "$row_s_low" != "NA" && "$row_s_high" != "NA" ]]; then
  GAMMA_S_LOW="$row_s_low"; GAMMA_S_HIGH="$row_s_high"
  abs_lo=""; abs_hi=""
elif [[ "$row_gamma_lo" != "NA" && "$row_gamma_hi" != "NA" ]]; then
  abs_lo="$row_gamma_lo"; abs_hi="$row_gamma_hi"
else
  echo "run_cache: marker '${CACHE_MARKER}' in ${CACHE_MARKER_TABLE} sets" >&2
  echo "      neither gamma_s_low/gamma_s_high (bounds as 2*s*Na) nor" >&2
  echo "      gamma_lo/gamma_hi (absolute bounds). It must set one pair;" >&2
  echo "      the row reads: ${row}" >&2
  exit 1
fi
[[ "$row_gamma_pts" != "NA" ]] && GAMMA_PTS="$row_gamma_pts"
[[ "$row_grids"     != "NA" ]] && GRIDS="$row_grids"
# From here on the resolved values are used under their lowercase names; the
# uppercase ones are what the checkpoint is told about, so the two cannot drift.
gamma_pts="$GAMMA_PTS"
grids="$GRIDS"

# ------------------------------------------------------------------ preflight --
# Everything resolvable is resolved before any GenerateCache starts: a typo in
# the second bed must not surface ten hours in.
run_dir_name="${DEMO_RUN}-${CACHE_MARKER}-${CACHE_RUN}"
bases=(); demo_dirs=(); bestfits=(); nas=(); cache_dirs=()
for bed in "$@"; do
  base="$(basename "$bed" .bed)"
  d="${OUT}/demo/${base}/${FS_RUN}-${DEMO_MODEL}-${DEMO_RUN}"
  best="${d}/best_inference.txt"
  [[ -s "$best" ]] || {
    echo "run_cache: no best_inference.txt for '${base}': ${best}" >&2
    echo "      Run scripts/demo/run_bestdemo.sh first, or check FS_RUN /" >&2
    echo "      DEMO_MODEL / DEMO_RUN." >&2
    exit 1; }

  # The two keys the original grepped for (submit_dadi-cli_cache_simp.sh:71-73).
  key="$(awk -F: '/^max_ll_key:/ {print $2; exit}' "$best")"
  na="$(awk -F: '/^Na:/ {print $2; exit}' "$best")"
  [[ -n "$key" ]] || { echo "run_cache: ${best} has no max_ll_key:" >&2; exit 1; }
  [[ -n "$na" ]] || { echo "run_cache: ${best} has no Na:" >&2; exit 1; }
  if ! awk -v n="$na" 'BEGIN {exit !(n + 0 > 0)}'; then
    echo "run_cache: ${base}: Na='${na}' is not positive, so 2*s*Na is meaningless" >&2
    exit 1
  fi

  fit="$(ls -1 "${d}"/*."${key}".InferDM.bestfits 2>/dev/null | head -1 || true)"
  [[ -n "$fit" ]] || {
    echo "run_cache: ${base}: no *.${key}.InferDM.bestfits in ${d}" >&2
    echo "      best_inference.txt names max_ll_key=${key}; the fit it refers" >&2
    echo "      to is gone. Re-run the demo stage." >&2
    exit 1; }

  bases+=("$base"); demo_dirs+=("$d"); bestfits+=("$fit"); nas+=("$na")
  cache_dirs+=("${OUT}/cache/${base}/${run_dir_name}")
done

[[ "${DRYRUN:-0}" == "1" ]] || activate_conda

# ---------------------------------------------------------------- the caches --
failed=0
for i in "${!bases[@]}"; do
  base="${bases[$i]}"
  na="${nas[$i]}"
  fit="${bestfits[$i]}"
  d="${cache_dirs[$i]}"
  out_bpkl="${d}/cache_${CACHE_MARKER}.spectra.bpkl"
  mkdir -p "$d"

  if [[ -s "$out_bpkl" && "${FORCE:-0}" != "1" ]]; then
    echo "run_cache: ${base}: cache already present, skipping: ${out_bpkl}"
    echo "      FORCE=1 rebuilds it."
    continue
  fi

  # Bounds: absolute if the marker gave them, else 2*s*Na from this bed's Na,
  # with awk "%.10f" exactly as dadi-cli_cache_simp.sh:79-81.
  if [[ -n "$abs_lo" ]]; then
    gamma_lo="$abs_lo"; gamma_hi="$abs_hi"
  else
    gamma_lo="$(awk "BEGIN {printf \"%.10f\", 2 * ${GAMMA_S_LOW} * ${na}}")"
    gamma_hi="$(awk "BEGIN {printf \"%.10f\", 2 * ${GAMMA_S_HIGH} * ${na}}")"
  fi
  # dadi-cli requires both bounds strictly positive (_check_positive_num), and
  # "%.10f" rounds a very small lower bound to 0.0000000000.
  if ! awk -v l="$gamma_lo" -v h="$gamma_hi" 'BEGIN {exit !(l + 0 > 0 && h + 0 > l + 0)}'; then
    echo "run_cache: ${base}: gamma bounds came out as ${gamma_lo} ${gamma_hi}," \
         "which dadi-cli will reject. Raise GAMMA_S_LOW." >&2
    failed=1
    continue
  fi

  { echo "CACHE_RUN=${CACHE_RUN}"
    echo "CACHE_MARKER=${CACHE_MARKER}"
    echo "FS_RUN=${FS_RUN}"
    echo "DEMO_RUN=${DEMO_RUN}"
    echo "DEMO_MODEL=${DEMO_MODEL}"
    echo "model_sel=${model_sel}"
    echo "best_inference=${demo_dirs[$i]}/best_inference.txt"
    echo "demo_popt=${fit}"
    echo "Na=${na}"
    echo "gamma_lo=${gamma_lo}"
    echo "gamma_hi=${gamma_hi}"
    echo "gamma_pts=${gamma_pts}"
    echo "grids=${grids}"
    echo "projection=${PROJECTION}"
    echo "cpus=${SGE_CACHE_PE}"
    echo "marker_table=${CACHE_MARKER_TABLE}"
    echo "marker_note=${row_note}"
    echo "host=$(hostname)"
    echo "date=$(date -Iseconds 2>/dev/null || date)"
  } > "${d}/spec.txt"

  # $GRIDS is three words and $gamma_lo/$gamma_hi two: both are meant to split.
  # shellcheck disable=SC2206
  cmd=(dadi-cli GenerateCache
       --model "$model_sel"
       --demo-popt "$fit"
       --sample-sizes "$PROJECTION"
       --grids ${grids}
       --gamma-pts "$gamma_pts"
       --gamma-bounds "$gamma_lo" "$gamma_hi"
       --output "$out_bpkl"
       --cpus "$SGE_CACHE_PE")

  # The command that actually runs, beside the output, runnable as it stands.
  { echo "#!/bin/bash"
    echo "# The dadi-cli GenerateCache command this run executed."
    echo "# ${base}, marker ${CACHE_MARKER}, $(date -Iseconds 2>/dev/null || date)"
    printf '%q ' "${cmd[@]}"
    echo
  } > "${d}/command_cache.sh"

  echo "run_cache: ${base}: Na=${na} gamma ${gamma_lo}..${gamma_hi}" \
       "gamma_pts=${gamma_pts} grids='${grids}'"
  printf '%s\n' "${cmd[*]}"
  if [[ "${DRYRUN:-0}" == "1" ]]; then
    echo "run_cache: DRYRUN=1, nothing was run; see ${d}/command_cache.sh"
    continue
  fi

  if "${cmd[@]}" 2>&1 | tee "${d}/cache.log"; then
    echo "run_cache: ${base}: wrote ${out_bpkl}" \
         "($(wc -c < "$out_bpkl" | tr -d ' ') bytes)"
  else
    echo "run_cache: ${base}: dadi-cli GenerateCache failed; see ${d}/cache.log" >&2
    # Whatever it left behind is not a cache, and must not look like one to the
    # skip check above on the next run.
    rm -f "$out_bpkl"
    failed=1
  fi
done

# ---------------------------------------------------------------- checkpoint --
# Diagnostic only: a bug in the check must never discard a ten-hour cache, so its
# failure warns and the stage's exit status is left alone.
if [[ "${DRYRUN:-0}" != "1" ]]; then
  if ! python3 "${REPO_DIR}/scripts/checks/check_cache.py" \
         --out "${OUT}/cache" --run "${run_dir_name}" \
         --marker "${CACHE_MARKER}" \
         --min-size-frac "${CACHE_MIN_SIZE_FRAC}" \
         --gamma-pts "${GAMMA_PTS}" --projection "${PROJECTION}" \
         --min-usable "${CACHE_MIN_USABLE_FRAC}" \
         --gamma-s-high "${GAMMA_S_HIGH}" "$@"; then
    echo "WARNING: check_cache.py failed; the cache above is unaffected." >&2
    echo "         Inspect the .bpkl by hand before trusting this run." >&2
  fi
fi

exit "$failed"
