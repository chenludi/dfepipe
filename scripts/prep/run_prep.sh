#!/bin/bash
# prep: filter a bed of candidate functional regions down to the sites the SFS
# will be built from.
#
#   bash scripts/prep/run_prep.sh <bed> [bed...]
#
# The filter chain is ported verbatim from the two scripts that produced the
# published results, in the order they run them:
#
#   submit_exHETquiencmask.sh:35  bedtools subtract  HET_BED                (heterochromatin)
#   submit_exHETquiencmask.sh:37  bedtools subtract  CANDIDATE_NEUTRAL_BED  (quiescent)
#   filt_nc-unmask.sh:38          bedtools subtract  CDS_BED                (coding)
#   filt_nc-unmask.sh:46-47       bedtools intersect MASK_BED, then sort
#
# The quiescent subtraction is not optional bookkeeping: those sites become the
# neutral CONTROL set in the sfs stage, so a test site must never also be one.
# The output name keeps the historical convention s.nc_exHET_noq.<base>.bed,
# which is what filt_nc-unmask.sh produced ("s.nc_" prefixed onto the name of
# its input, itself already "exHET_noq.<base>.bed").
#
# Unlike the originals this also merges the masked result, so overlapping input
# intervals cannot double-count downstream, and it records base pairs after every
# step in lengths.tsv so a surprising result localises to one filter.
#
# Results are keyed by PREP_RUN: output goes to $OUT/prep/<base>/$PREP_RUN/ and
# each run writes spec.txt recording the config it used. An existing non-empty
# run directory is never overwritten implicitly; use a new PREP_RUN, or FORCE=1.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# shellcheck source=../../config/env.sh
source "${REPO_DIR}/config/env.sh"

: "${OUT:?OUT not set}"
: "${PREP_RUN:?PREP_RUN not set}"
: "${HET_BED:?HET_BED not set}"
: "${CANDIDATE_NEUTRAL_BED:?CANDIDATE_NEUTRAL_BED not set}"
: "${CDS_BED:?CDS_BED not set}"
: "${MASK_BED:?MASK_BED not set}"
: "${BEDTOOLS_MOD:?BEDTOOLS_MOD not set}"

if [[ $# -eq 0 ]]; then
  echo "usage: run_prep.sh <bed> [bed...]" >&2
  exit 2
fi

# bp of a bed, the same count the originals used (count_lengthbedfile.sh), with
# an empty file counting 0 rather than printing nothing.
bp() { awk '{s += ($3 > $2) ? $3 - $2 : $2 - $3} END {printf "%d\n", s+0}' "$1"; }

# --------------------------------------------------------------- preflight ----
# Everything that can be checked before any work is checked first, so a typo in
# one bed of twenty does not surface after an hour of bedtools.
for f in "$HET_BED" "$CANDIDATE_NEUTRAL_BED" "$CDS_BED" "$MASK_BED"; do
  [[ -s "$f" ]] || { echo "prep: reference bed missing or empty: $f" >&2; exit 1; }
done
for bed in "$@"; do
  [[ -s "$bed" ]] || { echo "prep: input bed missing or empty: $bed" >&2; exit 1; }
done
for bed in "$@"; do
  base="$(basename "$bed" .bed)"
  d="${OUT}/prep/${base}/${PREP_RUN}"
  if [[ -d "$d" && -n "$(ls -A "$d" 2>/dev/null)" ]]; then
    if [[ "${FORCE:-0}" == "1" ]]; then
      echo "prep: FORCE=1 - overwriting the existing run in ${d}" >&2
    else
      echo "prep: run directory exists and is not empty: ${d}" >&2
      echo "      Results are keyed by PREP_RUN, so an existing run is never" >&2
      echo "      overwritten implicitly. Pick a new PREP_RUN, or re-run with" >&2
      echo "      FORCE=1 to replace this one." >&2
      exit 1
    fi
  fi
done

load_modules "$BEDTOOLS_MOD"

# ------------------------------------------------------------------- filter ----
for bed in "$@"; do
  base="$(basename "$bed" .bed)"
  d="${OUT}/prep/${base}/${PREP_RUN}"
  out_bed="${d}/s.nc_exHET_noq.${base}.bed"
  mkdir -p "$d"

  tmp="$(mktemp -d "${d}/.prep.XXXXXX")"
  # shellcheck disable=SC2064  # $tmp is expanded now on purpose, per iteration
  trap "rm -rf '${tmp}'" EXIT

  lengths="${tmp}/lengths.tsv"
  printf 'stage\tbp\n' > "$lengths"
  printf 'input\t%s\n' "$(bp "$bed")" >> "$lengths"

  echo "prep ${base}: 1/4 subtract heterochromatin"
  bedtools subtract -a "$bed" -b "$HET_BED" > "${tmp}/exHET.bed"
  printf 'minus_HET\t%s\n' "$(bp "${tmp}/exHET.bed")" >> "$lengths"

  echo "prep ${base}: 2/4 subtract quiescent (the control set)"
  bedtools subtract -a "${tmp}/exHET.bed" -b "$CANDIDATE_NEUTRAL_BED" \
    > "${tmp}/exHET_noq.bed"
  printf 'minus_quiescent\t%s\n' "$(bp "${tmp}/exHET_noq.bed")" >> "$lengths"

  echo "prep ${base}: 3/4 subtract CDS"
  bedtools subtract -a "${tmp}/exHET_noq.bed" -b "$CDS_BED" > "${tmp}/nc.bed"
  printf 'minus_CDS\t%s\n' "$(bp "${tmp}/nc.bed")" >> "$lengths"

  echo "prep ${base}: 4/4 intersect the accessibility mask, sort, merge"
  bedtools intersect -a "${tmp}/nc.bed" -b "$MASK_BED" \
    | sort -T "$tmp" -k1,1 -k2,2n \
    | bedtools merge -i - > "${tmp}/final.bed"
  printf 'after_mask\t%s\n' "$(bp "${tmp}/final.bed")" >> "$lengths"

  # Move into place only once every step has succeeded, so an interrupted run
  # cannot leave a half-filtered bed looking like a finished one.
  mv "${tmp}/final.bed" "$out_bed"
  mv "$lengths" "${d}/lengths.tsv"

  { echo "PREP_RUN=${PREP_RUN}"
    echo "input_bed=${bed}"
    echo "output_bed=${out_bed}"
    echo "HET_BED=${HET_BED}"
    echo "CANDIDATE_NEUTRAL_BED=${CANDIDATE_NEUTRAL_BED}"
    echo "CDS_BED=${CDS_BED}"
    echo "MASK_BED=${MASK_BED}"
    echo "filter_order=subtract_HET,subtract_quiescent,subtract_CDS,intersect_mask+sort+merge"
    echo "BEDTOOLS_MOD=${BEDTOOLS_MOD}"
    echo "bedtools_version=$(bedtools --version 2>/dev/null || echo unknown)"
    echo "host=$(hostname)"
    echo "date=$(date -Iseconds 2>/dev/null || date)"
  } > "${d}/spec.txt"

  rm -rf "$tmp"
  trap - EXIT
  echo "prep ${base}: wrote ${out_bed}"
  cat "${d}/lengths.tsv"
done

# ---------------------------------------------------------------- checkpoint ----
# Diagnostic only: a bug in the check must never discard a finished stage, so its
# failure warns and the stage's exit status is left alone.
if ! python3 "${REPO_DIR}/scripts/checks/check_prep.py" \
       --out "${OUT}/prep" --run "${PREP_RUN}" "$@"; then
  echo "WARNING: check_prep.py failed; the prep output above is unaffected." >&2
  echo "         Inspect lengths.tsv by hand before trusting this run." >&2
fi

exit 0
