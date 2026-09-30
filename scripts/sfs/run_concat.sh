#!/bin/bash
# sfs, step 7: concatenate the per-chromosome spectra and fsmut tables into the
# genome-wide files, then run the checkpoint.
#
#   bash scripts/sfs/run_concat.sh <bed> [bed...]
#
# <bed> is the ORIGINAL bed, the same argument prep and run_sfs_chrom.sh took.
# Ported from piperef_stddfe_2_concatefs.sh; the arithmetic stays in the two
# legacy R scripts, which are called, not rewritten:
#
#   concatefs.R     adds the spectra bin by bin, after checking that every file
#                   shares the same header and mask lines (same projection, same
#                   population) - concatefs.R:57
#   concate_mut.R   sums the site and length columns and recomputes the
#                   length-weighted average mutation rate - concate_mut.R:27,71
#
# Differences from the original, all deliberate:
#   * `stat -c %s` is GNU-only; `wc -c` is portable and is what the size sanity
#     check actually needs. The thresholds (661 and 500 BYTES - the original's
#     message says KB, but the numbers are bytes) are unchanged.
#   * a missing chromosome is NAMED. The original demanded exactly 22 files and
#     printed only "There are not 22 chromosomes", which is how chr4 went
#     missing from one bin unnoticed (record_250624_243-40-45.sh).
#   * the per-chromosome files are left in place under bychr/ instead of being
#     moved to an archive directory: results are keyed by FS_RUN, so nothing has
#     to be moved out of the way to keep runs apart.
#   * R writes into a temp directory and the result is renamed into place, so an
#     interrupted job cannot leave a truncated .fs under its final name.
#
# The genome-wide names are the historical ones and are what every later stage
# and every published result keys on:
#   has_qui.<shortbase>_<CTRL_TAG>.chrs_1t22.fs
#   quiecent_ctrl.<shortbase>_<CTRL_TAG>.chrs_1t22.fs
#   fsmut.has_qui.<shortbase>_<CTRL_TAG>.chrs_1t22.tab
#   fsmut.quiecent_ctrl.<shortbase>_<CTRL_TAG>.chrs_1t22.tab
# "chrs_1t22" is a literal, as in the original, and does not track CHROMS.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# shellcheck source=../../config/env.sh
source "${REPO_DIR}/config/env.sh"

: "${OUT:?OUT not set}"
: "${REFS_DIR:?REFS_DIR not set}"
: "${FS_RUN:?FS_RUN not set}"
: "${PREP_RUN:?PREP_RUN not set}"
: "${CHROMS:?CHROMS not set}"
: "${CTRL_TAG:?CTRL_TAG not set}"
: "${MIN_SNPS:?MIN_SNPS not set}"
: "${R_MOD:?R_MOD not set}"

# The two legacy R scripts live with the reference data, beside the rest of the
# original pipeline. Override LEGACY_SCRIPT_DIR if they sit elsewhere.
LEGACY_SCRIPT_DIR="${LEGACY_SCRIPT_DIR:-${REFS_DIR}/script}"
CONCAT_FS_R="${LEGACY_SCRIPT_DIR}/concatefs.R"
CONCAT_MUT_R="${LEGACY_SCRIPT_DIR}/concate_mut.R"

: "${MIN_FS_BYTES:?MIN_FS_BYTES not set}"
: "${MIN_FSMUT_BYTES:?MIN_FSMUT_BYTES not set}"

if [[ $# -eq 0 ]]; then
  echo "usage: run_concat.sh <bed> [bed...]" >&2
  exit 2
fi

for f in "$CONCAT_FS_R" "$CONCAT_MUT_R"; do
  [[ -s "$f" ]] || {
    echo "concat: legacy R script not found: $f" >&2
    echo "        Set LEGACY_SCRIPT_DIR to the directory holding concatefs.R" >&2
    echo "        and concate_mut.R." >&2
    exit 1; }
done

# CHROMS is bare numbers in pipeline.conf; everything else is chr-prefixed.
chroms=()
for c in $CHROMS; do
  case "$c" in chr*) chroms+=("$c") ;; *) chroms+=("chr${c}") ;; esac
done

load_modules "$R_MOD"

failed=0

for bed in "$@"; do
  base="$(basename "$bed" .bed)"
  bedstem="s.nc_exHET_noq.${base}"
  # The genome-wide name drops prep's prefix, exactly as
  # piperef_stddfe_2_concatefs.sh:12-14 did: first "s.nc_", then "nc_". Both are
  # plain first-occurrence substitutions, quirks included, because the published
  # file names were produced by these two lines.
  shortbase="${bedstem/s.nc_/}"
  shortbase="${shortbase/nc_/}"
  outbase="${shortbase}_${CTRL_TAG}"

  d="${OUT}/sfs/${base}/${FS_RUN}"
  bychr="${d}/bychr"
  if [[ ! -d "$bychr" ]]; then
    echo "concat: no per-chromosome output for '${base}': ${bychr}" >&2
    echo "        Run scripts/sfs/run_sfs_chrom.sh first, or check FS_RUN." >&2
    failed=1
    continue
  fi

  for set_name in has_qui quiecent_ctrl; do
    out_fs="${d}/${set_name}.${outbase}.chrs_1t22.fs"
    out_tab="${d}/fsmut.${set_name}.${outbase}.chrs_1t22.tab"

    if [[ -s "$out_fs" || -s "$out_tab" ]]; then
      if [[ "${FORCE:-0}" == "1" ]]; then
        echo "concat: FORCE=1 - replacing ${set_name} for ${base}" >&2
      else
        echo "concat: ${base} ${set_name} is already concatenated: ${out_fs}" >&2
        echo "        Results are keyed by FS_RUN, so an existing run is never" >&2
        echo "        overwritten implicitly. Pick a new FS_RUN, or FORCE=1." >&2
        failed=1
        continue
      fi
    fi

    fs_files=(); tab_files=(); missing_fs=(); missing_tab=()
    for chrom in "${chroms[@]}"; do
      # The .fs name ends in the VCF's own name, so it is globbed rather than
      # rebuilt from VCF_TMPL: the file on disk is the authority.
      f="$(ls -1 "${bychr}/smer.${set_name}.${chrom}.${bedstem}."*.fs 2>/dev/null | head -1 || true)"
      t="${bychr}/fsmut.smer.${set_name}.${chrom}.${bedstem}.${chrom}.tab"

      if [[ -n "$f" && -f "$f" && "$(wc -c < "$f")" -gt "$MIN_FS_BYTES" ]]; then
        fs_files+=("$f")
      else
        missing_fs+=("$chrom")
        if [[ -n "$f" && -f "$f" ]]; then
          echo "concat: ${base} ${set_name} ${chrom}: .fs is only $(wc -c < "$f") bytes" \
               "(under ${MIN_FS_BYTES}); it holds no counts" >&2
        fi
      fi

      if [[ -f "$t" && "$(wc -c < "$t")" -gt "$MIN_FSMUT_BYTES" ]]; then
        tab_files+=("$t")
      else
        missing_tab+=("$chrom")
      fi
    done

    if [[ "${#missing_fs[@]}" -ne 0 || "${#missing_tab[@]}" -ne 0 ]]; then
      echo "concat: ${base} ${set_name}: cannot concatenate, the genome-wide file" >&2
      echo "        would silently be missing chromosomes." >&2
      if [[ "${#missing_fs[@]}" -ne 0 ]]; then
        echo "        no usable .fs for:    ${missing_fs[*]}" >&2
      fi
      if [[ "${#missing_tab[@]}" -ne 0 ]]; then
        echo "        no usable fsmut for:  ${missing_tab[*]}" >&2
      fi
      echo "        check_sfs below says, per chromosome, whether that is expected." >&2
      failed=1
      continue
    fi

    tmp="$(mktemp -d "${d}/.concat.XXXXXX")"
    echo "concat: ${base} ${set_name}: ${#fs_files[@]} chromosomes"
    Rscript "$CONCAT_FS_R"  "$tmp" "$(basename "$out_fs")"  "$bychr" "${fs_files[@]}"
    Rscript "$CONCAT_MUT_R" "$tmp" "$(basename "$out_tab")" "$bychr" "${tab_files[@]}"

    if [[ ! -s "${tmp}/$(basename "$out_fs")" || ! -s "${tmp}/$(basename "$out_tab")" ]]; then
      echo "concat: ${base} ${set_name}: the R scripts produced nothing usable" >&2
      rm -rf "$tmp"
      failed=1
      continue
    fi
    mv "${tmp}/$(basename "$out_fs")"  "$out_fs"
    mv "${tmp}/$(basename "$out_tab")" "$out_tab"
    rm -rf "$tmp"
    echo "concat: wrote ${out_fs}"
    echo "concat: wrote ${out_tab}"
  done
done

# ---------------------------------------------------------------- checkpoint --
# Diagnostic only: a bug in the check must never discard finished output, so its
# failure warns and the stage's own exit status is left alone.
log_dir="${OUT}/sfs/log"
if ! python3 "${REPO_DIR}/scripts/checks/check_sfs.py" \
       --out "${OUT}/sfs" --run "${FS_RUN}" \
       --prep-out "${OUT}/prep" --prep-run "${PREP_RUN}" \
       --chroms "${chroms[*]}" --log-dir "$log_dir" \
       --min-snps "${MIN_SNPS}" "$@"; then
  echo "WARNING: check_sfs.py failed; the sfs output above is unaffected." >&2
  echo "         Inspect the per-chromosome files by hand before trusting this run." >&2
fi

exit "$failed"
