#!/bin/bash
# sfs, one chromosome: pair each test region with nearby candidate-neutral
# controls, then build the SFS and the fsmut table for both sets.
#
#   bash scripts/sfs/run_sfs_chrom.sh <chrom> <bed> [bed...]
#
# <chrom> is "1" or "chr1" (pipeline.conf's CHROMS is bare numbers; every bed,
# VCF and output name is chr-prefixed, so the two meet here). <bed> is the
# ORIGINAL bed, the same argument prep took: this stage reads prep's output,
# $OUT/prep/<base>/$PREP_RUN/s.nc_exHET_noq.<base>.bed.
#
# The recipe is ported from windctrl_splitchr.sh (control pairing) and
# get_fs_mut_bychr.sh (SFS + mutation rate), in the order windctrl_fsmut_chrs.sh
# runs them, with the file:line of each step:
#
#   1.  grep the chromosome out of the prep bed        windctrl_splitchr.sh:37
#   2.  sort it, and merge only if MERGE_BEFORE_WINDOW=1
#       (the original always merged)                   windctrl_splitchr.sh:41-42
#       NOTE this is a no-op either way: run_prep.sh:112 already merges its
#       output, so the bed arriving here is merged whatever this is set to.
#       See todo/tocheck.md item 1.
#   3.  quiescent_extest = the chromosome's candidate-
#       neutral bed MINUS the test bed, so a control
#       can never also be a test site                  windctrl_splitchr.sh:58
#   4.  one window per merged region,
#       [mid-CTRL_FLANK, mid+CTRL_FLANK] clipped at 0,
#       carrying the region's own coordinates in 4-5   windctrl_splitchr.sh:68-76
#   5.  windows overlapping any quiescent_extest site
#       -> columns 4-5 -> dedup, sort, merge
#          = smer.has_qui.*, THE test set downstream
#            uses, and the l reported in the summary   windctrl_splitchr.sh:90-105
#   6.  quiescent_extest within the windows, sorted
#       and merged = smer.quiecent_ctrl.*              windctrl_splitchr.sh:96,114-116
#   7.  for each of those two beds: intersect the
#       chromosome's VCF, project to PROJECTION, write
#       the .fs and the fsmut table                    get_fs_mut_bychr.sh:57-91
#
# The window is centred on each merged region's MIDPOINT and its width does not
# depend on the region's length, so a region longer than 2*CTRL_FLANK can fall
# entirely outside its own window and end up with no control. That is the
# published behaviour; check_sfs reports how much of the test set survives it.
#
# Output names are the historical ones verbatim - run_concat.sh, every later
# stage and every published result key on them:
#   bychr/smer.has_qui.<chrom>.<bedstem>.bed
#   bychr/smer.quiecent_ctrl.<chrom>.<bedstem>.bed
#   bychr/smer.<set>.<chrom>.<bedstem>.<vcfstem>.fs
#   bychr/fsmut.smer.<set>.<chrom>.<bedstem>.<chrom>.tab
# where <bedstem> is the prep bed's name without .bed and <vcfstem> is the
# chromosome's VCF name without .vcf.gz.
#
# Results are keyed by FS_RUN: $OUT/sfs/<base>/$FS_RUN/. Because this stage runs
# as an array job the guard is per chromosome, not per directory: a chromosome
# that FINISHED (it has a lengths.<chrom>.tsv) is never silently redone - use a
# new FS_RUN, or FORCE=1 - while one that was interrupted part-way simply runs
# again. Everything is built in a temp directory and renamed into place, so a
# killed task cannot leave a truncated file under a final name.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# shellcheck source=../../config/env.sh
source "${REPO_DIR}/config/env.sh"

: "${OUT:?OUT not set}"
: "${FS_RUN:?FS_RUN not set}"
: "${PREP_RUN:?PREP_RUN not set: the sfs stage reads the prep output}"
: "${PROJECTION:?PROJECTION not set}"
: "${CTRL_FLANK:?CTRL_FLANK not set}"
: "${MERGE_BEFORE_WINDOW:?MERGE_BEFORE_WINDOW not set}"
: "${CTRL_TAG:?CTRL_TAG not set}"
# The neutral controls come in one of two forms, and exactly one key is set.
# Two keys rather than one that is type-sniffed: a mistyped directory path would
# otherwise be silently taken for a single bed.
#   NEUTRAL_CTRL_DIR  one bed per chromosome, <CTRL_TAG>_<chrom>.bed. What this
#                     project has, and it saves all 22 tasks reading a 700 Mb
#                     file.
#   NEUTRAL_CTRL_BED  one genome-wide bed, sliced per chromosome below. Simpler
#                     to supply.
if [[ -n "${NEUTRAL_CTRL_DIR:-}" && -n "${NEUTRAL_CTRL_BED:-}" ]]; then
  echo "sfs: NEUTRAL_CTRL_DIR and NEUTRAL_CTRL_BED are both set; pick one." >&2
  echo "      NEUTRAL_CTRL_DIR=${NEUTRAL_CTRL_DIR}" >&2
  echo "      NEUTRAL_CTRL_BED=${NEUTRAL_CTRL_BED}" >&2
  exit 2
fi
if [[ -z "${NEUTRAL_CTRL_DIR:-}" && -z "${NEUTRAL_CTRL_BED:-}" ]]; then
  echo "sfs: the neutral controls are not set. Set NEUTRAL_CTRL_DIR to a" >&2
  echo "      directory of ${CTRL_TAG}_<chrom>.bed files, or NEUTRAL_CTRL_BED" >&2
  echo "      to one genome-wide bed (see config/profiles/)." >&2
  exit 2
fi
: "${VCF_TMPL:?VCF_TMPL not set: the per-chromosome VCF, with {chrom} in it}"
: "${MUTRATE_TAB:?MUTRATE_TAB not set: the per-site mutation-rate directory}"
: "${BEDTOOLS_MOD:?BEDTOOLS_MOD not set}"
# get_fs_mut_bychr.sh:26 and dadi_sfs_YRI.py:29 hardcoded popfile_2504_YRI108.txt
# and "YRI". Both are data choices, so pipeline.conf leaves them empty and the
# site profile sets them.
: "${POPFILE:?POPFILE not set: the dadi population panel, e.g. popfile_2504_YRI108.txt (add it to config/profiles/)}"
: "${POP:?POP not set: which population in POPFILE to build the SFS for, e.g. YRI (add it to config/profiles/)}"

# MUTRATE_PATTERN (pipeline.conf, from get_fs_mut_bychr.sh:70) names the files
# in MUTRATE_TAB; {chrnum} is the chromosome without the chr prefix. Held in a
# variable first: a literal '}' inside a ${VAR:-default} would close the
# expansion early and silently truncate it.
: "${MUTRATE_PATTERN:?MUTRATE_PATTERN not set}"
_mutrate_pattern="$MUTRATE_PATTERN"
MUTRATE_TMPL="${MUTRATE_TMPL:-${MUTRATE_TAB}/${_mutrate_pattern}}"

if [[ $# -lt 2 ]]; then
  echo "usage: run_sfs_chrom.sh <chrom> <bed> [bed...]" >&2
  exit 2
fi

raw_chrom="$1"; shift
case "$raw_chrom" in
  chr*) chrom="$raw_chrom" ;;
  *)    chrom="chr${raw_chrom}" ;;
esac
chrnum="${chrom#chr}"
beds=("$@")

bp() { awk '{s += ($3 > $2) ? $3 - $2 : $2 - $3} END {printf "%d\n", s+0}' "$1"; }

vcf="${VCF_TMPL//\{chrom\}/${chrom}}"
vcfname="$(basename "$vcf")"
vcfstem="${vcfname%.gz}"; vcfstem="${vcfstem%.vcf}"
if [[ -n "${NEUTRAL_CTRL_DIR:-}" ]]; then
  qui_bed="${NEUTRAL_CTRL_DIR}/${CTRL_TAG}_${chrom}.bed"; qui_whole=""
else
  qui_bed=""; qui_whole="$NEUTRAL_CTRL_BED"
fi
mutbed="${MUTRATE_TMPL//\{chrnum\}/${chrnum}}"

# --------------------------------------------------------------- preflight ----
# Everything checkable before any work is checked first: 22 array tasks should
# not each spend an hour before hitting the same typo.
for f in "$vcf" "${qui_bed:-$qui_whole}" "$mutbed" "$POPFILE"; do
  [[ -s "$f" ]] || { echo "sfs ${chrom}: missing or empty input: $f" >&2; exit 1; }
done

prep_beds=()
for bed in "${beds[@]}"; do
  base="$(basename "$bed" .bed)"
  prep_bed="${OUT}/prep/${base}/${PREP_RUN}/s.nc_exHET_noq.${base}.bed"
  [[ -s "$prep_bed" ]] || {
    echo "sfs ${chrom}: no prep output for '${base}': ${prep_bed}" >&2
    echo "      Run the prep stage first, or check PREP_RUN." >&2
    exit 1; }
  prep_beds+=("$prep_bed")
done

for bed in "${beds[@]}"; do
  base="$(basename "$bed" .bed)"
  done_marker="${OUT}/sfs/${base}/${FS_RUN}/bychr/lengths.${chrom}.tsv"
  if [[ -f "$done_marker" ]]; then
    if [[ "${FORCE:-0}" == "1" ]]; then
      echo "sfs ${chrom}: FORCE=1 - replacing the finished ${chrom} output for ${base}" >&2
    else
      echo "sfs ${chrom}: ${base} already has output for ${chrom}: ${done_marker}" >&2
      echo "      Results are keyed by FS_RUN, so a finished chromosome is never" >&2
      echo "      redone implicitly. Pick a new FS_RUN, or re-run with FORCE=1." >&2
      exit 1
    fi
  fi
done

load_modules "$BEDTOOLS_MOD"
activate_conda

# --------------------------------------------------------------------- SFS ----
# Build the spectrum and the sites table for one bed, in $work. This is
# dadi_sfs_YRI.py inlined: the same make_data_dict_vcf / from_data_dict(
# projections=[PROJECTION], polarized=False) call and the same sitescounts
# columns, so the fsmut tables stay readable by concate_mut.R. The `fsfile`
# column records the FINAL path, not the temporary one written to.
#
# Progress goes to stderr: stdout is the path of the sites table, which the
# caller captures.
build_fs() {
  local in_bed="$1" work="$2" dest="$3"
  local stem out_fs counts maskvcf
  stem="$(basename "$in_bed" .bed)"
  maskvcf="${work}/${stem}.${vcfname}"
  out_fs="${work}/${stem}.${vcfstem}.fs"
  counts="${work}/${stem}.sitescounts.tab"

  echo "sfs ${chrom}: intersect ${vcfname} with $(basename "$in_bed")" >&2
  bedtools intersect -a "$vcf" -b "$in_bed" -header | gzip > "$maskvcf"

  echo "sfs ${chrom}: project to ${PROJECTION} and write the spectrum" >&2
  python3 - "$PROJECTION" "$POPFILE" "$maskvcf" "$out_fs" "$in_bed" "$POP" \
            "$counts" "${dest}/${stem}.${vcfstem}.fs" >&2 <<'PY'
import gzip
import sys

import dadi

projnum = int(sys.argv[1])
popfile, invcf, out_fs, bed, pop, counts, label = sys.argv[2:9]

dd = dadi.Misc.make_data_dict_vcf(invcf, popfile)
fs = dadi.Spectrum.from_data_dict(dd, [pop], projections=[projnum], polarized=False)
fs.to_file(out_fs)

# sites_fsfile: the sum of the spectrum's counts row, exactly as
# dadi_sfs_YRI.py measured it and as concatefs.R later adds it up.
with open(out_fs) as fh:
    lines = fh.readlines()
sites_fs = sum(map(float, lines[1].split())) if len(lines) >= 2 else 0.0

# Polymorphic sites in the intersected VCF, i.e. before projection. Sites
# genotyped in fewer than `projnum` chromosomes are dropped by the projection,
# so sites_fs/sites_vcf is the fraction kept, and it scales the callable length.
sites_vcf = 0
with gzip.open(invcf, "rt", encoding="utf-8") as fh:
    for line in fh:
        if not line.startswith("#"):
            sites_vcf += 1

length_inputbed = 0
with open(bed) as fh:
    for line in fh:
        line = line.strip()
        if line and not line.startswith("#"):
            f = line.split("\t")
            if len(f) >= 3:
                length_inputbed += abs(int(f[2]) - int(f[1]))

if sites_vcf:
    ratio = sites_fs / sites_vcf
    length_adj = length_inputbed * sites_fs / sites_vcf
else:
    ratio = 0
    length_adj = 0

with open(counts, "w") as fh:
    fh.write("fsfile\tsites_fsfile\tsites_vcf\tratio_keptsites\t"
             "length_inputbed\tlength_adj_inputbed\tprojection_num\n")
    fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
        label, sites_fs, sites_vcf, ratio, length_inputbed, length_adj, projnum))
PY
  printf '%s\n' "$counts"
}

# The mutation-rate half of the fsmut table (get_fs_mut_bychr.sh:75-88).
# `averagemut` is the unweighted mean of the rate column over intersected ROWS,
# which is what the published tables hold; it is not bp-weighted. Unlike the
# original this writes its intermediates to the temp dir rather than into the
# mutation-rate reference directory, which must stay read-only.
build_mut() {
  local in_bed="$1" work="$2"
  local stem hit out l_bed l_hit coverage averagemut
  stem="$(basename "$in_bed" .bed)"
  hit="${work}/mut.${stem}.${chrom}.bed"
  out="${work}/mut.${stem}.${chrom}.tab"

  echo "sfs ${chrom}: intersect the mutation-rate table" >&2
  bedtools intersect -a "$mutbed" -b "$in_bed" > "$hit"
  l_bed="$(bp "$in_bed")"
  l_hit="$(bp "$hit")"
  if [[ "$l_bed" -gt 0 ]]; then
    coverage="$(awk -v a="$l_hit" -v b="$l_bed" 'BEGIN{printf "%.5f", a/b}')"
  else
    coverage="0.00000"
  fi
  averagemut="$(awk -F'\t' '{s += $NF} END {if (NR > 0) printf "%.9f", s/NR; else printf "0"}' "$hit")"

  printf 'mutbed\tlength_withmutinfo\tinputbed\tlength_inputbed\tratio_withmutinfo\tavepersite_diploi_Roulette_mutrate\n' > "$out"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(basename "$hit")" "$l_hit" "$(basename "$in_bed")" "$l_bed" "$coverage" "$averagemut" >> "$out"
  printf '%s\n' "$out"
}

# ------------------------------------------------------------------- per bed --
for i in "${!beds[@]}"; do
  bed="${beds[$i]}"
  prep_bed="${prep_beds[$i]}"
  base="$(basename "$bed" .bed)"
  bedstem="$(basename "$prep_bed" .bed)"      # s.nc_exHET_noq.<base>
  d="${OUT}/sfs/${base}/${FS_RUN}"
  bychr="${d}/bychr"
  mkdir -p "$bychr"

  tmp="$(mktemp -d "${bychr}/.sfs.${chrom}.XXXXXX")"
  # shellcheck disable=SC2064  # $tmp is expanded now on purpose, per iteration
  trap "rm -rf '${tmp}'" EXIT

  chrbed="${tmp}/chr.bed"
  # grep exits 1 on no match, which `set -e` would treat as a crash; a bed with
  # no intervals on this chromosome is legitimate, and check_sfs reports SKIP.
  grep "^${chrom}[[:space:]]" "$prep_bed" > "$chrbed" || true
  if [[ ! -s "$chrbed" ]]; then
    echo "sfs ${chrom}: ${base} has no intervals on ${chrom}; nothing to do"
    rm -rf "$tmp"; trap - EXIT
    continue
  fi

  lengths="${tmp}/lengths.tsv"
  printf 'stage\tbp\n' > "$lengths"
  printf 'prep_chr\t%s\n' "$(bp "$chrbed")" >> "$lengths"

  # MERGE_BEFORE_WINDOW decides what the control window is centred on. See the
  # config: merging was added to windctrl_splitchr.sh on 2024/11/04, so results
  # built before that date need MERGE_BEFORE_WINDOW=0 to be reproduced.
  if [[ "${MERGE_BEFORE_WINDOW}" == "1" ]]; then
    echo "sfs ${chrom}: ${base} 1/5 sort and merge the test bed"
    sort -T "$tmp" -k1,1 -k2,2n "$chrbed" > "${tmp}/sorted.bed"
    bedtools merge -i "${tmp}/sorted.bed" > "${tmp}/merged.bed"
  else
    echo "sfs ${chrom}: ${base} 1/5 sort the test bed (MERGE_BEFORE_WINDOW=0: pre-2024/11/04 behaviour)"
    sort -T "$tmp" -k1,1 -k2,2n "$chrbed" > "${tmp}/merged.bed"
  fi
  printf 'merged_chr\t%s\n' "$(bp "${tmp}/merged.bed")" >> "$lengths"

  echo "sfs ${chrom}: ${base} 2/5 subtract the test sites from the control candidates"
  if [[ -n "$qui_whole" ]]; then
    awk -v c="$chrom" '$1 == c' "$qui_whole" > "${tmp}/ctrl.bed"
    qui_bed="${tmp}/ctrl.bed"
    [[ -s "$qui_bed" ]] || {
      echo "sfs ${chrom}: no control intervals on ${chrom} in ${qui_whole}" >&2
      echo "      Controls are needed for every chromosome with test sites;" >&2
      echo "      check the chromosome naming matches CHROMS." >&2
      exit 1; }
  fi
  bedtools subtract -a "$qui_bed" -b "${tmp}/merged.bed" > "${tmp}/extest.bed"

  echo "sfs ${chrom}: ${base} 3/5 window each merged region by +/- ${CTRL_FLANK}"
  awk -v win="$CTRL_FLANK" '{
      mid = ($2 + $3) / 2
      start = int(mid - win); if (start < 0) start = 0
      end = int(mid + win)
      printf("%s\t%d\t%d\t%d\t%d\n", $1, start, end, $2, $3)
  }' "${tmp}/merged.bed" > "${tmp}/window.bed"

  echo "sfs ${chrom}: ${base} 4/5 keep the test regions that have a control"
  bedtools intersect -wa -a "${tmp}/window.bed" -b "${tmp}/extest.bed" \
    > "${tmp}/q.window.bed"
  awk '{printf("%s\t%d\t%d\n", $1, $4, $5)}' "${tmp}/q.window.bed" \
    | sort -T "$tmp" | uniq > "${tmp}/has_qui.raw.bed"
  sort -T "$tmp" -k1,1 -k2,2n "${tmp}/has_qui.raw.bed" \
    | bedtools merge -i - | sort -T "$tmp" -k1,1 -k2,2n > "${tmp}/has_qui.bed"
  printf 'has_qui\t%s\n' "$(bp "${tmp}/has_qui.bed")" >> "$lengths"

  echo "sfs ${chrom}: ${base} 5/5 collect the controls inside those windows"
  bedtools intersect -a "${tmp}/extest.bed" -b "${tmp}/window.bed" \
    > "${tmp}/ctrl.raw.bed"
  sort -T "$tmp" -k1,1 -k2,2n "${tmp}/ctrl.raw.bed" \
    | bedtools merge -i - | sort -T "$tmp" -k1,1 -k2,2n > "${tmp}/ctrl.bed"
  printf 'quiecent_ctrl\t%s\n' "$(bp "${tmp}/ctrl.bed")" >> "$lengths"

  if [[ ! -s "${tmp}/has_qui.bed" ]]; then
    echo "sfs ${chrom}: ${base} has intervals on ${chrom} but not one of them has a" >&2
    echo "      candidate-neutral site within ${CTRL_FLANK} bp of its midpoint, so no" >&2
    echo "      spectrum can be built. Consider a larger CTRL_FLANK." >&2
    rm -rf "$tmp"; trap - EXIT
    exit 1
  fi

  # The pairing is finished and self-consistent, so the two beds are named now;
  # the spectra are the slow part and are moved into place as each completes.
  has_qui="${bychr}/smer.has_qui.${chrom}.${bedstem}.bed"
  ctrl="${bychr}/smer.quiecent_ctrl.${chrom}.${bedstem}.bed"
  mv "${tmp}/has_qui.bed" "$has_qui"
  mv "${tmp}/ctrl.bed" "$ctrl"

  for in_bed in "$has_qui" "$ctrl"; do
    stem="$(basename "$in_bed" .bed)"
    counts="$(build_fs "$in_bed" "$tmp" "$bychr")"
    mut="$(build_mut "$in_bed" "$tmp")"
    # The fsmut table is the two halves pasted side by side, in this order:
    # concate_mut.R reads it back both by column name and by position.
    paste "$counts" "$mut" > "${tmp}/fsmut.tab"
    mv "${tmp}/${stem}.${vcfstem}.fs" "${bychr}/${stem}.${vcfstem}.fs"
    mv "${tmp}/fsmut.tab" "${bychr}/fsmut.${stem}.${chrom}.tab"
  done

  # lengths.<chrom>.tsv lands last: its presence is what marks the chromosome
  # finished, both for the guard above and for check_sfs.
  mv "$lengths" "${bychr}/lengths.${chrom}.tsv"

  # spec.txt is per bed and per run, and 22 array tasks race to write it, so it
  # is written to the temp dir and renamed: a rename is atomic, a redirect is not.
  { echo "FS_RUN=${FS_RUN}"
    echo "PREP_RUN=${PREP_RUN}"
    echo "prep_bed=${prep_bed}"
    echo "CTRL_FLANK=${CTRL_FLANK}"
    echo "CTRL_TAG=${CTRL_TAG}"
    echo "PROJECTION=${PROJECTION}"
    echo "POP=${POP}"
    echo "POPFILE=${POPFILE}"
    echo "VCF_TMPL=${VCF_TMPL}"
    echo "NEUTRAL_CTRL_DIR=${NEUTRAL_CTRL_DIR:-}"
    echo "NEUTRAL_CTRL_BED=${NEUTRAL_CTRL_BED:-}"
    echo "neutral_ctrl_used=${qui_bed:-$qui_whole}"
    echo "MUTRATE_TMPL=${MUTRATE_TMPL}"
    echo "BEDTOOLS_MOD=${BEDTOOLS_MOD}"
    echo "bedtools_version=$(bedtools --version 2>/dev/null || echo unknown)"
  } > "${tmp}/spec.txt"
  mv "${tmp}/spec.txt" "${d}/spec.txt"

  rm -rf "$tmp"
  trap - EXIT
  echo "sfs ${chrom}: ${base} done"
  cat "${bychr}/lengths.${chrom}.tsv"
done

# No checkpoint here: the verdict is per bed across ALL chromosomes, so it can
# only be reached once the array has finished. run_concat.sh runs check_sfs.py.
exit 0
