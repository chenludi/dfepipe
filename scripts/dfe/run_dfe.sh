#!/bin/bash
# dfe, the fits: infer the DFE of the TEST spectrum with the demography held
# fixed, from many starting points.
#
#   bash scripts/dfe/run_dfe.sh <pdf> <bed> [bed...]
#
# <pdf> is one of PDFS other than `neu` (gamma, neugamma, neugammazero); `neu`
# has no parameters to fit and is scripts/dfe/run_fitneu.sh. <bed> is the
# ORIGINAL bed, the same argument every earlier stage took.
#
# Ported from piperef_stddfe_7_fitgamneug.sh + make_p0_dfe_0430.py (which
# starting points, and the InferDFE command line) and
# submit_dadi-cli_dfebash.sh (how a pack of command lines is run).
#
# The starting points, unchanged (make_p0_dfe_0430.py:144-168, 266-306):
#
#   shape0 over num_bins log-spaced values from <PDF>_LOWER[0] to <PDF>_UPPER[0]
#   scale0 likewise from <PDF>_LOWER[1] to <PDF>_UPPER[1]
#   pneu   over num_bins LINEARLY spaced values (neugamma/neugammazero only)
#   keep a combination when 2e-6*Na <= shape0*scale0 <= 2*0.5*Na, so the mean
#   |gamma| of the starting DFE lies inside the cache's gamma range
#   the surviving combinations are numbered 0..N-1 in that order, and every
#   PACKNUMBER of them go into one pack file dadi-cli_dfe.<first index>.sh
#
# The same bounds serve as the p0 range and as the optimiser's bounds, exactly
# as the original passed -lr/-ur and -p0lr/-p0ur the same two strings.
#
# Modes:
#   (default)              generate the packs, then run all of them in order
#   DFE_GENERATE_ONLY=1    generate the packs and stop (what submit_dfe.sh does)
#   DFE_PACK=<pack name>   run that one pack; generates nothing, guards nothing
#                          (what worker_dfe.sh does, one SGE task per pack)
#
# Differences from the original, all deliberate:
#   * no `rm -r $dirout` (piperef_stddfe_7_fitgamneug.sh:46) and no rename to
#     `<dirout>.bak` (make_p0_dfe_0430.py:284-295). Results are keyed by
#     DFE_RUN; a non-empty run directory stops the stage unless FORCE=1.
#   * bounds and bin counts come from the config, not from the per-pdf branches
#     hardcoded at piperef_stddfe_7_fitgamneug.sh:63-86.
#   * `--cpus` is spelled in full; the original wrote `--cpu`, which worked only
#     as an argparse abbreviation. Same choice as run_demo.sh.
#   * a fit whose .InferDFE.bestfits already exists is skipped, so RESUME=1
#     finishes a timed-out run without repeating it; FORCE=1 redoes every fit.
#   * the malformed makefile_p0_dfeinput.tab (nine header columns, eight data
#     fields) is replaced by spec.txt, which every stage in this package writes.
#   * every starting point's outcome is appended to fit_outcomes.tsv
#     (marker / outcome / seconds / note) as its fit finishes. check_dfe.py
#     reads it; the fourth outcome, not_run, is derived there.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# shellcheck source=../../config/env.sh
source "${REPO_DIR}/config/env.sh"

: "${OUT:?OUT not set}"
: "${FS_RUN:?FS_RUN not set: the dfe stage reads the sfs output}"
: "${DEMO_RUN:?DEMO_RUN not set}"
: "${DEMO_MODEL:?DEMO_MODEL not set}"
: "${CACHE_RUN:?CACHE_RUN not set: it keys the dfe run directory}"
: "${CACHE_MARKER:?CACHE_MARKER not set}"
: "${DFE_RUN:?DFE_RUN not set}"
: "${CTRL_TAG:?CTRL_TAG not set}"
: "${PDFS:?PDFS not set}"
: "${OPTCOUNT:?OPTCOUNT not set}"
: "${MAXEVAL:?MAXEVAL not set}"
: "${NUM_CONVERGENCE:?NUM_CONVERGENCE not set}"
: "${PACKNUMBER:?PACKNUMBER not set}"
: "${SGE_DFE_PE:?SGE_DFE_PE not set: the core count InferDFE is given}"

# NOT in pipeline.conf. Both are overridable from the environment and default to
# what ran; they belong in pipeline.conf alongside the other pdfs' bounds.
#   * neugammazero's bounds. piperef_stddfe_7_fitgamneug.sh:71-78 gave it a
#     shape upper bound of 3, not the 2 neugamma gets; the rest is identical.
#   * DFE_TIMEOUT, the per-command wall clock from submit_dadi-cli_dfebash.sh:66.
NEUGAMMAZERO_LOWER="${NEUGAMMAZERO_LOWER:-1e-10 1e-10 0}"
NEUGAMMAZERO_UPPER="${NEUGAMMAZERO_UPPER:-3 1e8 1}"
NEUGAMMAZERO_BINS="${NEUGAMMAZERO_BINS:-10}"
DFE_TIMEOUT="${DFE_TIMEOUT:-60m}"

if [[ $# -lt 2 ]]; then
  echo "usage: run_dfe.sh <pdf> <bed> [bed...]" >&2
  exit 2
fi
pdf="$1"; shift

# Every pdf in PDFS must be handled or named in the refusal. The original fell
# through to a bare `exit 1` for neugammazero (piperef_stddfe_8_plotdfe.sh:106).
case "$pdf" in
  neu)
    echo "run_dfe: pdf 'neu' has no parameters to fit." >&2
    echo "      Use scripts/dfe/run_fitneu.sh for the neutral null." >&2
    exit 2 ;;
  gamma)         lower="$GAMMA_LOWER";         upper="$GAMMA_UPPER";         nbins="$GAMMA_BINS" ;;
  neugamma)      lower="$NEUGAMMA_LOWER";      upper="$NEUGAMMA_UPPER";      nbins="$NEUGAMMA_BINS" ;;
  neugammazero)  lower="$NEUGAMMAZERO_LOWER";  upper="$NEUGAMMAZERO_UPPER";  nbins="$NEUGAMMAZERO_BINS" ;;
  *)
    echo "run_dfe: pdf '${pdf}' is not supported." >&2
    echo "      Supported: gamma, neugamma, neugammazero (and neu, via" >&2
    echo "      run_fitneu.sh). PDFS is currently '${PDFS}'. To add a pdf," >&2
    echo "      give it a <PDF>_LOWER / <PDF>_UPPER / <PDF>_BINS triple and a" >&2
    echo "      branch here, and patch dadi/dadi-cli to know the name." >&2
    exit 2 ;;
esac

# --------------------------------------------------------------- preflight ----
# Every input for every bed is located before any fitting starts.
beds=(); test_fs=(); mut_t=(); mut_c=(); popts=(); nas=(); caches=(); dirs=()
for bed in "$@"; do
  base="$(basename "$bed" .bed)"
  # The genome-wide spectrum names, as run_concat.sh builds them.
  shortbase="s.nc_exHET_noq.${base}"
  shortbase="${shortbase/s.nc_/}"
  shortbase="${shortbase/nc_/}"
  fsbase2="${shortbase}_${CTRL_TAG}.chrs_1t22"
  sfs_dir="${OUT}/sfs/${base}/${FS_RUN}"
  demo_dir="${OUT}/demo/${base}/${FS_RUN}-${DEMO_MODEL}-${DEMO_RUN}"
  best="${demo_dir}/best_inference.txt"

  t_fs="${sfs_dir}/has_qui.${fsbase2}.fs"
  t_mut="${sfs_dir}/fsmut.has_qui.${fsbase2}.tab"
  c_mut="${sfs_dir}/fsmut.quiecent_ctrl.${fsbase2}.tab"
  for f in "$t_fs" "$t_mut" "$c_mut"; do
    [[ -s "$f" ]] || {
      echo "run_dfe: missing sfs output for '${base}': ${f}" >&2
      echo "      Run scripts/sfs/run_concat.sh first, or check FS_RUN." >&2
      exit 1; }
  done
  [[ -s "$best" ]] || {
    echo "run_dfe: no demographic best fit for '${base}': ${best}" >&2
    echo "      Run scripts/demo/run_bestfit.sh first, or check DEMO_RUN." >&2
    exit 1; }

  key="$(awk -F: '/^max_ll_key:/ {print $2}' "$best")"
  popt="$(ls -1 "${demo_dir}"/*."${key}".InferDM.bestfits 2>/dev/null | head -1 || true)"
  [[ -n "$popt" ]] || {
    echo "run_dfe: no .bestfits for max_ll_key=${key} in ${demo_dir}" >&2
    exit 1; }
  na="$(awk -F: '/^Na:/ {print $2}' "$best")"
  [[ -n "$na" ]] || { echo "run_dfe: no Na: line in ${best}" >&2; exit 1; }

  # DFE_CACHE_FILE overrides the cache path for one bed only; it exists so an
  # existing cache built outside the pipeline can be reused (the validation run
  # does exactly that). NOT in pipeline.conf; see the report.
  cache="${DFE_CACHE_FILE:-${OUT}/cache/${base}/${DEMO_RUN}-${CACHE_MARKER}-${CACHE_RUN}/cache_${CACHE_MARKER}.spectra.bpkl}"
  [[ -s "$cache" ]] || {
    echo "run_dfe: no cache for '${base}': ${cache}" >&2
    echo "      Run scripts/cache/run_cache.sh first, check CACHE_RUN /" >&2
    echo "      CACHE_MARKER, or set DFE_CACHE_FILE to an existing cache." >&2
    exit 1; }

  beds+=("$base"); test_fs+=("$t_fs"); mut_t+=("$t_mut"); mut_c+=("$c_mut")
  popts+=("$popt"); nas+=("$na"); caches+=("$cache")
  dirs+=("${OUT}/dfe/${base}/${CACHE_RUN}-${CACHE_MARKER}-${pdf}-${DFE_RUN}")
done

pack_mode="${DFE_PACK:-}"

if [[ -z "$pack_mode" ]]; then
  for d in "${dirs[@]}"; do
    if [[ -d "$d" && -n "$(ls -A "$d" 2>/dev/null)" ]]; then
      if [[ "${FORCE:-0}" == "1" ]]; then
        echo "run_dfe: FORCE=1 - overwriting the existing run in ${d}" >&2
      elif [[ "${RESUME:-0}" == "1" ]]; then
        # A timed-out array leaves a part-filled run; resuming keeps the fits
        # that finished and runs only the missing starting points.
        echo "run_dfe: RESUME=1 - keeping the fits already in ${d}" >&2
      else
        echo "run_dfe: run directory exists and is not empty: ${d}" >&2
        echo "      Results are keyed by DFE_RUN, so an existing run is never" >&2
        echo "      overwritten implicitly. RESUME=1 keeps what finished and" >&2
        echo "      runs only the missing starting points (use this after an" >&2
        echo "      h_rt kill); FORCE=1 redoes every fit; or pick a new DFE_RUN." >&2
        exit 1
      fi
    fi
  done
fi

activate_conda

# ---------------------------------------------------------- the command packs --
if [[ -z "$pack_mode" ]]; then
  for i in "${!beds[@]}"; do
    d="${dirs[$i]}"
    mkdir -p "$d"
    echo "run_dfe: ${beds[$i]} ${pdf}: generating the starting points"
    python3 - "$d" "$pdf" "${test_fs[$i]}" "${caches[$i]}" "${popts[$i]}" \
                 "${mut_t[$i]}" "${mut_c[$i]}" "${nas[$i]}" \
                 "$lower" "$upper" "$nbins" "$OPTCOUNT" "$MAXEVAL" \
                 "$SGE_DFE_PE" "$NUM_CONVERGENCE" "$PACKNUMBER" <<'PY'
"""Write the InferDFE command packs. A faithful port of make_p0_dfe_0430.py.

Differences are listed in run_dfe.sh's header.
"""
import csv
import os
import sys
from itertools import product

import numpy as np

(dirout, pdf, testfs, cachefile, demopopt, fsmutt, fsmutc, na_arg,
 lower_s, upper_s, nbins_s, optcount, maxeval, num_cpu, num_convergence,
 packnumber_s) = sys.argv[1:17]

ROULETTE_SCALE = 1.015e-7
nbins = int(nbins_s)
packnumber = int(packnumber_s)
Na = float(na_arg)


def die(msg):
    sys.stderr.write("run_dfe: %s\n" % msg)
    sys.exit(1)


def getul(fsmutfile):
    """2*u*L from the genome-wide row of an fsmut table (make_p0...:170-188).

    The factor 2 is the original's; it cancels in the ratio below but is kept so
    the ratio's floating-point value is bit-for-bit the published one.
    """
    with open(fsmutfile) as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    if "etal" not in row.get("fsfile", ""):
        die("%s has no genome-wide ('etal') row first" % fsmutfile)
    length = float(row["length_adj_inputbed"])
    u = ROULETTE_SCALE * float(row["avepersite_diploi_Roulette_mutrate"]) / 2
    return 2 * float(u * length)


def valid_p0s(pdf, lo, hi, gamma0_range, nbins):
    """The starting grid (make_p0_dfe_0430.py:144-168).

    shape and scale are log-spaced, pneu LINEARLY spaced, and a combination
    survives when shape*scale - the gamma distribution's mean - falls inside the
    cache's gamma range. The iteration order fixes the index of every starting
    point, and those indices name the output files, so it must not change.
    """
    shape0s = np.logspace(np.log10(lo[0]), np.log10(hi[0]), num=nbins)
    scale0s = np.logspace(np.log10(lo[1]), np.log10(hi[1]), num=nbins)
    if pdf in ("neugamma", "neugammazero"):
        pneus = np.linspace(lo[2], hi[2], num=nbins)
        return [(a, b, p) for a, b, p in product(shape0s, scale0s, pneus)
                if gamma0_range[0] <= a * b <= gamma0_range[1]]
    if pdf == "gamma":
        return [(a, b) for a in shape0s for b in scale0s
                if gamma0_range[0] <= a * b <= gamma0_range[1]]
    die("pdf '%s' has no starting grid" % pdf)


lower = [float(v) for v in lower_s.split()]
upper = [float(v) for v in upper_s.split()]
want = 3 if pdf in ("neugamma", "neugammazero") else 2
if len(lower) != want or len(upper) != want:
    die("%s needs %d bounds, got %d lower and %d upper"
        % (pdf, want, len(lower), len(upper)))

mut_ctrl = getul(fsmutc)
mut_test = getul(fsmutt)
if mut_ctrl <= 0:
    die("the control mutation count came out as %r; check %s" % (mut_ctrl, fsmutc))
mutratio = mut_test / mut_ctrl

# The cache's gamma range, 2*s*Na with s from 1e-6 to 0.5 (make_p0...:266).
# Hardcoded there and matching config's GAMMA_S_LOW / GAMMA_S_HIGH defaults; the
# cache stage owns those keys, so they are not re-read here.
gamma0_range = [2 * 0.000001 * Na, 2 * 0.5 * Na]

pairs = valid_p0s(pdf, lower, upper, gamma0_range, nbins)
if not pairs:
    die("not one of the %d^%d starting combinations has shape*scale inside "
        "[%g, %g]: widen the bounds, or check Na=%g"
        % (nbins, want, gamma0_range[0], gamma0_range[1], Na))

# str() of a float is its shortest round-tripping form, which is what the
# published command lines carry ('2.0', '100000000.0', '1e-10').
lstr = " ".join(str(v) for v in lower)
ustr = " ".join(str(v) for v in upper)

packs, output = [], None
for i, p0 in enumerate(pairs):
    line = (
        "dadi-cli InferDFE"
        " --fs " + testfs +
        " --cache1d " + cachefile +
        " --pdf1d " + pdf +
        " --p0 " + " ".join(str(float(v)) for v in p0) +
        " --lbounds " + lstr +
        " --ubounds " + ustr +
        " --demo-popt " + demopopt +
        " --ratio " + str(mutratio) +
        " --output-prefix " + os.path.join(dirout, str(i)) +
        " --optimizations " + str(optcount) +
        " --maxeval " + str(maxeval) +
        " --nomisid" +
        " --cpus " + str(num_cpu) +
        " --check-convergence " + str(num_convergence)
    )
    if i % packnumber == 0:
        output = os.path.join(dirout, "dadi-cli_dfe.%d.sh" % i)
        packs.append(output)
        open(output, "w").close()
    with open(output, "a") as fh:
        fh.write(line + "\n")

with open(os.path.join(dirout, "spec.txt"), "w") as fh:
    for key, value in (
            ("pdf", pdf), ("n_starting_points", len(pairs)),
            ("n_packs", len(packs)), ("packnumber", packnumber),
            ("num_bins", nbins), ("lbounds", lstr), ("ubounds", ustr),
            ("p0_lbounds", lstr), ("p0_ubounds", ustr),
            ("Na", Na), ("gamma0_low", gamma0_range[0]),
            ("gamma0_high", gamma0_range[1]),
            ("mutation_ctrl_diploid", mut_ctrl),
            ("mutation_test_diploid", mut_test), ("mutratio", mutratio),
            ("test_fs", testfs), ("cachefile", cachefile),
            ("demopopt", demopopt), ("fsmut_test", fsmutt),
            ("fsmut_ctrl", fsmutc), ("optimizations", optcount),
            ("maxeval", maxeval), ("cpus", num_cpu),
            ("check_convergence", num_convergence)):
        fh.write("%s=%s\n" % (key, value))

print("run_dfe: %d starting points in %d pack(s), mutratio=%s"
      % (len(pairs), len(packs), mutratio))
PY
  done
  if [[ "${DFE_GENERATE_ONLY:-0}" != "0" ]]; then
    echo "run_dfe: DFE_GENERATE_ONLY set; the packs are written, nothing was run"
    exit 0
  fi
fi

# ------------------------------------------------------- the outcome log ----
# fit_outcomes.tsv records what became of every starting point this job
# attempts, appended as the fit finishes, so a task killed by h_rt still leaves
# a record of what it got through. It is APPEND-ONLY and never rewritten: a
# RESUME re-run adds rows, and check_dfe.py takes the LAST row for a marker.
# Concurrent array tasks share the file, so the append is serialised with flock
# where flock exists (one write of one short line otherwise).
log_outcome() {  # <logfile> <marker> <outcome> <seconds> <note>
  local f="$1"; shift
  local line; line="$(printf '%s\t%s\t%s\t%s' "$@")"
  if command -v flock >/dev/null 2>&1; then
    printf '%s\n' "$line" | flock "${f}.lock" bash -c 'cat >> "$1"' _ "$f"
  else
    printf '%s\n' "$line" >> "$f"
  fi
}

# ------------------------------------------------------------------- the fits --
# Each command line is run on its own with a timeout, as
# submit_dadi-cli_dfebash.sh:66-81 does: a starting point that wanders into a
# numerically impossible region must not take the rest of the pack down with it.
n_done=0; n_skipped=0; n_failed=0
for i in "${!beds[@]}"; do
  d="${dirs[$i]}"
  outlog="${d}/fit_outcomes.tsv"
  # noclobber makes the header write O_EXCL, so two array tasks racing to
  # create the log cannot truncate each other's rows.
  ( set -o noclobber; printf 'marker\toutcome\tseconds\tnote\n' > "$outlog" ) 2>/dev/null || true
  if [[ -n "$pack_mode" ]]; then
    packs=("${d}/${pack_mode}")
    [[ -s "${packs[0]}" ]] || {
      echo "run_dfe: no such pack: ${packs[0]}" >&2
      echo "      Run run_dfe.sh without DFE_PACK first to generate them." >&2
      exit 1; }
  else
    packs=()
    while IFS= read -r p; do packs+=("$p"); done < <(ls -1 "${d}"/dadi-cli_dfe.*.sh)
  fi

  for pack in "${packs[@]}"; do
    echo "run_dfe: ${beds[$i]} ${pdf}: $(basename "$pack")"
    while IFS= read -r cmd || [[ -n "$cmd" ]]; do
      [[ -z "${cmd// /}" ]] && continue
      case "$cmd" in \#*) continue ;; esac

      # The output prefix is the last field of --output-prefix; its .bestfits is
      # what marks this starting point finished.
      prefix="$(awk '{for (k = 1; k < NF; k++) if ($k == "--output-prefix") print $(k+1)}' <<< "$cmd")"
      marker="$(basename "${prefix:-}")"
      if [[ -n "$prefix" && -s "${prefix}.InferDFE.bestfits" && "${FORCE:-0}" != "1" ]]; then
        n_skipped=$((n_skipped + 1))
        log_outcome "$outlog" "$marker" converged 0 "already present, not re-run"
        continue
      fi

      t0=$SECONDS
      if timeout "$DFE_TIMEOUT" bash -c "$cmd"; then
        rc=0
        n_done=$((n_done + 1))
      else
        rc=$?
        if [[ "$rc" -eq 124 ]]; then
          echo "run_dfe: timed out after ${DFE_TIMEOUT}: ${cmd}" >&2
        else
          echo "run_dfe: exit ${rc}: ${cmd}" >&2
        fi
        n_failed=$((n_failed + 1))
      fi

      # A non-empty .bestfits is the fit; rc 124 is the timeout kill; anything
      # else that left no .bestfits ran its optimisations without ever getting
      # a likelihood. `rc` in the note keeps a real crash distinguishable.
      if [[ -n "$marker" ]]; then
        if [[ -s "${prefix}.InferDFE.bestfits" ]]; then
          log_outcome "$outlog" "$marker" converged "$((SECONDS - t0))" "rc=${rc}"
        elif [[ "$rc" -eq 124 ]]; then
          log_outcome "$outlog" "$marker" timeout "$((SECONDS - t0))" \
                      "killed after ${DFE_TIMEOUT}"
        else
          log_outcome "$outlog" "$marker" no_likelihood "$((SECONDS - t0))" \
                      "rc=${rc}, no .bestfits"
        fi
      fi
    done < "$pack"
  done
done

echo "run_dfe: ${pdf}: ${n_done} fitted, ${n_skipped} already present," \
     "${n_failed} failed or timed out"

# A pack in which nothing at all worked is a failure; individual starting points
# failing is normal and expected, which is why the original let them pass.
if [[ "$n_done" -eq 0 && "$n_skipped" -eq 0 ]]; then
  echo "run_dfe: not one fit succeeded" >&2
  exit 1
fi
exit 0
