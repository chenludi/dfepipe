#!/bin/bash
# summary: one row per bed, gathering everything the DFE run produced into the
# table you actually read.
#
#   bash scripts/summary/run_summary.sh <bed> [bed...]
#
# <bed> is the ORIGINAL bed, the same argument every earlier stage took.
#
# Ported from piperef_stddfe_9_sumdfe.sh + summary_getul.py. Per bed it collects:
#
#   ll_neu            the `ll: ` line of run_fitneu.sh's
#                     <test fs base>.dadi_cal_ll_neu_0723.log
#   per pdf           the winning row of large_ll.txt - filename, log-likelihood
#                     and every fitted parameter
#   per pdf           the six discretised bins of dfe_<...>_<pdf>.csv
#   u, l, l_inputbed  from the genome-wide row of the test fsmut table, as
#                     summary_getul.py:19-33 computes them:
#                       u          = ROULETTE_SCALE * avepersite_diploi_Roulette_mutrate / 2
#                                    (ROULETTE_SCALE defaults to the original's
#                                    hardcoded 1.015e-7, summary_getul.py:18)
#                       l          = length_adj_inputbed
#                       l_inputbed = length_inputbed
#                     and ul = u * l, which several downstream scripts expect.
#
# The genome-wide row is the FIRST data row of the fsmut table, which is what
# `next(reader)` picks up. Checked rather than assumed: its length_inputbed
# (77,600,048 for phastCons17way 0-5) equals the sum of the 22 per-chromosome
# rows below it exactly.
#
# Differences from the original, all deliberate:
#   * the join is keyed on the bed, not positional. The original builds a
#     two-line block per field group and `paste`s them side by side, so any
#     block with the wrong number of lines - a missing log, an absent
#     large_ll.txt - shifts every later column by one and still produces a
#     well-formed row. Here a missing value is NA in its own column, and the
#     column count is fixed by the header.
#   * column names are prefixed by pdf: `gamma_shape`, `neugamma_shape`. The
#     original repeats the bare names once per pdf, so a reader cannot address
#     them and pandas silently renames the duplicates.
#   * no `awk '!seen[$0]++]'` pass. That existed to drop the header line the
#     original re-emitted for every bed; it also silently collapsed two beds
#     whose values happened to be identical.
#   * the output is keyed by SUMMARY_RUN, not deleted and rewritten. The
#     original `rm`s its output file first, so a failed run destroys the
#     previous table.
#
# Not ported: the best-pdf model selection of summary_getbestpdf.py, which turns
# this table into `bestpdf_sumdfe_*.tab`. That is the next script in the stage;
# MODELSEL_TABLE and friends in pipeline.conf are already its configuration.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# shellcheck source=../../config/env.sh
source "${REPO_DIR}/config/env.sh"

: "${OUT:?OUT not set}"
: "${FS_RUN:?FS_RUN not set}"
: "${DEMO_RUN:?DEMO_RUN not set}"
: "${CACHE_RUN:?CACHE_RUN not set}"
: "${CACHE_MARKER:?CACHE_MARKER not set}"
: "${DFE_RUN:?DFE_RUN not set}"
: "${SUMMARY_RUN:?SUMMARY_RUN not set: it keys the summary directory}"
: "${CTRL_TAG:?CTRL_TAG not set}"
: "${PDFS:?PDFS not set}"
: "${ROULETTE_SCALE:?ROULETTE_SCALE not set}"

if [[ $# -eq 0 ]]; then
  echo "usage: run_summary.sh <bed> [bed...]" >&2
  exit 2
fi

d="${OUT}/summary/${SUMMARY_RUN}"
out_tab="${d}/sumdfe.tsv"
if [[ -s "$out_tab" && "${FORCE:-0}" != "1" ]]; then
  echo "summary: ${out_tab} already exists." >&2
  echo "      Results are keyed by SUMMARY_RUN, so an existing table is never" >&2
  echo "      overwritten implicitly. Pick a new SUMMARY_RUN, or FORCE=1." >&2
  exit 1
fi
mkdir -p "$d"

# No activate_conda: this stage reads text files with the Python standard
# library only. Nothing here imports dadi, so it does not need the env.
tmp="$(mktemp -d "${d}/.summary.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT

python3 - "$tmp" "$OUT" "$FS_RUN" "$DEMO_RUN" "$CACHE_RUN" "$CACHE_MARKER" \
         "$DFE_RUN" "$CTRL_TAG" "$PDFS" "$ROULETTE_SCALE" "$DFE_BINS" "$@" <<'PY'
import csv
import glob
import os
import sys

(tmpdir, out, fs_run, demo_run, cache_run, cache_marker, dfe_run, ctrl_tag,
 pdfs_arg, roulette_arg, bins_arg) = sys.argv[1:12]
beds = sys.argv[12:]
pdfs = pdfs_arg.split()
roulette_scale = float(roulette_arg)

NEU_LOG_GLOB = "*.dadi_cal_ll_neu_0723.log"
BIN_PREFIX = "bin_"

# The columns are fixed by the config, never by what happened to be on disk. A
# header that grew or shrank with the files found would be exactly the problem
# this stage exists to avoid.
#
# What each pdf's large_ll.txt holds: dadi-cli's own parameters, then theta and
# the three run_bestdfe.sh derives (:309-313 and the header it writes).
PDF_PARAMS = {
    "gamma":    ["shape", "scale", "theta", "scale_div2Na", "sdadi", "Ne"],
    "neugamma": ["shape", "scale", "pneu", "theta", "scale_div2Na", "sdadi", "Ne"],
}

# The bin labels, formatted exactly as run_bestdfe.sh:420 does - the edges go
# through float() first, so "0" becomes "0.0" and "inf" stays "inf".
edges = [float(b) for b in bins_arg.split()]
BIN_LABELS = ["%s <= |s| < %s" % (edges[j], edges[j + 1])
              for j in range(len(edges) - 1)]


def die(msg):
    sys.stderr.write("summary: %s\n" % msg)
    sys.exit(1)


def short_base(base):
    """The published naming: prep's prefix dropped, first `s.nc_` then `nc_`.

    piperef_stddfe_2_concatefs.sh:12-14, quirks included - both are plain
    first-occurrence substitutions, which is why `nc_` can bite twice.
    """
    stem = "s.nc_exHET_noq.%s" % base
    stem = stem.replace("s.nc_", "", 1)
    stem = stem.replace("nc_", "", 1)
    return "%s_%s" % (stem, ctrl_tag)


def read_neu_ll(directory):
    """The `ll: ` line of run_fitneu.sh's log, or None."""
    hits = sorted(glob.glob(os.path.join(directory, NEU_LOG_GLOB)))
    if not hits:
        return None
    try:
        with open(hits[0]) as fh:
            for line in fh:
                if line.startswith("ll:"):
                    return float(line.split(":", 1)[1].strip())
    except (OSError, ValueError):
        return None
    return None


def read_best_fit(path):
    """The winning row of large_ll.txt as a dict, or None.

    run_bestdfe.sh writes it sorted by descending log-likelihood, so the first
    data row is the best fit.
    """
    try:
        with open(path) as fh:
            rd = csv.reader(fh, delimiter="\t")
            header = next(rd)
            for raw in rd:
                if raw:
                    return dict(zip(header, raw))
    except (OSError, StopIteration):
        return None
    return None


def read_dfe_bins(path):
    """{bin label: proportion} from the discretised DFE csv, or None.

    The csv is two lines: a header whose first cell is empty, and one data row
    whose first cell is the run name.
    """
    try:
        with open(path) as fh:
            rd = csv.reader(fh)
            header = next(rd)
            row = next(rd)
    except (OSError, StopIteration):
        return None
    return dict(zip(header[1:], row[1:]))


def read_ul(path):
    """(u, l, l_inputbed) from the genome-wide row of a fsmut table.

    summary_getul.py:19-33. The genome-wide row is the first data row, and its
    `fsfile` carries the `.etal` marker concate_mut.R:40 adds; the original treats a
    row without it as a format error, and so do we.
    """
    try:
        with open(path) as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            row = next(rd)
    except (OSError, StopIteration):
        return None
    if "etal" not in (row.get("fsfile") or ""):
        die("%s: the first row's fsfile has no '.etal' marker, so it is not the "
            "genome-wide row. The table may be truncated." % path)
    try:
        u = roulette_scale * float(row["avepersite_diploi_Roulette_mutrate"]) / 2
        l = float(row["length_adj_inputbed"])
        l_inputbed = float(row["length_inputbed"])
    except (KeyError, TypeError, ValueError) as exc:
        die("%s: cannot read u/l: %s" % (path, exc))
    return u, l, l_inputbed


# ----------------------------------------------------------------- the rows --
fit_pdfs = [p for p in pdfs if p != "neu"]
for pdf in fit_pdfs:
    if pdf not in PDF_PARAMS:
        die("pdf %r has no entry in PDF_PARAMS, so its columns are unknown" % pdf)

# Per pdf: filename, likelihood, the fitted parameters, then the bins - the
# published block order.
columns = ["outputbase", "ll_neu"]
for pdf in fit_pdfs:
    columns.append("%s_filename" % pdf)
    columns.append("ll_%s" % pdf)
    columns += ["%s_%s" % (pdf, p) for p in PDF_PARAMS[pdf]]
    columns += ["%s_%s%s" % (pdf, BIN_PREFIX, b) for b in BIN_LABELS]
columns += ["u", "l", "l_inputbed", "ul"]

rows = []
for bed in beds:
    base = os.path.basename(bed)
    if base.endswith(".bed"):
        base = base[:-4]
    outputbase = short_base(base)
    fsbase = "%s.chrs_1t22" % outputbase

    row = {"outputbase": outputbase}

    dfe_root = os.path.join(out, "dfe", base)
    run_key = "%s-%s-%%s-%s" % (cache_run, cache_marker, dfe_run)

    row["ll_neu"] = read_neu_ll(os.path.join(dfe_root, run_key % "neu"))

    for pdf in fit_pdfs:
        d = os.path.join(dfe_root, run_key % pdf)
        best = read_best_fit(os.path.join(d, "large_ll.txt"))
        if best is None:
            row["%s_filename" % pdf] = None
            row["ll_%s" % pdf] = None
        else:
            row["%s_filename" % pdf] = best.get("filename")
            row["ll_%s" % pdf] = best.get("Log(likelihood)")
            for key, value in best.items():
                if key in ("filename", "Log(likelihood)"):
                    continue
                row["%s_%s" % (pdf, key)] = value
        csvname = "dfe_%s_%s.csv" % (fsbase, pdf)
        bins = read_dfe_bins(os.path.join(d, csvname))
        if bins is not None:
            for label, value in bins.items():
                row["%s_%s%s" % (pdf, BIN_PREFIX, label)] = value

    fsmut = os.path.join(out, "sfs", base, fs_run,
                         "fsmut.has_qui.%s.tab" % fsbase)
    ul = read_ul(fsmut) if os.path.exists(fsmut) else None
    if ul is None:
        sys.stderr.write("summary: %s: no usable fsmut table at %s\n"
                         % (base, fsmut))
        row.update({"u": None, "l": None, "l_inputbed": None, "ul": None})
    else:
        u, l, l_inputbed = ul
        row.update({"u": u, "l": l, "l_inputbed": l_inputbed, "ul": u * l})

    rows.append(row)

# An unexpected parameter means run_bestdfe.sh changed and PDF_PARAMS did not.
# Absorbing it silently would make the header vary again, so stop instead.
for r in rows:
    for key in r:
        if key not in columns and key not in ("outputbase", "ll_neu",
                                              "u", "l", "l_inputbed", "ul"):
            die("%s is not a known column. If run_bestdfe.sh now writes it, add "
                "it to PDF_PARAMS in this script." % key)

path = os.path.join(tmpdir, "sumdfe.tsv")
with open(path, "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n",
                   quoting=csv.QUOTE_NONE)
    w.writerow(columns)
    for r in rows:
        line = ["NA" if r.get(c) is None else r[c] for c in columns]
        if len(line) != len(columns):
            die("row for %s has %d fields, header has %d"
                % (r.get("outputbase"), len(line), len(columns)))
        w.writerow(line)

print("summary: %d bed(s), %d columns -> sumdfe.tsv" % (len(rows), len(columns)))
PY

mv "${tmp}/sumdfe.tsv" "$out_tab"
echo "summary: wrote ${out_tab}"

# ---------------------------------------------------------------- checkpoint --
# Diagnostic only: a bug in the check must never discard finished output, so its
# failure warns and the stage's own exit status is left alone.
if ! python3 "${REPO_DIR}/scripts/checks/check_summary.py" \
       --out "${OUT}/summary" --run "${SUMMARY_RUN}" --pdfs "$PDFS" "$@"; then
  echo "WARNING: check_summary.py failed; the table above is unaffected." >&2
fi
exit 0
