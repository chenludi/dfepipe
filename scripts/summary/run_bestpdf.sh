#!/bin/bash
# summary, model selection: which pdf does each bed's spectrum actually support,
# and what is its DFE?
#
#   bash scripts/summary/run_bestpdf.sh
#
# Takes no arguments: it reads the sumdfe.tsv that run_summary.sh wrote for this
# SUMMARY_RUN and writes bestpdf_sumdfe.tsv beside it.
#
# Ported from summary_getbestpdf.py. The comparison, unchanged (:44-57):
#
#   delta_gamma_neu      = ll_gamma    - ll_neu
#   delta_neugamma_gamma = ll_neugamma - ll_gamma
#   delta_neugamma_neu   = ll_neugamma - ll_neu
#
#   gamma beats neutral  when delta_gamma_neu      > cut[gamma-neu]
#   then, from whichever won:
#     from gamma    neugamma wins when delta_neugamma_gamma > cut[neugamma-gamma]
#     from neutral  neugamma wins when delta_neugamma_neu   > cut[neugamma-neu]
#
# The cutoffs are one row per confidence level in MODELSEL_TABLE, so every bed
# gets a `pdf_best_<confidence>` column, plus `pdf_best_allconfidence` which is
# filled only when every level agrees. Then the winning pdf's six bins are
# copied into `<primary>_bin_*`, with the neutral case being [1, 0, 0, 0, 0, 0]
# (:124-126), and num_nonneutral = (1 - first bin) * ul (:170-176).
#
# MODELSEL_PRIMARY is `delta_50`: an empirical choice of Chenlu's, not derived
# from a test statistic. The chi_5% / chi_0.10% rows are the nominal
# likelihood-ratio cutoffs; the delta_N rows are flat log-likelihood margins,
# which is what the published figures use.
#
# Differences from the original, all deliberate:
#   * the DFE bins are looked up by name. The original walks columns to the
#     right of a `pdf` heading collecting anything matching /\|\s*s\s*\|/, and
#     tells the two pdfs apart by `pdf` vs `pdf.1` - the names pandas invents
#     for the duplicate headings (:85-123). run_summary.sh prefixes its columns
#     by pdf, so that machinery is not needed.
#   * a small mutational target uses MODELSEL_SMALL_PRIMARY instead. The
#     original takes one --confidence for every bed. A flat margin of 50
#     log-likelihood units is a different test on 1,500 SNPs than on 500,000,
#     so below MODELSEL_SMALL_L the nominal chi-square cutoff is used and the
#     row records which applied. Set MODELSEL_SMALL_L=0 for the original
#     behaviour.
#   * `ul` comes from run_summary.sh rather than being recomputed here.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-$(dirname "$(dirname "$here")")}"
# shellcheck source=../../config/pipeline.conf
source "${REPO_DIR}/config/pipeline.conf"
# shellcheck source=../../config/env.sh
source "${REPO_DIR}/config/env.sh"

: "${OUT:?OUT not set}"
: "${SUMMARY_RUN:?SUMMARY_RUN not set}"
: "${MODELSEL_TABLE:?MODELSEL_TABLE not set: the confidence/cutoff csv}"
: "${MODELSEL_PRIMARY:?MODELSEL_PRIMARY not set: which confidence level the DFE columns use}"
: "${MODELSEL_SMALL_PRIMARY:?MODELSEL_SMALL_PRIMARY not set}"
: "${MODELSEL_SMALL_L:?MODELSEL_SMALL_L not set: the mutational target below which the small-target cutoff applies}"
: "${DFE_BINS:?DFE_BINS not set}"

d="${OUT}/summary/${SUMMARY_RUN}"
in_tab="${d}/sumdfe.tsv"
out_tab="${d}/bestpdf_sumdfe.tsv"

[[ -s "$in_tab" ]] || {
  echo "bestpdf: no summary table at ${in_tab}" >&2
  echo "      Run scripts/summary/run_summary.sh first, or check SUMMARY_RUN." >&2
  exit 1; }
[[ -s "$MODELSEL_TABLE" ]] || {
  echo "bestpdf: no model-selection table at ${MODELSEL_TABLE}" >&2
  exit 1; }
if [[ -s "$out_tab" && "${FORCE:-0}" != "1" ]]; then
  echo "bestpdf: ${out_tab} already exists; FORCE=1 to replace it." >&2
  exit 1
fi

tmp="$(mktemp -d "${d}/.bestpdf.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT

# No activate_conda: the standard library only, as run_summary.sh.
python3 - "$tmp" "$in_tab" "$MODELSEL_TABLE" "$MODELSEL_PRIMARY" \
         "$MODELSEL_SMALL_PRIMARY" "$MODELSEL_SMALL_L" "$DFE_BINS" <<'PY'
import csv
import os
import sys

(tmpdir, in_tab, cut_table, primary, small_primary, small_l_arg,
 bins_arg) = sys.argv[1:8]
small_l = float(small_l_arg)

edges = [float(b) for b in bins_arg.split()]
BIN_LABELS = ["%s <= |s| < %s" % (edges[j], edges[j + 1])
              for j in range(len(edges) - 1)]
NEUTRAL_DFE = [1.0] + [0.0] * (len(BIN_LABELS) - 1)


def die(msg):
    sys.stderr.write("bestpdf: %s\n" % msg)
    sys.exit(1)


def num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read_cutoffs(path):
    """[(confidence, {comparison: cutoff})] in file order."""
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        die("%s has no rows" % path)
    out = []
    for r in rows:
        try:
            out.append((r["confidence"], {
                "gamma-neu": float(r["gamma-neu"]),
                "neugamma-gamma": float(r["neugamma-gamma"]),
                "neugamma-neu": float(r["neugamma-neu"]),
            }))
        except (KeyError, TypeError, ValueError) as exc:
            die("%s: bad row %r (%s). Needs confidence, gamma-neu, "
                "neugamma-gamma, neugamma-neu" % (path, r, exc))
    return out


def best_pdf(deltas, cut):
    """summary_getbestpdf.py:44-57, unchanged.

    `deltas` is None for any comparison whose likelihood is missing, in which
    case no claim can be made and the answer is None rather than `neutral` -
    the original would have raised on the subtraction.
    """
    d_gn, d_ngg, d_ngn = deltas
    if d_gn is None:
        return None
    if d_gn > cut["gamma-neu"]:
        best = "gamma"
        if d_ngg is not None and d_ngg > cut["neugamma-gamma"]:
            best = "neugamma"
    else:
        best = "neutral"
        if d_ngn is not None and d_ngn > cut["neugamma-neu"]:
            best = "neugamma"
    return best


cutoffs = read_cutoffs(cut_table)
levels = [c for c, _ in cutoffs]
if primary not in levels:
    die("MODELSEL_PRIMARY=%r is not a confidence level in %s. It has: %s"
        % (primary, cut_table, ", ".join(levels)))
if small_primary not in levels:
    die("MODELSEL_SMALL_PRIMARY=%r is not a confidence level in %s. It has: %s"
        % (small_primary, cut_table, ", ".join(levels)))

with open(in_tab) as fh:
    rd = csv.reader(fh, delimiter="\t")
    columns = next(rd)
    rows = [dict(zip(columns, r)) for r in rd if r]
if not rows:
    die("%s has a header but no rows" % in_tab)

new_cols = (["ll_gamma-ll_neu", "ll_neugamma-ll_gamma", "ll_neugamma-ll_neu"]
            + ["pdf_best_%s" % c for c in levels]
            + ["pdf_best_allconfidence", "modelsel_level", "pdf_plotting"]
            + ["best_bin_%s" % b for b in BIN_LABELS]
            + ["num_nonneutral"])

for row in rows:
    ll_neu = num(row.get("ll_neu"))
    ll_gamma = num(row.get("ll_gamma"))
    ll_neugamma = num(row.get("ll_neugamma"))

    def delta(a, b):
        return None if a is None or b is None else a - b

    deltas = (delta(ll_gamma, ll_neu), delta(ll_neugamma, ll_gamma),
              delta(ll_neugamma, ll_neu))
    for name, value in zip(new_cols[:3], deltas):
        row[name] = "NA" if value is None else repr(value)

    picks = {}
    for confidence, cut in cutoffs:
        pick = best_pdf(deltas, cut)
        picks[confidence] = pick
        row["pdf_best_%s" % confidence] = pick or "NA"

    decided = [p for p in picks.values() if p is not None]
    row["pdf_best_allconfidence"] = (
        decided[0] if decided and len(set(decided)) == 1
        and len(decided) == len(picks) else "")

    # Which level decides the reported DFE. A flat margin is a different test
    # on a small target than on a large one, so below MODELSEL_SMALL_L the
    # nominal chi-square cutoff applies instead.
    l = num(row.get("l"))
    level = (small_primary if (l is not None and small_l > 0 and l < small_l)
             else primary)
    row["modelsel_level"] = level
    chosen = picks.get(level)
    row["pdf_plotting"] = chosen or "NA"

    if chosen is None:
        values = [None] * len(BIN_LABELS)
    elif chosen == "neutral":
        values = list(NEUTRAL_DFE)
    else:
        values = [num(row.get("%s_bin_%s" % (chosen, b))) for b in BIN_LABELS]
    for label, value in zip(BIN_LABELS, values):
        row["best_bin_%s" % label] = "NA" if value is None else repr(value)

    ul = num(row.get("ul"))
    first = values[0]
    if ul is None or first is None:
        row["num_nonneutral"] = "NA"
    elif first > 1:
        sys.stderr.write("bestpdf: %s: the neutral bin is %r, above 1\n"
                         % (row.get("outputbase"), first))
        row["num_nonneutral"] = "NA"
    else:
        row["num_nonneutral"] = repr((1 - first) * ul)

out_columns = columns + [c for c in new_cols if c not in columns]
path = os.path.join(tmpdir, "bestpdf_sumdfe.tsv")
with open(path, "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n",
                   quoting=csv.QUOTE_NONE)
    w.writerow(out_columns)
    for row in rows:
        w.writerow([row.get(c, "NA") for c in out_columns])

for row in rows:
    print("bestpdf: %s: %s at %s (all levels: %s)"
          % (row.get("outputbase"), row.get("pdf_plotting"),
             row.get("modelsel_level"),
             row.get("pdf_best_allconfidence") or "disagree"))
print("bestpdf: %d bed(s), %d columns -> bestpdf_sumdfe.tsv"
      % (len(rows), len(out_columns)))
PY

mv "${tmp}/bestpdf_sumdfe.tsv" "$out_tab"
echo "bestpdf: wrote ${out_tab}"
exit 0
