#!/bin/bash
# dfe, the pick: collect the fits, choose the best one, discretise its DFE and
# draw the diagnostics.
#
#   bash scripts/dfe/run_bestdfe.sh <pdf> <bed> [bed...]
#
# <pdf> is one of PDFS other than `neu`; <bed> is the ORIGINAL bed. Run it once
# every pack of run_dfe.sh has finished.
#
# Ported from piperef_stddfe_8_plotdfe.sh, plot_converge_dfe_0430.py and
# plot_DFE_single.py. The rules, unchanged:
#
#   1. every *.InferDFE.bestfits in the run directory contributes ONE row: the
#      BEST data row of its "# Converged results" block. A file with no such
#      block contributes nothing.   (plot_converge_dfe_0430.py:282-321)
#   2. per row, scale_div2Na = scale/(2*Na), sdadi = shape*scale_div2Na and
#      Ne = theta/(2*u_test*L_test*2)                                 (:324-333)
#   3. the cutoff is the k-th largest log-likelihood with
#      k = max(round(n*0.05), 2), and rows with ll >= it are kept     (:369-374)
#   4. large_ll.txt holds the kept rows, sorted by descending ll, with the
#      columns filename, Log(likelihood), <the pdf's parameters>, theta,
#      scale_div2Na, sdadi, Ne                                        (:412-430)
#   5. the top rows get a model-vs-data figure and a folded model spectrum, from
#      cache1d.integrate(params, None, PDFs.<pdf>, theta, None)       (:450-516)
#   6. the discretised DFE is the FIRST row of large_ll.txt through the gamma
#      CDF over DFE_BINS, scaled by (1-pneu), with pneu added to the first bin
#      (plot_DFE_single.py:74-126). pneu is 0 when the pdf has no such column.
#
# Note on rule 1: the original's parser overwrites its per-file entry on every
# converged row and does not break, so the row that survived was the LAST of the
# converged block - the WORST converged fit in that file. We take the best;
# see docs/tocheck.md for what that changes against the published values.
# That is what produced every published number (the reference gamma fit's file
# holds -610.1733440069229 first and the published -610.1733440139265 last), so
# it is ported exactly. It is in the report.
#
# Differences from the original, all deliberate:
#   * `neugammazero` is handled. The original hit a bare `exit 1`
#     (piperef_stddfe_8_plotdfe.sh:105-107), and so did any other pdf.
#   * the mutation-load calculation (plot_mutload_single.py, :112-119) is NOT
#     ported. It is out of scope and the original's own header flags it as
#     wrong; u and L still reach the summary through the fsmut tables.
#   * nothing is copied to a dated ~/mv_<date> directory or to plots_dfe_mutload.
#   * figure names: `ll_<param>_all.pdf` (every converged fit) and
#     `ll_<param>_top.pdf` (the fits above the cutoff). The original put the
#     cutoff value in the file name; it is in spec.txt instead. The per-fit
#     figures are `sfs_rank<k>.pdf` / `sfs_rank<k>_poisson.pdf` with the model
#     spectrum `model_sfs_rank<k>.fs`, so rank 1 has a stable name for
#     check_dfe.py; the original named them after rounded parameter values,
#     which collided and silently overwrote several of the top ten.
#   * the cutoff index is clamped to the last available fit. The original
#     indexes llsort[max(round(n*0.05), 2)] unguarded, so a run with two or
#     fewer converged fits dies with an IndexError.
#   * the `ll ~ param` figure over ALL converged fits is drawn with the points
#     correctly paired. The original sorts `lls` in place while it is still the
#     y-vector of that scatter (plot_converge_dfe_0430.py:337-341, 353), so
#     every point in those six figures pairs a parameter with another fit's
#     likelihood. It affects no number, only the figures, and reproducing a
#     knowingly scrambled figure would be worse than fixing it - so this one is
#     corrected, and it is in the report.
#   * PDF only; PNG=1 additionally writes PNG.
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
: "${DEMO_MODEL:?DEMO_MODEL not set}"
: "${CACHE_RUN:?CACHE_RUN not set}"
: "${CACHE_MARKER:?CACHE_MARKER not set}"
: "${DFE_RUN:?DFE_RUN not set}"
: "${CTRL_TAG:?CTRL_TAG not set}"
: "${DFE_BINS:?DFE_BINS not set}"
: "${NUM_CONVERGENCE:?NUM_CONVERGENCE not set}"
: "${DFE_MAX_NOTRUN_FRAC:?DFE_MAX_NOTRUN_FRAC not set}"
: "${DFE_MIN_CONVERGED_FRAC:?DFE_MIN_CONVERGED_FRAC not set}"
: "${DFE_MIN_CONVERGED_N:?DFE_MIN_CONVERGED_N not set}"
: "${RESID_RARE_MAX_COUNT:?RESID_RARE_MAX_COUNT not set}"
: "${RESID_MAX_RARE_WARN:?RESID_MAX_RARE_WARN not set}"
: "${RESID_MAX_RARE_FAIL:?RESID_MAX_RARE_FAIL not set}"
: "${RESID_BULK_BIAS_WARN:?RESID_BULK_BIAS_WARN not set}"
: "${RESID_REST_OVER_RARE_MIN:?RESID_REST_OVER_RARE_MIN not set}"
: "${G2_PER_SNP_WARN:?G2_PER_SNP_WARN not set}"
: "${G2_PER_SNP_FAIL:?G2_PER_SNP_FAIL not set}"

: "${DFE_CUTOFF_RANK:?DFE_CUTOFF_RANK not set}"
: "${DFE_PLOT_TOP:?DFE_PLOT_TOP not set}"

if [[ $# -lt 2 ]]; then
  echo "usage: run_bestdfe.sh <pdf> <bed> [bed...]" >&2
  exit 2
fi
pdf="$1"; shift

# Every pdf in PDFS is handled or named in the refusal; nothing falls through.
case "$pdf" in
  neu)
    echo "run_bestdfe: pdf 'neu' has nothing to pick: its single likelihood is" >&2
    echo "      written by scripts/dfe/run_fitneu.sh, which runs its own check." >&2
    exit 2 ;;
  gamma|neugamma|neugammazero) ;;
  *)
    echo "run_bestdfe: pdf '${pdf}' is not supported." >&2
    echo "      Supported: gamma, neugamma, neugammazero (and neu, via" >&2
    echo "      run_fitneu.sh). PDFS is currently '${PDFS:-unset}'." >&2
    exit 2 ;;
esac

# --------------------------------------------------------------- preflight ----
beds=(); test_fs=(); mut_t=(); mut_c=(); nas=(); caches=(); dirs=(); csvnames=()
for bed in "$@"; do
  base="$(basename "$bed" .bed)"
  shortbase="s.nc_exHET_noq.${base}"
  shortbase="${shortbase/s.nc_/}"
  shortbase="${shortbase/nc_/}"
  fsbase2="${shortbase}_${CTRL_TAG}.chrs_1t22"
  sfs_dir="${OUT}/sfs/${base}/${FS_RUN}"
  demo_dir="${OUT}/demo/${base}/${FS_RUN}-${DEMO_MODEL}-${DEMO_RUN}"
  best="${demo_dir}/best_inference.txt"
  d="${OUT}/dfe/${base}/${CACHE_RUN}-${CACHE_MARKER}-${pdf}-${DFE_RUN}"

  t_fs="${sfs_dir}/has_qui.${fsbase2}.fs"
  t_mut="${sfs_dir}/fsmut.has_qui.${fsbase2}.tab"
  c_mut="${sfs_dir}/fsmut.quiecent_ctrl.${fsbase2}.tab"
  for f in "$t_fs" "$t_mut" "$c_mut" "$best"; do
    [[ -s "$f" ]] || {
      echo "run_bestdfe: missing upstream output for '${base}': ${f}" >&2
      exit 1; }
  done
  na="$(awk -F: '/^Na:/ {print $2}' "$best")"
  [[ -n "$na" ]] || { echo "run_bestdfe: no Na: line in ${best}" >&2; exit 1; }

  cache="${DFE_CACHE_FILE:-${OUT}/cache/${base}/${DEMO_RUN}-${CACHE_MARKER}-${CACHE_RUN}/cache_${CACHE_MARKER}.spectra.bpkl}"
  [[ -s "$cache" ]] || {
    echo "run_bestdfe: no cache for '${base}': ${cache}" >&2
    echo "      Check CACHE_RUN / CACHE_MARKER, or set DFE_CACHE_FILE." >&2
    exit 1; }

  if ! compgen -G "${d}/*.InferDFE.bestfits" > /dev/null; then
    echo "run_bestdfe: no fits for '${base}' ${pdf} in ${d}" >&2
    echo "      Run scripts/dfe/run_dfe.sh (or scripts/sge/submit_dfe.sh)" >&2
    echo "      first, or check DFE_RUN." >&2
    exit 1
  fi

  beds+=("$base"); test_fs+=("$t_fs"); mut_t+=("$t_mut"); mut_c+=("$c_mut")
  nas+=("$na"); caches+=("$cache"); dirs+=("$d")
  # The historical csv name, which the summary stage reads.
  csvnames+=("dfe_${fsbase2}_${pdf}.csv")
done

activate_conda

# ------------------------------------------------------------------- the pick --
failed=0
for i in "${!beds[@]}"; do
  base="${beds[$i]}"
  d="${dirs[$i]}"
  echo "run_bestdfe: ${base} ${pdf}: collecting the fits in ${d}"
  if ! python3 - "$d" "$pdf" "${test_fs[$i]}" "${caches[$i]}" "${mut_t[$i]}" \
                "${nas[$i]}" "${csvnames[$i]}" "$DFE_CUTOFF_RANK" \
                "$DFE_PLOT_TOP" "${PNG:-0}" $DFE_BINS <<'PY'
"""Collect the DFE fits, pick the best, discretise it, plot.

A faithful port of plot_converge_dfe_0430.py and plot_DFE_single.py; the
differences are listed in run_bestdfe.sh's header.
"""
import csv
import glob
import os
import pickle
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")                     # no display on a compute node
import matplotlib.pyplot as plt           # noqa: E402
import pandas as pd                       # noqa: E402
from scipy.stats import gamma as gamma_dist  # noqa: E402

import dadi                               # noqa: E402
import dadi.DFE as DFE                    # noqa: E402

(dfe_dir, pdf, testfs, cachefile, fsmutt, na_arg, csvname, cutoff_rank_arg,
 plot_top_arg) = sys.argv[1:10]
want_png = sys.argv[10] not in ("", "0", "no", "false")
bins = [float(b) for b in sys.argv[11:]]

ROULETTE_SCALE = 1.015e-7
Na = float(na_arg)
cutoff_rank = int(cutoff_rank_arg)
plot_top = int(plot_top_arg)


def die(msg):
    sys.stderr.write("run_bestdfe: %s\n" % msg)
    sys.exit(1)


def save(fig, stem):
    """PDF always, PNG as well under PNG=1."""
    fig.savefig(os.path.join(dfe_dir, stem + ".pdf"))
    if want_png:
        fig.savefig(os.path.join(dfe_dir, stem + ".png"), dpi=200)
    plt.close(fig)


def getul2(fsmutfile):
    """2*u*L from the genome-wide row of an fsmut table (plot_converge...:100)."""
    with open(fsmutfile) as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    if "etal" not in row.get("fsfile", ""):
        die("%s has no genome-wide ('etal') row first" % fsmutfile)
    length = float(row["length_adj_inputbed"])
    u = ROULETTE_SCALE * float(row["avepersite_diploi_Roulette_mutrate"]) / 2
    return 2 * float(u * length)


if len(bins) < 2:
    die("DFE_BINS needs at least two edges, got %r" % bins)

PDF_FUNCS = {
    "gamma": getattr(DFE.PDFs, "gamma", None),
    "neugamma": getattr(DFE.PDFs, "neugamma", None),
    "neugammazero": getattr(DFE.PDFs, "neugammazero", None),
}
pdf_func = PDF_FUNCS.get(pdf)
if pdf_func is None:
    die("dadi.DFE.PDFs has no '%s'. gamma is stock; neugamma and neugammazero "
        "come from setup/patches/dadi_DFE_PDFs.patch - run setup/setup.sh" % pdf)

mutation_test_diploid = getul2(fsmutt)
if mutation_test_diploid <= 0:
    die("the test mutation count came out as %r; check %s"
        % (mutation_test_diploid, fsmutt))
scale_scaleparam = 1.0 / (2 * Na)

# --------------------------------------------------------- collect the fits --
# One row per .bestfits file: the BEST (highest ll) row of its converged block.
# The original kept the LAST row instead - its dict was keyed by filename and
# updated once per row inside the convergence-tracking loop, so each row
# overwrote the previous. That was incidental, not intended (2026-09-30), and
# the difference from the published values is recorded in docs/tocheck.md.
files = sorted(glob.glob(os.path.join(dfe_dir, "*.InferDFE.bestfits")))
data_dict = {}
names = None
for path in files:
    with open(path) as fh:
        lines = fh.read().splitlines()
    found = False
    for line in lines:
        stripped = line.strip()
        if stripped == "# Converged results":
            found = True
            continue
        if not found:
            continue
        if stripped == "#":                    # end of the converged block
            break
        if stripped.startswith("# Log(likelihood)"):
            # dadi-cli names the columns; that is the authority on this pdf's
            # parameters, so no per-pdf table is needed here.
            names = [n.strip() for n in stripped.lstrip("# ").split("\t")]
            continue
        if stripped.startswith("#") or not stripped:
            continue
        try:
            values = [float(v) for v in stripped.split("\t")]
        except ValueError:
            break
        if names is None or len(values) != len(names):
            sys.stderr.write(
                "run_bestdfe: %s has %d values but %s columns; skipping\n"
                % (os.path.basename(path), len(values),
                   len(names) if names else "no"))
            break
        row = dict(zip(names, values))
        key = os.path.basename(path)
        prev = data_dict.get(key)
        if prev is None or row[names[0]] > prev[names[0]]:
            data_dict[key] = row

if not data_dict:
    die("none of the %d .InferDFE.bestfits file(s) in %s holds a converged "
        "result" % (len(files), dfe_dir))
for required in ("Log(likelihood)", "theta", "shape", "scale"):
    if required not in names:
        die("the .bestfits header has no '%s' column: %s" % (required, names))

ll_name = names[0]
for entry in data_dict.values():
    entry["scale_div2Na"] = entry["scale"] * scale_scaleparam
    entry["sdadi"] = entry["shape"] * entry["scale_div2Na"]
    # Na again, from the DFE fit's own theta; the original called it Ne.
    entry["Ne"] = entry["theta"] / (2 * mutation_test_diploid)
columns = list(names) + ["scale_div2Na", "sdadi", "Ne"]
params_plot = columns[1:]

# ------------------------------------------------------------------- cutoff --
# The k-th largest ll, k = max(round(n*0.05), cutoff_rank), clamped to the last
# fit so a run with few converged fits does not die with an IndexError.
lls_sorted = sorted((e[ll_name] for e in data_dict.values()), reverse=True)
k = max(int(round(len(lls_sorted) * 0.05)), cutoff_rank)
if k >= len(lls_sorted):
    k = len(lls_sorted) - 1
cutoff = lls_sorted[k]

best_key = max(data_dict, key=lambda key: data_dict[key][ll_name])
best = data_dict[best_key]
kept = {key: value for key, value in data_dict.items()
        if value[ll_name] >= cutoff}
second = None
if len(kept) > 1:
    ordered = sorted(kept, key=lambda key: kept[key][ll_name], reverse=True)
    second = kept[ordered[1]]


def circle(entry, param, label):
    plt.scatter(entry[param], entry[ll_name], edgecolors="red", marker="o",
                linewidths=1, facecolors="none", label=label)


def ll_scatter(entries, param, stem, title):
    fig = plt.figure()
    fig.clear()
    plt.scatter([e[param] for e in entries.values()],
                [e[ll_name] for e in entries.values()])
    if second is not None:
        circle(second, param, "2nd best ll")
    circle(best, param, "best ll")
    plt.xlabel(param)
    plt.ylabel(ll_name)
    plt.title(title)
    save(fig, stem)


for param in params_plot:
    ll_scatter(data_dict, param,
               "ll_%s_all" % param,
               "%s: ll ~ %s (%d converged fits)"
               % (pdf, param, len(data_dict)))
    ll_scatter(kept, param,
               "ll_%s_top" % param,
               "%s: ll ~ %s (ll >= %.4f, %d fits)" % (pdf, param, cutoff, len(kept)))

# -------------------------------------------------------------- large_ll.txt --
# The kept fits, sorted by descending ll. Column names and order are historical:
# the summary stage and plot_DFE_single's logic both key on them.
sorted_rows = sorted(kept.items(), key=lambda kv: kv[1][ll_name], reverse=True)
large_ll = os.path.join(dfe_dir, "large_ll.txt")
with open(large_ll, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=["filename"] + columns,
                            delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for filename, values in sorted_rows:
        row = {"filename": filename}
        row.update(values)
        writer.writerow(row)
print("run_bestdfe: %d of %d converged fits kept (ll >= %.6f) in %s"
      % (len(kept), len(data_dict), cutoff, large_ll))

# ------------------------------------------------- model spectra and figures --
cache1d = pickle.load(open(cachefile, "rb"))
sfs_test = dadi.Spectrum.from_file(testfs)
param_names = [n for n in names if n not in (ll_name, "theta")]

for rank, (filename, values) in enumerate(sorted_rows[:plot_top], start=1):
    popt = [values[n] for n in param_names]
    theta_test = values["theta"]
    model_sfs = cache1d.integrate(popt, None, pdf_func, theta_test, None)
    label = "%s rank %d  ll:%s  %s" % (
        pdf, rank, values[ll_name],
        "; ".join("%s:%s" % (n, v) for n, v in zip(param_names, popt)))

    fig = plt.figure()
    fig.clear()
    dadi.Plotting.plot_1d_comp_multinom(model_sfs, sfs_test)
    plt.title(label, fontsize=7)
    save(fig, "sfs_rank%d" % rank)

    fig = plt.figure()
    fig.clear()
    dadi.Plotting.plot_1d_comp_Poisson(model_sfs, sfs_test)
    plt.title(label, fontsize=7)
    save(fig, "sfs_rank%d_poisson" % rank)

    model_sfs.fold().to_file(
        os.path.join(dfe_dir, "model_sfs_rank%d.fs" % rank))

# --------------------------------------------------------- the discretised DFE --
# The FIRST row of large_ll.txt through the gamma CDF over the bins, scaled by
# (1-pneu), with pneu added to the first bin (plot_DFE_single.py:119-126). The
# scale is scale_div2Na, so the bins are in units of |s|, not of gamma.
top = sorted_rows[0][1]
shape = top["shape"]
scale_div2Na = top["scale_div2Na"]
pneu = float(top.get("pneu", 0.0))


def cdf(x):
    return gamma_dist.cdf(x, shape, 0, scale=scale_div2Na)


values = [cdf(bins[j + 1]) - cdf(bins[j]) for j in range(len(bins) - 1)]
values = [v * (1 - pneu) for v in values]
values[0] = values[0] + pneu

csvpath = os.path.join(dfe_dir, csvname)
labels = ["%s <= |s| < %s" % (bins[j], bins[j + 1]) for j in range(len(bins) - 1)]
rowname = os.path.splitext(csvname)[0]
frame = pd.DataFrame([values], columns=labels, index=[rowname])
frame.to_csv(csvpath)
print("run_bestdfe: wrote %s" % csvpath)
print("run_bestdfe: best ll=%s from %s; %s"
      % (top[ll_name], sorted_rows[0][0],
         "; ".join("%s=%s" % (n, top[n]) for n in param_names)))

fig = plt.figure(figsize=(8, 4))
plt.bar(range(len(labels)), values)
plt.xticks(range(len(labels)), labels, fontsize=8)
plt.xlabel("selection coefficient |s|")
plt.ylabel("Proportion of mutations")
plt.title(rowname, fontsize=9)
plt.tight_layout()
save(fig, os.path.splitext(csvname)[0])

with open(os.path.join(dfe_dir, "spec.bestdfe.txt"), "w") as fh:
    for key, value in (
            ("pdf", pdf), ("Na", Na),
            ("mutation_test_diploid", mutation_test_diploid),
            ("n_bestfits_files", len(files)),
            ("n_converged", len(data_dict)),
            ("cutoff_rank", k), ("likelihood_cutoff", cutoff),
            ("n_kept", len(kept)), ("best_file", sorted_rows[0][0]),
            ("best_ll", top[ll_name]), ("cachefile", cachefile),
            ("test_fs", testfs), ("fsmut_test", fsmutt),
            ("dfe_bins", " ".join(str(b) for b in bins)),
            ("dfe_csv", csvname), ("plot_top", plot_top)):
        fh.write("%s=%s\n" % (key, value))
PY
  then
    echo "run_bestdfe: ${base} ${pdf}: the pick failed" >&2
    failed=1
  fi
done

# ---------------------------------------------------------------- checkpoint --
# Diagnostic only: a bug in the check must never discard a finished stage, so its
# failure warns and the stage's own exit status is left alone.
if ! python3 "${REPO_DIR}/scripts/checks/check_dfe.py" \
       --out "${OUT}/dfe" --run "${CACHE_RUN}-${CACHE_MARKER}-${pdf}-${DFE_RUN}" --pdf "${pdf}" \
       --neu-run "${CACHE_RUN}-${CACHE_MARKER}-neu-${DFE_RUN}" \
       --sfs-out "${OUT}/sfs" --fs-run "${FS_RUN}" \
       --rare-max-count "${RESID_RARE_MAX_COUNT}" \
       --max-rare-warn "${RESID_MAX_RARE_WARN}" \
       --max-rare-fail "${RESID_MAX_RARE_FAIL}" \
       --bulk-bias-warn "${RESID_BULK_BIAS_WARN}" \
       --rest-over-rare-min "${RESID_REST_OVER_RARE_MIN}" \
       --g2-warn "${G2_PER_SNP_WARN}" --g2-fail "${G2_PER_SNP_FAIL}" \
       --num-convergence "${NUM_CONVERGENCE}" \
       --max-notrun-frac "${DFE_MAX_NOTRUN_FRAC}" \
       --min-converged-frac "${DFE_MIN_CONVERGED_FRAC}" \
       --min-converged-n "${DFE_MIN_CONVERGED_N}" "$@"; then
  echo "WARNING: check_dfe.py failed; the dfe output above is unaffected." >&2
  echo "         Inspect large_ll.txt and sfs_rank1.pdf by hand before" >&2
  echo "         trusting this run." >&2
fi

exit "$failed"
