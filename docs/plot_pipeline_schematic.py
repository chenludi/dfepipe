#!/usr/bin/env python3
"""Draw the pipeline as the scripts actually submit it, into docs/.

    python3 docs/plot_pipeline_schematic.py        # writes the png and the pdf

Every box is one qsub, the launcher jobs included. The label beside a box is the
script its worker runs, and the second line inside says how the job is sized.
Those come from the submitters, not from memory:

    prep          submit_prep.sh   one job, all beds     prep/run_prep.sh
    sfs           submit_sfs.sh    1-(chroms x beds)     sfs/run_sfs_chrom.sh
    concat        step             one job               sfs/run_concat.sh
    demo          submit_demo.sh   1-(p0 files x beds)   demo/run_demo.sh
    bestdemo      step             one job               demo/run_bestdemo.sh
    cache         submit_cache.sh  1-(beds)              cache/run_cache.sh
    fitneu        step             one job               dfe/run_fitneu.sh
    dfe           submit_dfe.sh    1-(packs), per pdf    dfe/run_dfe.sh
    bestdfe       step, per pdf    one job               dfe/run_bestdfe.sh
    summary       step             one job               summary/run_summary.sh
    bestpdf       step             one job               summary/run_bestpdf.sh
    the launchers step             one job               chains/run_all.sh FROM=...

The three phases are real, not presentation. run_all.sh exits after submitting a
launcher because submit_cache.sh reads best_inference.txt (:58) and
submit_dfe.sh generates the command packs (:144), neither of which exists at
phase 1 qsub time - and an array's size must be known when it is submitted.

Arrows are -hold_jid. The amber line inside each box is its gate: every worker
sources gate_upstream.sh and calls dfepipe_gate_or_exit before doing any work,
refusing to start if that check file holds a FAIL. The gate names the nearest
stage that has actually written a check file, which is why concat gates on prep
rather than sfs, and bestdfe on the cache rather than on the dfe array.
"""
import pathlib
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = pathlib.Path(__file__).resolve().parent

INK, GREY, EDGE, CHECK = "#1a1a1a", "#8a8a8a", "#4a6d8c", "#8a6d28"
BOX, ARRAY, LAUNCH = "#ffffff", "#e8f0f7", "#faf6ec"
BAND = ["#f7f5f0", "#eef2f4", "#f4f0f5"]

BOXW, BOXH = 2.9, 0.86
LW = 2.45          # half-width of a wide launcher box

#  name, x, y, sizing, script, gate, writes, label side
NODES = [
    ("prep",     0.0,  8.40, "one job, all beds",     "prep/run_prep.sh",      None,      "check_bed, check_prep", "r"),
    ("sfs",      0.0,  7.25, "array   chroms x beds", "sfs/run_sfs_chrom.sh",  "prep",    None,                    "r"),
    ("concat",   0.0,  6.10, "one job",               "sfs/run_concat.sh",     "prep",    "check_sfs",             "r"),
    ("demo",     0.0,  4.95, "array   p0 x beds",     "demo/run_demo.sh",      "sfs",     None,                    "r"),
    ("bestdemo", 0.0,  3.80, "one job",               "demo/run_bestdemo.sh",  "sfs",     "check_demo",            "r"),
    ("launch:cache", 0.0, 2.70, "one job",            "chains/run_all.sh FROM=cache", None, None,                  "r"),
    ("cache",    0.0,  1.30, "array   beds",          "cache/run_cache.sh",    "demo",    "check_cache",           "r"),
    ("launch:dfe",   0.0, 0.20, "one job",            "chains/run_all.sh FROM=dfe",   None, None,                  "r"),
    ("fitneu",  -2.6, -1.45, "one job",               "dfe/run_fitneu.sh",     "cache",   "check_dfe (neu)",       "l"),
    ("dfe",      2.2, -1.45, "array   packs, per pdf","dfe/run_dfe.sh",        "cache",   None,                    "r"),
    ("bestdfe",  2.2, -2.80, "one job per pdf",       "dfe/run_bestdfe.sh",    "cache",   "check_dfe",             "r"),
    ("summary", -0.2, -4.15, "one job",               "summary/run_summary.sh","dfe",     "check_summary",         "r"),
    ("bestpdf", -0.2, -5.40, "one job",               "summary/run_bestpdf.sh","summary", None,                    "r"),
]

HOLDS = [("prep", "sfs"), ("sfs", "concat"), ("concat", "demo"),
         ("demo", "bestdemo"), ("bestdemo", "launch:cache"),
         ("launch:cache", "cache"), ("cache", "launch:dfe"),
         ("launch:dfe", "fitneu"), ("launch:dfe", "dfe"), ("dfe", "bestdfe"),
         ("fitneu", "summary"), ("bestdfe", "summary"), ("summary", "bestpdf")]

BANDS = [(9.05, 2.20, "phase 1    FROM=prep",
          "ends with a launcher, because submit_cache.sh reads "
          "best_inference.txt and bestdemo has not written it yet"),
         (2.10, -0.33, "phase 2    FROM=cache",
          "ends with a launcher, because submit_dfe.sh generates the "
          "command packs and their number sizes the dfe array"),
         (-0.43, -6.00, "phase 3    FROM=dfe",
          "runs through to the table")]

fig, ax = plt.subplots(figsize=(15.0, 14.6))
pos = {}

for i, (top, bot, name, why) in enumerate(BANDS):
    ax.add_patch(FancyBboxPatch((-7.15, bot), 12.9, top - bot,
                                boxstyle="round,pad=0.10", linewidth=0,
                                facecolor=BAND[i], zorder=0))
    ax.text(-7.00, top - 0.24, name, fontsize=12, fontweight="bold",
            color=INK, va="top", ha="left", zorder=1)
    ax.text(-7.00, top - 0.58, "\n".join(textwrap.wrap(why, 44)),
            fontsize=7.9, color=GREY, va="top", ha="left",
            style="italic", linespacing=1.45, zorder=1)

for name, x, y, sizing, script, gate, writes, side in NODES:
    pos[name] = (x, y)
    launcher = name.startswith("launch")
    w = LW if launcher else BOXW / 2
    face = LAUNCH if launcher else (ARRAY if sizing.startswith("array") else BOX)
    ax.add_patch(FancyBboxPatch((x - w, y - BOXH / 2), 2 * w, BOXH,
                                boxstyle="round,pad=0.05", facecolor=face,
                                edgecolor=EDGE, linewidth=1.6,
                                linestyle="--" if launcher else "-", zorder=3))
    dy = 0.20 if gate else 0.11
    ax.text(x, y + dy, name, fontsize=12, fontweight="bold", color=INK,
            ha="center", va="center", zorder=4)
    ax.text(x, y + dy - 0.29, sizing, fontsize=8.6,
            color=EDGE if sizing.startswith("array") else GREY,
            ha="center", va="center", zorder=4)
    if gate:
        ax.text(x, y - 0.27, "gate: %s" % gate, fontsize=7.8, color=CHECK,
                ha="center", va="center", zorder=4)

    if side == "b":
        ax.text(x, y - BOXH / 2 - 0.17, script, fontsize=8.6, color=INK,
                ha="center", va="center", family="monospace", zorder=4)
        if writes:
            ax.text(x, y - BOXH / 2 - 0.42, "writes " + writes, fontsize=7.8,
                    color=CHECK, ha="center", va="center", zorder=4)
    else:
        tx = x + w + 0.18 if side == "r" else x - w - 0.18
        ha = "left" if side == "r" else "right"
        ax.text(tx, y + 0.13, script, fontsize=8.6, color=INK, ha=ha,
                va="center", family="monospace", zorder=4)
        if writes:
            ax.text(tx, y - 0.15, "writes " + writes, fontsize=7.8,
                    color=CHECK, ha=ha, va="center", zorder=4)

for a, b in HOLDS:
    (x1, y1), (x2, y2) = pos[a], pos[b]
    ax.add_patch(FancyArrowPatch(
        (x1, y1 - BOXH / 2), (x2, y2 + BOXH / 2),
        arrowstyle="-|>", mutation_scale=14, linewidth=1.6, color=EDGE,
        zorder=2, connectionstyle="arc3,rad=%.2f" % (0.0 if x1 == x2 else 0.16)))

ax.text(-7.00, 9.85, "dfepipe: bed files in, DFE table out",
        fontsize=17, fontweight="bold", color=INK, ha="left")
ax.text(-7.00, 9.50, "bash scripts/chains/run_all.sh <bed> [bed...]",
        fontsize=9.6, color=GREY, ha="left", family="monospace")
ax.text(-7.00, 9.24,
        "Every box is one qsub. Blue fill = array job, dashed = a launcher that "
        "re-runs run_all.sh for the next phase. Arrows are -hold_jid.",
        fontsize=8.6, color=GREY, ha="left")

ax.text(-7.00, -6.40,
        "Each bed keeps its own control set, demography, cache and DFE "
        "throughout. The beds never pool; only the job packing groups them.",
        fontsize=9.0, color=INK, ha="left")
ax.text(-7.00, -6.67,
        "prep is deliberately not an array: run_prep.sh writes the stage's one "
        "check file itself, so parallel tasks would race on it.",
        fontsize=8.6, color=GREY, ha="left")
ax.text(-7.00, -6.94,
        "PDFS = neu gamma neugamma, so phase 3 runs a dfe array and a bestdfe "
        "for gamma and for neugamma; neu is fitted once, by fitneu.",
        fontsize=8.6, color=GREY, ha="left")

ax.set_xlim(-7.3, 8.0)
ax.set_ylim(-7.15, 10.2)
ax.axis("off")
fig.tight_layout()
for ext in ("png", "pdf"):
    p = OUT / ("pipeline_schematic.%s" % ext)
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    print("wrote %s" % p)
