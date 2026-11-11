#!/usr/bin/env python3
"""Draw what each stage reads and what it writes, into docs/.

    python3 docs/plot_io_schematic.py       # writes the png and the pdf

A companion to plot_pipeline_schematic.py, which draws the same pipeline as a
job graph - one box per qsub. This one ignores jobs entirely and draws the files:
the inputs on the left, each stage's output directory in the middle, and on the
right the handful of files a reader actually opens.

Two things it exists to make obvious, both of which are otherwise only learnable
by reading the submitters:

  * Every input except the bed is config. The bed is the single positional
    argument to run_all.sh; everything else arrives through a site profile
    (config/profiles/template.conf) or config/pipeline.conf.
  * A stage's directory name carries its upstream run key as well as its own, so
    a path states its whole provenance. demo is <FS_RUN>-<DEMO_MODEL>-<DEMO_RUN>,
    cache is <DEMO_RUN>-<CACHE_MARKER>-<CACHE_RUN>, dfe is
    <CACHE_RUN>-<CACHE_MARKER>-<pdf>-<DFE_RUN>. That chaining is why two runs
    sharing a bed do not collide, and why a path is enough to say which cache a
    fit came from.

The path templates are taken from the scripts, not from memory: run_prep.sh,
run_sfs_chrom.sh, run_demo.sh, run_cache.sh, run_dfe.sh, run_bestdfe.sh,
run_summary.sh and run_bestpdf.sh.
"""
import pathlib
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt           # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parent

INK, GREY, EDGE, ACCENT = "#1a1a1a", "#8a8a8a", "#4a6d8c", "#8a6d28"
STAGE, INBOX, OUTBOX, CHECKC = "#ffffff", "#f2f6f9", "#faf6ec", "#f4f0f5"

XIN, XST, XOUT = -8.2, -1.35, 5.6      # column left edges
WIN, WST, WOUT = 5.9, 6.1, 5.6         # column widths

# Boxes are stacked from their own heights, not from a fixed step: the stage
# boxes differ by several lines and a fixed step overlapped them.
GAP_IN, GAP_OUT, GAP_STAGE = 0.26, 0.30, 0.62
LINE = 0.30

#  stage, directory template, files it writes. y is computed below.
STAGES = [
    ("prep",    "prep/<base>/<PREP_RUN>/",
     ["s.nc_exHET_noq.<base>.bed", "lengths.tsv"]),
    ("sfs",     "sfs/<base>/<FS_RUN>/",
     ["<fsbase>.fs   test and control spectra",
      "fsmut.<fsbase>.tab", "bychr/lengths.<chrom>.tsv"]),
    ("demo",    "demo/<base>/<FS_RUN>-<DEMO_MODEL>-<DEMO_RUN>/",
     ["best_inference.txt", "sfs_bestfit.pdf", "ll_*.pdf, Na_*.pdf"]),
    ("cache",   "cache/<base>/<DEMO_RUN>-<CACHE_MARKER>-<CACHE_RUN>/",
     ["cache_<CACHE_MARKER>.spectra.bpkl"]),
    ("dfe",     "dfe/<base>/<CACHE_RUN>-<CACHE_MARKER>-<pdf>-<DFE_RUN>/",
     ["<k>.InferDFE.bestfits   one per fit",
      "dfe_<group>_<pdf>.pdf and .csv",
      "sfs_rank<k>_poisson.pdf", "spec.bestdfe.txt, large_ll.txt"]),
    ("summary", "summary/<SUMMARY_RUN>/",
     ["sumdfe.tsv", "bestpdf_sumdfe.tsv"]),
]

#  into which stage, title, the variables
INPUTS = [
    ("prep", "the test bed", ["the one positional argument", "to run_all.sh"]),
    ("prep", "the four filters prep applies",
     ["MASK_BED      bases to keep", "CDS_BED, HET_BED, QUI_BED  removed"]),
    ("sfs", "the variants",
     ["VCF_TMPL   literal {chrom}", "POPFILE, POP"]),
    ("sfs", "the neutral controls",
     ["NEUTRAL_CTRL_DIR  <CTRL_TAG>_<chrom>.bed", "or NEUTRAL_CTRL_BED"]),
    ("sfs", "the mutation rates",
     ["MUTRATE_TAB", "MUTRATE_PATTERN, ROULETTE_SCALE"]),
    ("demo", "the starting-point grid",
     ["$REFS_DIR/dadi_initialparams/", "<DEMO_MODEL>_<DEMO_P0_SET>/"]),
    ("dfe", "what to fit",
     ["PDFS   gamma neugamma", "DFE_BINS, PACKNUMBER, DFE_TIMEOUT"]),
]

#  from which stage, what a reader opens it for
OUTPUTS = [
    ("demo", "the demography", ["best_inference.txt", "sfs_bestfit.pdf"]),
    ("dfe", "the discrete DFE, from the best fit",
     ["dfe_<group>_<pdf>.pdf", "dfe_<group>_<pdf>.csv", "spec.bestdfe.txt"]),
    ("dfe", "did the fit describe the data",
     ["sfs_rank1_poisson.pdf"]),
    ("summary", "the table the figures are built from",
     ["bestpdf_sumdfe.tsv", "   delta_50_* = the chosen pdf's bins"]),
]

CHECKS = ("every stage also writes its own verdict",
          ["prep/check_bed.<PREP_RUN>.tsv",
           "prep/check_prep.<PREP_RUN>.tsv",
           "sfs/check_sfs.<FS_RUN>.tsv",
           "demo/check_demo.<FS_RUN>-<DEMO_MODEL>-<DEMO_RUN>.tsv",
           "cache/<base>/<...>/check_cache.tsv",
           "dfe/check_dfe.<CACHE_RUN>-<CACHE_MARKER>-<pdf>-<DFE_RUN>.tsv",
           "summary/check_summary.<SUMMARY_RUN>.tsv",
           "",
           "A stage that cannot write a PASS row stops the chain,",
           "so these are what to read first when a run stops."])


def box(ax, x, y, w, h, face, **kw):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.07",
                                facecolor=face, edgecolor=kw.pop("edge", EDGE),
                                linewidth=kw.pop("lw", 1.5), zorder=3, **kw))


def arrow(ax, x0, y0, x1, y1, colour=EDGE, style="-|>", lw=1.5, ls="-"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style,
                                 mutation_scale=15, linewidth=lw,
                                 linestyle=ls, color=colour, zorder=2,
                                 shrinkA=2, shrinkB=2))


fig, ax = plt.subplots(figsize=(18.5, 13.6))
ypos, yh = {}, {}


def stage_height(writes):
    return 0.46 + LINE * (len(writes) + 1)      # +1 for the path template line


# ------------------------------------------------------------- the stages ----
y = 8.9
for name, template, writes in STAGES:
    h = stage_height(writes)
    y -= h / 2
    ypos[name], yh[name] = y, h
    box(ax, XST, y - h / 2, WST, h, STAGE, lw=1.9)
    ax.text(XST + 0.22, y + h / 2 - 0.26, name, fontsize=14,
            fontweight="bold", color=INK, va="center", zorder=4)
    ax.text(XST + 0.24, y + h / 2 - 0.26 - LINE, "$OUT/" + template,
            fontsize=8.8, color=EDGE, family="monospace", va="center", zorder=4)
    for i, f in enumerate(writes):
        ax.text(XST + 0.42, y + h / 2 - 0.26 - LINE * (i + 2), f, fontsize=8.6,
                color=GREY, family="monospace", va="center", zorder=4)
    y -= h / 2 + GAP_STAGE

for a, b in zip([s_[0] for s_ in STAGES], [s_[0] for s_ in STAGES][1:]):
    arrow(ax, XST + WST / 2, ypos[a] - yh[a] / 2, XST + WST / 2,
          ypos[b] + yh[b] / 2, colour=INK, lw=2.0)


def stack(items, xleft, width, face, edge, textcol, gap, fs_title, fs_body,
          ceiling):
    """Lay the side boxes out as one top-down flow.

    Each group would like to sit centred on its own stage, but a group that is
    taller than its stage's share would then run into the one above. So a box is
    placed at whichever is lower: where its stage wants it, or just under the
    previous box. Centring alone overlapped prep's inputs with sfs's.
    """
    bystage = {}
    order = []
    for stage, title, lines in items:
        if stage not in bystage:
            order.append(stage)
        bystage.setdefault(stage, []).append((title, lines))
    cursor = ceiling
    for stage in order:
        group = bystage[stage]
        hs = [0.40 + 0.27 * len(l) for _, l in group]
        total = sum(hs) + gap * (len(hs) - 1)
        top = min(ypos[stage] + total / 2, cursor)
        for (title, lines), h in zip(group, hs):
            cy = top - h / 2
            box(ax, xleft, cy - h / 2, width, h, face, lw=1.2, edge=edge)
            ax.text(xleft + 0.20, cy + h / 2 - 0.22, title, fontsize=fs_title,
                    fontweight="bold", color=INK, va="center", zorder=4)
            for i, t in enumerate(lines):
                ax.text(xleft + 0.36, cy + h / 2 - 0.50 - 0.27 * i, t,
                        fontsize=fs_body, color=textcol, family="monospace",
                        va="center", zorder=4)
            yield stage, cy
            top -= h + gap
        cursor = top


for stage, cy in stack(INPUTS, XIN, WIN, INBOX, EDGE, GREY, GAP_IN, 10, 8.5,
                       ceiling=8.75):
    arrow(ax, XIN + WIN, cy, XST, ypos[stage], colour=EDGE, lw=1.3)

for stage, cy in stack(OUTPUTS, XOUT, WOUT, OUTBOX, ACCENT, ACCENT, GAP_OUT,
                       10, 8.6, ceiling=8.75):
    arrow(ax, XST + WST, ypos[stage], XOUT, cy, colour=ACCENT, lw=1.3, ls="--")

# ------------------------------------------------------------- the checks ----
title, lines = CHECKS
h = 0.46 + 0.265 * len(lines)
# box() draws upward from its y, so the bottom has to allow for the height or
# the panel climbs back into the summary box.
ybot = min(ypos.values()) - yh[STAGES[-1][0]] / 2 - 0.75 - h
box(ax, XIN, ybot, WIN + 1.9, h, CHECKC, lw=1.2, edge=GREY)
ax.text(XIN + 0.20, ybot + h - 0.26, title, fontsize=10, fontweight="bold",
        color=INK, va="center", zorder=4)
for i, t in enumerate(lines):
    ax.text(XIN + 0.36, ybot + h - 0.58 - 0.265 * i, t, fontsize=8.3,
            color=GREY,
            family="monospace" if "/" in t else "sans-serif",
            style="italic" if t.startswith(("A stage", "so these")) else "normal",
            va="center", zorder=4)

ax.text(XOUT, ybot + h - 0.30,
        "\n".join(textwrap.wrap(
            "Run keys are per stage, and each directory name carries its "
            "upstream key as well as its own. So two runs that share a bed never "
            "collide, and a path on its own says which spectra, which "
            "demography and which cache a fit came from.", 44)),
        fontsize=9.2, color=INK, va="top", ha="left", linespacing=1.55,
        zorder=4)

ax.text(XIN, 9.75, "dfepipe: what each stage reads and writes", fontsize=17,
        fontweight="bold", color=INK, va="center", ha="left")
ax.text(XIN, 9.32,
        "inputs are config, except the bed  |  $OUT holds one directory per "
        "stage  |  companion to pipeline_schematic, which draws the jobs",
        fontsize=9.6, color=GREY, va="center", ha="left", style="italic")

ax.set_xlim(XIN - 0.5, XOUT + WOUT + 0.6)
ax.set_ylim(ybot - 0.5, 10.15)
ax.axis("off")
fig.tight_layout()
fig.savefig(OUT / "io_schematic.pdf")
fig.savefig(OUT / "io_schematic.png", dpi=170)
print("wrote %s and %s" % (OUT / "io_schematic.pdf", OUT / "io_schematic.png"))
