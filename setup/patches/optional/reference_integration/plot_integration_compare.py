#!/usr/bin/env python3
"""EnhA neugamma DFE under the two cache-integration methods.

quad = stock dadi: the two exterior tail weights come from scipy.integrate.quad
CDF  = the Topatch/integration patch: those weights come from the gamma CDF

Drawn in the Fig. 2A style (plots/dfe_bars.py, plots/fig_style.py) so the bars are
directly comparable with the manuscript figure.

Input:  results_quad/dfe_EnhA_qui2.chrs_1t22_neugamma.csv
        results_cdf/dfe_EnhA_qui2.chrs_1t22_neugamma.csv
Output: EnhA_integration_compare.pdf (PNG=1 also writes .png)
        EnhA_integration_compare.csv  (the plotted numbers plus their difference)
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

# The Fig. 2A drawing code lives in the manuscript tree, outside this package.
# Walk up until we find it rather than counting directory levels, so the script
# keeps working wherever this directory is placed.
REL = os.path.join("manuscript", "manuscript_science_260821", "plots")
_d = HERE
while _d != os.path.dirname(_d):
    if os.path.isdir(os.path.join(_d, REL)):
        break
    _d = os.path.dirname(_d)
PLOTS = os.path.join(_d, REL)
if not os.path.isdir(PLOTS):
    sys.exit(f"cannot find the manuscript plotting helpers ({REL}) above {HERE}")
sys.path.insert(0, PLOTS)

import matplotlib                                     # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
import dfe_bars as db                                 # noqa: E402
import fig_style as fs                                # noqa: E402

COLORS = {"quad (stock dadi)": "blue", "CDF (patched)": "darkorange"}
SIZE = dict(w=4.6, h=2.8, left=0.62, right=0.08, bottom=0.62, top=0.22)
BAR_W = 0.3


def load(sub, label):
    """Read one dfe_*.csv into the five Fig. 2A classes.

    The csv splits the strongest class into 0.01-0.5 and >=0.5; Fig. 2A shows them
    summed, as plots/figS2/make_source_data.py does.
    """
    f = os.path.join(HERE, sub, "dfe_EnhA_qui2.chrs_1t22_neugamma.csv")
    row = pd.read_csv(f).iloc[0]
    out = {c: float(row[c]) for c in db.CLASSES[:-1]}
    out[db.CLASSES[-1]] = float(row["0.01 <= |s| < 0.5"]) + float(row["0.5 <= |s| < inf"])
    out["annotation"] = label
    return out


def main():
    rows = [load("results_quad", "quad (stock dadi)"), load("results_cdf", "CDF (patched)")]
    t = pd.DataFrame(rows)

    fig, ax = fs.new_panel(SIZE)
    db.grouped_bars(ax, t, COLORS, BAR_W)
    db.class_axis(ax)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Probability mass")
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=COLORS[c]) for c in t.annotation],
              labels=list(t.annotation), frameon=False, fontsize=7,
              loc="upper right", title="Cache integration", title_fontsize=7)
    ax.set_title("EnhA, neugamma: exterior-tail integration", fontsize=9, loc="left", pad=8)
    fs.save_panel(fig, os.path.join(HERE, "EnhA_integration_compare"))

    diff = t.set_index("annotation")[db.CLASSES].T
    diff["difference"] = diff["CDF (patched)"] - diff["quad (stock dadi)"]
    diff.to_csv(os.path.join(HERE, "EnhA_integration_compare.csv"))
    print(diff.to_string())
    print(f"\nlargest absolute difference in any |s| class: {diff['difference'].abs().max():.2e}")


if __name__ == "__main__":
    main()
