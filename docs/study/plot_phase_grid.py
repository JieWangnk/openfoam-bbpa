#!/usr/bin/env python3
"""plot_phase_grid.py -- the chapter-06 schematic: left vs centred windows.

A drawing, no case data: two periods of a 4-bin phase grid, the same
bin labels on both rows, with the integration windows placed left
(default) and centred. The same label names a different window in each
row; a centred bin 0 straddles the cycle boundary and belongs to the
cycle containing its label.

Usage:  python3 plot_phase_grid.py [--out figures/phase-grid.png]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))

BLUE = "#2a78d6"      # the highlighted window (same bin in both rows)
ORANGE = "#eb6834"    # bin-label markers
INK = "#0b0b0b"
INK2 = "#52514e"
GRAY1, GRAY2 = "#e6e5e1", "#d5d4cf"   # alternating window fills
SURFACE = "#fcfcfb"

T, NBINS = 1.0, 4
D = T / NBINS         # bin width, Delta


def row(ax, y0, h, offset, title, highlight):
    """One alignment row: windows as rectangles, labels as diamonds."""
    for c in range(2):
        for j in range(NBINS):
            x = c*T + j*D + offset
            hl = (c, j) == highlight
            ax.add_patch(Rectangle(
                (x, y0), D, h,
                facecolor=BLUE if hl else (GRAY1 if j % 2 else GRAY2),
                alpha=0.85 if hl else 1.0, edgecolor=SURFACE, lw=1.2))
            ax.text(c*T + j*D + 0.5*D + offset, y0 + 0.5*h, f"bin {j}",
                    ha="center", va="center", fontsize=8,
                    color="white" if hl else INK2)
    for c in range(2):
        for j in range(NBINS):
            ax.plot(c*T + j*D, y0, "D", color=ORANGE, ms=5, zorder=5,
                    clip_on=False)
    ax.text(-0.20, y0 + 0.5*h, title, ha="right", va="center",
            fontsize=9.5, color=INK)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "figures",
                                                  "phase-grid.png"))
    a = ap.parse_args()

    fig, ax = plt.subplots(figsize=(9.2, 4.0), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    # Rows: left alignment on top, centred below. Highlight bin 1 of
    # cycle 1 in both, so the half-bin shift of the SAME label is visible.
    row(ax, 1.30, 0.42, 0.0, "binAlignment left\n(default)", (1, 1))
    row(ax, 0.45, 0.42, -0.5*D, "binAlignment centred", (1, 1))

    # Cycle boundaries and captions.
    for x in (0, 1, 2):
        ax.plot([x, x], [0.18, 2.02], color=INK, lw=1.1, zorder=4)
    for c, x in [(0, 0.5), (1, 1.5)]:
        ax.text(x, 2.06, f"cycle {c}", ha="center", va="bottom",
                fontsize=9.5, color=INK2)

    # The label of the highlighted bin (t = 1 + 1*D on both rows).
    xl = 1 + D
    ax.text(xl + 0.30, 2.28, "the same label,  "
            "t = phaseOrigin + cycle·period + j·Δ",
            fontsize=8.5, color=INK, ha="left", va="bottom")
    for ytip in (1.28, 0.47):
        ax.annotate("", xy=(xl, ytip), xytext=(xl + 0.31, 2.26),
                    arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))

    # Bin width Delta over one left-aligned window.
    ax.add_patch(FancyArrowPatch((0.25, 1.82), (0.5, 1.82),
                                 arrowstyle="<|-|>", mutation_scale=8,
                                 color=INK2, lw=0.9))
    ax.text(0.375, 1.87, "Δ = period / nBins", ha="center",
            va="bottom", fontsize=8.5, color=INK)

    # phaseOrigin.
    ax.annotate("phaseOrigin", xy=(0, 0.18), xytext=(-0.72, -0.06),
                fontsize=8.5, color=INK, va="top",
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))

    # The straddling centred bin 0 of cycle 1: window (1 - D/2, 1 + D/2).
    ax.annotate("a centred bin 0 straddles the cycle boundary;\n"
                "it belongs to the cycle containing its label",
                xy=(0.97, 0.45), xytext=(0.42, -0.06), fontsize=8.5,
                color=INK, ha="center", va="top",
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8,
                                shrinkB=6))

    # Orange diamond key.
    ax.plot(1.42, -0.30, "D", color=ORANGE, ms=5, clip_on=False)
    ax.text(1.46, -0.30, "bin labels (the output time directories)",
            fontsize=8.5, color=INK, va="center")

    ax.set_xlim(-0.85, 2.12)
    ax.set_ylim(-0.44, 2.5)
    ax.axis("off")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.tight_layout()
    fig.savefig(a.out, facecolor=SURFACE, bbox_inches="tight")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
