#!/usr/bin/env python3
"""plot_first_run.py -- the chapter-02 figure: U vs U_Bin vs U_PA.

Reads the tutorials/bbpaVerify case AFTER it has been run (chapter 01:
./Allclean && ./Allrun) and plots, at the cell where the velocity moves
the most, the three quantities that share a time directory:

  U      the instantaneous solution at each solver write (context, gray)
  U_Bin  each retained cycle's own bin mean, drawn over its centred
         window (blue segments)
  U_PA   the cross-cycle mean as of each cycle's dump (orange diamonds)

The picture shows the whole v2 story at a glance: the discarded startup
cycle leaves nothing, the first retained cycle's U_PA sits exactly on
its U_Bin, and the second cycle's U_PA moves to the mean of the two.

Usage:
  python3 plot_first_run.py [--case ../../tutorials/bbpaVerify]
                            [--out figures/u-bin-pa.png]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CASE = os.path.join(HERE, "..", "..", "tutorials", "bbpaVerify")

# Palette: categorical slots 1-2 of the validated default palette
# (docs: dataviz reference instance); gray ink for context and text.
BLUE = "#2a78d6"      # U_Bin
ORANGE = "#eb6834"    # U_PA
INK = "#0b0b0b"
INK2 = "#52514e"      # U (instantaneous) + secondary text
SURFACE = "#fcfcfb"

PERIOD, NBINS = 1.0, 10
BIN_DT = PERIOD / NBINS


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=DEFAULT_CASE)
    ap.add_argument("--out", default=os.path.join(HERE, "figures",
                                                  "u-bin-pa.png"))
    a = ap.parse_args()
    case = os.path.abspath(a.case)

    sys.path.insert(0, case)
    from check_bbpa import read_internal, time_dirs

    tds = time_dirs(case)
    inst = [(t, read_internal(os.path.join(case, d, "U")))
            for t, d in tds if os.path.isfile(os.path.join(case, d, "U"))
            and os.path.isfile(os.path.join(case, d, "p"))]
    if len(inst) < 10:
        sys.exit(f"{case}: no per-write U fields -- run the tutorial "
                 "first (chapter 01: ./Allclean && ./Allrun)")

    # Probe: the cell whose Ux moves the most over the run. A uniform
    # field (the 0/ initial condition) parses as one row; broadcast it.
    ncell = max(U.shape[0] for _, U in inst)
    ux = np.array([np.repeat(U[:, 0], ncell) if U.shape[0] == 1
                   else U[:, 0] for _, U in inst])     # (ntime, ncell)
    cell = int(np.argmax(np.std(ux, axis=0)))
    times = np.array([t for t, _ in inst])
    u = ux[:, cell]

    def at(t, fname):
        for tv, d in tds:
            if abs(tv - t) < 1e-7:
                p = os.path.join(case, d, fname)
                if os.path.isfile(p):
                    return read_internal(p)[cell, 0]
        return None

    # Retained-cycle windows that hold products (cycle 0 is discarded).
    labels, ubin, upa = [], [], []
    for c in range(0, 4):
        for j in range(NBINS):
            t = c*PERIOD + j*BIN_DT
            b, p = at(t, "U_Bin"), at(t, "U_PA")
            if b is not None or p is not None:
                labels.append(t)
                ubin.append(b)
                upa.append(p)
    if not labels:
        sys.exit(f"{case}: no U_Bin/U_PA files -- run the tutorial first")

    fig, ax = plt.subplots(figsize=(9.2, 4.6), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    # Cycle bands and their captions.
    tmax = times.max()
    for x in np.arange(0, tmax + 1e-9, PERIOD):
        ax.axvline(x, color=INK2, lw=0.6, alpha=0.35, zorder=0)
    ymin, ymax = u.min(), u.max()
    pad = 0.18*(ymax - ymin)
    ytop = ymax + pad
    for c, txt in [(0, "cycle 0 — discarded,\nno products"),
                   (1, "cycle 1"), (2, "cycle 2")]:
        ax.text(c + 0.5, ytop, txt, ha="center", va="top",
                fontsize=8.5, color=INK2)

    # U: the instantaneous context trace.
    ax.plot(times, u, "-", color=INK2, lw=1.2, alpha=0.75,
            marker="o", ms=2.6, label="U (instantaneous, per write)")

    # U_Bin: each retained cycle's own bin mean over its centred window.
    first = True
    for t, b in zip(labels, ubin):
        if b is None:
            continue
        ax.plot([t - BIN_DT/2, t + BIN_DT/2], [b, b], color=BLUE,
                lw=2.6, solid_capstyle="butt",
                label="U_Bin (this cycle's window mean)" if first else None)
        first = False

    # U_PA: the ensemble mean as of that cycle's dump.
    tpa = [t for t, p in zip(labels, upa) if p is not None]
    vpa = [p for p in upa if p is not None]
    ax.plot(tpa, vpa, "D", color=ORANGE, ms=4.6, mew=0, linestyle="none",
            label="U_PA (cross-cycle mean at the dump)")

    # Two teaching notes, placed in the empty pocket above each trough.
    span = ymax - ymin
    y1 = at(1.6, "U_PA")
    if y1 is not None:
        ax.annotate("one cycle folded:\nU_PA sits exactly on U_Bin",
                    xy=(1.6, y1 + 0.02*span), xytext=(1.72, ymin + 0.62*span),
                    ha="center", fontsize=8.5, color=INK,
                    arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    y2 = at(2.7, "U_PA")
    if y2 is not None:
        ax.annotate("the flow nearly repeats, so the\ntwo-cycle mean barely "
                    "moves; the tiny\nshift is what _PAVarianceMag stores",
                    xy=(2.7, y2 + 0.02*span), xytext=(2.72, ymin + 0.56*span),
                    ha="center", fontsize=8.5, color=INK,
                    arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))

    ax.set_xlabel("time [s]", color=INK)
    ax.set_ylabel(f"Ux at probe cell {cell} [m/s]", color=INK)
    ax.set_title("tutorials/bbpaVerify: the three quantities that share "
                 "a time directory", color=INK, fontsize=11, pad=34)
    ax.set_xlim(-0.05, tmax + 0.05)
    ax.set_ylim(ymin - 0.4*pad, ytop + 0.6*pad)
    ax.grid(axis="y", color=INK2, alpha=0.15, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8.5)
    leg = ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0),
                    ncol=3, fontsize=8.5, frameon=False,
                    columnspacing=1.4, handletextpad=0.6)
    for t in leg.get_texts():
        t.set_color(INK)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.tight_layout()
    fig.savefig(a.out, facecolor=SURFACE, bbox_inches="tight")
    print(f"wrote {a.out}  (probe cell {cell})")


if __name__ == "__main__":
    main()
