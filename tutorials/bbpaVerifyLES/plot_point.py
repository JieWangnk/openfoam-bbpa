#!/usr/bin/env python3
"""Follow ONE cell through time and show what BBPA actually does to it.

Three curves, one point in the field:

  U          the raw instantaneous velocity the solver computed
  U_Bin      that cycle's bin value -- U time-averaged over the bin containing
             each write time. Piecewise constant by construction, and DIFFERENT
             every cycle, because the turbulence is different every cycle.
  U_BinPA    the phase average -- the mean of the bin values across cycles.
             Phase-locked: the same curve repeats every period.

The band around U_BinPA is +/- one standard deviation across cycles, i.e.
sqrt(U_BinPAVariance). It is the honest error bar on the phase average: it says
how much the individual cycles disagree at that phase, which is exactly the
quantity that is identically zero in a clean periodic case (and therefore
untestable there).

Usage: python3 plot_point.py [--cell N] [--comp 0|1|2] [--out FILE]
"""
import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_bbpa import read_internal, time_dirs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=".")
    ap.add_argument("--cell", type=int, default=560)
    ap.add_argument("--comp", type=int, default=0, help="0=Ux 1=Uy 2=Uz")
    ap.add_argument("--ncells", type=int, default=1440)
    ap.add_argument("--period", type=float, default=1.0)
    ap.add_argument("--nbins", type=int, default=25)
    ap.add_argument("--out", default="bbpa_point.png")
    ap.add_argument("--style", choices=["centre", "step"], default="centre",
                    help="centre: plot each bin at its CENTRE time and connect "
                         "(the paper's convention). step: staircase at the "
                         "left-edge label the solver emits.")
    a = ap.parse_args()
    C, K = a.cell, a.comp
    comp_name = "xyz"[K]

    def val(path):
        arr = read_internal(path)
        if arr.shape[0] == 1:          # 'uniform' -> one row; broadcast
            arr = np.repeat(arr, a.ncells, axis=0)
        return arr[C, K]

    tds = time_dirs(a.case)

    inst = [(t, d) for t, d in tds
            if os.path.exists(os.path.join(a.case, d, "U"))
            and os.path.exists(os.path.join(a.case, d, "p"))]
    t_i = np.array([t for t, _ in inst])
    u_i = np.array([val(os.path.join(a.case, d, "U")) for _, d in inst])

    def series(suffix):
        rows = [(t, d) for t, d in tds
                if os.path.exists(os.path.join(a.case, d, f"U{suffix}"))]
        if not rows:
            return np.array([]), np.array([])
        return (np.array([t for t, _ in rows]),
                np.array([val(os.path.join(a.case, d, f"U{suffix}"))
                          for _, d in rows]))

    t_b, u_b = series("_Bin")
    t_p, u_p = series("_BinPA")
    t_v, u_v = series("_BinPAVariance")

    # A bin is emitted under its LEFT-EDGE label t_j, but a rectangular-window
    # average represents the signal at the window's CENTRE t_j + dt_bin/2: for a
    # locally linear signal the top-hat average equals the midpoint value
    # exactly. Labelling at the left edge therefore imposes a half-bin PHASE LAG
    # -- a first-order error that is pure bookkeeping, not physics. main.tex uses
    # the centre convention for exactly this reason ("the representative phase of
    # a rectangular-window average ... so the half-bin labelling offset is
    # removed rather than charged to the bin-width error"), which also makes the
    # residual curvature- rather than slope-limited.
    half = 0.5*a.period/a.nbins
    shift = half if a.style == "centre" else 0.0

    # The variance is written only for bins with MORE THAN ONE cycle
    # (binCounts_[binI] > 1 -- a single cycle has no cross-cycle spread yet),
    # so there are fewer _BinPAVariance files than _BinPA files and the two
    # series must be matched on TIME, not by position. Assuming equal lengths
    # silently drops the band.
    sd_at = {round(t, 9): np.sqrt(max(v, 0.0)) for t, v in zip(t_v, u_v)}

    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(11, 7.2), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    for k in range(1, int(t_i[-1]/a.period) + 1):
        for x in (ax, ax2):
            x.axvline(k*a.period, color="0.85", lw=1, zorder=0)

    ax.plot(t_i, u_i, color="0.55", lw=1.0, label="$U$  instantaneous (what the solver did)")
    if t_b.size:
        if a.style == "centre":
            ax.plot(t_b + shift, u_b, "o-", color="#1f77b4", lw=1.6, ms=4.5,
                    label="$U_{Bin}$  this cycle's bin (differs every cycle)")
        else:
            ax.step(t_b, u_b, where="post", color="#1f77b4", lw=1.8,
                    label="$U_{Bin}$  this cycle's bin (differs every cycle)")
            ax.plot(t_b, u_b, ".", color="#1f77b4", ms=5)
    if t_p.size:
        band = np.array([sd_at.get(round(t, 9), np.nan) for t in t_p])
        if a.style == "centre":
            ax.plot(t_p + shift, u_p, "o-", color="#d62728", lw=2.2, ms=4.5,
                    label="$U_{BinPA}$  phase average (repeats every cycle)")
            if np.isfinite(band).any():
                ax.fill_between(t_p + shift, u_p - band, u_p + band,
                                color="#d62728", alpha=0.22, lw=0,
                                label=r"$\pm\sqrt{U_{BinPAVariance}}$  cycle-to-cycle spread")
        else:
            ax.step(t_p, u_p, where="post", color="#d62728", lw=2.4,
                    label="$U_{BinPA}$  phase average (repeats every cycle)")
            if np.isfinite(band).any():
                ax.fill_between(t_p, u_p - band, u_p + band, step="post",
                                color="#d62728", alpha=0.25, lw=0,
                                label=r"$\pm\sqrt{U_{BinPAVariance}}$  cycle-to-cycle spread")

    ax.set_ylabel(f"$U_{comp_name}$  [m/s]")
    ax.set_title(f"BBPA at a single point (cell {C}) — period {a.period:g}, "
                 f"bin width {a.period/a.nbins:g}"
                 + (f"  [bins at their CENTRE phase, $t_j+\\Delta t_{{bin}}/2$]"
                    if a.style == "centre" else "  [staircase at left-edge label]")
                 + "\nthe phase average is the SAME curve each cycle; "
                 f"the bin values are not")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.95)
    ax.grid(alpha=0.25)

    # lower panel: how far each cycle's bin sits from the phase average
    if t_b.size and t_p.size:
        pa_at = {round(t, 9): v for t, v in zip(t_p, u_p)}
        d_t, d_v = [], []
        for t, v in zip(t_b, u_b):
            p = pa_at.get(round(t, 9))
            if p is not None:
                d_t.append(t); d_v.append(v - p)
        ax2.axhline(0, color="#d62728", lw=1.5)
        ax2.plot(np.array(d_t) + shift, d_v, "o-", color="#1f77b4", ms=4, lw=1.2)
        ax2.set_ylabel("$U_{Bin}-U_{BinPA}$")
        ax2.grid(alpha=0.25)
        ax2.set_title("this is the turbulent fluctuation the phase average removes "
                      "(exactly 0 in cycle 0: with one cycle, the average IS the bin)",
                      fontsize=9, loc="left")
    ax2.set_xlabel("time  [s]")

    fig.tight_layout()
    fig.savefig(a.out, dpi=140)
    print(f"wrote {a.out}")

    if t_b.size and t_p.size:
        print(f"\ncell {C}, component U{comp_name}:")
        print(f"  instantaneous  range  {u_i.min():+.4f} .. {u_i.max():+.4f}")
        print(f"  U_Bin          range  {u_b.min():+.4f} .. {u_b.max():+.4f}")
        print(f"  U_BinPA        range  {u_p.min():+.4f} .. {u_p.max():+.4f}")
        if sd_at:
            sv = np.array(list(sd_at.values()))
            print(f"  cycle-to-cycle sd    mean {sv.mean():.4f}, max {sv.max():.4f}")
        if d_v:
            print(f"  |U_Bin - U_BinPA|    mean {np.mean(np.abs(d_v)):.4f}, "
                  f"max {np.max(np.abs(d_v)):.4f}")


if __name__ == "__main__":
    sys.exit(main())
