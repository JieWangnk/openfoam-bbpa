#!/usr/bin/env python3
"""Three runs, one accident, two ways of restarting from it.

All three see the SAME SIGKILL, which left the last complete write at t=3.58 --
phase 0.58, strictly inside bin 14 (bin width 0.04, so 0.58/0.04 = 14.5).

  reference   never killed at all                       -- the truth
  latestTime  restarted from 3.58, the naive choice     -- bin 14 is straddled:
              it can only cover 0.58-0.60 of its window, so the guard drops it
              and that bin averages 4 cycles instead of 5
  binBoundary restarted from 3.56, the bin's own edge    -- bin 14 covers
              0.56-0.60 in full, nothing is straddled, nothing is dropped

The point: rewinding to the bin boundary costs 0.02 s of re-run (one bin width,
20 timesteps) and buys back a whole cycle at that phase. Co-location guarantees
this is always possible: if binWidth is a multiple of writeInterval then every
bin boundary IS a solver write time, so a complete state always exists there.

Usage: python3 plot_restart_compare.py [--cell N]
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

NC = 1440


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default="../bbpaVerifyLES_ref")
    ap.add_argument("--mid", default="../bbpaVerifyLES_mid")
    ap.add_argument("--edge", default=".")
    ap.add_argument("--cell", type=int, default=560)
    ap.add_argument("--comp", type=int, default=0)
    ap.add_argument("--period", type=float, default=1.0)
    ap.add_argument("--nbins", type=int, default=25)
    ap.add_argument("--wlo", type=float, default=4.0)
    ap.add_argument("--strad", type=int, default=14)
    ap.add_argument("--out", default="bbpa_restart_compare.png")
    a = ap.parse_args()
    C, K = a.cell, a.comp
    bw = a.period/a.nbins

    def val(p):
        arr = read_internal(p)
        if arr.shape[0] == 1:
            arr = np.repeat(arr, NC, axis=0)
        return arr[C, K]

    def window(case, suf):
        d = {}
        for t, dd in time_dirs(case):
            p = os.path.join(case, dd, f"U{suf}")
            if os.path.exists(p) and a.wlo - 1e-9 <= t < a.wlo + a.period - 1e-9:
                d[round(t, 9)] = p
        ts = sorted(d)
        return np.array(ts) - a.wlo, np.array([val(d[t]) for t in ts])

    def counts(case):
        txt = open(os.path.join(case, "5.04", "uniform", "bbpaVProperties")).read()
        m = re.search(r"\nU\s*\n\{(.*?)\n\}", txt, re.S)
        raw = re.search(r"binCounts\s+(.*?);", m.group(1), re.S).group(1).strip()
        u = re.match(r"(\d+)\s*\{\s*([0-9.]+)\s*\}", raw)
        if u:
            return np.array([float(u.group(2))]*int(u.group(1)))
        return np.array([float(x) for x in re.findall(r"[0-9.]+", raw)[1:]])

    cases = {}
    for key, path in (("ref", a.ref), ("mid", a.mid), ("edge", a.edge)):
        if not os.path.isdir(path):
            print(f"[skip] {key}: {path} missing")
            continue
        ph, pa = window(path, "_BinPA")
        if not len(ph):
            print(f"[skip] {key}: no _BinPA in the final window")
            continue
        cases[key] = dict(ph=ph, pa=pa, cnt=counts(path))
    _, vr = window(a.ref, "_BinPAVariance")
    sd = np.sqrt(np.maximum(vr, 0))

    if "ref" not in cases:
        print("no reference; nothing to compare"); return 1

    x = cases["ref"]["ph"] + 0.5*bw          # bins at their centre phase
    xs = a.strad*bw + 0.5*bw                 # the straddled bin's centre

    fig, (ax, ax2, ax3) = plt.subplots(
        3, 1, figsize=(11.5, 9.4), sharex=True,
        gridspec_kw={"height_ratios": [2.6, 1.9, 1.1]},
    )
    for axx in (ax, ax2, ax3):
        axx.axvline(xs, color="k", ls=":", lw=1.6)
        axx.grid(alpha=0.25)

    ax.plot(x, cases["ref"]["pa"], "o-", color="#2ca02c", lw=2.2, ms=4,
            label="never killed  (the truth)")
    if "mid" in cases:
        ax.plot(x, cases["mid"]["pa"], "s--", color="#d62728", lw=1.5, ms=3.5,
                label="restart from latestTime 3.58  (mid-bin: bin 14 dropped)")
    if "edge" in cases:
        ax.plot(x, cases["edge"]["pa"], "^:", color="#1f77b4", lw=1.5, ms=3.5,
                label="restart from bin boundary 3.56  (nothing lost)")
    ax.fill_between(x, cases["ref"]["pa"] - sd, cases["ref"]["pa"] + sd,
                    color="#2ca02c", alpha=0.13, lw=0,
                    label=r"$\pm\sqrt{U_{BinPAVariance}}$")
    ax.set_ylabel(f"$U_{'xyz'[K]}$  [m/s]")
    ax.set_title("Same SIGKILL, two ways to restart from it "
                 f"(cell {C}, final phase average)\n"
                 "rewinding to the bin boundary costs 0.02 s of re-run "
                 "and loses nothing")
    ax.legend(loc="upper right", fontsize=8.5)

    ax2.axhline(0, color="0.4", lw=1)
    if "mid" in cases:
        ax2.plot(x, cases["mid"]["pa"] - cases["ref"]["pa"], "s--",
                 color="#d62728", lw=1.5, ms=4,
                 label="latestTime $-$ truth")
    if "edge" in cases:
        ax2.plot(x, cases["edge"]["pa"] - cases["ref"]["pa"], "^:",
                 color="#1f77b4", lw=1.5, ms=4,
                 label="bin boundary $-$ truth  (flat zero)")
    ax2.set_ylabel("difference from truth")
    ax2.legend(loc="upper right", fontsize=8)
    ax2.set_title("the naive restart deviates at exactly one phase; "
                  "the rewound restart is identically zero everywhere",
                  fontsize=9, loc="left")

    ax3.plot(x, cases["ref"]["cnt"], "o-", color="#2ca02c", lw=2, ms=3)
    if "mid" in cases:
        ax3.plot(x, cases["mid"]["cnt"], "s--", color="#d62728", lw=1.5, ms=4)
    if "edge" in cases:
        ax3.plot(x, cases["edge"]["cnt"], "^:", color="#1f77b4", lw=1.5, ms=4)
    ax3.set_ylabel("cycles\nfolded")
    ax3.set_ylim(3.5, 5.5)
    ax3.set_xlabel("phase within the cycle  [s]")
    ax3.set_title("and why: the naive restart credits bin 14 only 4 cycles",
                  fontsize=9, loc="left")
    ax3.annotate("bin 14\n(straddled)", xy=(xs, 4.0), xytext=(xs + 0.09, 3.75),
                 arrowprops=dict(arrowstyle="->", lw=1.2), fontsize=8.5)

    fig.tight_layout()
    fig.savefig(a.out, dpi=140)
    print(f"wrote {a.out}\n")

    print(f"cell {C}, final phase average vs the never-killed truth:")
    for key, lab in (("mid", "restart from latestTime 3.58 (mid-bin)"),
                     ("edge", "restart from bin boundary 3.56   ")):
        if key not in cases:
            continue
        d = np.abs(cases[key]["pa"] - cases["ref"]["pa"])
        print(f"  {lab}: max |diff| = {d.max():.3e}  at bin {int(d.argmax())}")
        print(f"      {'':38s} bins losing a cycle: "
              f"{int((cases[key]['cnt'] < cases['ref']['cnt']).sum())}")


if __name__ == "__main__":
    sys.exit(main())
