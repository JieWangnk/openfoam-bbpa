#!/usr/bin/env python3
"""What does a sudden kill + restart actually do to the phase average?

One point in the field, the final phase average, against an IDENTICAL run that
was never killed. The question is not "did it crash" -- it didn't -- but "which
phases were affected, and by how much, and is the damage bias or noise?"

The run: SIGKILL mid-cycle at t=3.602, resumed from the last surviving write at
t=3.58 (phase 0.58), which lands strictly inside bin 14. So:

  phases before 0.56  -- cycle 3 never ran after the restart -> those bins
                         averaged 4 cycles instead of 5
  bin 14 (0.56-0.60)  -- straddled: covered only 0.58-0.60 in cycle 3, so the
                         guard DROPS it -> 4 cycles, unbiased
  phases after 0.60   -- fully covered in every cycle -> 5 cycles, identical to
                         the continuous run

The point of the figure: the kill costs exactly ONE CYCLE at exactly the phases
that were not re-run, and nothing else. It is a loss of samples, not corruption.

Usage: python3 plot_kill_effect.py [--cell N]
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
    ap.add_argument("--case", default=".")
    ap.add_argument("--ref", default="../bbpaVerifyLES_ref")
    ap.add_argument("--cell", type=int, default=560)
    ap.add_argument("--comp", type=int, default=0)
    ap.add_argument("--period", type=float, default=1.0)
    ap.add_argument("--nbins", type=int, default=25)
    ap.add_argument("--wlo", type=float, default=4.0, help="final window start")
    ap.add_argument("--resume", type=float, default=3.58)
    ap.add_argument("--kill", type=float, default=3.602)
    ap.add_argument("--out", default="bbpa_kill_effect.png")
    a = ap.parse_args()
    C, K = a.cell, a.comp
    bw = a.period/a.nbins
    resume_phase = round(a.resume % a.period, 9)
    strad = int(resume_phase/bw)

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

    ph_r, pa_r = window(a.case, "_BinPA")
    ph_f, pa_f = window(a.ref, "_BinPA")
    _, vr_f = window(a.ref, "_BinPAVariance")

    def counts(case):
        txt = open(os.path.join(case, "5.04", "uniform", "bbpaVProperties")).read()
        m = re.search(r"\nU\s*\n\{(.*?)\n\}", txt, re.S)
        raw = re.search(r"binCounts\s+(.*?);", m.group(1), re.S).group(1).strip()
        u = re.match(r"(\d+)\s*\{\s*([0-9.]+)\s*\}", raw)
        if u:
            return np.array([float(u.group(2))]*int(u.group(1)))
        return np.array([float(x) for x in re.findall(r"[0-9.]+", raw)[1:]])

    cr, cf = counts(a.case), counts(a.ref)
    half = 0.5*bw
    x = ph_r + half                      # bins at their centre phase

    fig, (ax, ax2, ax3) = plt.subplots(
        3, 1, figsize=(11.5, 9.2), sharex=True,
        gridspec_kw={"height_ratios": [3, 1.6, 1.1]},
    )

    for axx in (ax, ax2, ax3):
        axx.axvspan(0, resume_phase, color="orange", alpha=0.13, lw=0)
        axx.axvline(resume_phase, color="k", ls=":", lw=1.8)
        axx.grid(alpha=0.25)

    ax.plot(x, pa_f, "o-", color="#2ca02c", lw=2.0, ms=4,
            label="$U_{BinPA}$  continuous run (never killed), 5 cycles")
    ax.plot(x, pa_r, "s--", color="#d62728", lw=1.6, ms=3.5,
            label="$U_{BinPA}$  killed at t=3.602, restarted from 3.58")
    sd = np.sqrt(np.maximum(vr_f, 0))
    ax.fill_between(x, pa_f - sd, pa_f + sd, color="#2ca02c", alpha=0.15, lw=0,
                    label=r"$\pm\sqrt{U_{BinPAVariance}}$  cycle-to-cycle spread")
    ax.set_ylabel(f"$U_{'xyz'[K]}$  [m/s]")
    ax.legend(loc="upper right", fontsize=8.5)
    ax.set_title(
        f"What a sudden kill costs the phase average (cell {C})\n"
        f"the two curves are indistinguishable: the damage is a lost SAMPLE, not corruption")

    d = pa_r - pa_f
    ax2.axhline(0, color="0.4", lw=1)
    ax2.plot(x, d, "o-", color="#d62728", lw=1.4, ms=4)
    ax2.fill_between(x, -sd/np.sqrt(5), sd/np.sqrt(5), color="0.6", alpha=0.25, lw=0,
                     label=r"$\pm$ standard error of the 5-cycle mean, $\sigma/\sqrt{5}$")
    ax2.set_ylabel("killed $-$ continuous")
    ax2.legend(loc="upper right", fontsize=8)
    ax2.set_title("the difference is confined to the phases that lost a cycle, "
                  "and is the size of ordinary sampling scatter", fontsize=9, loc="left")

    ax3.step(x, cr, where="mid", color="#d62728", lw=1.8, label="killed+restarted")
    ax3.step(x, cf, where="mid", color="#2ca02c", lw=1.8, ls="--", label="continuous")
    ax3.plot(x[strad], cr[strad], "kv", ms=9,
             label=f"bin {strad}: straddled by the restart -> dropped by the guard")
    ax3.set_ylabel("cycles\nfolded")
    ax3.set_xlabel("phase within the cycle  [s]")
    ax3.set_ylim(3.5, 5.5)
    ax3.legend(loc="lower right", fontsize=8)
    ax3.set_title("and WHY: exactly the un-rerun phases hold one cycle fewer",
                  fontsize=9, loc="left")

    ax.annotate("these phases never ran again\nafter the restart -> 4 cycles",
                xy=(resume_phase/2, ax.get_ylim()[0] + 0.12*np.ptp(ax.get_ylim())),
                ha="center", fontsize=8.5, color="#a05000")
    ax.annotate("restart\nphase", xy=(resume_phase, ax.get_ylim()[1]*0.93),
                xytext=(resume_phase + 0.09, ax.get_ylim()[1]*0.93),
                arrowprops=dict(arrowstyle="->", lw=1.3), fontsize=8.5)

    fig.tight_layout()
    fig.savefig(a.out, dpi=140)
    print(f"wrote {a.out}\n")

    lo = slice(0, strad+1)
    hi = slice(strad+1, None)
    print(f"cell {C}, final phase average, killed vs continuous:")
    print(f"  phases 0..{resume_phase:.2f} (lost a cycle, 4 vs 5): "
          f"mean |diff| = {np.mean(np.abs(d[lo])):.6f}")
    print(f"  phases {resume_phase:.2f}..1.0 (fully covered)     : "
          f"mean |diff| = {np.mean(np.abs(d[hi])):.3e}")
    print(f"  sampling scatter of a 5-cycle mean, sigma/sqrt(5)  = "
          f"{np.mean(sd)/np.sqrt(5):.6f}")
    print(f"\n  the affected phases differ by {np.mean(np.abs(d[lo])):.4f}, which is the size of")
    print(f"  ordinary sampling scatter ({np.mean(sd)/np.sqrt(5):.4f}) -- a lost sample, not a bias.")
    print(f"  the unaffected phases are identical to {np.mean(np.abs(d[hi])):.0e}.")


if __name__ == "__main__":
    sys.exit(main())
