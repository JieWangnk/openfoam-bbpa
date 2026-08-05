#!/usr/bin/env python3
"""Did the phase average SURVIVE a restart, or silently start over?

A run was stopped mid-study and continued from its last write (foamRun with
startFrom latestTime, or a recovery after a kill). This tracks ONE point through
the whole history and asks two questions:

  1. CONTINUITY -- does the phase-mean field after the restart carry the
     pre-restart cycles, or did it reset to a single-cycle value? The failure
     mode is not a crash: a reset accumulator produces a plausible-looking field
     that is silently the average of far fewer cycles than the counters claim.
     That is exactly the class of bug that has bitten this library before
     (bin-window corruption, WSS cold-start, boundary drop on restore), and
     every one of them was invisible without a check like this.

     The tell is cross-cycle variance. It is undefined for N<2 cycles, so it can
     only appear once >=2 cycles are folded in. If it is present at the first
     cycle completed AFTER the restart, a pre-restart cycle must have been
     carried across.

  2. COST -- what did the restart lose? Compared against an identical CONTINUOUS
     run (no interruption, pass with --ref), the restarted case must differ,
     because the partial cycle in progress at the stop was never checkpointed and
     is genuinely gone. The honest question is not "are they identical" (they
     cannot be) but "does the restarted case average exactly the cycles it saw".

Field names: this OF14 port writes the phase mean as U_PA, the cross-cycle
variance magnitude as U_PAVarianceMag, and the current-cycle field as U_Bin.
Older openfoam-bbpa cases used U_BinPA / U_BinPAVariance -- point --pa-suffix and
--var-suffix at those to run this check there.

Usage: python3 check_restart_continuity.py [--cell N] [--ref ../bbpaVerifyLES_ref]
       python3 check_restart_continuity.py --pa-suffix _BinPA --var-suffix _BinPAVariance
"""
import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_bbpa import read_internal, time_dirs

NC = 1440


def val(path, cell, comp):
    a = read_internal(path)
    if a.shape[0] == 1:
        a = np.repeat(a, NC, axis=0)
    comp = min(comp, a.shape[1] - 1)     # scalar magnitude fields have one component
    return a[cell, comp]


def windows(case, period, suffix):
    """-> {cycle_window: {write_time: (t, dir)}} for files named U<suffix>."""
    out = {}
    for t, d in time_dirs(case):
        p = os.path.join(case, d, f"U{suffix}")
        if not os.path.exists(p):
            continue
        w = int(np.floor(t/period + 1e-9))
        out.setdefault(w, {})[round(t, 9)] = (t, d)
    return out


def counts_from_props(case):
    """binCounts for U from the newest Properties dict."""
    for t, d in reversed(time_dirs(case)):
        p = os.path.join(case, d, "uniform", "bbpaVProperties")
        if not os.path.exists(p):
            continue
        txt = open(p).read()
        m = re.search(r"\nU\s*\n\{(.*?)\n\}", txt, re.S)
        if not m:
            continue
        body = m.group(1)
        mc = re.search(r"binCounts\s+(.*?);", body, re.S)
        raw = mc.group(1).strip() if mc else ""
        mu = re.match(r"(\d+)\s*\{\s*([0-9.]+)\s*\}", raw)
        if mu:                       # uniform form: "25 { 5 }"
            return d, [float(mu.group(2))]*int(mu.group(1))
        nums = re.findall(r"[0-9.]+", raw)
        if nums:
            return d, [float(x) for x in nums[1:]] if len(nums) > 1 else []
    return None, []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=".")
    ap.add_argument("--ref", default="../bbpaVerifyLES_ref")
    ap.add_argument("--cell", type=int, default=560)
    ap.add_argument("--comp", type=int, default=0)
    ap.add_argument("--period", type=float, default=1.0)
    ap.add_argument("--pa-suffix", default="_PA",
                    help="phase-mean field suffix (U_PA here; U_BinPA on old cases)")
    ap.add_argument("--var-suffix", default="_PAVarianceMag",
                    help="cross-cycle variance suffix (U_PAVarianceMag here)")
    a = ap.parse_args()
    C, K = a.cell, a.comp
    PA, VAR = a.pa_suffix, a.var_suffix

    print("=" * 74)
    print(f"RESTART CONTINUITY at cell {C}, component U{'xyz'[K]}")
    print("=" * 74)

    for label, case in (("RESTARTED", a.case), ("REFERENCE (continuous)", a.ref)):
        if not os.path.isdir(case):
            print(f"\n{label}: MISSING ({case})")
            continue
        print(f"\n--- {label} ---")
        pa = windows(case, a.period, PA)
        vr = windows(case, a.period, VAR)
        print(f"  cycle windows holding U{PA}: {sorted(pa)}")
        print(f"  writes per window          : "
              f"{ {w: len(v) for w, v in sorted(pa.items())} }")
        print(f"  windows holding U{VAR}: {sorted(vr)}")
        d, cnts = counts_from_props(case)
        if cnts:
            u = sorted(set(cnts))
            print(f"  final binCounts (from {d}): "
                  f"{'all ' + str(u[0]) if len(u) == 1 else f'{min(cnts)}..{max(cnts)} (mixed)'}")

    # ---- the continuity question, at one point, window by window ----
    print("\n" + "=" * 74)
    print(f"Does U{PA} carry its history across the restart?")
    print("=" * 74)
    print(f"{'window':>8} {'t':>8} {f'U{PA}[{chr(120+K)}]':>12} {'sqrt(varMag)':>13}   interpretation")
    pa = windows(a.case, a.period, PA)
    vr = windows(a.case, a.period, VAR)
    probe = 18                       # a fixed mid-cycle phase, same index each window
    for w in sorted(pa):
        bt = sorted(pa[w])
        if probe >= len(bt):
            continue
        key = bt[probe]
        t, d = pa[w][key]
        v = val(os.path.join(a.case, d, f"U{PA}"), C, K)
        has_var = w in vr and key in vr[w]
        if has_var:
            tv, dv = vr[w][key]
            sd = f"{np.sqrt(max(val(os.path.join(a.case, dv, f'U{VAR}'), C, 0), 0)):13.6f}"
            note = ">=2 cycles folded"
        else:
            sd = f"{'-':>13}"
            note = "variance absent => only 1 cycle folded"
        print(f"{w:>8} {t:>8.2f} {v:>12.6f} {sd}   {note}")

    print(f"\nIf the restart had RESET the accumulator, the window written after it")
    print(f"would show a variance of '-' (one cycle) and a U{PA} equal to that single")
    print(f"cycle's phase value. If it CONTINUED, cross-cycle variance is present from")
    print(f"the first post-restart cycle and the phase mean evolves smoothly from the")
    print(f"pre-restart history.")


if __name__ == "__main__":
    sys.exit(main())
