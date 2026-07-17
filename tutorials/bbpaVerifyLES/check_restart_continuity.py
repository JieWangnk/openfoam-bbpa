#!/usr/bin/env python3
"""Did the phase average SURVIVE an accidental kill, or silently start over?

The run was SIGKILLed mid-cycle at t~3.5 and restarted from its last write.
This tracks ONE point through the whole history and asks two questions:

  1. CONTINUITY -- does U_BinPA after the restart carry the pre-kill cycles, or
     did it reset to a single-cycle value? The failure mode is not a crash: a
     reset accumulator produces a plausible-looking field that is silently the
     average of far fewer cycles than the counters claim. That is exactly the
     class of bug that has bitten this library before (bin-window corruption,
     WSS cold-start, boundary drop on restore), and every one of them was
     invisible without a check like this.

  2. COST -- what did the kill actually lose? Compared against an identical
     CONTINUOUS run (no kill), the restarted case must differ, because the
     3.0->3.5 partial cycle was never checkpointed and is genuinely gone. The
     honest question is not "are they identical" (they cannot be) but "does the
     restarted case average exactly the cycles it actually saw".

Usage: python3 check_restart_continuity.py [--cell N] [--ref ../bbpaVerifyLES_ref]
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
    return a[cell, comp]


def windows(case, period, nbins, suffix):
    """-> {cycle_window_start: {bin_time: value}} for files named U<suffix>."""
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
    ap.add_argument("--nbins", type=int, default=25)
    a = ap.parse_args()
    C, K = a.cell, a.comp

    print("=" * 74)
    print(f"RESTART CONTINUITY at cell {C}, component U{'xyz'[K]}")
    print("=" * 74)

    for label, case in (("RESTARTED (killed at ~3.5)", a.case), ("REFERENCE (continuous)", a.ref)):
        if not os.path.isdir(case):
            print(f"\n{label}: MISSING ({case})")
            continue
        print(f"\n--- {label} ---")
        pa = windows(case, a.period, a.nbins, "_BinPA")
        vr = windows(case, a.period, a.nbins, "_BinPAVariance")
        bn = windows(case, a.period, a.nbins, "_Bin")
        print(f"  cycle windows holding _BinPA: {sorted(pa)}")
        print(f"  bins per window             : "
              f"{ {w: len(v) for w, v in sorted(pa.items())} }")
        print(f"  windows holding _BinPAVariance: {sorted(vr)}")
        d, cnts = counts_from_props(case)
        if cnts:
            u = sorted(set(cnts))
            print(f"  final binCounts (from {d}): "
                  f"{'all ' + str(u[0]) if len(u) == 1 else f'{min(cnts)}..{max(cnts)} (mixed)'}")

    # ---- the continuity question, at one point, window by window ----
    print("\n" + "=" * 74)
    print("Does U_BinPA carry its history across the kill?")
    print("=" * 74)
    print(f"{'window':>8} {'bin t':>8} {'U_BinPA':>12} {'sqrt(var)':>12}   interpretation")
    pa = windows(a.case, a.period, a.nbins, "_BinPA")
    vr = windows(a.case, a.period, a.nbins, "_BinPAVariance")
    probe_bin = 18                       # a mid-cycle bin
    for w in sorted(pa):
        bt = sorted(pa[w])
        if probe_bin >= len(bt):
            continue
        key = bt[probe_bin]
        t, d = pa[w][key]
        v = val(os.path.join(a.case, d, "U_BinPA"), C, K)
        sd = ""
        if w in vr and key in vr[w]:
            tv, dv = vr[w][key]
            sd = f"{np.sqrt(max(val(os.path.join(a.case, dv, 'U_BinPAVariance'), C, K), 0)):12.6f}"
        else:
            sd = f"{'-':>12}"
        note = "variance absent => only 1 cycle folded" if not (w in vr and key in vr[w]) else ""
        print(f"{w:>8} {t:>8.2f} {v:>12.6f} {sd}   {note}")

    print("\nIf the restart had RESET the accumulator, the window written after the")
    print("kill would show a variance of '-' (one cycle) and a U_BinPA equal to that")
    print("single cycle's bin. If it CONTINUED, the variance is present and the")
    print("phase average keeps evolving smoothly from the pre-kill history.")


if __name__ == "__main__":
    sys.exit(main())
