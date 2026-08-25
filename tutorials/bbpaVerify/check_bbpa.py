#!/usr/bin/env python3
"""Independent verification of the BBPA function object, from its own output.

The case keeps the per-cycle capture on (fields ( U { writePerCycle on; } )),
so every folded cycle leaves its OWN bin means on disk: U_Bin and U_Bin_UU at
that cycle's label times. The cross-cycle files (U_PA, U_PA_UU,
U_PAVarianceMag) must then be exact functions of those per-cycle files:

  1. U_PA at the final window  == the plain mean of U_Bin over the folded
     cycles (checked per bin, per cell, at round-off).
  2. U_PA_UU                   == the plain mean of U_Bin_UU likewise.
  3. U_PAVarianceMag           == the cross-cycle sample variance of the
     per-cycle bin means: sum over components of
     sum_c (U_Bin_c - mean)^2 / (nCycles - 1).   (needs >= 2 cycles)
  4. At the FIRST folded cycle's window, U_PA == U_Bin identically (an
     ensemble of one).

These are re-derived here in numpy, independently of the C++ Welford chain,
so a folding, variance or naming bug shows up as a mismatch beyond round-off
(the tolerance is set by writePrecision in controlDict).

They cannot see a mis-placed bin WINDOW (the per-cycle files would be wrong
consistently), so a fifth check is physical:

  5. U_PA at each bin label is close to the instantaneous U saved at that
     same label time, averaged over the folded cycles. The bins are centred
     on their labels, so the two agree to O((2*pi*binDeltaT/period)^2) for a
     smooth flow -- percent level here, against ~30% if the window sat one
     half-bin off.

The machine-precision check of the time integration itself (every step, exact
overlap weights, variable dt) lives in ../../tests/ (./Alltest): it needs the
full per-step trajectory, which a tutorial case does not keep.

Usage:  python3 check_bbpa.py [--case .] [--period 1.0] [--nbins 10]
                              [--tol 1e-8] [--place-tol 0.05]
"""
import argparse
import os
import re
import sys

import numpy as np


def read_internal(path):
    """Parse an OpenFOAM ascii volField internalField -> (n,c) array.
    Handles vector (3), symmTensor (6) and scalar (1), uniform or nonuniform."""
    txt = open(path).read()
    m = re.search(r"internalField\s+uniform\s*\(([^)]*)\)", txt)
    if m:
        return np.array([[float(x) for x in m.group(1).split()]])
    m = re.search(r"internalField\s+uniform\s+([-0-9.eE+]+)\s*;", txt)
    if m:
        return np.array([[float(m.group(1))]])
    m = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s*\n?\s*(\d+)\s*\(", txt)
    if not m:
        raise RuntimeError(f"cannot parse internalField in {path}")
    kind, n = m.group(1), int(m.group(2))
    body = txt[m.end():]
    chunk = body[: body.index("\n)")+2] if "\n)" in body else body
    if kind == "scalar":
        vals = re.findall(r"[-0-9.eE+]+", chunk)
        arr = np.array([[float(v)] for v in vals[:n]])
    else:
        vals = re.findall(r"\(([^)]*)\)", chunk)
        arr = np.array([[float(x) for x in v.split()] for v in vals[:n]])
    if arr.shape[0] != n:
        raise RuntimeError(f"{path}: expected {n} entries, parsed {arr.shape[0]}")
    return arr


def time_dirs(case):
    out = []
    for d in os.listdir(case):
        try:
            t = float(d)
        except ValueError:
            continue
        if os.path.isdir(os.path.join(case, d)):
            out.append((t, d))
    return sorted(out)


def find_dir(tds, t):
    for tv, d in tds:
        if abs(tv - t) < 1e-7:
            return d
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=".")
    ap.add_argument("--period", type=float, default=1.0)
    ap.add_argument("--nbins", type=int, default=10)
    ap.add_argument("--field", default="U")
    ap.add_argument("--tol", type=float, default=1e-8,
                    help="tolerance for the exact identities (round-off, "
                         "limited by writePrecision)")
    ap.add_argument("--place-tol", type=float, default=0.05,
                    help="window-placement tolerance, relative to the "
                         "window's peak |U_PA|")
    a = ap.parse_args()
    case = os.path.abspath(a.case)
    F = a.field

    tds = time_dirs(case)
    if not tds:
        print("[err] no time dirs"); return 1

    # ---- the FO's own bookkeeping: phase origin, folded-cycle count --------
    label_origin, ncycles = None, None
    for t, d in reversed(tds):
        pd = os.path.join(case, d, "uniform")
        if not os.path.isdir(pd):
            continue
        for pf in os.listdir(pd):
            if pf.endswith("Properties"):
                txt = open(os.path.join(pd, pf)).read()
                m = re.search(r"startTime\s+([0-9.eE+-]+);", txt)
                if m and label_origin is None:
                    label_origin = float(m.group(1))
                m = re.search(r"nCycles\s+(\d+);", txt)
                if m and ncycles is None:
                    ncycles = int(m.group(1))
        if label_origin is not None:
            break
    if label_origin is None:
        print("[err] no *Properties found; cannot get the phase origin")
        return 1
    print(f"[fo]  label origin {label_origin}  folded cycles {ncycles}")

    dT = a.period / a.nbins

    # ---- which cycles left a full per-cycle window on disk -----------------
    # A folded cycle c has U_Bin in every label dir  c*period + j*dT.
    # The discarded startup cycle(s) and the still-open cycle leave nothing
    # (finalised-only products), so these windows ARE the folded ensemble.
    kmax = int(max(t for t, _ in tds) / a.period) + 1
    folded = []
    for c in range(kmax + 1):
        dirs = [find_dir(tds, label_origin + c*a.period + j*dT)
                for j in range(a.nbins)]
        if all(d is not None and
               os.path.isfile(os.path.join(case, d, F + "_Bin"))
               for d in dirs):
            folded.append((c, dirs))
    if not folded:
        print(f"[err] no complete {F}_Bin window on disk"); return 1
    print(f"[truth] per-cycle windows on disk: cycles "
          f"{[c for c, _ in folded]}")
    if ncycles is not None and len(folded) != ncycles:
        print(f"[err] {len(folded)} per-cycle windows but the FO folded "
              f"{ncycles} cycles"); return 1

    first_dirs = folded[0][1]
    final_dirs = folded[-1][1]
    k = len(folded)

    # ---- the five checks, bin by bin ---------------------------------------
    print(f"\n{'bin':>3} {'t_dir':>8} {'_PA':>12} {'_PA_UU':>12} "
          f"{'_PAVarMag':>12} {'firstPA':>12} {'placement':>12}")
    worst = {"_PA": 0.0, "_PA_UU": 0.0, "_PAVarMag": 0.0,
             "firstPA": 0.0, "placement": 0.0}
    nbad = 0

    # global scale for the placement check: peak |U_PA| over the window
    pa_all = [read_internal(os.path.join(case, d, F + "_PA"))
              for d in final_dirs]
    scale = max(np.max(np.abs(p)) for p in pa_all)

    for j in range(a.nbins):
        row = {}
        bins = [read_internal(os.path.join(case, dirs[j], F + "_Bin"))
                for _, dirs in folded]
        mean = np.mean(bins, axis=0)

        # 1. cross-cycle mean
        row["_PA"] = np.max(np.abs(pa_all[j] - mean))

        # 2. cross-cycle mean of the per-cycle second moments
        p_uu = os.path.join(case, final_dirs[j], F + "_PA_UU")
        if os.path.isfile(p_uu):
            uu_bins = [read_internal(os.path.join(case, dirs[j], F + "_Bin_UU"))
                       for _, dirs in folded]
            row["_PA_UU"] = np.max(np.abs(read_internal(p_uu)
                                          - np.mean(uu_bins, axis=0)))

        # 3. cross-cycle sample variance of the per-cycle bin means
        p_v = os.path.join(case, final_dirs[j], F + "_PAVarianceMag")
        if k > 1 and os.path.isfile(p_v):
            m2 = np.sum([np.sum((b - mean)**2, axis=1) for b in bins], axis=0)
            row["_PAVarMag"] = np.max(np.abs(
                read_internal(p_v)[:, 0] - m2 / (k - 1.0)))
        elif k > 1:
            print(f"[err] {p_v} missing"); return 1

        # 4. ensemble of one: first folded window's _PA == its _Bin
        row["firstPA"] = np.max(np.abs(
            read_internal(os.path.join(case, first_dirs[j], F + "_PA"))
            - bins[0]))

        # 5. placement: _PA vs the label-time instantaneous samples
        inst = [read_internal(os.path.join(case, dirs[j], F))
                for _, dirs in folded]
        row["placement"] = np.max(np.abs(pa_all[j] - np.mean(inst, axis=0))) \
            / scale

        bad = (any(row.get(c, 0.0) > a.tol
                   for c in ("_PA", "_PA_UU", "_PAVarMag", "firstPA"))
               or row["placement"] > a.place_tol)
        nbad += 1 if bad else 0
        for c in worst:
            worst[c] = max(worst[c], row.get(c, 0.0))
        cells = " ".join(f"{row.get(c, float('nan')):12.3e}"
                         for c in ("_PA", "_PA_UU", "_PAVarMag",
                                   "firstPA", "placement"))
        print(f"{j:>3} {final_dirs[j]:>8} {cells}"
              f"{'   <-- MISMATCH' if bad else ''}")

    print()
    for c in ("_PA", "_PA_UU", "_PAVarMag", "firstPA"):
        print(f"  worst |FO - independent| {c:<10} = {worst[c]:.3e}"
              f"  (tol {a.tol:g})")
    print(f"  worst placement (rel. to peak)      = {worst['placement']:.3e}"
          f"  (tol {a.place_tol:g})")
    if nbad == 0:
        print("PASS: the cross-cycle mean, second moment and variance all\n"
              "      reproduce an independent re-derivation from the\n"
              "      per-cycle _Bin files, and every bin window sits where\n"
              "      its label says.")
        return 0
    print(f"FAIL: {nbad} bin(s) disagree beyond tolerance.")
    return 2


if __name__ == "__main__":
    sys.exit(main())
