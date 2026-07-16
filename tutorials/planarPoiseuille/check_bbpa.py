#!/usr/bin/env python3
"""Independent verification of the BBPA function object.

Strategy: the case saves the FULL instantaneous U at EVERY timestep, so we have
the exact trajectory the FO saw.  We re-derive the bin averages here from that
trajectory, with logic written independently of the C++, and compare against the
FO's own U_PA output.  A bug in the C++ binning (wrong window, wrong time
weight, wrong cycle folding, a bin silently left at zero) shows up as a mismatch.

This is a same-trajectory identity check, so agreement should be to round-off
(~1e-10), not "close enough".  Anything larger is a real defect.

BBPA definition being checked (per the source):
  bin j of cycle n covers phase [j*dT_bin, (j+1)*dT_bin), dT_bin = period/nBins,
  phase measured from startTime_.  Within a cycle the bin value is the
  TIME-WEIGHTED mean  sum(U*dt)/sum(dt)  over the steps landing in that bin.
  Across cycles the reported U_PA[j] is the plain mean of the per-cycle bin
  values (Welford), i.e. mean over the nCycles that contributed.

Usage:  python3 check_bbpa.py [--case .] [--period 1.0] [--nbins 10]
                             [--start-time auto] [--cycles all] [--tol 1e-8]
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


# back-compat alias
read_internal_vector = read_internal


def sqr_sym(U):
    """<U_i U_j> as OpenFOAM's symmTensor component order (xx xy xz yy yz zz)."""
    x, y, z = U[:, 0], U[:, 1], U[:, 2]
    return np.stack([x*x, x*y, x*z, y*y, y*z, z*z], axis=1)


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=".")
    ap.add_argument("--period", type=float, default=1.0)
    ap.add_argument("--nbins", type=int, default=10)
    ap.add_argument("--start-time", default="auto",
                    help="BBPA phase origin; 'auto' reads it from *Properties")
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--start-cycle", type=int, default=0,
                    help="the FO's startCycle: cycles below this are not folded")
    ap.add_argument("--include-t0", action="store_true",
                    help="fold the initial condition with a full dt, as the FO does "
                         "if its execute() is called at the first time")
    a = ap.parse_args()
    case = os.path.abspath(a.case)

    tds = time_dirs(case)
    if not tds:
        print("[err] no time dirs"); return 1

    # ---- the FO's own bookkeeping (phase origin + how many cycles it folded)
    start_time, ncycles = None, None
    for t, d in reversed(tds):
        for f in ("uniform",):
            pd = os.path.join(case, d, f)
            if not os.path.isdir(pd):
                continue
            for pf in os.listdir(pd):
                if pf.endswith("Properties"):
                    txt = open(os.path.join(pd, pf)).read()
                    m = re.search(r"startTime\s+([0-9.eE+-]+);", txt)
                    if m and start_time is None:
                        start_time = float(m.group(1))
                    m = re.search(r"nCycles\s+(\d+);", txt)
                    if m and ncycles is None:
                        ncycles = int(m.group(1))
        if start_time is not None:
            break
    if a.start_time != "auto":
        start_time = float(a.start_time)
    if start_time is None:
        print("[err] no *Properties found; cannot get the phase origin"); return 1
    print(f"[fo]  phase origin startTime={start_time}  nCycles={ncycles}")

    dT = a.period / a.nbins

    # ---- ground truth: every saved instantaneous U (dirs holding a real 'p') --
    traj = [(t, d) for (t, d) in tds
            if os.path.exists(os.path.join(case, d, "U"))
            and os.path.exists(os.path.join(case, d, "p"))]
    if len(traj) < 10:
        print(f"[err] only {len(traj)} instantaneous states; need the per-step save"); return 1
    print(f"[truth] {len(traj)} per-step states from t={traj[0][0]} to {traj[-1][0]}")

    # dt for each state = t_k - t_{k-1}  (matches the FO: it adds dt of the step
    # that has just been advanced TO this time)
    times = np.array([t for t, _ in traj])
    # accumulate: for each state k>0 the FO accumulated U(t_k)*dt_k into the bin
    # of phase(t_k)
    sums, sqsums, wts = {}, {}, {}
    # The FO's execute() runs at every time it is called and does
    #   accumulate(phaseBin(), field, time_.deltaT())
    # If it is called at the very first time (the initial condition), it folds
    # U(t0) in with a full dt even though no time has elapsed there. Whether
    # that happens is a property of the solver's FO hook, so make it switchable
    # and let the data decide rather than guessing.
    k0 = 0 if a.include_t0 else 1
    for k in range(k0, len(traj)):
        t = times[k]
        dt = (times[k] - times[k-1]) if k > 0 else (times[1] - times[0])
        if t < start_time - 1e-12:
            continue
        ph = (t - start_time) % a.period
        j = min(int(ph / dT), a.nbins - 1)
        n = int((t - start_time) // a.period)          # cycle index
        U = read_internal(os.path.join(case, traj[k][1], "U"))
        key = (n, j)
        sums[key] = sums.get(key, 0.0) + U * dt
        # raw second moment <U_i U_j>: the FO accumulates sqr(U)*dt in the SAME
        # bin, so re-derive it the same way. R_ij^bin and k_bin come from this.
        sqsums[key] = sqsums.get(key, 0.0) + sqr_sym(U) * dt
        wts[key] = wts.get(key, 0.0) + dt

    # Per-cycle bin means, then the cross-cycle mean = what U_PA should hold.
    #
    # CRITICAL: compare against exactly the cycles the FO actually FOLDED, not
    # every cycle present on disk.  A cycle is folded only when the solver
    # CROSSES its upper boundary (finalizeAllForCycle runs on the cycle
    # transition), so the final, still-open cycle is NOT in U_PA -- and with
    # floating-point drift the run may stop a hair below the boundary and leave
    # the last cycle unfolded too.  Averaging cycles the FO never folded makes
    # a correct FO look broken.
    present = sorted({n for (n, j) in sums})
    if ncycles is None:
        print("[err] no nCycles in Properties; cannot know which cycles were folded")
        return 1
    present = [n for n in present if n >= a.start_cycle]
    cycles = present[:ncycles]          # folded = the first nCycles from startCycle
    print(f"[truth] cycles on disk: {present} ; FO folded {ncycles} -> comparing {cycles}")
    expect, expect_uu, expect_var = {}, {}, {}
    for j in range(a.nbins):
        ok = [n for n in cycles if (n, j) in wts and wts[(n, j)] > 0]
        if not ok:
            continue
        # per-cycle bin means / bin second moments
        cyc_mean = [sums[(n, j)] / wts[(n, j)] for n in ok]
        cyc_uu   = [sqsums[(n, j)] / wts[(n, j)] for n in ok]
        expect[j]    = (np.mean(cyc_mean, axis=0), len(ok))   # <name>_PA
        expect_uu[j] = (np.mean(cyc_uu, axis=0), len(ok))     # <name>_PA_UU  = mean over cycles of <U_iU_j>_bin
        # <name>_PAVariance = M2/(n-1) where M2 = sum |x - mean|^2 of the per-cycle
        # BIN MEANS (Welford over cycles) -> the cross-cycle sample variance
        if len(ok) > 1:
            m = np.mean(cyc_mean, axis=0)
            m2 = np.sum([np.sum((v - m)**2, axis=1) for v in cyc_mean], axis=0)
            expect_var[j] = (m2 / (len(ok) - 1.0), len(ok))

    # ---- the FO's answer: U_PA in the phase-aligned dirs of the LAST written window
    # find the newest dir set holding U_PA
    pa = [(t, d) for (t, d) in tds if os.path.exists(os.path.join(case, d, "U_PA"))]
    if not pa:
        print("[err] no U_PA written"); return 1
    # bins of the most recent window: the newest nbins U_PA dirs
    latest = pa[-a.nbins:]
    print(f"[fo]  U_PA window: t={latest[0][0]} .. {latest[-1][0]}  ({len(latest)} bins)")

    print(f"\n{'bin':>3} {'t_dir':>16} {'_PA':>12} {'_PA_UU':>12} {'_PAVariance':>13}")
    worst = {"_PA": 0.0, "_PA_UU": 0.0, "_PAVariance": 0.0}
    nbad = 0
    for i, (t, d) in enumerate(latest):
        # bin index = POSITION in the window (see note above)
        j = i
        row = {}
        # --- first moment
        fo = read_internal(os.path.join(case, d, "U_PA"))
        ex = expect[j][0]
        if ex.shape[0] == 1 and fo.shape[0] > 1:
            ex = np.repeat(ex, fo.shape[0], axis=0)
        row["_PA"] = np.max(np.abs(fo - ex))
        # --- raw second moment <U_i U_j>  (R_ij and k_bin come from this)
        p_uu = os.path.join(case, d, "U_PA_UU")
        if os.path.exists(p_uu) and j in expect_uu:
            fo_uu = read_internal(p_uu)
            ex_uu = expect_uu[j][0]
            row["_PA_UU"] = np.max(np.abs(fo_uu - ex_uu))
        # --- cross-cycle variance of the bin mean
        p_v = os.path.join(case, d, "U_PAVariance")
        if os.path.exists(p_v) and j in expect_var:
            fo_v = read_internal(p_v)
            ex_v = expect_var[j][0].reshape(-1, 1)
            row["_PAVariance"] = np.max(np.abs(fo_v - ex_v))
        bad = any(v > a.tol for v in row.values())
        nbad += 1 if bad else 0
        for k2 in worst:
            worst[k2] = max(worst[k2], row.get(k2, 0.0))
        cells = " ".join(f"{row.get(k2, float('nan')):12.3e}" for k2 in ("_PA", "_PA_UU", "_PAVariance"))
        print(f"{j:>3} {d:>16} {cells}{'   <-- MISMATCH' if bad else ''}")

    print()
    for k2 in ("_PA", "_PA_UU", "_PAVariance"):
        print(f"  worst |FO - independent| {k2:<12} = {worst[k2]:.3e}")
    print(f"  (tol {a.tol:g})")
    if nbad == 0 and max(worst.values()) <= a.tol:
        print("PASS: _PA, _PA_UU and _PAVariance all reproduce an independent\n      re-derivation of the same trajectory to round-off.")
        return 0
    print(f"FAIL: {nbad} bin(s) disagree beyond tolerance.")
    return 2


if __name__ == "__main__":
    sys.exit(main())
