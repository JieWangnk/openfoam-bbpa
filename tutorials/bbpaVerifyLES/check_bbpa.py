#!/usr/bin/env python3
"""Independent verification of the BBPA function object (LES / turbulent case).

Same strategy as tutorials/bbpaVerify: the case saves the FULL instantaneous
field at EVERY timestep, so we hold the exact trajectory the FO saw.  We
re-derive the bin averages here with logic written independently of the C++ and
compare against the FO's own <field>_PA / _PA_UU / _PAVariance output.  This is
a same-trajectory identity check, so agreement should be to round-off (~1e-10),
not "close enough".  Anything larger is a real defect.

Two things this adds over the bbpaVerify checker:

  1. --field:  the FO tracks more than U.  wallShearStress is checked too.

  2. --patches: BOUNDARY fields are compared, not just the internalField.
     This is essential and not cosmetic: wallShearStress has a UNIFORM (0 0 0)
     internalField -- all of its data lives on the 'walls' patch.  A checker
     that only looks at internalField "passes" the WSS test by comparing zeros
     against zeros, which verifies nothing.  Likewise U's own boundary carries
     the mapped turbulent inlet record.

Why this case exists (the other two cannot do this job):
  - bbpaVerify (clean sine inlet) is deterministic and periodic: once converged
    every cycle is identical, so the cross-cycle variance is ~0 by construction
    and cannot catch a variance bug.
  - planarPoiseuille is a developing flow: it exercises the variance but is not
    periodic-with-turbulence.
  - here the inlet boundaryData is a periodic base profile (period 1.0 exactly)
    plus a ~0.07 rms perturbation that genuinely differs every cycle, so the
    mean, the second moment AND the variance are all non-trivial at once --
    the regime the paper actually runs in.

Usage:
  python3 check_bbpa.py --period 0.9997 --nbins 7 --start-cycle 0 --include-t0
  python3 check_bbpa.py --field wallShearStress --patches walls ...
"""
import argparse
import math
import os
import re
import sys

import numpy as np


# ---------------------------------------------------------------- parsing ----

def _parse_list_body(chunk, kind, n):
    """Parse the '( ... )' payload of an OpenFOAM List<kind> with n entries."""
    if kind == "scalar":
        vals = re.findall(r"[-0-9.eE+]+", chunk)
        arr = np.array([[float(v)] for v in vals[:n]])
    else:
        vals = re.findall(r"\(([^)]*)\)", chunk)
        arr = np.array([[float(x) for x in v.split()] for v in vals[:n]])
    if arr.shape[0] != n:
        raise RuntimeError(f"expected {n} entries, parsed {arr.shape[0]}")
    return arr


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
    return _parse_list_body(chunk, kind, n)


def read_patch(path, patch):
    """Parse boundaryField/<patch>/value -> (n,c) array, or None if absent.

    This is where wallShearStress actually lives (its internalField is all
    zeros), so a WSS check that skips this is vacuous.
    """
    txt = open(path).read()
    bf = txt.find("boundaryField")
    if bf < 0:
        return None
    m = re.search(r"\n\s*" + re.escape(patch) + r"\s*\n\s*\{", txt[bf:])
    if not m:
        return None
    seg = txt[bf + m.end():]
    depth = 1
    for i, ch in enumerate(seg):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                seg = seg[:i]
                break
    mv = re.search(r"value\s+nonuniform\s+List<(\w+)>\s*\n?\s*(\d+)\s*\(", seg)
    if mv:
        kind, n = mv.group(1), int(mv.group(2))
        body = seg[mv.end():]
        end = body.index("\n)") + 2 if "\n)" in body else len(body)
        return _parse_list_body(body[:end], kind, n)
    mv = re.search(r"value\s+uniform\s*\(([^)]*)\)", seg)
    if mv:
        return np.array([[float(x) for x in mv.group(1).split()]])
    mv = re.search(r"value\s+uniform\s+([-0-9.eE+]+)\s*;", seg)
    if mv:
        return np.array([[float(mv.group(1))]])
    if re.search(r"type\s+(empty|nonuniformTransform|symmetry)", seg):
        return None
    # calculated/zeroGradient patches with no explicit value carry no data
    return None


def read_all(path, patches):
    """-> dict {'internal': arr, '<patch>': arr, ...} (missing entries dropped)."""
    out = {"internal": read_internal(path)}
    for p in patches:
        a = read_patch(path, p)
        if a is not None:
            out[p] = a
    return out


def phase_bin(t, start_time, period, dT, nbins):
    """Bin index for time t -- the same snap-then-truncate the FO uses.

    A plain int(fmod(t-t0, period)/dT) is exact only when the result is not
    sitting on an integer boundary, and in an ON-GRID case it always is:
    fmod(2.9, 1.0) returns 0.8999999999999999, so 0.8999.../0.05 =
    17.999999999999996 truncates to 17 when the step belongs in bin 18. The
    error grows with t, which is why the naive version gives 51 steps in the
    first cycle's bin and 49 in the third instead of a constant 50.

    This is NOT a fudge to make the test pass: it is the same arithmetic the FO
    performs, applied independently to a time this script derived on its own.
    The check remains an identity check -- if the FO folded a step into a
    different bin than this does, the bin means still disagree by ~1e-2 and the
    test still fails loudly.
    """
    phase = math.fmod(t - start_time, period)
    if phase > period*(1.0 - 1e-12):
        phase = 0.0
    x = phase/dT
    r = math.floor(x)
    j = int(r + 1.0 if (x - r) > (1.0 - 1e-9) else r)
    return min(j, nbins - 1)


def cycle_number(t, start_time, period):
    """Cycle index for t, snapped the same way as phase_bin."""
    x = (t - start_time)/period
    r = math.floor(x)
    return int(r + 1.0 if (x - r) > (1.0 - 1e-9) else r)


def sqr_sym(U):
    """<U_i U_j> in OpenFOAM's symmTensor component order (xx xy xz yy yz zz)."""
    if U.shape[1] == 1:
        return U * U
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


# ------------------------------------------------------------------- main ----

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=".")
    ap.add_argument("--field", default="U", help="tracked field to verify")
    ap.add_argument("--patches", nargs="*", default=None,
                    help="boundary patches to verify (default: auto by field)")
    ap.add_argument("--period", type=float, default=0.9997)
    ap.add_argument("--nbins", type=int, default=7)
    ap.add_argument("--start-time", default="auto",
                    help="BBPA phase origin; 'auto' reads it from *Properties")
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--start-cycle", type=int, default=0)
    ap.add_argument("--include-t0", action="store_true",
                    help="fold the initial condition with a full dt, as the FO "
                         "does if its execute() is called at the first time")
    a = ap.parse_args()
    case = os.path.abspath(a.case)
    fld = a.field
    if a.patches is None:
        a.patches = ["walls"] if fld == "wallShearStress" else ["inlet", "walls"]

    tds = time_dirs(case)
    if not tds:
        print("[err] no time dirs"); return 1

    # ---- the FO's own bookkeeping (phase origin + how many cycles it folded)
    start_time, ncycles = None, None
    for t, d in reversed(tds):
        pd = os.path.join(case, d, "uniform")
        if not os.path.isdir(pd):
            continue
        for pf in os.listdir(pd):
            if pf.endswith("Properties"):
                txt = open(os.path.join(pd, pf)).read()
                m = re.search(r"startTime\s+([0-9.eE+-]+);", txt)
                if m and start_time is None:
                    start_time = float(m.group(1))
                m = re.search(rf"{fld}\b.*?nCycles\s+(\d+);", txt, re.S)
                if m is None:
                    m = re.search(r"nCycles\s+(\d+);", txt)
                if m and ncycles is None:
                    ncycles = int(m.group(1))
        if start_time is not None:
            break
    if a.start_time != "auto":
        start_time = float(a.start_time)
    if start_time is None:
        print("[err] no *Properties found; cannot get the phase origin"); return 1
    print(f"[fo]  field={fld}  phase origin startTime={start_time}  nCycles={ncycles}")

    dT = a.period / a.nbins

    # ---- ground truth: every saved instantaneous state -----------------------
    traj = [(t, d) for (t, d) in tds
            if os.path.exists(os.path.join(case, d, fld))
            and os.path.exists(os.path.join(case, d, "p"))]
    if len(traj) < 10:
        print(f"[err] only {len(traj)} instantaneous states; need the per-step save")
        return 1
    print(f"[truth] {len(traj)} per-step states from t={traj[0][0]} to {traj[-1][0]}")

    probe = read_all(os.path.join(case, traj[len(traj)//2][1], fld), a.patches)
    print(f"[truth] regions checked: " + ", ".join(
        f"{k}({v.shape[0]}x{v.shape[1]}, max|.|={np.abs(v).max():.3e})"
        for k, v in probe.items()))

    times = np.array([t for t, _ in traj])
    sums, sqsums, wts = {}, {}, {}
    k0 = 0 if a.include_t0 else 1
    for k in range(k0, len(traj)):
        t = times[k]
        dt = (times[k] - times[k-1]) if k > 0 else (times[1] - times[0])
        if t < start_time - 1e-12:
            continue
        j = phase_bin(t, start_time, a.period, dT, a.nbins)
        n = cycle_number(t, start_time, a.period)
        F = read_all(os.path.join(case, traj[k][1], fld), a.patches)
        key = (n, j)
        if key not in wts:
            sums[key] = {r: v * dt for r, v in F.items()}
            sqsums[key] = {r: sqr_sym(v) * dt for r, v in F.items()}
            wts[key] = dt
        else:
            for r, v in F.items():
                sums[key][r] = sums[key][r] + v * dt
                sqsums[key][r] = sqsums[key][r] + sqr_sym(v) * dt
            wts[key] += dt

    present = sorted({n for (n, j) in sums})
    if ncycles is None:
        print("[err] no nCycles in Properties; cannot know which cycles were folded")
        return 1
    present = [n for n in present if n >= a.start_cycle]
    cycles = present[:ncycles]
    print(f"[truth] cycles on disk: {present} ; FO folded {ncycles} -> comparing {cycles}")

    expect, expect_uu, expect_var = {}, {}, {}
    for j in range(a.nbins):
        ok = [n for n in cycles if (n, j) in wts and wts[(n, j)] > 0]
        if not ok:
            continue
        regions = sums[(ok[0], j)].keys()
        cyc_mean = {r: [sums[(n, j)][r] / wts[(n, j)] for n in ok] for r in regions}
        cyc_uu = {r: [sqsums[(n, j)][r] / wts[(n, j)] for n in ok] for r in regions}
        expect[j] = {r: np.mean(v, axis=0) for r, v in cyc_mean.items()}
        expect_uu[j] = {r: np.mean(v, axis=0) for r, v in cyc_uu.items()}
        if len(ok) > 1:
            ev = {}
            for r, v in cyc_mean.items():
                m = np.mean(v, axis=0)
                m2 = np.sum([np.sum((x - m)**2, axis=1) for x in v], axis=0)
                ev[r] = (m2 / (len(ok) - 1.0)).reshape(-1, 1)
            expect_var[j] = ev

    # ---- the FO's answer ----------------------------------------------------
    pa = [(t, d) for (t, d) in tds
          if os.path.exists(os.path.join(case, d, f"{fld}_PA"))]
    if not pa:
        print(f"[err] no {fld}_PA written"); return 1
    latest = pa[-a.nbins:]
    print(f"[fo]  {fld}_PA window: t={latest[0][0]} .. {latest[-1][0]} ({len(latest)} bins)")

    suffixes = ("_PA", "_PA_UU", "_PAVariance")
    print(f"\n{'bin':>3} {'region':>10} {'t_dir':>16} "
          f"{'_PA':>12} {'_PA_UU':>12} {'_PAVariance':>13}")
    worst = {s: 0.0 for s in suffixes}
    nbad = 0
    checked_nontrivial = {s: False for s in suffixes}
    for i, (t, d) in enumerate(latest):
        j = i                      # bin index = POSITION in the window
        if j not in expect:
            continue
        for r in expect[j]:
            row = {}
            p_pa = os.path.join(case, d, f"{fld}_PA")
            fo = read_all(p_pa, a.patches).get(r)
            if fo is None:
                continue
            ex = expect[j][r]
            if ex.shape[0] == 1 and fo.shape[0] > 1:
                ex = np.repeat(ex, fo.shape[0], axis=0)
            row["_PA"] = np.max(np.abs(fo - ex))
            if np.abs(ex).max() > 1e-12:
                checked_nontrivial["_PA"] = True

            p_uu = os.path.join(case, d, f"{fld}_PA_UU")
            if os.path.exists(p_uu) and j in expect_uu:
                fo_uu = read_all(p_uu, a.patches).get(r)
                if fo_uu is not None:
                    ex_uu = expect_uu[j][r]
                    if ex_uu.shape[0] == 1 and fo_uu.shape[0] > 1:
                        ex_uu = np.repeat(ex_uu, fo_uu.shape[0], axis=0)
                    row["_PA_UU"] = np.max(np.abs(fo_uu - ex_uu))
                    if np.abs(ex_uu).max() > 1e-12:
                        checked_nontrivial["_PA_UU"] = True

            p_v = os.path.join(case, d, f"{fld}_PAVariance")
            if os.path.exists(p_v) and j in expect_var:
                fo_v = read_all(p_v, a.patches).get(r)
                if fo_v is not None:
                    ex_v = expect_var[j][r]
                    if ex_v.shape[0] == 1 and fo_v.shape[0] > 1:
                        ex_v = np.repeat(ex_v, fo_v.shape[0], axis=0)
                    row["_PAVariance"] = np.max(np.abs(fo_v - ex_v))
                    if np.abs(ex_v).max() > 1e-12:
                        checked_nontrivial["_PAVariance"] = True

            bad = any(v > a.tol for v in row.values())
            nbad += 1 if bad else 0
            for k2 in worst:
                worst[k2] = max(worst[k2], row.get(k2, 0.0))
            cells = " ".join(f"{row.get(k2, float('nan')):12.3e}" for k2 in suffixes)
            print(f"{j:>3} {r:>10} {d:>16} {cells}"
                  f"{'   <-- MISMATCH' if bad else ''}")

    print()
    for k2 in suffixes:
        tag = "" if checked_nontrivial[k2] else "   [TRIVIAL: expected values are all ~0]"
        print(f"  worst |FO - independent| {k2:<12} = {worst[k2]:.3e}{tag}")
    print(f"  (tol {a.tol:g})")

    # A check that only ever compared zeros is not a passing check.
    vacuous = [k2 for k2 in suffixes if not checked_nontrivial[k2]]
    if nbad == 0 and max(worst.values()) <= a.tol:
        if vacuous:
            print(f"INCONCLUSIVE: {', '.join(vacuous)} compared only ~zero values; "
                  f"that verifies nothing.")
            return 3
        print(f"PASS: {fld}_PA, _PA_UU and _PAVariance all reproduce an independent\n"
              f"      re-derivation of the same trajectory to round-off.")
        return 0
    print(f"FAIL: {nbad} comparison(s) disagree beyond tolerance.")
    return 2


if __name__ == "__main__":
    sys.exit(main())
