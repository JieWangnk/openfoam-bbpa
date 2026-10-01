#!/usr/bin/env python3
"""
test_centred.py -- the centred-window ownership contract (review item 5).

Contract: a centred bin belongs to the cycle containing its LABEL; its
integration window [t_j - D/2, t_j + D/2) may overlap the neighbouring
nominal cycles; discardCycles applies to the labelled cycle; accumulation
begins half a bin before the first retained cycle, so the first retained
centred window is complete.

Rows (bin-constant trajectory, binAlignment centred, discardCycles 1):

  weights   every retained bin -- including bin 0, whose window straddles
            the nominal cycle origin -- carries the full weight D per
            cycle: the checkpoint's per-cycle counting is equal-weight
            (asserted via the signature's binWeights/binCounts).
  restart   checkpoint at t = 3.0, which under centred alignment lies
            INSIDE bin 0's window: bins 1..N-1 reproduce the continuous
            run exactly; the straddled bin's partially-covered cycle is
            DROPPED (the documented unbiased choice), so it re-converges
            with one fewer cycle rather than a biased one.

Runnable standalone (exit 1 on failure, for Alltest) and via pytest.
"""
from __future__ import annotations
import os
import re
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bbpa_verify_lib as rv

TOL = 1e-12
DT = rv.BIN_DT / 16
T_CKPT = 3.0
T_END = 6.0
FINAL_CS = 5.0

BODY = """        fields          (U);
        phaseGrid
        {
            period          1.0;
            nBins           8;
            phaseOrigin     0;
            binAlignment    centred;
        }
        discardCycles   1;
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
"""


def _err(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _run(name, end_time):
    rundir, logpath = rv.make_case(name, "bin_constant", DT, end_time, 3.0,
                                   bbpa_body=BODY)
    assert rv.block_mesh(rundir).returncode == 0, name + ": blockMesh failed"
    assert rv.run_serial(rundir).returncode == 0, name + ": foamRun failed"
    return rundir, logpath


def _bin_counts(rundir, t):
    """Parse U's binCounts from the checkpoint Properties at time t."""
    p = os.path.join(rundir, t, "uniform")
    f = [x for x in os.listdir(p) if x.endswith("Properties")][0]
    txt = open(os.path.join(p, f)).read()
    m = re.search(r"binCounts\s+(\d+)\s*\{\s*([0-9eE.+-]+)\s*\}", txt)
    if m:
        return [float(m.group(2))] * int(m.group(1))
    m = re.search(r"binCounts[^(]*\(([^)]*)\)", txt, re.S)
    return [float(x) for x in m.group(1).split()]


def run_rows():
    errors = []
    os.makedirs(rv.RUNROOT, exist_ok=True)

    # -- continuous reference -------------------------------------------------
    c1, _ = _run("cen_cont", T_END)
    counts = _bin_counts(c1, "6")
    # equal-weight counting: every retained bin, including the
    # boundary-straddling bin 0, folded the same number of cycles
    if len(set(counts)) != 1:
        errors.append(f"weights: unequal cycle counts across bins: {counts}")
    b1 = rv.read_bbpa_output(c1, FINAL_CS)
    if len(b1) != rv.I:
        errors.append(f"continuous: expected {rv.I} bins, got {len(b1)}")

    # -- restart INSIDE the centred boundary bin ------------------------------
    c2, l2 = _run("cen_restart", T_CKPT)
    with open(os.path.join(c2, "system", "controlDict"), "w") as f:
        f.write(rv.control_dict("bin_constant", DT, T_END, 3.0, "startTime",
                                l2, start_time=T_CKPT, bbpa_body=BODY))
    assert rv.sh("foamRun > log.foamRun2 2>&1", c2).returncode == 0
    log2 = open(os.path.join(c2, "log.foamRun2")).read()
    if "straddled" not in log2:
        errors.append("restart: the boundary-bin straddle was not detected")
    b2 = rv.read_bbpa_output(c2, FINAL_CS)
    for j in range(1, rv.I):
        if j not in b1 or j not in b2:
            errors.append(f"restart: bin {j} missing")
            continue
        for fld in ("PA", "PA_UU", "PAVariance"):
            e = _err(b1[j][fld], b2[j][fld])
            if e != 0.0:
                errors.append(f"restart: bin {j} {fld} continuous vs "
                              f"restarted differ by {e:.3e}")
    # bin 0: one cycle dropped (partial window), never biased
    cc = _bin_counts(c2, "6")
    if not (cc[0] == counts[0] - 1):
        errors.append(f"restart: straddled bin 0 should have exactly one "
                      f"fewer cycle ({cc[0]} vs {counts[0]})")
    if not all(cc[j] == counts[j] for j in range(1, rv.I)):
        errors.append(f"restart: non-straddled counts changed: {cc}")
    return errors


def test_centred():
    errors = run_rows()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    errors = run_rows()
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("centred rows: PASS -- equal-weight counting across the boundary")
    print("  bin; a restart inside the centred boundary window reproduces")
    print("  bins 1..N-1 exactly and drops exactly one (partial) cycle from")
    print("  the straddled bin, unbiased")
    sys.exit(0)
