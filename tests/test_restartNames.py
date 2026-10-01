#!/usr/bin/env python3
"""
test_restartNames.py -- restart across a change of time-name spelling.

The defect this guards against (production run, 2026-09): the restore looked a
bin's checkpoint directory up by the TEXT of its time name. OpenFOAM adapts
its time precision to the directory it restarts from, so the spelling the
restore computed ("134.04127679981") differed in the last digit from the one
the previous leg had written ("134.041276799808"), the guard that refuses an
unmatched checkpoint fired, and 44 of 48 bins were reset to a cold start.

Reproduction: a period with 14 significant digits (bin times need 15 digits to
spell), a non-binary dt (the solver time drifts), leg 1 written at
timePrecision 15, leg 2 restarted from latestTime with timePrecision 12.
Assertions, continuous vs restarted at the last full cycle written:

  1. _PA / _PA_UU / _PAVariance agree on ALL bins to 1e-10;
  2. the restart log carries no "Checkpoint mean field missing" warning and
     no "reset to a clean cold start".

Before the value-matched lookup (timeDirName.H) this test fails on both
counts; the directory names on disk and the names the restore asks for can be
inspected in _runs/RN_restart.
"""
from __future__ import annotations
import os
import re
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bbpa_verify_lib as rv

# A period whose multiples need 15 significant digits: every bin label
# cycleStart + j*binDeltaT then prints differently at 12 and at 15 digits.
rv.T = 1.0000000000001
rv.BIN_DT = rv.T / rv.I
DT = rv.BIN_DT / 16              # T/128: not a binary fraction, so t drifts
T_CKPT = 3 * rv.T
T_END = 6 * rv.T
TOL_RESTART = 1e-10


def _err(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def max_field_err(a, b, bins, fields=("PA", "PA_UU", "PAVariance")):
    out = {}
    for fld in fields:
        e = 0.0
        for j in bins:
            va, vb = a.get(j, {}).get(fld), b.get(j, {}).get(fld)
            if va is None or vb is None:
                continue
            e = max(e, _err(va, vb))
        out[fld] = e
    return out


def run_pair():
    os.makedirs(rv.RUNROOT, exist_ok=True)
    # adjustableRunTime lands the write on 3T exactly, so the restart directory
    # is spelt with 14 digits ("3.0000000000003") while the bin labels inside
    # a cycle need 15: OpenFOAM then re-derives its precision as 14 from that
    # directory, which is the production situation.
    kw = dict(time_precision=15, write_control="adjustableRunTime", adjust_time_step="yes")
    c1, l1 = rv.make_case("RN_cont", "bin_constant", DT, T_END, T_CKPT, **kw)
    assert rv.block_mesh(c1).returncode == 0
    assert rv.run_serial(c1).returncode == 0
    c2, l2 = rv.make_case("RN_restart", "bin_constant", DT, T_CKPT, T_CKPT, **kw)
    assert rv.block_mesh(c2).returncode == 0
    assert rv.run_serial(c2).returncode == 0
    # leg 2: the spelling changes (precision 12, and OpenFOAM re-derives its
    # precision from the directory it resumes from)
    with open(os.path.join(c2, "system", "controlDict"), "w") as f:
        f.write(rv.control_dict("bin_constant", DT, T_END, T_CKPT, "latestTime",
                                l2, time_precision=12,
                                write_control="adjustableRunTime",
                                adjust_time_step="yes"))
    r = rv.sh("foamRun > log.foamRun2 2>&1", c2)
    assert r.returncode == 0, "restart leg failed; see log.foamRun2"
    return c1, c2


def compute(c1, c2):
    cs = rv.latest_cycle_start(c1)
    assert cs is not None, "no full set of bins written by the continuous run"
    b1 = rv.read_bbpa_output(c1, cs)
    b2 = rv.read_bbpa_output(c2, cs)
    allb = [j for j in range(rv.I) if j in b1 and j in b2]
    assert len(allb) == rv.I, f"expected {rv.I} bins at cycle start {cs}, got {len(allb)}"
    with open(os.path.join(c2, "log.foamRun2")) as f:
        log2 = f.read()
    resets = len(re.findall(r"Checkpoint mean field missing|reset to a clean cold start", log2))
    return dict(cycle_start=cs, cont_vs_restart=max_field_err(b1, b2, allb), resets=resets)


def test_restart_names():
    c1, c2 = run_pair()
    r = compute(c1, c2)
    assert r["resets"] == 0, f"{r['resets']} bins reset at the restart"
    assert max(r["cont_vs_restart"].values()) < TOL_RESTART, r["cont_vs_restart"]


if __name__ == "__main__":
    c1, c2 = run_pair()
    r = compute(c1, c2)
    print(f"restart across a time-name spelling change: T={rv.T!r}, dt=T/128, "
          f"leg 1 timePrecision 15, leg 2 timePrecision 12, compared at cycle start {r['cycle_start']!r}")
    print(f"  bins reset at the restart (log warnings): {r['resets']}")
    print("  continuous vs restarted (all bins):")
    for k, v in r["cont_vs_restart"].items():
        print(f"      {k:<11} max|err| = {v:.3e}")
    ok = r["resets"] == 0 and max(r["cont_vs_restart"].values()) < TOL_RESTART
    print("  RESULT:", "PASS" if ok else "FAIL", f"(tol {TOL_RESTART:.0e})")
    sys.exit(0 if ok else 1)
