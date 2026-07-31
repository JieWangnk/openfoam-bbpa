#!/usr/bin/env python3
"""
test_restart.py -- restart-vs-uninterrupted correctness (roadmap step 4).

Continuous run to t=6 vs the same trajectory killed at the t=3 checkpoint and
restarted to t=6 (serial, default phase-aligned writeMode, cross-cycle VARYING
field so the variance restore path is exercised).  Asserts:

  1. continuous vs restarted on-disk _PA/_PA_UU/_PAVariance at the final dump:
     ALL bins, <= 1e-10 (roadmap tolerance; with the no-blend checkpoint the
     restored state is the pure finalised ensemble, so this lands at ~1e-16).
  2. restarted BBPA vs the offline restore-replica (BBPAReplica.restore fed
     segment 1's public outputs, then replaying segment 2): <= 1e-12 -- proves
     the restore CODE matches its spec, independently of (1).

Known gap (documented, reported, NOT asserted): writeMode companion writes no
phase-aligned checkpoint fields, so a companion-only run's restart restores
counters without state; strict restart metadata is roadmap step 8.

Runnable standalone (exit 1 on failure, for Alltest) and via pytest.
"""
from __future__ import annotations
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bbpa_verify_lib as rv
import bbpa_offline as bo

TOL_RESTART = 1e-10               # roadmap acceptance
TOL_IMPL = 1e-12
DT = rv.BIN_DT / 16               # exact 2**-7
T_CKPT = 3.0
T_END = 6.0
FINAL_CS = 5.0


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
    # continuous reference
    c1, l1 = rv.make_case("R_cont", "bin_constant", DT, T_END, 3.0)
    assert rv.block_mesh(c1).returncode == 0
    assert rv.run_serial(c1).returncode == 0
    # checkpointed + restarted (same trajectory)
    c2, l2 = rv.make_case("R_restart", "bin_constant", DT, T_CKPT, 3.0)
    assert rv.block_mesh(c2).returncode == 0
    assert rv.run_serial(c2).returncode == 0
    with open(os.path.join(c2, "system", "controlDict"), "w") as f:
        f.write(rv.control_dict("bin_constant", DT, T_END, 3.0, "startTime",
                                l2, start_time=T_CKPT))
    assert rv.sh("foamRun > log.foamRun2 2>&1", c2).returncode == 0
    return c1, c2, l2


def compute(c1, c2, l2):
    b1 = rv.read_bbpa_output(c1, FINAL_CS)
    b2 = rv.read_bbpa_output(c2, FINAL_CS)
    allb = [j for j in range(rv.I) if j in b1 and j in b2]
    assert len(allb) == rv.I, f"expected {rv.I} bins, got {len(allb)}"

    # offline restore-replica: seg1 -> restore from its public outputs -> seg2
    samples = rv.read_log(l2)
    seg1, seg2 = rv.split_restart_log(samples, T_CKPT)
    r1 = bo.BBPAReplica(rv.T, rv.I).replay(seg1)
    rr = bo.BBPAReplica.restore(r1.write_fields(), r1.binCounts, r1.nCycles,
                                r1.currentCycle, r1.previousCycle,
                                rv.T, rv.I).replay(seg2)
    wrest = rr.write_fields()

    return dict(
        cont_vs_restart=max_field_err(b1, b2, allb),
        restart_vs_replica=max_field_err(b2, wrest, allb),
    )


def test_restart():
    c1, c2, l2 = run_pair()
    r = compute(c1, c2, l2)
    assert max(r["cont_vs_restart"].values()) < TOL_RESTART, r["cont_vs_restart"]
    assert max(r["restart_vs_replica"].values()) < TOL_IMPL, r["restart_vs_replica"]


if __name__ == "__main__":
    c1, c2, l2 = run_pair()
    r = compute(c1, c2, l2)
    print(f"restart: ckpt at t={T_CKPT}, end t={T_END}, varying c_n (variance restored)")
    print("  continuous vs restarted (all bins):")
    for k, v in r["cont_vs_restart"].items():
        print(f"      {k:<11} max|err| = {v:.3e}")
    print("  restarted vs offline restore-replica (all bins):")
    for k, v in r["restart_vs_replica"].items():
        print(f"      {k:<11} max|err| = {v:.3e}")
    ok = (max(r["cont_vs_restart"].values()) < TOL_RESTART
          and max(r["restart_vs_replica"].values()) < TOL_IMPL)
    print("  RESULT:", "PASS" if ok else "FAIL",
          f"(restart tol {TOL_RESTART:.0e}, impl tol {TOL_IMPL:.0e})")
    sys.exit(0 if ok else 1)
