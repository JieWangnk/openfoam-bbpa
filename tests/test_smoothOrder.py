#!/usr/bin/env python3
"""
test_smoothOrder.py -- convergence order of the bin integral in the time step.

Smooth prescribed field U(t) = amp (cos, sin, camp cos 2)(2 pi phase/T) + c_n,
whose bin means have a closed form (bbpa_offline.closed_form_smooth_binmean).
Six cycles, dt = binDeltaT / 2^k for k = 0..5, under the two sampling
conventions the harness can prescribe:
  end  : the sample is the field at the END of the step (t), which is what a
         flow solver hands the function object -- the production convention.
         The accumulation is then the right-endpoint rectangle rule, first
         order in dt.
  mid  : the sample is the field at the step MIDPOINT (t - dt/2), the
         harness's default for the bin-constant tests (it names the bin the
         interval lies in). The accumulation is then the midpoint rule,
         second order in dt.
Asserts:
  1. implementation: on-disk _PA == BBPAReplica bit-faithfully at every step
     and under both conventions (the accumulator sees only (t, dt, U));
  2. discretisation: the max relative error of _PA over the finalised bins
     follows a log-log slope within 0.1 of 1 (end) and of 2 (mid).
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

T_END = 6.0
FINAL_CS = 5.0
KS = range(0, 6)
CONV = {"end": 0.0, "mid": 0.5}
EXPECTED = {"end": 1.0, "mid": 2.0}


def _err(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def run_sweep(conv):
    os.makedirs(rv.RUNROOT, exist_ok=True)
    rows = []
    for k in KS:
        dt = rv.BIN_DT / 2 ** k
        c, log = rv.make_case(f"S_smooth_{conv}_k{k}", "smooth", dt, T_END, 3.0,
                              sample_offset=CONV[conv])
        assert rv.block_mesh(c).returncode == 0
        assert rv.run_serial(c).returncode == 0
        disk = rv.read_bbpa_output(c, FINAL_CS)
        rep = bo.BBPAReplica(rv.T, rv.I).replay(rv.read_log(log))
        wrep = rep.write_fields()
        closed = bo.closed_form_smooth_binmean(rv.T, rv.I, rv.AMP, rv.CAMP, rv.c_of_n, int(rep.nCycles))
        impl = max(_err(disk[j]["PA"], wrep[j]["PA"]) for j in range(rv.I))
        scale = max(np.max(np.abs(np.asarray(closed[j], float))) for j in range(rv.I))
        disc = max(_err(disk[j]["PA"], closed[j]) for j in range(rv.I)) / scale
        rows.append((dt / rv.BIN_DT, impl, disc))
    return rows


def slope(rows):
    x = np.log([r[0] for r in rows]); y = np.log([r[2] for r in rows])
    return float(np.polyfit(x, y, 1)[0])


def test_smooth_order():
    for conv in CONV:
        rows = run_sweep(conv)
        assert max(r[1] for r in rows) == 0.0, f"{conv}: implementation differs from the replica"
        s = slope(rows)
        assert abs(s - EXPECTED[conv]) < 0.1, f"{conv}: log-log slope {s:.3f}"


if __name__ == "__main__":
    ok = True
    for conv in CONV:
        rows = run_sweep(conv)
        print(f"smooth field, sample at step {conv}: relative error of _PA vs the closed form, 6 cycles")
        print("  dt/binDeltaT   impl(on-disk vs replica)   rel. disc. error   ratio")
        prev = None
        for r in rows:
            ratio = "   ---" if prev is None else f"{prev / r[2]:6.2f}"
            print(f"  {r[0]:<13.6g} {r[1]:<25.3e} {r[2]:<18.3e} {ratio}")
            prev = r[2]
        s = slope(rows)
        good = max(r[1] for r in rows) == 0.0 and abs(s - EXPECTED[conv]) < 0.1
        ok = ok and good
        print(f"  log-log slope {s:.3f} (expected {EXPECTED[conv]:.0f}): {'PASS' if good else 'FAIL'}")
    print("  RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
