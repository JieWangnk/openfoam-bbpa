#!/usr/bin/env python3
"""
test_parallel.py -- parallel operation with mixed per-field products
(review-matrix row: MPI run, exact manifest and numerical equivalence).

The identical bin-constant prescribe trajectory run serial and decomposed
(np = 2), with MIXED products -- U carries everything, p only its mean.
The prescribed field is spatially uniform, so every cell (and hence every
rank) must hold the same value, making the serial and parallel outputs
directly comparable. Asserts:

  1. serial == parallel bit-for-bit for U's _PA/_PA_UU/_PAVariance at the
     final dump (the per-bin patch/processor reductions and the
     processor-local write paths agree);
  2. the parallel manifest matches the product configuration exactly
     (p_PA present on the ranks, p_PA_UU / p_PAVariance nowhere);
  3. serial == closed form (the numbers are RIGHT, both agreeing on a
     wrong answer would still fail).

Runnable standalone (exit 1 on failure, for Alltest) and via pytest.
"""
from __future__ import annotations
import glob
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bbpa_verify_lib as rv
import bbpa_offline as bo

TOL = 1e-12
DT = rv.BIN_DT / 16
END_TIME = 6.0
FINAL_CS = 5.0

BODY = """        fields          ( U  p { secondMoment off; crossCycleVariance off; } );
        period          1.0;
        nBins           8;
        discardCycles   0;
        phaseOrigin     0;
        output { layout phaseAligned; cycleWriteInterval 1; }
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
"""


def _err(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _make(name):
    rundir, _ = rv.make_case(name, "bin_constant", DT, END_TIME, 1.0,
                             bbpa_body=BODY)
    assert rv.block_mesh(rundir).returncode == 0, name + ": blockMesh failed"
    return rundir


def run_check():
    errors = []
    os.makedirs(rv.RUNROOT, exist_ok=True)
    closed = bo.closed_form_bin_constant(rv.A, rv.c_of_n, 6)

    rd_s = _make("par_serial")
    assert rv.run_serial(rd_s).returncode == 0, "serial run failed"
    got_s = rv.read_bbpa_output(rd_s, FINAL_CS)

    rd_p = _make("par_parallel")
    r = rv.run_parallel(rd_p, 2)
    if r.returncode != 0:
        return ["parallel run failed -- see " + rd_p + "/log.foamRun"]
    got_p = rv.read_bbpa_output(rd_p, FINAL_CS, parallel=True)

    for j in range(rv.I):
        if j not in got_s or j not in got_p:
            errors.append(f"bin {j} missing (serial {j in got_s}, "
                          f"parallel {j in got_p})")
            continue
        for fld in ("PA", "PA_UU", "PAVariance"):
            a, b = got_s[j].get(fld), got_p[j].get(fld)
            if a is None or b is None:
                errors.append(f"bin {j} {fld} missing")
                continue
            if _err(a, b) != 0.0:
                errors.append(f"bin {j} {fld}: serial vs parallel differ "
                              f"by {_err(a, b):.3e}")
            e = _err(a, closed[j][fld])
            if e >= TOL:
                errors.append(f"bin {j} {fld}: serial vs closed form "
                              f"{e:.3e}")

    # manifest on the ranks: p mean-only
    if not glob.glob(os.path.join(rd_p, "processor0", "[0-9]*", "p_PA")):
        errors.append("parallel: p_PA missing on processor0")
    for bad in ("p_PA_UU", "p_PAVariance"):
        hits = glob.glob(os.path.join(rd_p, "processor*", "[0-9]*", bad))
        if hits:
            errors.append(f"parallel: {bad} written though off: {hits[:2]}")
    return errors


def test_parallel():
    errors = run_check()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    errors = run_check()
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("parallel: PASS -- serial == parallel bit-for-bit, both == closed")
    print(f"  form to {TOL:.0e}, mixed-product manifest exact on the ranks")
    sys.exit(0)
