#!/usr/bin/env python3
"""
test_dictNames.py -- the v2 dictionary renames and their one-release
compatibility (roadmap step 6).

Two runs of the identical prescribe trajectory, configured once with the NEW
names and once with the DEPRECATED spellings:

  new:  discardCycles 1; maxCycles 4; phaseOrigin 0;
        meanWriteIntervalCycles 1; binAlignment left;
        phaseSelection { mode list; bins (1 4 6); }
  old:  startCycle 1;    cycles 4;    startTime 0;
        PAwriteInterval 1;           binCentred no;
        binsOfInterest (1 4 6);

Asserts:
  1. both runs produce IDENTICAL numbers (bit-for-bit) -- the rename is
     spelling only;
  2. the new-names run == the closed form (ensemble = cycles 1..4 on the
     listed bins);
  3. the old-names log carries a deprecation warning per old entry, and the
     new-names log carries none.

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

TOL = 1e-12
DT = rv.BIN_DT / 16
END_TIME = 6.0
BINS = (1, 4, 6)
N_CYCLES = 4                     # discard 1, max 4 -> ensemble = cycles 1..4

NEW_BODY = """        fields          (U);
        period          1.0;
        nBins           8;
        discardCycles   1;
        maxCycles       4;
        phaseOrigin     0;
        meanWriteIntervalCycles 1;
        binAlignment    left;
        phaseSelection { mode list; bins (1 4 6); }
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
"""

OLD_BODY = """        fields          (U);
        period          1.0;
        nBins           8;
        startCycle      1;
        cycles          4;
        startTime       0;
        PAwriteInterval 1;
        binCentred      no;
        binsOfInterest  (1 4 6);
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
"""

DEPRECATED_WORDS = ["startCycle", "cycles", "startTime",
                    "PAwriteInterval", "binCentred", "binsOfInterest"]


def _err(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _run(name, body):
    os.makedirs(rv.RUNROOT, exist_ok=True)
    rundir, _ = rv.make_case(name, "bin_constant", DT, END_TIME, 1.0,
                             bbpa_body=body)
    assert rv.block_mesh(rundir).returncode == 0, name + ": blockMesh failed"
    assert rv.run_serial(rundir).returncode == 0, name + ": foamRun failed"
    log = open(os.path.join(rundir, "log.foamRun")).read()
    return rundir, log


FINAL_CS = 5.0                   # final dump labels the previous cycle


def _read(rundir):
    # latest_cycle_start() demands a FULL bin set, so it cannot be used on a
    # sparse (phaseSelection list) run; the final dump location is known.
    return rv.read_bbpa_output(rundir, FINAL_CS)


def run_check():
    errors = []
    rd_new, log_new = _run("names_new", NEW_BODY)
    rd_old, log_old = _run("names_old", OLD_BODY)
    got_new, got_old = _read(rd_new), _read(rd_old)

    # 1. bit-identical between spellings
    for j in BINS:
        if j not in got_new or j not in got_old:
            errors.append(f"bin {j} missing (new: {j in got_new}, "
                          f"old: {j in got_old})")
            continue
        for fld in ("PA", "PA_UU", "PAVariance"):
            a, b = got_new[j].get(fld), got_old[j].get(fld)
            if a is None or b is None:
                errors.append(f"bin {j} {fld} missing")
            elif _err(a, b) != 0.0:
                errors.append(f"bin {j} {fld}: new vs old spelling differ "
                              f"by {_err(a, b):.3e}")

    # 2. new-names run == closed form (ensemble = cycles 1..4)
    closed = bo.closed_form_bin_constant(rv.A, lambda n: rv.c_of_n(n + 1),
                                         N_CYCLES)
    for j in BINS:
        if j not in got_new:
            continue
        for fld in ("PA", "PA_UU", "PAVariance"):
            e = _err(got_new[j][fld], closed[j][fld])
            if e >= TOL:
                errors.append(f"new-names bin {j} {fld} vs closed form: {e:.3e}")

    # 3. deprecation warnings present for old, absent for new
    for w in DEPRECATED_WORDS:
        if f"'{w}' is deprecated" not in log_old:
            errors.append(f"old-names log lacks the deprecation warning "
                          f"for '{w}'")
    if "is deprecated" in log_new:
        errors.append("new-names log unexpectedly carries a deprecation "
                      "warning")
    return errors


def test_dictNames():
    errors = run_check()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    errors = run_check()
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("dict names: PASS -- new spellings == old spellings bit-for-bit,")
    print("  new == closed form, deprecation warnings exactly where expected")
    sys.exit(0)
