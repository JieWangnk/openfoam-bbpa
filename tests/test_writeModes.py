#!/usr/bin/env python3
"""
test_writeModes.py -- the output-layout invariance matrix (roadmap Tutorial-1).

THE invariant under test: writeMode changes WHERE and WHEN results are written;
it must never change the NUMBERS. Every mode runs the identical bin-constant
prescribe trajectory; each mode's finalised cross-cycle products, read from
wherever that mode puts them, must equal the algorithm-independent closed form
to machine precision -- and modes that finalise the same number of cycles must
agree with each other bit-for-bit.

Mode -> where the finalised numbers land -> expected cycles at the read point:
  phaseAlignedDirs  phase-aligned dirs, every write            6 (final dump)
  perCycle          phase-aligned dirs, once per cycle         6 (final dump)
  bin               phase-aligned dirs, _BinPA* naming         6 (final dump)
                    + per-cycle U_Bin (LAST completed cycle's own bin mean,
                      NO cross-cycle averaging) = A[j] + c(5) exactly
  both              companion each write + final dump          6 (final dump)
  companion         solver write dirs (needs writeInterval=BIN_DT); read at
                    dirs cycleStart+j*BIN_DT DURING the final cycle -> the
                    pure finalised mean over the 5 complete cycles          5

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
DT = rv.BIN_DT / 16               # exact 2**-7: no float drift
END_TIME = 6.0
FINAL_CS = 5.0                    # final dump labels the previous cycle
N_FINAL = 6                       # cycles finalised at the final dump
N_COMPANION = 5                   # complete cycles when reading mid final cycle


def _err(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _run(name, write_interval, bbpa_extra):
    os.makedirs(rv.RUNROOT, exist_ok=True)
    rundir, logpath = rv.make_case(name, "bin_constant", DT, END_TIME,
                                   write_interval, bbpa_extra=bbpa_extra)
    assert rv.block_mesh(rundir).returncode == 0, name + ": blockMesh failed"
    assert rv.run_serial(rundir).returncode == 0, name + ": foamRun failed"
    return rundir


def _check_vs_closed_form(mode, got, closed, fields, errors):
    for j in range(rv.I):
        if j not in got:
            errors.append(f"{mode}: bin {j} missing on disk")
            continue
        for fld in fields:
            ref = closed[j][fld]
            val = got[j].get(fld)
            if val is None:
                errors.append(f"{mode}: bin {j} {fld} missing")
                continue
            e = _err(val, ref)
            if e >= TOL:
                errors.append(f"{mode}: bin {j} {fld} err {e:.3e}")


def run_matrix():
    closed6 = bo.closed_form_bin_constant(rv.A, rv.c_of_n, N_FINAL)
    closed5 = bo.closed_form_bin_constant(rv.A, rv.c_of_n, N_COMPANION)
    errors = []
    finalised = {}                # mode -> {bin: PA vec} for cross-mode bitwise check

    # -- phaseAlignedDirs (explicit default), perCycle: 3-field check at final dump
    for mode in ("phaseAlignedDirs", "perCycle"):
        rd = _run(mode, 3.0, f"writeMode       {mode};")
        cs = rv.latest_cycle_start(rd)
        if cs != FINAL_CS:
            errors.append(f"{mode}: expected final dump at {FINAL_CS}, got {cs}")
            continue
        got = rv.read_bbpa_output(rd, cs)
        _check_vs_closed_form(mode, got, closed6,
                              ("PA", "PA_UU", "PAVariance"), errors)
        finalised[mode] = {j: got[j]["PA"] for j in got}

    # -- bin: _BinPA* naming at the final dump + per-cycle U_Bin
    rd = _run("bin", 3.0, "writeMode       bin;\n        PAwriteInterval 1;")
    cs = rv.latest_cycle_start(rd, fname="U_BinPA")
    if cs != FINAL_CS:
        errors.append(f"bin: expected final dump at {FINAL_CS}, got {cs}")
    else:
        got = rv.read_bbpa_output(rd, cs, suffix="_BinPA")
        _check_vs_closed_form("bin", got, closed6,
                              ("PA", "PA_UU", "PAVariance"), errors)
        finalised["bin"] = {j: got[j]["PA"] for j in got}
        # U_Bin = the LAST completed cycle's own bin mean, exactly A[j]+c(5)
        bins = rv.read_bin_files(rd, cs)
        for j in range(rv.I):
            if j not in bins:
                errors.append(f"bin: U_Bin bin {j} missing")
                continue
            e = _err(bins[j], rv.A[j] + rv.c_of_n(5))
            if e >= TOL:
                errors.append(f"bin: U_Bin bin {j} err {e:.3e}")

    # -- both: companion every write + final phase-aligned dump (3 fields)
    rd = _run("both", rv.BIN_DT, "writeMode       both;")
    cs = rv.latest_cycle_start(rd)
    if cs != FINAL_CS:
        errors.append(f"both: expected final dump at {FINAL_CS}, got {cs}")
    else:
        got = rv.read_bbpa_output(rd, cs)
        _check_vs_closed_form("both", got, closed6,
                              ("PA", "PA_UU", "PAVariance"), errors)
        finalised["both"] = {j: got[j]["PA"] for j in got}

    # -- companion: U_PA beside the solver fields; 5 complete cycles mid-cycle
    rd = _run("companion", rv.BIN_DT, "writeMode       companion;")
    got = rv.read_companion_output(rd, FINAL_CS)
    _check_vs_closed_form("companion", got, closed5, ("PA",), errors)

    # -- cross-mode bitwise identity for the modes that finalise 6 cycles
    ref_mode = "phaseAlignedDirs"
    if ref_mode in finalised:
        for mode, vals in finalised.items():
            if mode == ref_mode:
                continue
            for j, ref in finalised[ref_mode].items():
                if j not in vals:
                    errors.append(f"cross: {mode} bin {j} missing")
                    continue
                e = _err(vals[j], ref)
                if e != 0.0:
                    errors.append(f"cross: {mode} bin {j} differs from "
                                  f"{ref_mode} by {e:.3e}")
    return errors


def test_writeModes():
    errors = run_matrix()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    errors = run_matrix()
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("writeModes matrix: PASS")
    print("  phaseAlignedDirs / perCycle / bin / both: finalised _PA, _PA_UU,")
    print(f"  _PAVariance == closed form (6 cycles) to {TOL:.0e}; bitwise")
    print("  identical across modes. bin: U_Bin == A[j]+c(5) exactly.")
    print("  companion: solver-dir U_PA == closed form (5 cycles).")
    sys.exit(0)
