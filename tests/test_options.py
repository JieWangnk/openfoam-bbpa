#!/usr/bin/env python3
"""
test_options.py -- the remaining Tutorial-1 dictionary-option rows.

Same bin-constant prescribe trajectory as the other tests (dt = binDeltaT/16,
exact); each option row asserts its documented contract to machine precision:

  startCycle 2       discarded cycles contribute nothing: the finalised
                     ensemble is cycles 2..5, i.e. closed form over
                     c(n+2), n=0..3.
  cycles 3           accumulation stops after exactly 3 completed cycles
                     (0,1,2) and the values FREEZE: the final dump equals the
                     3-cycle closed form.  (Locks in the off-by-one fix: the
                     limit cycle used to stay stranded unfinalised, silently
                     delivering an (N-1)-cycle ensemble.)
  PAwriteInterval 2  (writeMode bin) U_Bin appears every cycle; U_BinPA is
                     absent at an odd-cycle dump and present at the final one.
  PAwriteInterval 0  U_BinPA appears ONLY at the final write.
  binsOfInterest     (2 5): only the listed bins are accumulated and written;
                     their values equal the full-run closed form (sparse
                     tracking must not change the numbers); unlisted bins
                     produce no output.

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
FINAL_CS = 5.0


def _err(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _run(name, bbpa_extra, write_interval=1.0):
    os.makedirs(rv.RUNROOT, exist_ok=True)
    rundir, logpath = rv.make_case(name, "bin_constant", DT, END_TIME,
                                   write_interval, bbpa_extra=bbpa_extra)
    assert rv.block_mesh(rundir).returncode == 0, name + ": blockMesh failed"
    assert rv.run_serial(rundir).returncode == 0, name + ": foamRun failed"
    return rundir, logpath


def _bin_dir(rundir, cs, j):
    for tv, full in rv._time_dirs(rundir).items():
        if abs(tv - (cs + j * rv.BIN_DT)) < 1e-7:
            return full
    return None


def _check_ensemble(mode, got, closed, errors,
                    fields=("PA", "PA_UU", "PAVariance"), bins=None):
    for j in (bins if bins is not None else range(rv.I)):
        if j not in got:
            errors.append(f"{mode}: bin {j} missing")
            continue
        for fld in fields:
            ref, val = closed[j][fld], got[j].get(fld)
            if val is None:
                errors.append(f"{mode}: bin {j} {fld} missing")
            elif _err(val, ref) >= TOL:
                errors.append(f"{mode}: bin {j} {fld} err {_err(val, ref):.3e}")


def run_rows():
    errors = []

    # -- startCycle 2: ensemble = cycles 2..5 --------------------------------
    rd, _ = _run("opt_startCycle", "discardCycles   2;")
    got = rv.read_bbpa_output(rd, rv.latest_cycle_start(rd))
    closed = bo.closed_form_bin_constant(rv.A, lambda n: rv.c_of_n(n + 2), 4)
    _check_ensemble("startCycle", got, closed, errors)

    # -- cycles 3: exactly cycles 0,1,2; frozen thereafter -------------------
    rd, _ = _run("opt_cycles", "maxCycles       3;")
    cs = rv.latest_cycle_start(rd)
    got = rv.read_bbpa_output(rd, cs)
    closed = bo.closed_form_bin_constant(rv.A, rv.c_of_n, 3)
    _check_ensemble("cycles", got, closed, errors)
    if cs != FINAL_CS:
        errors.append(f"cycles: frozen values must still be dumped at the "
                      f"final write (expected cs {FINAL_CS}, got {cs})")

    # -- PAwriteInterval 2 (bin mode): cadence of _BinPA vs _Bin -------------
    rd, _ = _run("opt_pawi2",
                 "writeMode       bin;\n        meanWriteIntervalCycles 2;")
    # final dump (currentCycle 6, even; also last write): full _BinPA set
    got = rv.read_bbpa_output(rd, FINAL_CS, suffix="_BinPA")
    closed = bo.closed_form_bin_constant(rv.A, rv.c_of_n, 6)
    _check_ensemble("PAwriteInterval2-final", got, closed, errors)
    # odd-cycle dump at cs=4.0 (written when currentCycle 5): U_Bin yes, U_BinPA no
    d = _bin_dir(rd, 4.0, 0)
    if d is None:
        errors.append("PAwriteInterval2: cs=4.0 bin-0 dir missing")
    else:
        if not os.path.isfile(os.path.join(d, "U_Bin")):
            errors.append("PAwriteInterval2: U_Bin missing at odd-cycle dump")
        if os.path.isfile(os.path.join(d, "U_BinPA")):
            errors.append("PAwriteInterval2: U_BinPA written at odd-cycle dump")

    # -- PAwriteInterval 0: _BinPA only at the final write -------------------
    rd, _ = _run("opt_pawi0",
                 "writeMode       bin;\n        meanWriteIntervalCycles 0;")
    got = rv.read_bbpa_output(rd, FINAL_CS, suffix="_BinPA")
    _check_ensemble("PAwriteInterval0-final", got, closed, errors)
    for cs_early in (3.0, 4.0):
        d = _bin_dir(rd, cs_early, 0)
        if d and os.path.isfile(os.path.join(d, "U_BinPA")):
            errors.append(f"PAwriteInterval0: U_BinPA written early (cs {cs_early})")

    # -- sparse binsOfInterest (2 5): only those bins, same numbers ----------
    rd, _ = _run("opt_sparse",
                 "phaseSelection { mode list; bins (2 5); }")
    got = rv.read_bbpa_output(rd, FINAL_CS)
    _check_ensemble("sparse", got, closed, errors, bins=(2, 5))
    for j in range(rv.I):
        if j in (2, 5):
            continue
        d = _bin_dir(rd, FINAL_CS, j)
        if d and os.path.isfile(os.path.join(d, "U_PA")):
            errors.append(f"sparse: unlisted bin {j} was written")

    return errors


def test_options():
    errors = run_rows()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    errors = run_rows()
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("options rows: PASS")
    print("  startCycle 2 (ensemble = cycles 2..5), cycles 3 (exactly 3, frozen),")
    print("  PAwriteInterval 2/0 cadence, sparse binsOfInterest (2 5) --")
    print(f"  all values == closed form to {TOL:.0e}")
    sys.exit(0)
