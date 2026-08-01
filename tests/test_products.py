#!/usr/bin/env python3
"""
test_products.py -- per-field product controls (roadmap step 7 / v2 fields()).

Same bin-constant prescribe trajectory; each row runs one per-field
configuration and asserts BOTH the file manifest (a product that is off must
produce no file anywhere) AND that the products still present equal the
closed form to machine precision -- turning a product off must not perturb
the ones that remain.

  plain      fields (U);                                back-compat: all on
  meanOnly   U { secondMoment off; crossCycleVariance off; }
                -> U_PA only; no U_PA_UU, no U_PAVariance anywhere
  sqrOnly    U { phaseMean off; crossCycleVariance off; }
                -> U_PA_UU only; no U_PA anywhere (the mean chain still
                   accumulates internally; only its FILE is off)
  restart    meanOnly, checkpointed at t=3 and restarted to t=6:
                continuous == restarted exactly (reduced products must not
                break the checkpoint round-trip)

The invalid combination (crossCycleVariance without phaseMean) aborts --
asserted in test_invalidInput.

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

MEAN_ONLY = ("fields          ( U { secondMoment off; "
             "crossCycleVariance off; } );")
SQR_ONLY = ("fields          ( U { phaseMean off; "
            "crossCycleVariance off; } );")


def _err(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _body(fields_line):
    return f"""        {fields_line}
        period          1.0;
        nBins           8;
        discardCycles   0;
        phaseOrigin     0;
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
"""


def _run(name, fields_line, end_time=END_TIME):
    os.makedirs(rv.RUNROOT, exist_ok=True)
    rundir, _ = rv.make_case(name, "bin_constant", DT, end_time, 1.0,
                             bbpa_body=_body(fields_line))
    assert rv.block_mesh(rundir).returncode == 0, name + ": blockMesh failed"
    assert rv.run_serial(rundir).returncode == 0, name + ": foamRun failed"
    return rundir


def _manifest(rundir):
    """Which product files exist anywhere in the run's time dirs."""
    found = {"U_PA": 0, "U_PA_UU": 0, "U_PAVariance": 0}
    for tv, full in rv._time_dirs(rundir).items():
        for f in found:
            if os.path.isfile(os.path.join(full, f)):
                found[f] += 1
    return found


def _check_closed(name, got, closed, fields, errors):
    for j in range(rv.I):
        if j not in got:
            errors.append(f"{name}: bin {j} missing")
            continue
        for fld in fields:
            val = got[j].get(fld)
            if val is None:
                errors.append(f"{name}: bin {j} {fld} missing")
            elif _err(val, closed[j][fld]) >= TOL:
                errors.append(f"{name}: bin {j} {fld} err "
                              f"{_err(val, closed[j][fld]):.3e}")


def run_rows():
    errors = []
    closed = bo.closed_form_bin_constant(rv.A, rv.c_of_n, 6)

    # -- plain wordList: back-compat, all products ---------------------------
    rd = _run("prod_plain", "fields          (U);")
    m = _manifest(rd)
    if not (m["U_PA"] and m["U_PA_UU"] and m["U_PAVariance"]):
        errors.append(f"plain: expected all products, manifest {m}")
    _check_closed("plain", rv.read_bbpa_output(rd, FINAL_CS), closed,
                  ("PA", "PA_UU", "PAVariance"), errors)

    # -- mean only ----------------------------------------------------------
    rd = _run("prod_meanOnly", MEAN_ONLY)
    m = _manifest(rd)
    if not m["U_PA"]:
        errors.append("meanOnly: U_PA missing")
    if m["U_PA_UU"] or m["U_PAVariance"]:
        errors.append(f"meanOnly: off products written, manifest {m}")
    _check_closed("meanOnly", rv.read_bbpa_output(rd, FINAL_CS), closed,
                  ("PA",), errors)

    # -- second moment only --------------------------------------------------
    rd = _run("prod_sqrOnly", SQR_ONLY)
    m = _manifest(rd)
    if not m["U_PA_UU"]:
        errors.append("sqrOnly: U_PA_UU missing")
    if m["U_PA"] or m["U_PAVariance"]:
        errors.append(f"sqrOnly: off products written, manifest {m}")
    _check_closed("sqrOnly", rv.read_bbpa_output(rd, FINAL_CS), closed,
                  ("PA_UU",), errors)

    # -- restart with reduced products ---------------------------------------
    c2, l2 = rv.make_case("prod_restart", "bin_constant", DT, 3.0, 3.0,
                          bbpa_body=_body(MEAN_ONLY))
    assert rv.block_mesh(c2).returncode == 0
    assert rv.run_serial(c2).returncode == 0
    with open(os.path.join(c2, "system", "controlDict"), "w") as f:
        f.write(rv.control_dict("bin_constant", DT, END_TIME, 3.0,
                                "startTime", l2, start_time=3.0,
                                bbpa_body=_body(MEAN_ONLY)))
    assert rv.sh("foamRun > log.foamRun2 2>&1", c2).returncode == 0
    cont = rv.read_bbpa_output(_run("prod_restart_cont", MEAN_ONLY), FINAL_CS)
    rest = rv.read_bbpa_output(c2, FINAL_CS)
    for j in range(rv.I):
        if j not in cont or j not in rest:
            errors.append(f"restart: bin {j} missing")
            continue
        e = _err(cont[j]["PA"], rest[j]["PA"])
        if e != 0.0:
            errors.append(f"restart: bin {j} PA continuous vs restarted "
                          f"differ by {e:.3e}")
    return errors


def test_products():
    errors = run_rows()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    errors = run_rows()
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("products rows: PASS -- plain/meanOnly/sqrOnly manifests exact,")
    print(f"  surviving products == closed form to {TOL:.0e}, and a meanOnly")
    print("  restart reproduces the continuous run exactly")
    sys.exit(0)
