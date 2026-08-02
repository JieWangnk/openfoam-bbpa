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
  sqrOnly    U { writePhaseMean off; crossCycleVariance off; }
                -> U_PA_UU only among products; no U_PA anywhere. The mean
                   chain still persists as the U_PAmeanState checkpoint
                   (every finalised accumulated bin is restart-complete
                   regardless of which products are published)
  varOnly    U { writePhaseMean off; secondMoment off; }
                -> U_PAVariance == closed form with NO U_PA file; the mean
                   chain persists as the checkpoint-only U_PAmeanState
                   (Choice A: calculation dependency and file output are
                   separate concerns)
  mixed      fields ( U  p { secondMoment off; crossCycleVariance off; } )
                -> per-field allocation independent: U carries all three
                   products, p only p_PA (p is not prescribed, so only its
                   manifest is asserted)
  perCycOff  writeMode bin with U { writePerCycle off; }
                -> no U_Bin anywhere; the cross-cycle _BinPA* unaffected
  restart    meanOnly, checkpointed at t=3 and restarted to t=6:
                continuous == restarted exactly; ALSO varOnly restarted ->
                variance survives via the meanState checkpoint

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
SQR_ONLY = ("fields          ( U { writePhaseMean off; "
            "crossCycleVariance off; } );")
VAR_ONLY = ("fields          ( U { writePhaseMean off; "
            "secondMoment off; } );")
MIXED = ("fields          ( U  p { secondMoment off; "
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


def _manifest(rundir, names=("U_PA", "U_PA_UU", "U_PAVarianceMag",
                             "U_PAmeanState", "U_Bin", "p_PA", "p_PA_UU")):
    """How often each product file appears across the run's time dirs."""
    found = {f: 0 for f in names}
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
    if not (m["U_PA"] and m["U_PA_UU"] and m["U_PAVarianceMag"]):
        errors.append(f"plain: expected all products, manifest {m}")
    _check_closed("plain", rv.read_bbpa_output(rd, FINAL_CS), closed,
                  ("PA", "PA_UU", "PAVariance"), errors)

    # -- mean only ----------------------------------------------------------
    rd = _run("prod_meanOnly", MEAN_ONLY)
    m = _manifest(rd)
    if not m["U_PA"]:
        errors.append("meanOnly: U_PA missing")
    if m["U_PA_UU"] or m["U_PAVarianceMag"]:
        errors.append(f"meanOnly: off products written, manifest {m}")
    _check_closed("meanOnly", rv.read_bbpa_output(rd, FINAL_CS), closed,
                  ("PA",), errors)

    # -- second moment only --------------------------------------------------
    rd = _run("prod_sqrOnly", SQR_ONLY)
    m = _manifest(rd)
    if not m["U_PA_UU"]:
        errors.append("sqrOnly: U_PA_UU missing")
    if m["U_PA"] or m["U_PAVarianceMag"]:
        errors.append(f"sqrOnly: off products written, manifest {m}")
    if not m["U_PAmeanState"]:
        errors.append("sqrOnly: U_PAmeanState checkpoint missing")
    _check_closed("sqrOnly", rv.read_bbpa_output(rd, FINAL_CS), closed,
                  ("PA_UU",), errors)

    # -- variance only (Choice A): no _PA file, meanState checkpoint ---------
    rd = _run("prod_varOnly", VAR_ONLY)
    m = _manifest(rd)
    if not m["U_PAVarianceMag"]:
        errors.append("varOnly: U_PAVariance missing")
    if not m["U_PAmeanState"]:
        errors.append("varOnly: U_PAmeanState checkpoint missing")
    if m["U_PA"] or m["U_PA_UU"]:
        errors.append(f"varOnly: off products written, manifest {m}")
    _check_closed("varOnly", rv.read_bbpa_output(rd, FINAL_CS), closed,
                  ("PAVariance",), errors)

    # -- mixed fields: per-field allocation independent ----------------------
    rd = _run("prod_mixed", MIXED)
    m = _manifest(rd)
    if not (m["U_PA"] and m["U_PA_UU"] and m["U_PAVarianceMag"] and m["p_PA"]):
        errors.append(f"mixed: expected U all + p_PA, manifest {m}")
    if m["p_PA_UU"]:
        errors.append(f"mixed: p_PA_UU written though off, manifest {m}")
    _check_closed("mixed-U", rv.read_bbpa_output(rd, FINAL_CS), closed,
                  ("PA", "PA_UU", "PAVariance"), errors)

    # -- writePerCycle off in bin mode: no U_Bin, _BinPA unaffected ----------
    rd = _run("prod_perCycOff",
              "fields          ( U { writePerCycle off; } );\n"
              "        writeMode       bin;")
    m = _manifest(rd)
    if m["U_Bin"]:
        errors.append(f"perCycOff: U_Bin written though off ({m['U_Bin']}x)")
    got = rv.read_bbpa_output(rd, FINAL_CS, suffix="_BinPA")
    _check_closed("perCycOff", got, closed, ("PA", "PA_UU", "PAVariance"),
                  errors)

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

    # -- varOnly restart: the variance survives via the meanState checkpoint
    c3, l3 = rv.make_case("prod_restart_var", "bin_constant", DT, 3.0, 3.0,
                          bbpa_body=_body(VAR_ONLY))
    assert rv.block_mesh(c3).returncode == 0
    assert rv.run_serial(c3).returncode == 0
    with open(os.path.join(c3, "system", "controlDict"), "w") as f:
        f.write(rv.control_dict("bin_constant", DT, END_TIME, 3.0,
                                "startTime", l3, start_time=3.0,
                                bbpa_body=_body(VAR_ONLY)))
    assert rv.sh("foamRun > log.foamRun2 2>&1", c3).returncode == 0
    rest_v = rv.read_bbpa_output(c3, FINAL_CS)
    for j in range(rv.I):
        if j not in rest_v:
            errors.append(f"varOnly restart: bin {j} missing")
            continue
        e = _err(rest_v[j]["PAVariance"], closed[j]["PAVariance"])
        if e >= TOL:
            errors.append(f"varOnly restart: bin {j} PAVariance err {e:.3e}")
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
    print("products rows: PASS -- plain/meanOnly/sqrOnly/varOnly/mixed/")
    print(f"  perCycOff manifests exact, surviving products == closed form to")
    print(f"  {TOL:.0e}; meanOnly restart == continuous exactly; varOnly")
    print("  variance survives restart via the meanState checkpoint")
    sys.exit(0)
