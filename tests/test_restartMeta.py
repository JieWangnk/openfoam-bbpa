#!/usr/bin/env python3
"""
test_restartMeta.py -- strict restart metadata (roadmap step 8).

Every checkpoint now carries a configuration signature (formatVersion,
period, phase origin, nBins, binAlignment, selected bins, per-field type +
product settings, counters, bin weights, partial-cycle presence). On
restart, any change that alters accumulator SHAPE or MEANING must abort
with a message naming the entry; changes that only discard state warn.

Each row: run segment 1 (t 0..3) with config A, rewrite controlDict with
config B (+ startTime 3), rerun -- assert the expected abort/warn.

  ABORT rows: nBins changed, period changed, binAlignment changed,
              phase selection changed, explicit phaseOrigin changed,
              field removed, product turned ON (secondMoment off -> on),
              tampered/newer formatVersion.
  WARN rows:  product turned OFF (runs to completion, warns, surviving
              _PA still == closed form).

Runnable standalone (exit 1 on failure, for Alltest) and via pytest.
"""
from __future__ import annotations
import glob
import os
import re
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bbpa_verify_lib as rv
import bbpa_offline as bo

TOL = 1e-12
DT = rv.BIN_DT / 16
T_CKPT = 3.0
T_END = 6.0
FINAL_CS = 5.0


def _body(extra=""):
    return f"""        fields          (U);
        period          1.0;
        nBins           8;
        discardCycles   0;
        phaseOrigin     0;
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
{extra}"""


ABORT_ROWS = [
    ("nBins",       _body(),                     _body("        nBins           10;\n").replace("nBins           8;\n", ""),  "'nBins' changed"),
    ("period",      _body(),                     _body().replace("period          1.0;", "period          0.5;"),             "'period' changed"),
    ("alignment",   _body(),                     _body("        binAlignment    centred;\n"),                                 "'binAlignment' changed"),
    ("selection",   _body(),                     _body("        accumulationSelection { mode list; bins (2 5); }\n"),                "phase selection changed"),
    ("phaseOrigin", _body(),                     _body().replace("phaseOrigin     0;", "phaseOrigin     0.25;"),              "phase origin changed"),
    ("fieldGone",   _body().replace("fields          (U);", "fields          (U p);"),
                    _body(),                                                                                                   "no longer in the 'fields' list"),
    ("prodOn",      _body().replace("fields          (U);",
                                    "fields          ( U { secondMoment off; crossCycleVariance off; } );"),
                    _body(),                                                                                                   "turned ON across the restart"),
]


def _err(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _seg1(name, bodyA):
    rundir, logpath = rv.make_case(name, "bin_constant", DT, T_CKPT, 3.0,
                                   bbpa_body=bodyA)
    assert rv.block_mesh(rundir).returncode == 0, name + ": blockMesh failed"
    assert rv.run_serial(rundir).returncode == 0, name + ": segment 1 failed"
    return rundir, logpath


def _restart(rundir, logpath, bodyB):
    with open(os.path.join(rundir, "system", "controlDict"), "w") as f:
        f.write(rv.control_dict("bin_constant", DT, T_END, 3.0, "startTime",
                                logpath, start_time=T_CKPT, bbpa_body=bodyB))
    r = rv.sh("foamRun > log.foamRun2 2>&1", rundir)
    log = open(os.path.join(rundir, "log.foamRun2")).read()
    return r.returncode, log


def run_rows():
    errors = []
    os.makedirs(rv.RUNROOT, exist_ok=True)

    for name, bodyA, bodyB, fragment in ABORT_ROWS:
        rundir, logpath = _seg1("meta_" + name, bodyA)
        rc, log = _restart(rundir, logpath, bodyB)
        if rc == 0:
            errors.append(f"{name}: restart ran; expected an abort")
        elif fragment not in log:
            tail = "\n".join(log.splitlines()[-12:])
            errors.append(f"{name}: aborted without '{fragment}'. "
                          f"Tail:\n{tail}")

    # -- tampered / newer formatVersion --------------------------------------
    rundir, logpath = _seg1("meta_format", _body())
    props = glob.glob(os.path.join(rundir, "3", "uniform", "*Properties"))
    if not props:
        errors.append("format: checkpoint Properties file not found")
    else:
        s = open(props[0]).read()
        s2, n = re.subn(r"formatVersion\s+2;", "formatVersion   99;", s)
        if n != 1:
            errors.append("format: could not tamper formatVersion")
        else:
            open(props[0], "w").write(s2)
            rc, log = _restart(rundir, logpath, _body())
            if rc == 0:
                errors.append("format: restart ran; expected an abort")
            elif "formatVersion" not in log:
                errors.append("format: aborted without naming formatVersion")

    # -- tampered timeQuadrature: part of the numerical meaning ---------------
    rundir, logpath = _seg1("meta_quad", _body())
    props = glob.glob(os.path.join(rundir, "3", "uniform", "*Properties"))
    s0 = open(props[0]).read()
    open(props[0], "w").write(
        s0.replace("timeQuadrature  endPoint;", "timeQuadrature  trapezoidal;"))
    rc, log = _restart(rundir, logpath, _body())
    if rc == 0:
        errors.append("quad: restart ran; expected an abort")
    elif "timeQuadrature" not in log:
        errors.append("quad: aborted without naming timeQuadrature")

    # -- WARN row: product turned OFF runs, warns, survivors exact -----------
    bodyA = _body()
    bodyB = _body().replace(
        "fields          (U);",
        "fields          ( U { secondMoment off; crossCycleVariance off; } );")
    rundir, logpath = _seg1("meta_prodOff", bodyA)
    rc, log = _restart(rundir, logpath, bodyB)
    if rc != 0:
        errors.append("prodOff: expected the restart to run (warn only)")
    else:
        if "turned off across the restart" not in log:
            errors.append("prodOff: no product-discard warning in the log")
        closed = bo.closed_form_bin_constant(rv.A, rv.c_of_n, 6)
        got = rv.read_bbpa_output(rundir, FINAL_CS)
        for j in range(rv.I):
            if j not in got:
                errors.append(f"prodOff: bin {j} missing")
                continue
            e = _err(got[j]["PA"], closed[j]["PA"])
            if e >= TOL:
                errors.append(f"prodOff: bin {j} PA err {e:.3e}")

    return errors


def test_restartMeta():
    errors = run_rows()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    errors = run_rows()
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("restart metadata: PASS -- 8 shape/meaning changes abort naming the")
    print("  entry; a product turned off warns, runs, and the surviving _PA")
    print(f"  still equals the closed form to {TOL:.0e}")
    sys.exit(0)
