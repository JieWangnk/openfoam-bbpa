#!/usr/bin/env python3
"""
test_output.py -- the v2 output{} block (roadmap step 7b / review item 3).

output{} decides WHERE and WHEN; per-field options decide WHAT; the
accumulator has one numerical meaning throughout. No legacy writeMode word
survives inside it, and the per-cycle _Bin capture is NOT a layout (it is
the per-field writePerCycle product).

  defaults    output { }                -> layout phaseAligned, dump once
              per completed cycle (cs grids at 4.0 AND 5.0), numbers ==
              closed form
  companion   output { layout companion; }  -> solver-dir U_PA == closed
              form over the 5 complete cycles; no phase-aligned dirs
  both        output { layout both; }   -> phase-aligned == closed form
              AND the companion files beside the solver's own fields
  cadence     output { cycleWriteInterval 2; }  -> dumps only after even
              cycles (cs 3.0 and 5.0 present; cs 4.0 absent)
  meanCadence output { meanWriteIntervalCycles 0; } -> _PA at the final
              dump only (absent at cs 4.0)
  noFinal     output { cycleWriteInterval 4; writeFinal off; } -> dump at
              cycle 4 only (cs 3.0); NO final dump (cs 5.0 absent)
  noProps     output { writeProperties off; } -> no Properties checkpoint
              anywhere (and the run is not restartable, warned at read)
  capture     fields ( U { writePerCycle on; } ) + output { } -> v2 naming:
              U_Bin beside U_PA in the phase-aligned dirs; no U_BinPA
              spelling anywhere; numbers == closed form
  deprecated  writeMode perCycle still works and warns

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


def _err(a, b):
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


def _body(fields="fields          (U);", output="output { }"):
    return f"""        {fields}
        period          1.0;
        nBins           8;
        discardCycles   0;
        phaseOrigin     0;
        {output}
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
"""


def _run(name, body, write_interval=1.0):
    os.makedirs(rv.RUNROOT, exist_ok=True)
    rundir, _ = rv.make_case(name, "bin_constant", DT, END_TIME,
                             write_interval, bbpa_body=body)
    assert rv.block_mesh(rundir).returncode == 0, name + ": blockMesh failed"
    assert rv.run_serial(rundir).returncode == 0, name + ": foamRun failed"
    return rundir


def _file_at(rundir, t, fname):
    for tv, full in rv._time_dirs(rundir).items():
        if abs(tv - t) < 1e-7:
            return os.path.isfile(os.path.join(full, fname))
    return False


def _check_closed(name, got, closed, errors,
                  fields=("PA", "PA_UU", "PAVariance")):
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
    closed6 = bo.closed_form_bin_constant(rv.A, rv.c_of_n, 6)
    closed5 = bo.closed_form_bin_constant(rv.A, rv.c_of_n, 5)

    # -- defaults: per-cycle cadence, numbers exact --------------------------
    rd = _run("out_defaults", _body())
    _check_closed("defaults", rv.read_bbpa_output(rd, FINAL_CS), closed6,
                  errors)
    if not _file_at(rd, 4.0, "U_PA"):
        errors.append("defaults: per-cycle cadence missing the cs=4.0 dump")

    # -- companion -----------------------------------------------------------
    rd = _run("out_companion", _body(output="output { layout companion; }"),
              write_interval=rv.BIN_DT)
    got = rv.read_companion_output(rd, FINAL_CS)
    for j in range(rv.I):
        if j not in got:
            errors.append(f"companion: bin {j} missing")
        elif _err(got[j]["PA"], closed5[j]["PA"]) >= TOL:
            errors.append(f"companion: bin {j} err "
                          f"{_err(got[j]['PA'], closed5[j]['PA']):.3e}")

    # -- both ----------------------------------------------------------------
    rd = _run("out_both", _body(output="output { layout both; }"),
              write_interval=rv.BIN_DT)
    _check_closed("both", rv.read_bbpa_output(rd, FINAL_CS), closed6, errors)

    # -- cycleWriteInterval 2 ------------------------------------------------
    rd = _run("out_cadence", _body(output="output { cycleWriteInterval 2; }"))
    for cs, expect in ((3.0, True), (4.0, False), (5.0, True)):
        if _file_at(rd, cs, "U_PA") != expect:
            errors.append(f"cadence: U_PA at cs={cs} expected {expect}")
    _check_closed("cadence", rv.read_bbpa_output(rd, FINAL_CS), closed6,
                  errors)

    # -- meanWriteIntervalCycles 0: _PA at the final dump only ---------------
    rd = _run("out_meanCadence",
              _body(output="output { meanWriteIntervalCycles 0; }"))
    if _file_at(rd, 4.0, "U_PA"):
        errors.append("meanCadence: U_PA written before the final dump")
    _check_closed("meanCadence", rv.read_bbpa_output(rd, FINAL_CS), closed6,
                  errors)

    # -- writeFinal off ------------------------------------------------------
    rd = _run("out_noFinal",
              _body(output="output { cycleWriteInterval 4; writeFinal off; }"))
    if not _file_at(rd, 3.0, "U_PA"):
        errors.append("noFinal: the cycle-4 dump (cs=3.0) is missing")
    if _file_at(rd, 5.0, "U_PA"):
        errors.append("noFinal: a final dump appeared despite writeFinal off")

    # -- writeProperties off -------------------------------------------------
    rd = _run("out_noProps", _body(output="output { writeProperties off; }"))
    props = glob.glob(os.path.join(rd, "[0-9]*", "uniform", "*Properties"))
    if props:
        errors.append(f"noProps: checkpoint written anyway: {props[:2]}")

    # -- per-field capture under the v2 naming -------------------------------
    rd = _run("out_capture",
              _body(fields="fields          ( U { writePerCycle on; } );"))
    if not _file_at(rd, FINAL_CS, "U_Bin"):
        errors.append("capture: U_Bin missing at the final cycle's grid")
    binpa = glob.glob(os.path.join(rd, "[0-9]*", "U_BinPA*"))
    if binpa:
        errors.append(f"capture: legacy _BinPA spelling appeared: {binpa[:2]}")
    _check_closed("capture", rv.read_bbpa_output(rd, FINAL_CS), closed6,
                  errors)

    # -- deprecated writeMode still works and warns --------------------------
    rd = _run("out_legacy", _body(output="writeMode       perCycle;"))
    log = open(os.path.join(rd, "log.foamRun")).read()
    if "'writeMode' is deprecated" not in log:
        errors.append("legacy: no deprecation warning for writeMode")
    _check_closed("legacy", rv.read_bbpa_output(rd, FINAL_CS), closed6,
                  errors)

    return errors


def test_output():
    errors = run_rows()
    assert not errors, "\n".join(errors)


if __name__ == "__main__":
    errors = run_rows()
    if errors:
        print("FAIL:")
        for e in errors:
            print("  ", e)
        sys.exit(1)
    print("output{} rows: PASS -- defaults/companion/both/cadence/meanCadence/")
    print("  noFinal/noProps/capture manifests and cadences exact, all")
    print(f"  numbers == closed form to {TOL:.0e}; writeMode warns deprecated")
    sys.exit(0)
